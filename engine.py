"""SI internally: K, Pa, mol/s, J/mol. Steady state, one liquid + vapor."""
from functools import lru_cache
import math
import numpy as np
from reference import load_reference
from dwsim_model import model

NAMES = {'Water': 'Eau', 'Methanol': 'Méthanol', 'Propane': 'Propane',
         'N-butane': 'n-Butane', 'N-pentane': 'n-Pentane', 'N-hexane': 'n-Hexane'}

@lru_cache(maxsize=4)
def make_model(filename):
    engine = model(filename)
    return engine, engine, {}

def validate(T, P, F, zs):
    if not all(math.isfinite(v) for v in [T, P, F, *zs]):
        raise ValueError('Toutes les valeurs doivent être finies.')
    if T <= 0 or P <= 0 or F <= 0:
        raise ValueError('Température absolue, pression et débit doivent être positifs.')
    if min(zs) < 0 or abs(sum(zs)-1) > 1e-7:
        raise ValueError('Les fractions molaires doivent être positives et leur somme égale à 1.')

def simulate(filename, Tin, Tout, P, F, zs):
    validate(Tin, P, F, zs)
    validate(Tout, P, F, zs)
    flash, constants, _ = make_model(filename)
    if len(zs) != len(constants.CASs):
        raise ValueError('Nombre de fractions incorrect.')
    feed = flash.flash(T=Tin, P=P, zs=zs)
    out = flash.flash(T=Tout, P=P, zs=zs)
    beta = out.VF
    if not math.isfinite(beta) or not -1e-9 <= beta <= 1+1e-9:
        raise ValueError('Le calcul de phases n’a pas convergé.')
    beta = min(1., max(0., beta))
    vap = out.gas if beta > 1e-10 else None
    liq = out.liquid0 if beta < 1-1e-10 else None
    x = np.array(liq.zs) if liq is not None else np.zeros(len(zs))
    y = np.array(vap.zs) if vap is not None else np.zeros(len(zs))
    L, V = F*(1-beta), F*beta
    component_error = F*np.array(zs) - L*x - V*y
    hin = F*feed.H()/1000
    hl = L*liq.H()/1000 if liq is not None else 0.
    hv = V*vap.H()/1000 if vap is not None else 0.
    duty = F*(out.H()-feed.H())/1000
    if np.max(np.abs(component_error))/F > 1e-6:
        raise ValueError('Bilan matière non fermé : résultats refusés.')
    if liq is not None and vap is not None:
        fug_error = max(abs(np.log((x[i]*liq.phis()[i])/(y[i]*vap.phis()[i])))
                        for i, z in enumerate(zs) if z > 1e-12)
        if fug_error > 1e-5:
            raise ValueError('Équilibre de fugacité non convergé : résultats refusés.')
    else:
        fug_error = 0.
    return dict(beta=beta, L=L, V=V, x=x, y=y, hin=hin, hl=hl, hv=hv, Q=duty,
                energy_error=hin+duty-hl-hv, component_error=component_error,
                fugacity_error=fug_error, feed=feed, outlet=out,
                mass_in=F*sum(z*m for z,m in zip(zs,constants.MWs))/1000,
                mass_liq=L*sum(z*m for z,m in zip(x,constants.MWs))/1000,
                mass_vap=V*sum(z*m for z,m in zip(y,constants.MWs))/1000)

def binary_diagram(P, points=51):
    flash, _, _ = make_model('eau_methanol.dwxml')
    rows, failures = [], []
    for z in np.linspace(0., 1., points):
        try:
            bubble = flash.flash(P=P, VF=0, zs=[float(z), float(1-z)])
            dew = flash.flash(P=P, VF=1, zs=[float(z), float(1-z)])
            if not all(math.isfinite(t) for t in [bubble.T, dew.T]):
                raise ValueError('Température non finie')
            rows.append(dict(fraction_eau=float(z), bulle_K=bubble.T, rosee_K=dew.T))
        except Exception as exc:
            failures.append((float(z), str(exc)))
    return rows, failures
