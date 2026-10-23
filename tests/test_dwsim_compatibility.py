"""Independent reference values produced by the official DWSIM 10.2.11 engine."""
import json
from pathlib import Path
import numpy as np
import pytest
from dwsim_model import model
from engine import binary_diagram

ROWS=json.loads((Path(__file__).parents[1]/'data'/'dwsim_validation.json').read_text(encoding='utf-8'))['records']

@pytest.mark.parametrize('r',[r for r in ROWS if 'phase' in r])
def test_official_phase_properties(r):
    m=model(r['file']); z=np.array(r['z']); kind=r['phase'][0]
    np.testing.assert_allclose(np.exp(m.lnphi(r['T'],r['P'],z,kind)),r['phi'],rtol=1e-10,atol=1e-12)
    assert m.enthalpy(r['T'],r['P'],z,kind)==pytest.approx(r['H'],abs=.01)

@pytest.mark.parametrize('r',[r for r in ROWS if r.get('kind')=='flash_tight'])
def test_official_flash(r):
    s=model(r['file']).flash(T=r['T'],P=r['P'],zs=r['z'])
    # DWSIM's internal flash stops with small residuals even when the requested
    # outer tolerance is tightened further. The reference grid bounds this
    # numerical solver difference at <1e-6 in beta, <0.05 J/mol in mixture H.
    assert s.VF==pytest.approx(r['beta'],abs=1e-6)
    assert s.H()==pytest.approx(r['H'],abs=.05)
    if s.gas:
        np.testing.assert_allclose(s.gas.zs,r['y'],atol=1e-6,rtol=0)
    if s.liquid0:
        np.testing.assert_allclose(s.liquid0.zs,r['x'],atol=1e-6,rtol=0)

@pytest.mark.parametrize('r',[r for r in ROWS if r.get('kind') in ['bulle','rosee']])
def test_official_txy(r):
    vf=0 if r['kind']=='bulle' else 1
    t=model(r['file']).saturation(101325,[r['z0'],1-r['z0']],vf).T
    assert t==pytest.approx(r['T'],abs=1e-5)

@pytest.mark.parametrize('pressure',[5e4,2e5,5e5,1e6])
def test_txy_other_pressures(pressure):
    rows,errors=binary_diagram(pressure,11)
    assert len(rows)==11 and not errors
    assert all(r['rosee_K']>=r['bulle_K']-1e-6 for r in rows)
