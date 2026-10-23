"""PRSV2-VL compatibility equations and PT flash.

Equations adapted from DWSIM, Copyright 2012-2024 Daniel Wagner O. de
Medeiros, PRSV2-VL.vb and PengRobinsonStryjekVera2VL.vb (GPL-3.0-or-later).
This module is distributed under the same license; see LICENSE.
Units: K, Pa, mol, J. Saved chemical names are intentionally preserved.
"""
from functools import lru_cache
from types import SimpleNamespace
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares
from reference import load_reference

R = 8.314
SQRT2 = np.sqrt(2.)
PARAMS = {'Water':(-.06635,.0199,.443), 'Methanol':(-.16816,-1.34,.588),
          'Propane':(.03136,.2757,.447)}
# DWSIM looks up saved names literally. Its table has Butane/Pentane/Hexane,
# not N-butane/N-pentane/N-hexane; these therefore use its PR fallback.

class Phase:
    def __init__(self, model, T, P, zs, kind):
        self.model, self.T, self.P, self.kind = model,T,P,kind
        self.zs = np.array(zs,dtype=float)
    def phis(self):
        return np.exp(self.model.lnphi(self.T,self.P,self.zs,self.kind))
    def H(self):
        return self.model.enthalpy(self.T,self.P,self.zs,self.kind)

class State:
    def __init__(self, model, T, P, zs, beta, x, y):
        self.T,self.P,self.zs,self.VF = T,P,list(zs),float(beta)
        self.gas = Phase(model,T,P,y,'V') if beta > 1e-12 else None
        self.liquid0 = Phase(model,T,P,x,'L') if beta < 1-1e-12 else None
    def H(self):
        return ((1-self.VF)*self.liquid0.H() if self.liquid0 else 0.) + (self.VF*self.gas.H() if self.gas else 0.)

class Model:
    def __init__(self, filename):
        self.filename = filename
        self.compounds = load_reference(filename)['compounds']
        self.names = [c['Name'] for c in self.compounds]
        self.Tcs = np.array([float(c['Critical_Temperature']) for c in self.compounds])
        self.Pcs = np.array([float(c['Critical_Pressure']) for c in self.compounds])
        self.omegas = np.array([float(c['Acentric_Factor']) for c in self.compounds])
        self.MWs = np.array([float(c['Molar_Weight']) for c in self.compounds])
        self.CASs = [c['CAS_Number'] for c in self.compounds]
        self.kappas = np.array([PARAMS.get(n,(0.,0.,0.)) for n in self.names])
        self.k1 = np.zeros((len(self.names),len(self.names)))
        self.k2 = np.zeros_like(self.k1)
        if 'Water' in self.names and 'Methanol' in self.names:
            i,j = self.names.index('Water'),self.names.index('Methanol')
            self.k1[i,j],self.k1[j,i] = -.0839,-.0968
            self.k2 = self.k1.T.copy()
        self.bi = .0778*R*self.Tcs/self.Pcs
        self.cp_coeffs = np.array([[float(c['Ideal_Gas_Heat_Capacity_Const_'+v]) for v in 'ABCDE'] for c in self.compounds])
        if any(c['IdealgasCpEquation'] != '16' for c in self.compounds):
            raise ValueError('Corrélation Cp non prise en charge dans ce fichier.')

    def attraction(self,T):
        tr = T/self.Tcs
        sq = np.sqrt(tr)
        w = self.omegas
        k1,k2,k3 = self.kappas.T
        prsv = .378893 + 1.4897153*w - .17131848*w*w + .0196544*w**3
        prsv += (k1+k2*(k3-tr)*(1-sq))*(1+sq)*(.7-tr)
        fallback = np.where(w<=.491,.37464+1.5422*w-.26992*w*w,
                            .379642+1.48503*w-.164423*w*w+.016666*w**3)
        ci = np.where(k1*k2*k3!=0,prsv,fallback)
        ai = .45724*(1+ci*(1-sq))**2*R*R*self.Tcs**2/self.Pcs
        return ai,ci

    def parameters(self,T,P,z):
        ai,ci = self.attraction(T)
        den = z[:,None]*self.k1 + z[None,:]*self.k2
        term = np.divide(self.k1*self.k2,den,out=np.zeros_like(den),where=den!=0)
        aij = np.sqrt(ai[:,None]*ai[None,:])*(1-term)
        am = float(z @ aij @ z)
        bm = float(z @ self.bi)
        A,B = am*P/(R*T)**2,bm*P/(R*T)
        roots = np.roots([1.,B-1.,A-3*B*B-2*B,-A*B+B*B+B**3])
        roots = sorted(float(v.real) for v in roots if abs(v.imag)<1e-8 and v.real>B)
        if not roots:
            raise ValueError('Aucune racine physique de l’équation PRSV2-VL.')
        return ai,ci,term,aij,am,bm,A,B,roots

    def lnphi(self,T,P,z,kind):
        z = np.asarray(z,dtype=float)
        ai,ci,term,aij,am,bm,A,B,roots = self.parameters(T,P,z)
        Z = roots[0] if kind=='L' else roots[-1]
        # DWSIM CalcLnFug: correction to the attraction partial derivatives.
        sums = aij @ z
        den = (z[:,None]*self.k1 + z[None,:]*self.k1.T)**2
        numer = z[:,None]*z[None,:]*np.sqrt(ai[:,None]*ai[None,:])*self.k1*self.k2
        numer *= ((1-z[:,None])*self.k1-z[None,:]*self.k2)
        sums += np.divide(numer,den,out=np.zeros_like(den),where=den!=0).sum(axis=1)
        logratio = np.log((Z+(1+SQRT2)*B)/(Z+(1-SQRT2)*B))
        return self.bi/bm*(Z-1)-np.log(Z-B) - A*(2*sums/am-self.bi/bm)*logratio/(2*SQRT2*B)

    @lru_cache(maxsize=512)
    def ideal_integrals(self,T):
        vals = []
        for a,b,c,d,e in self.cp_coeffs:
            # ChemSep equation 16 is J/(kmol K), converted here to J/(mol K).
            vals.append(quad(lambda t:(a+np.exp(b/t+c+d*t+e*t*t))/1000,298.15,T,epsabs=1e-7)[0])
        return np.array(vals)

    def enthalpy(self,T,P,z,kind):
        z = np.asarray(z,dtype=float)
        ai,ci,term,aij,am,bm,A,B,roots = self.parameters(T,P,z)
        Z = roots[0] if kind=='L' else roots[-1]
        # Historical DWSIM H_PR_MIX derivative, retained for saved-model parity.
        derivative = -R/2*np.sqrt(.45724/T)*np.sum(z[:,None]*z[None,:]*(1-term)*(
            ci[None,:]*np.sqrt(ai[:,None]*self.Tcs[None,:]/self.Pcs[None,:])
            +ci[:,None]*np.sqrt(ai[None,:]*self.Tcs[:,None]/self.Pcs[:,None])))
        logratio = np.log((Z+(1+SQRT2)*B)/(Z+(1-SQRT2)*B))
        residual = R*T*(Z-1)+(T*derivative-am)/(2*SQRT2*bm)*logratio
        return float(z @ self.ideal_integrals(float(T)) + residual)

    def wilson(self,T,P):
        return self.Pcs/P*np.exp(5.373*(1+self.omegas)*(1-self.Tcs/T))

    @staticmethod
    def rr(z,K,previous=.5):
        delta=K-1
        f0=float(z @ delta)
        f1=float(z @ (delta/K))
        if max(abs(delta))<1e-9:
            return previous
        if f0<=0: return 0.
        if f1>=0: return 1.
        return brentq(lambda b:float(z @ (delta/(1+b*delta))),0.,1.,xtol=1e-14)

    def flash(self,T=None,P=None,zs=None,VF=None):
        if VF is not None:
            return self.saturation(P,zs,VF)
        z=np.array(zs,dtype=float)
        K=self.wilson(T,P)
        beta=.5
        for iteration in range(500):
            beta=self.rr(z,K,beta)
            x=z/(1+beta*(K-1)); x/=sum(x)
            y=K*x; y/=sum(y)
            logK=self.lnphi(T,P,x,'L')-self.lnphi(T,P,y,'V')
            err=max(abs(logK-np.log(K)))
            if err<2e-10:
                if max(abs(logK))<1e-7:
                    roots=self.parameters(T,P,z)[-1]
                    beta=1. if roots[-1]>.5 else 0.
                    x=y=z.copy()
                return State(self,T,P,z,beta,x,y)
            damping=.5 if iteration>30 else 1.
            K=np.exp(np.clip(damping*logK+(1-damping)*np.log(K),-50,50))
        raise ValueError('Flash PT non convergé après 500 itérations.')

    def saturation(self,P,zs,VF):
        z=np.array(zs,dtype=float)
        active=np.flatnonzero(z>0)
        if len(active)==1:
            i=int(active[0])
            cp=self.compounds[i]
            if cp['VaporPressureEquation']!='101':
                raise ValueError('Corrélation de pression de vapeur non prise en charge.')
            a,b,c,d,e=[float(cp['Vapor_Pressure_Constant_'+v]) for v in 'ABCDE']
            # DWSIM's pure-component shortcut uses the saved vapor-pressure
            # correlation, not equality of EOS fugacities.
            def f(T): return a+b/T+c*np.log(T)+d*T**e-np.log(P)
            low=float(cp['Vapor_Pressure_TMIN']); high=min(float(cp['Vapor_Pressure_TMAX']),self.Tcs[i])
            T=brentq(f,low,high,xtol=1e-8)
            return SimpleNamespace(T=T)
        if len(z)!=2 or len(active)!=2:
            raise ValueError('Diagramme disponible uniquement pour le mélange binaire.')
        # Two independent fugacity equations; unknown T and other-phase fraction.
        def equations(u):
            t=u[0]; c=1/(1+np.exp(-u[1]))
            other=np.array([c,1-c])
            x,y=(z,other) if VF==0 else (other,z)
            return np.log(x)+self.lnphi(t,P,x,'L')-np.log(y)-self.lnphi(t,P,y,'V')
        pure_t=[self.saturation(P,np.eye(2)[i],VF).T for i in range(2)]
        guess=float(z @ pure_t)
        if VF==0:
            q=z*self.wilson(guess,P); q/=sum(q)
        else:
            q=z/self.wilson(guess,P); q/=sum(q)
        initial=[np.clip(guess,250,450),np.log(q[0]/q[1])]
        result=least_squares(equations,initial,bounds=([200.,-30.],[600.,30.]),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=150)
        if not result.success or max(abs(equations(result.x)))>1e-7:
            raise ValueError('Point de saturation non convergé.')
        c=1/(1+np.exp(-result.x[1])); other=np.array([c,1-c])
        x,y=(z,other) if VF==0 else (other,z)
        zl=self.parameters(result.x[0],P,x)[-1][0]
        zv=self.parameters(result.x[0],P,y)[-1][-1]
        if zv-zl<1e-4:
            raise ValueError('Solution de saturation triviale : phases indiscernables.')
        return SimpleNamespace(T=float(result.x[0]))

@lru_cache(maxsize=4)
def model(filename):
    return Model(filename)
