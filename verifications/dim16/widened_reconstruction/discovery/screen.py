"""Numerical sensitivity of a published spherical SDP; no new certificate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json, math, hashlib
import numpy as np
from scipy.optimize import minimize, minimize_scalar
D=Path(__file__).resolve().parent
SOURCE=D/'factors.json'

class Witness:
    def __init__(self,path=SOURCE):
        data=json.loads(Path(path).read_text())
        assert data['dimension']==16 and len(data['blocks'])==34
        self.margins=np.array([float(b['lambda_decimal']) for b in data['blocks']])
        self.L=[np.zeros((b['size'],b['size'])) for b in data['blocks']]
        for b in data['blocks']:
            for i,j,v in b['L_entries']:self.L[b['index']][i,j]=float(v)
        assert len(self.L)==34 and all(self.margins>0)
        self.a=np.array([x[0,0]**2+self.margins[k] for k,x in enumerate(self.L[:16])])
        self.B=self.L[16].T@self.L[16]+self.margins[16]*np.eye(2)
        self.U=1+sum(self.a)+self.B[0,0]+sum((self.L[17]@np.ones(17))**2)+17*self.margins[17]
    @staticmethod
    def gegen(x,n,d):
        x=np.asarray(x);out=[np.ones_like(x)]
        if d:out.append(x)
        for k in range(1,d):out.append(((2*k+n-2)*x*out[-1]-k*out[-2])/(k+n-2))
        return np.stack(out,axis=-1)
    def kernel(self,u,v,t):
        u,v,t=np.broadcast_arrays(u,v,t);g=t-u*v;z=(1-u*u)*(1-v*v)
        q=[np.ones_like(u),g]
        for k in range(1,16):q.append(((2*k+13)*g*q[-1]-k*z*q[-2])/(k+13))
        value=np.zeros_like(u)+self.B[1,1]
        for k in range(17):
            pu=self.gegen(u,16+2*k,16-k);pv=self.gegen(v,16+2*k,16-k)
            L=self.L[17+k]
            dot=np.sum((pu@L.T)*(pv@L.T),axis=-1)+self.margins[17+k]*np.sum(pu*pv,axis=-1)
            value+=q[k]*dot
        return value
    def cycle(self,u,v,t):
        return self.kernel(u,v,t)+self.kernel(u,t,v)+self.kernel(v,t,u)
    def pair(self,t):
        t=np.asarray(t);p=self.gegen(t,16,16)
        cross=(p@self.L[17].T)@(self.L[17]@np.ones(17))+self.margins[17]*p.sum(axis=-1)+self.B[0,1]
        return 1+p[:,1:]@self.a+2*cross+self.kernel(t,t,np.ones_like(t))

def triple_grid(c):
    aa=np.unique(np.r_[np.linspace(-1,c,55),[-.5,0,.5,c]])
    us=[];vs=[];ts=[]
    for u in aa:
        for v in aa[aa>=u]:
            r=np.sqrt(max(0,(1-u*u)*(1-v*v)))
            lo=max(v,u*v-r);hi=min(c,u*v+r)
            if lo<=hi:
                ts.extend(np.linspace(lo,hi,35));us.extend([u]*35);vs.extend([v]*35)
    return np.array(us),np.array(vs),np.array(ts)

def run():
    w=Witness();out=[]
    for c in [.5,.5005,.501,.503,.505,.51]:
        tx=np.linspace(-1,c,10001);px=w.pair(tx);ip=np.argmax(px);pairmax=float(px[ip]);pairpoint=float(tx[ip])
        peaks=np.r_[ip,np.where((px[1:-1]>=px[:-2])&(px[1:-1]>=px[2:]))[0]+1]
        for i in peaks:
            a,b=tx[max(0,i-2)],tx[min(len(tx)-1,i+2)]
            r=minimize_scalar(lambda t:-w.pair(np.array([t]))[0],bounds=(a,b),method='bounded',options={'xatol':1e-14})
            if -r.fun>pairmax:pairmax=float(-r.fun);pairpoint=float(r.x)
        u,v,t=triple_grid(c);vals=w.cycle(u,v,t);ids=np.argsort(vals)[-12:];triplemax=float(vals[ids[-1]]);point=[float(u[ids[-1]]),float(v[ids[-1]]),float(t[ids[-1]])]
        # The main grid is Gram-feasible. Local SLSQP retains determinant and cap.
        for i in ids:
            def fn(x):return -float(w.cycle(*map(np.asarray,x)))
            def det(x):a,b,d=x;return 1+2*a*b*d-a*a-b*b-d*d
            r=minimize(fn,[u[i],v[i],t[i]],method='SLSQP',bounds=[(-1,c)]*3,
                constraints=[{'type':'ineq','fun':det}],options={'maxiter':80,'ftol':1e-11})
            if det(r.x)>=-1e-9 and -r.fun>triplemax:triplemax=float(-r.fun);point=r.x.tolist()
        # Only a CONDITIONAL repair cost, using an illustrative LP guard B=9000.
        guardB=9000.;p=max(0,pairmax);tt=max(0,triplemax)
        row=dict(cap=c,radius_ratio=1/math.sqrt(2*(1-c)),contact_objective=w.U,
            pair_maximum=pairmax,pair_point=pairpoint,cyclic_maximum=triplemax,cyclic_point=point,
            sampled_triples=len(vals),hypothetical_guard_B=guardB,
            hypothetical_repaired_count=w.U+(guardB-1)*(p+(2*guardB-1)*tt/3),certified=False)
        out.append(row);print(row,flush=True)
        (D/'SCREEN.json').write_text(json.dumps(dict(scope='Numerical source sensitivity only. Contact identity not independently replayed; expanded maxima are not proved; B=9000 is illustrative, not a supplied certificate.',source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),results=out),indent=2)+'\n')
if __name__=='__main__':run()
