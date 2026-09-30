import numpy as np, math
from scipy.linalg import expm
from scipy.optimize import minimize
def sym(M): return 0.5*(M+M.T)
def logdet(M):
    s,l=np.linalg.slogdet(sym(M)); return l if s>0 else -np.inf
def psd_sqrt(M):
    w,V=np.linalg.eigh(sym(M)); return V@np.diag(np.sqrt(np.maximum(w,0)))@V.T
def eps_oracle(Sigma,Pi,C,n_starts=16,seed=0):
    """参数化保证 detΣ_e=detΣ·2^{-2C} 且 Σ_e⪯Σ 精确成立"""
    d=Sigma.shape[0]; sq=psd_sqrt(Sigma); W=sym(sq@Pi@sq)
    npair=d*(d-1)//2; trg=2.0*C*math.log(2.0)
    def build(p):
        v=p[:npair]; z=p[npair:]
        Sk=np.zeros((d,d)); idx=0
        for i in range(d):
            for j in range(i+1,d):
                Sk[i,j]=v[idx]; Sk[j,i]=-v[idx]; idx+=1
        U=expm(Sk)
        zz=np.concatenate([z,[0.0]]); zz=zz-zz.max()
        s=trg*np.exp(zz)/np.exp(zz).sum()
        S=-U@np.diag(s)@U.T
        return sym(sq@expm(S)@sq), S
    def obj(p):
        Se,_=build(p); return float(np.trace(Pi@Se))
    rng=np.random.default_rng(seed); best=None;bx=None
    starts=[np.zeros(npair+d-1)]
    for _ in range(n_starts):
        starts.append(rng.standard_normal(npair+d-1))
    for s0 in starts:
        r=minimize(obj,s0,method='Nelder-Mead',options=dict(maxiter=6000,maxfev=12000,xatol=1e-11,fatol=1e-15))
        if best is None or r.fun<best: best,bx=r.fun,r.x
    # 局部精修
    for _ in range(3):
        r=minimize(obj,bx,method='Nelder-Mead',options=dict(maxiter=20000,maxfev=40000,xatol=1e-13,fatol=1e-17))
        if r.fun<best: best,bx=r.fun,r.x
    Se,S=build(bx)
    return Se,best,float(np.linalg.eigvalsh(S).max())

for d,C in ((2,6.0),(2,1.5),(4,3.0)):
    rng=np.random.default_rng(3)
    A=rng.standard_normal((d,d)); Sigma=sym(A@A.T)+np.eye(d)
    B=rng.standard_normal((d,d)); Pi=sym(B@B.T)+0.3*np.eye(d)
    an=None
    if logdet(Pi)>-np.inf:
        c=math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d); an=c*np.linalg.inv(Pi)
        feas=float(np.linalg.eigvalsh(Sigma-an).min())>=-1e-9
    else: feas=False
    Se,val,evmax=eps_oracle(Sigma,Pi,C)
    bound=d*math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d)
    print(f"d={d} C={C}: oracle={val:.10f}  AM-GM下界={bound:.10f}  解析可行={feas}"
          + (f"  解析={float(np.trace(Pi@an)):.10f} 相对差={abs(val/float(np.trace(Pi@an))-1):.2e}" if feas else "")
          + f"  max_eig(S)={evmax:.2e} (应<=0)  Σ_e⪯Σ? {float(np.linalg.eigvalsh(Sigma-Se).min())>=-1e-8}"
          + f"  率误差={abs(0.5*(logdet(Sigma)-logdet(Se))/math.log(2)-C):.2e}")
