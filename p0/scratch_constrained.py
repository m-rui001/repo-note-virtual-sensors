import numpy as np, math
from scipy.linalg import expm
from scipy.optimize import minimize
np.set_printoptions(precision=6, suppress=True)

def sym(M): return 0.5*(M+M.T)
def logdet(M):
    s,l=np.linalg.slogdet(sym(M)); return l if s>0 else -np.inf
def psd_sqrt(M):
    w,V=np.linalg.eigh(sym(M)); return V@np.diag(np.sqrt(np.maximum(w,0)))@V.T
def min_eig(M): return float(np.linalg.eigvalsh(sym(M)).min())

def eps_constrained(Sigma, Pi, C, n_starts=24, seed=0, verbose=False):
    d=Sigma.shape[0]; sq=psd_sqrt(Sigma)
    W=sym(sq@Pi@sq)
    tri=np.triu_indices(d)
    trg=-2.0*C*math.log(2.0)
    def unpack(v):
        S=np.zeros((d,d)); S[tri]=v; S=sym(S)
        S=S-(np.trace(S)-trg)/d*np.eye(d)     # 精确满足 tr(S)=trg
        return S
    def obj(v,pen):
        S=unpack(v)
        ev=np.linalg.eigvalsh(S)
        f=float(np.trace(W@expm(S)))
        return f+pen*float(np.sum(np.maximum(ev,0.0)**2))
    def obj_only(v): return obj(v,0.0)
    rng=np.random.default_rng(seed)
    best=None;bestv=None
    # 起点：解析解（若可行）+ 各向同性 + 随机
    starts=[]
    if logdet(Pi)>-np.inf:
        try:
            c=math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d)
            Se=c*np.linalg.inv(Pi)
            M=np.linalg.solve(sq, np.linalg.solve(sq, Se).T).T if False else np.linalg.inv(sq)@Se@np.linalg.inv(sq)
            M=sym(M)
            w,V=np.linalg.eigh(M); w=np.maximum(w,1e-12)
            from scipy.linalg import logm
            S0=sym(logm(V@np.diag(w)@V.T).real)
            starts.append(S0[tri])
        except Exception as e: pass
    starts.append(np.zeros(d*(d+1)//2))
    for _ in range(n_starts):
        v=rng.standard_normal(d*(d+1)//2)
        starts.append(v/np.maximum(np.abs(v).max(),1e-9)*1.0)
    for pen in (1e2, 1e4, 1e6, 1e8):
        for s0 in starts:
            res=minimize(lambda v: obj(v,pen), s0, method='Nelder-Mead',
                         options=dict(maxiter=4000, maxfev=8000, xatol=1e-10, fatol=1e-14))
            val=obj_only(res.x)
            feasible = min_eig(unpack(res.x))<=1e-7
            if feasible and (best is None or val<best):
                best=val; bestv=res.x.copy()
    if bestv is None:
        return None, np.inf, False
    Se=sym(sq@expm(unpack(bestv))@sq)
    return Se, best, (min_eig(unpack(bestv))<=1e-6)

# --- 测试 1：可行情形应与解析解一致
d=2
Sigma=np.diag([4.0,1.0]); th=math.radians(30)
Rm=np.array([[math.cos(th),-math.sin(th)],[math.sin(th),math.cos(th)]])
Pi=Rm@np.diag([16.0,1.0])@Rm.T
C=6.0
c=math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d)
Se_an=c*np.linalg.inv(Pi)
dJ=lambda P,S: float(np.trace(P@S))
print("测试1 解析可行? ", min_eig(Sigma-Se_an)>=-1e-9, " 解析 ΔJ=", dJ(Pi,Se_an))
Se_n,val,feas=eps_constrained(Sigma,Pi,C)
print("     数值 ΔJ=", dJ(Pi,Se_n), " feasible=",feas, " 相对差=", abs(dJ(Pi,Se_n)/dJ(Pi,Se_an)-1))
print("     AM-GM 下界=", d*math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d))

# --- 测试 2：紧预算（解析不可行）
C=1.5
c=math.exp((logdet(Pi)+logdet(Sigma))/d)*2.0**(-2.0*C/d)
Se_an=c*np.linalg.inv(Pi)
print("\n测试2 (C=1.5) 解析可行? ", min_eig(Sigma-Se_an)>=-1e-9, " 解析ΔJ=",dJ(Pi,Se_an))
Se_n,val,feas=eps_constrained(Sigma,Pi,C)
print("     数值 ΔJ=", dJ(Pi,Se_n), " feasible=",feas)
print("     Σ_e=\n", Se_n, " ⪯Σ?", min_eig(Sigma-Se_n)>=-1e-9)
