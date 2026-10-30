import numpy as np,itertools
from scipy.stats import t as tdist
from interval_methods import wilson,region_interval
EXACT={k:np.array(list(itertools.product(range(k),repeat=k)),dtype=np.int16) for k in [3,5]}

def fieller_coefficients(y,n,M):
 k=len(y);scale=np.mean(n);yy=y/scale;xx=n/scale
 cov=np.cov(np.stack([yy,xx]),ddof=1)*(1-k/M)/k;crit=tdist.ppf(.975,k-1)**2
 return np.array([1-crit*cov[1,1],-2*(yy.mean()-crit*cov[0,1]),yy.mean()**2-crit*cov[0,0]])

def polynomial_set(coef):
 a,b,c=coef;tol=1e-12*max(1,np.max(abs(coef)));roots=[]
 if abs(a)<tol:
  if abs(b)>tol:roots=[-c/b]
 else:
  disc=b*b-4*a*c
  if disc>=-tol:roots=[(-b-np.sqrt(max(disc,0)))/(2*a),(-b+np.sqrt(max(disc,0)))/(2*a)]
 roots=[float(np.clip(x,0,1)) for x in roots if -1e-12<=x<=1+1e-12]
 cuts=sorted(set([0.,1.]+roots));segs=[]
 for l,h in zip(cuts[:-1],cuts[1:]):
  mid=(l+h)/2
  if np.polyval(coef,mid)<=tol:segs.append((l,h))
 for x in roots:
  if np.polyval(coef,x)<=tol and not any(l-1e-12<=x<=h+1e-12 for l,h in segs):segs.append((x,x))
 if not segs and abs(a)+abs(b)<tol and c<=tol:segs=[(0.,1.)]
 segs.sort();merged=[]
 for l,h in segs:
  if merged and l<=merged[-1][1]+1e-12:merged[-1]=(merged[-1][0],max(merged[-1][1],h))
  else:merged.append((l,h))
 return merged

def percentile_interval(y,n,rng):
 k=len(y);idx=EXACT[k] if k in EXACT else rng.integers(0,k,size=(999,k))
 vals=np.sort(y[idx].sum(1)/n[idx].sum(1));N=len(vals)
 return float(vals[int(np.ceil(.025*N))-1]),float(vals[int(np.ceil(.975*N))-1])

def one_draw(arr,ids,M,truth,outer,boot):
 ss=arr[ids];y=ss[:,0];n=ss.sum(1);k=len(ids);nobs=int(n.sum());phat=float(y.sum()/nobs)
 wl,wh=wilson(float(y.sum()),nobs);rl,rh=region_interval(ss,M)
 if k==M:
  fieller=[(phat,phat)];bl=bh=phat
 else:
  fieller=polynomial_set(fieller_coefficients(y,n,M));bl,bh=percentile_interval(y,n,boot)
 assert any(l-1e-10<=phat<=h+1e-10 for l,h in fieller),(ss,fieller)
 factor=np.sqrt(1-k/M);fl=max(0,phat+factor*(bl-phat));fh=min(1,phat+factor*(bh-phat))
 N=int(arr.sum());K=int(arr[:,0].sum());iid=outer.hypergeometric(K,N-K,nobs);il,ih=wilson(iid,nobs)
 methods={'cell_Wilson':[(wl,wh)],'region_ratio_t':[(rl[0],rh[0])],'Fieller':fieller,'region_percentile':[(bl,bh)],'percentile_FPC':[(fl,fh)],'random_cell_Wilson':[(il,ih)]}
 out=[]
 for method,segs in methods.items():
  width=sum(h-l for l,h in segs);assert len(segs) in [1,2]
  out.append(dict(method=method,estimate=phat if method!='random_cell_Wilson' else iid/nobs,lower1=segs[0][0],upper1=segs[0][1],lower2=segs[1][0] if len(segs)>1 else np.nan,upper2=segs[1][1] if len(segs)>1 else np.nan,components=len(segs),width=width,full_set=bool(width>=1-1e-12),covered=any(l-1e-12<=truth<=h+1e-12 for l,h in segs),n_cells=nobs))
 return out
