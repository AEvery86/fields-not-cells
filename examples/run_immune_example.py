"""Run the supplied count-only example without installing or reading numeric metadata."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import fovplan

def validate(counts):
    x=np.asarray(counts,dtype=float)
    if x.ndim!=2 or x.shape[0]<2 or x.shape[1]<2:
        raise ValueError('Provide at least two pilot fields and a focal category plus its complement.')
    if not np.isfinite(x).all() or (x<0).any() or (x!=np.floor(x)).any() or (x.sum(1)<=0).any():
        raise ValueError('Counts must be finite nonnegative integers, with positive totals in every field.')
    if ((x.sum(0)==0)|(x.sum(0)==x.sum())).any():
        raise ValueError('Boundary pilot proportion: no identifiable precision plan; expand the pilot.')
    p=x.sum(0)/x.sum(); resid=x-x.sum(1)[:,None]*p
    if (resid.var(0,ddof=1)==0).any():
        raise ValueError('Zero pilot residual variance cannot establish adequate precision; expand the pilot.')
    return counts

def main():
    counts=validate(pd.read_csv(ROOT/'examples/immune_pilot_G24.tsv',sep='\t'))
    for name,target,relative in [('relative',.25,True),('absolute',.02,False)]:
        result=fovplan.plan(counts,target,relative=relative,assurance=.9)
        expected=pd.read_csv(ROOT/f'examples/immune_pilot_G24_{name}_output.tsv',sep='\t')
        pd.testing.assert_frame_equal(result,expected,check_dtype=False,rtol=1e-10,atol=1e-12)
        print(name); print(result.to_string(index=False))
    print('Both expected-output checks passed.')
if __name__=='__main__': main()
