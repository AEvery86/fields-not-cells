"""G21-C3. Whole-section truth in the two mouse Xenium kidney sections (ShamR, Day14R).
Applies c2_whole_section.analyse with 500- and 250-unit virtual fields (stored coordinate units)."""
import pandas as pd
import c2_whole_section as c
cells = []
for s in ['ShamR', 'Day14R']:
    d = pd.read_csv(f'../data/mouse_example_{s}.tsv.gz', sep='\t')
    cells.append(pd.DataFrame(dict(section=s, donor=s, organ='mouse kidney', x_um=d.x, y_um=d.y, label=d.celltype)))
cells = pd.concat(cells)
for F, MIN in [(500, 100), (250, 30)]:
    at, sa, it, pl = c.analyse(cells, F, MIN, out=f'../results/C3_mouse_F{F}')
    print('F', F); print(at[at.p_true >= .01].groupby('section').agg(D=('D', 'median'), M=('M', 'first'), cover=('coverage_of_cells', 'first'), nbar=('nbar', 'first')))
    print(sa.round(4).to_string())
    print(it.groupby(['k', 'method']).coverage.median().unstack().round(3))
    print(pl[['cells_assur', 'fovplan_assur', 'fovplan80_assur']].mean().round(3).to_dict(), pl[['cells_k', 'fovplan_k', 'fovplan80_k']].median().to_dict())
