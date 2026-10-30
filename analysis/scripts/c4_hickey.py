"""G22-C4. Whole-section human references in two further organs: Hickey et al. 2023 CODEX atlas of the
human intestine (Dryad doi:10.5061/dryad.pk0p2ngrf; public domain). 8 donors x 8 intestinal sites,
small bowel (SB) and colon (CL); 25 annotated cell types.

Coordinates are image pixels; 0.37744 um per pixel is ASSUMED (20x CODEX on Keyence; not stated in README).
(1) c2_whole_section.analyse on every imaged region with 500-um virtual fields (whole-region truth);
(2) donor-level false-positive analysis: for each intestinal site, the 8 donors are split 4 vs 4 in all
    70 ways (no true group difference by construction of the random split); six tests per cell type,
    using 10 random fields per region (seed 20261061) and all fields.
usage: c4_hickey.py CELLS_TSV_GZ OUTPREFIX
"""
import sys, itertools, numpy as np, pandas as pd
sys.path.insert(0, '.')
import c2_whole_section as c2
from b2_simulation_da import tests, TESTS

PX = 0.37744
src, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(src, sep='\t', dtype={'tissue_unit': str, 'region': str}) if src.endswith('.gz') else pd.read_pickle(src)
d = d.dropna(subset=['cell_type', 'unique_region', 'x', 'y'])
cells = pd.DataFrame(dict(section=d.unique_region, donor=d.donor, organ=d.tissue.map({'SB': 'human small intestine', 'CL': 'human colon'}),
                          site=d.Tissue_location, x_um=d.x * PX, y_um=d.y * PX, label=d.cell_type))
print(cells.groupby('organ').agg(regions=('section', 'nunique'), donors=('donor', 'nunique'), cells=('label', 'size')))
if '--skip-c2' not in sys.argv:
    at, sa, it, pl = c2.analyse(cells.drop(columns='site'), F=500, MIN=100, out=out)
    ok = at[at.p_true >= 0.01]
    print(ok.groupby('organ').agg(D=('D', 'median'), q25=('D', lambda s: s.quantile(.25)), q75=('D', lambda s: s.quantile(.75)),
                                  combos=('D', 'size'), M=('M', 'median'), nbar=('nbar', 'median'), cover=('coverage_of_cells', 'median')))
    print(sa.groupby(['organ', 'k'])[['tv_fields', 'tv_cells']].median())
    print(it.groupby(['organ', 'k', 'method']).coverage.median().unstack())
    print(pl.groupby('organ')[['cells_assur', 'fovplan_assur', 'fovplan80_assur']].mean())

# donor-level false positives
rng = np.random.default_rng(20261061)
cells['field'] = pd.concat([c2.tile(g, 500) for _, g in cells.groupby('section')]).reindex(cells.index)
rows = []
for site, g in cells.groupby('site'):
    W = g.pivot_table(index=['donor', 'section', 'field'], columns='label', values='x_um', aggfunc='size', fill_value=0)
    W = W[W.sum(1) >= 100]
    tot = W.sum(1); prop = W.sum(0) / tot.sum()
    cats = list(prop[prop >= 0.01].index)
    donors = sorted(W.index.get_level_values('donor').unique())
    if len(donors) < 6:
        continue
    W = W.reset_index()
    labcols = [c for c in W.columns if c not in ('donor', 'section', 'field')]
    sub10 = W.iloc[rng.permutation(len(W))].groupby('section').head(10)
    for scen, F in [('all_fields', W), ('10_fields', sub10)]:
        nf = F[labcols].sum(axis=1).to_numpy(float)
        cl = F.donor.to_numpy()
        for si, s in enumerate(itertools.combinations(donors, len(donors) // 2)):
            if s[0] != donors[0]:
                continue  # each unordered split once
            grp = np.isin(cl, s).astype(int)
            for c in cats:
                ps = tests(F[c].to_numpy(float), nf, cl, grp)
                rows.append(dict(site=site, scenario=scen, split=si, category=c, **dict(zip(TESTS, ps))))
r = pd.DataFrame(rows)
r.to_csv(out + '_da_null.tsv.gz', sep='\t', index=False, compression='gzip')
f = r.groupby('scenario')[TESTS].apply(lambda x: (x < 0.05).mean())
print('donor-level null false-positive rates'); print(f.round(3))
f.to_csv(out + '_da_fpr.tsv', sep='\t')
print('sites', r.site.nunique(), 'splits per site', r.groupby('site').split.nunique().to_dict())
