"""Summary tables (Supplementary Tables S10-S13) and Supplementary Figure S13 for A9-A12."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = '../results/a9_12/'
SIX = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']
grp = lambda s: s.replace({'human small intestine': 'Human intestine', 'human colon': 'Human intestine', 'mouse kidney': 'Mouse kidney'})

# ---- A9 preferential selection (S10)
c = pd.read_csv(R + 'A9_composition.tsv.gz', sep='\t'); v = pd.read_csv(R + 'A9_coverage.tsv.gz', sep='\t'); p = pd.read_csv(R + 'A9_planning.tsv.gz', sep='\t')
SCH = ['random', 'systematic', 'interior', 'dense', 'enriched']
NAME = {'random': 'Simple random', 'systematic': 'Systematic (serpentine)', 'interior': 'Interior fields only', 'dense': 'Cell-dense fields preferred', 'enriched': 'Structure-enriched'}
s10 = pd.DataFrame(index=SCH)
s10['Median TV, 3 fields'] = c[c.k == 3].groupby('scheme').tv.median()
s10['Median TV, 10 fields'] = c[c.k == 10].groupby('scheme').tv.median()
v['ab'] = v.rel_bias.abs() * 100
s10['|Relative bias| %, other types (3 fields)'] = v[(v.k == 3) & ~v.focal].groupby('scheme').ab.median()
s10['|Relative bias| %, enriched type (3 fields)'] = v[(v.k == 3) & v.focal].groupby('scheme').ab.median()
s10['Coverage %, other types (3 fields)'] = 100 * v[(v.k == 3) & ~v.focal].groupby('scheme').coverage.median()
s10['Coverage %, enriched type (3 fields)'] = 100 * v[(v.k == 3) & v.focal].groupby('scheme').coverage.median()
s10['fovplan 90% within ±25%, other types'] = 100 * p[~p.focal].groupby('scheme').share.mean()
s10['fovplan 90% within ±25%, enriched type'] = 100 * p[p.focal].groupby('scheme').share.mean()
s10.index = [NAME[i] for i in s10.index]
s10 = s10.round(3); s10.to_csv(R + 'S10_preferential_selection.tsv', sep='\t')
print(s10.to_string())

# ---- A10 planner sensitivity (S11)
a = pd.read_csv(R + 'A10.tsv.gz', sep='\t')
M = ['cellcount', 'bootstrap', 'fovplan', 'fovplan80', 'fovplan90']
ORD = ['rel 0.1 / 5 pilot', 'rel 0.5 / 5 pilot', 'abs 0.02 / 5 pilot', 'rel 0.25 / 3 pilot', 'rel 0.25 / 8 pilot', 'rel 0.25 / 10 pilot']
LAB = {'rel 0.1 / 5 pilot': '±10%, 5 pilot fields', 'rel 0.5 / 5 pilot': '±50%, 5 pilot fields', 'abs 0.02 / 5 pilot': '±0.02 absolute, 5 pilot fields',
       'rel 0.25 / 3 pilot': '±25%, 3 pilot fields', 'rel 0.25 / 8 pilot': '±25%, 8 pilot fields', 'rel 0.25 / 10 pilot': '±25%, 10 pilot fields'}
mean = (100 * a.groupby(['setting', 'method']).share.mean()).unstack()[M].reindex(ORD)
rng_ = (100 * a.groupby(['setting', 'source', 'method']).share.mean()).unstack()
k = a.groupby(['setting', 'method']).k.median().unstack()[M].reindex(ORD)
s11 = pd.DataFrame(index=[LAB[o] for o in ORD])
for m in M:
    lo = rng_[m].groupby('setting').min().reindex(ORD); hi = rng_[m].groupby('setting').max().reindex(ORD)
    s11[m] = [f'{x:.1f} ({l:.1f}–{h:.1f}); k {int(kk)}' for x, l, h, kk in zip(mean[m], lo, hi, k[m])]
s11['Plans needing all remaining fields, fovplan 90%'] = (100 * (a[a.method == 'fovplan90'].k >= a[a.method == 'fovplan90'].pool).groupby(a.setting).mean()).reindex(ORD).round(1).to_numpy()
s11.to_csv(R + 'S11_planner_sensitivity.tsv', sep='\t'); print(s11.to_string())

# ---- A11 field size (S12)
pr = pd.read_csv(R + 'A11_precision.tsv.gz', sep='\t'); de = pd.read_csv(R + 'A11_design.tsv.gz', sep='\t')
pr['g'] = grp(pr.organ); de['g'] = grp(de.organ)
s12 = de.groupby(['g', 'F']).agg(D=('D', 'median'), cells_per_field=('nbar', 'median'), fields=('M', 'median'), area_required=('area_req_mm2', 'median'))
for b in [1, 2, 4]:
    s12[f'within ±25%, {b} mm2'] = 100 * pr[pr.budget_mm2 == b].groupby(['g', 'F']).share.mean()
s12 = s12.round(2); s12.to_csv(R + 'S12_field_size.tsv', sep='\t'); print(s12.to_string())

# ---- A12 annotation noise (S13)
dd = pd.read_csv(R + 'A12_design.tsv.gz', sep='\t'); fp = pd.read_csv(R + 'A12_fp.tsv.gz', sep='\t'); pl = pd.read_csv(R + 'A12_planning.tsv.gz', sep='\t')
PO = ['none', 'random 5%', 'random 10%', 'field-varying 0-20%']
dd['src'] = dd.source.replace({'human small intestine': 'Human small intestine', 'human colon': 'Human colon'})
s13 = dd.groupby(['perturbation', 'src']).D.median().unstack().reindex(PO)
f = (100 * fp.groupby(['perturbation'])[SIX].apply(lambda d: (d < 0.05).mean())).reindex(PO)
s13 = s13.join(f.rename(columns=lambda x: 'FPR % ' + x))
s13 = s13.join((100 * pl.groupby(['perturbation', 'method']).share.mean()).unstack().reindex(PO).rename(columns=lambda x: 'within ±25% ' + x))
s13 = s13.round(1); s13.to_csv(R + 'S13_annotation_noise.tsv', sep='\t'); print(s13.to_string())

# ---- Figure S13
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(2, 2, figsize=(10, 8), gridspec_kw=dict(hspace=0.75, wspace=0.35))
B, O, A, G = '#2a78d6', '#eb6834', '#1baf7a', '#8a63c9'
# a
x = np.arange(len(SCH)); w = 0.38
a1 = ax[0, 0]
a1.bar(x - w / 2, 100 * p[~p.focal].groupby('scheme').share.mean().reindex(SCH), w, color=B, label='Other cell types')
a1.bar(x + w / 2, 100 * p[p.focal].groupby('scheme').share.mean().reindex(SCH), w, color=O, label='Type used to enrich selection')
a1.axhline(95, ls='--', c='0.4', lw=0.8); a1.set_xticks(x); a1.set_xticklabels(['Random', 'Systematic', 'Interior', 'Cell-dense', 'Enriched'])
a1.set_ylabel('fovplan 90% plan within ±25% (%)'); a1.set_ylim(0, 105); a1.legend(frameon=False, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
a1.set_title('a  Preferential field selection (whole-tissue truth)', loc='left', fontweight='bold', fontsize=8.5)
# b
a2 = ax[0, 1]; xs = np.arange(len(ORD))
for m, col, mk in [('cellcount', O, 'o'), ('bootstrap', G, 's'), ('fovplan', B, '^'), ('fovplan90', A, 'D')]:
    a2.plot(xs, mean[m], marker=mk, color=col, lw=1.2, ms=4, label={'cellcount': 'Cell count', 'bootstrap': 'Bootstrap', 'fovplan': 'fovplan point', 'fovplan90': 'fovplan 90%'}[m])
a2.axhline(95, ls='--', c='0.4', lw=0.8); a2.set_xticks(xs); a2.set_xticklabels(['±10%\n5 pilot', '±50%\n5 pilot', '±0.02\n5 pilot', '±25%\n3 pilot', '±25%\n8 pilot', '±25%\n10 pilot'], fontsize=7)
a2.set_ylabel('Within target (%)'); a2.set_ylim(0, 105); a2.legend(frameon=False, fontsize=7, loc='lower right', ncol=2)
a2.set_title('b  Precision target and pilot size', loc='left', fontweight='bold', fontsize=8.5)
# c
a3 = ax[1, 0]
for gname, ls in [('Human intestine', '-'), ('Mouse kidney', '--')]:
    for F, col in [(250.0, A), (500.0, B), (1000.0, O)]:
        sub = pr[(pr.g == gname) & (pr.F == F)].groupby('budget_mm2').share.mean() * 100
        a3.plot(sub.index, sub.values, ls=ls, marker='o', ms=3.5, color=col, label=f'{gname}, {int(F)}' if gname == 'Human intestine' else None)
a3.set_xticks([1, 2, 4]); a3.set_xlabel('Imaged area (mm², or stored units² × 10⁶ for mouse)'); a3.set_ylabel('Within ±25% (%)'); a3.set_ylim(0, 100)
a3.legend(['250-μm fields', '500-μm fields', '1,000-μm fields'], frameon=False, fontsize=7, loc='lower right')
a3.text(1.05, 92, 'solid: human intestine; dashed: mouse kidney', fontsize=6.8, color='0.3')
a3.set_title('c  Field size at equal imaged area', loc='left', fontweight='bold', fontsize=8.5)
# d
a4 = ax[1, 1]; xp = np.arange(len(PO))
for src, col in [('GSE325587', B), ('Lupus', A), ('GSE211785', '0.55'), ('Human small intestine', O), ('Human colon', G)]:
    a4.plot(xp, s13[src], marker='o', ms=4, color=col, label='D, ' + {'GSE325587': 'DKD kidney', 'Lupus': 'lupus kidney', 'GSE211785': 'kidney atlas'}.get(src, src))
a4.set_ylabel('Median design effect D'); a4.set_xticks(xp); a4.set_xticklabels(['None', 'Random 5%', 'Random 10%', 'Field-varying\n0–20%'], fontsize=7)
a5 = a4.twinx(); a5.spines['top'].set_visible(False)
a5.plot(xp, s13['FPR % cell_chisq'], 'k--', marker='x', ms=4, label='FPR, cell-level test (kidney)')
a5.plot(xp, s13['FPR % indiv_quasibin'], 'k:', marker='+', ms=5, label='FPR, individual quasi-binomial')
a5.set_ylabel('False-positive rate (%)'); a5.set_ylim(0, 100)
h1, l1 = a4.get_legend_handles_labels(); h2, l2 = a5.get_legend_handles_labels()
a4.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=2); a4.set_ylim(0, 45)
a4.set_title('d  Cell-type annotation error', loc='left', fontweight='bold', fontsize=8.5)
fig.savefig(R + 'Supplementary_Figure_S13.png', dpi=300, bbox_inches='tight'); fig.savefig(R + 'Supplementary_Figure_S13.pdf', bbox_inches='tight')
print('figure saved')
