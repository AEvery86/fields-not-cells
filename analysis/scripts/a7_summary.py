"""Summaries and Supplementary Figure S12 for A7 (mixed models and DA tools)."""
import glob, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a7_glmm_da import ALL

R = '../results/a7/'
LAB = {'cell_chisq': 'Cell-level chi-square', 'field_wilcoxon': 'Field Wilcoxon', 'field_quasibin': 'Field quasi-binomial',
       'field_cluster': 'Field logistic, clustered SE', 'glmm_ri_lrt': 'GLMM (individual) LRT', 'glmm_ri_wald': 'GLMM (individual) Wald',
       'glmm_olre_lrt': 'GLMM (individual + field) LRT', 'glmm_olre_wald': 'GLMM (individual + field) Wald',
       'indiv_wilcoxon': 'Individual Wilcoxon', 'indiv_quasibin': 'Individual quasi-binomial', 'betabin': 'Beta-binomial (individual)',
       'nb_glm': 'Negative binomial (individual)', 'propeller_logit': 'propeller (logit)', 'propeller_asin': 'propeller (arcsine)'}
ORDER = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'glmm_ri_wald', 'glmm_ri_lrt', 'glmm_olre_wald',
         'glmm_olre_lrt', 'indiv_wilcoxon', 'indiv_quasibin', 'betabin', 'nb_glm', 'propeller_logit', 'propeller_asin']
UNIT = {**{k: 'cells' for k in ['cell_chisq']}, **{k: 'fields' for k in ['field_wilcoxon', 'field_quasibin']},
        **{k: 'mixed' for k in ['field_cluster', 'glmm_ri_wald', 'glmm_ri_lrt', 'glmm_olre_wald', 'glmm_olre_lrt']},
        **{k: 'individuals' for k in ['indiv_wilcoxon', 'indiv_quasibin', 'betabin', 'nb_glm', 'propeller_logit', 'propeller_asin']}}
COL = {'cells': '#eb6834', 'fields': '#2a78d6', 'mixed': '#8a63c9', 'individuals': '#1baf7a'}


def rej(df, keys):
    return df.groupby(keys)[ORDER].apply(lambda d: (d < 0.05).mean())


k = pd.read_csv(R + 'kidney.tsv.gz', sep='\t')
kn = k[k.scenario != 'real_contrast']
t_k = rej(kn, ['source', 'scenario']).T
real = k[k.scenario == 'real_contrast']


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


t_real = real.groupby('source')[ORDER].apply(lambda d: pd.Series({c: int((bh(d[c].fillna(1)) < 0.05).sum()) for c in ORDER})).T
h = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(R + 'hickey_*.tsv.gz'))])
t_h = rej(h, ['scenario']).T
s = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(R + 'sim_*.tsv.gz'))])
t_s = rej(s, ['delta', 'n_per_group', 'k']).T

tab = pd.DataFrame({'Unit': [UNIT[c] for c in ORDER], 'Test': [LAB[c] for c in ORDER]}, index=ORDER)
for (src, sc), col in t_k.items():
    tab[f'{src} null, {sc}'] = (100 * col).round(1)
for sc, col in t_h.items():
    tab[f'Intestine donor null, {sc}'] = (100 * col).round(1)
for (d, npg, kk), col in t_s.items():
    if npg == 6:
        tab[f'Sim {"null" if d == 0 else "power"}, 6/group, {kk} fields'] = (100 * col).round(1)
for src, col in t_real.items():
    tab[f'{src} real contrast, BH discoveries'] = col
tab.to_csv(R + 'A7_summary_table.tsv', sep='\t', index=False)
t_s.to_csv(R + 'A7_simulation_rejection.tsv', sep='\t')
print(tab.to_string())
print('\nmissing P values per test:\n', pd.concat([k, h, s])[ORDER].isna().mean().round(4))
print('\nsim power, all designs:\n', (100 * t_s.xs(0.4, level='delta', axis=1)).round(1).to_string())
print('\nsim null, all designs:\n', (100 * t_s.xs(0.0, level='delta', axis=1)).round(1).to_string())

# ---------------- figure S12
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(1, 4, figsize=(11, 4.2), sharey=True, gridspec_kw=dict(wspace=0.08))
y = np.arange(len(ORDER))[::-1]
panels = [('a  Kidney, random patient splits', [(t_k[('GSE325587', 'all_fields')], 'o', 'DKD (all fields)'), (t_k[('Lupus', 'all_fields')], 's', 'Lupus (all fields)')]),
          ('b  Intestine, random donor splits', [(t_h['all_fields'], 'o', 'All fields'), (t_h['10_fields'], 's', '10 fields/region')]),
          ('c  Simulation, no true difference', [(t_s[(0.0, 6, 3)], 'o', '6/group, 3 fields'), (t_s[(0.0, 6, 36)], 's', '6/group, 36 fields')]),
          ('d  Simulation, power (δ = 0.4)', [(t_s[(0.4, 6, 3)], 'o', '6/group, 3 fields'), (t_s[(0.4, 6, 36)], 's', '6/group, 36 fields')])]
from matplotlib.lines import Line2D
for a, (title, series) in zip(ax, panels):
    for j, (ser, m, lab) in enumerate(series):
        vals = 100 * ser.reindex(ORDER).to_numpy()
        a.scatter(vals, y + (0.15 if j == 0 else -0.15), marker=m, s=22, c=[COL[UNIT[c]] for c in ORDER],
                  edgecolors='k' if j == 0 else 'none', linewidths=0.4, zorder=3)
    if 'power' not in title:
        a.axvline(5, ls='--', c='0.4', lw=0.8)
        a.set_xscale('symlog', linthresh=10); a.set_xlim(0, 110); a.set_xticks([0, 5, 10, 25, 50, 100]); a.set_xticklabels(['0', '5', '10', '25', '50', '100'])
        a.set_xlabel('Rejection at P < 0.05 (%)'); loc = 'lower right'
    else:
        a.set_xlim(0, 105); a.set_xticks([0, 25, 50, 75, 100]); a.set_xlabel('Power at P < 0.05 (%)'); loc = 'lower left'
    a.set_title(title, loc='left', fontsize=8.5, fontweight='bold', pad=16)
    a.grid(axis='x', c='0.92', lw=0.6, zorder=0)
    hs = [Line2D([], [], ls='', marker=series[0][1], mfc='0.6', mec='k', mew=0.4, ms=5),
          Line2D([], [], ls='', marker=series[1][1], mfc='0.6', mec='none', ms=5)]
    a.legend(hs, [series[0][2], series[1][2]], frameon=False, fontsize=6.8, loc='lower left', bbox_to_anchor=(-0.02, 0.995), ncol=2, handletextpad=0.1, columnspacing=0.8, borderaxespad=0)
ax[0].set_yticks(y); ax[0].set_yticklabels([LAB[c] for c in ORDER])
hu = [Line2D([], [], ls='', marker='o', mfc=COL[c], mec='none', ms=6) for c in ['cells', 'fields', 'mixed', 'individuals']]
fig.legend(hu, ['Cells as replicates', 'Fields as replicates', 'Field data, individual-level model', 'Individual totals'],
           loc='lower center', ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(0.55, -0.01))
fig.subplots_adjust(left=0.2, right=0.99, top=0.86, bottom=0.2)
fig.savefig(R + 'Supplementary_Figure_S12.png', dpi=300); fig.savefig(R + 'Supplementary_Figure_S12.pdf')
print('figure saved')
