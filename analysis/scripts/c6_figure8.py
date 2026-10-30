"""G22 Figure 8: whole-section human references and cross-organ comparison."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a5_figures import BLUE, ORANGE, AQUA, INK, INK2, GRID, tag, save
YELLOW = '#eda100'; AQ_TXT = '#12855c'; G = '../results/'; R = '../results/'


def fig8():
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.5)
    # a TV vs k
    ax = fig.add_subplot(gs[0, 0]); tag(ax, 'a', -0.32)
    sa = pd.read_csv(G + 'C4_hickey_sampling.tsv', sep='\t')
    for organ, ls in [('human small intestine', '-'), ('human colon', '--')]:
        g = sa[sa.organ == organ].groupby('k')[['tv_fields', 'tv_cells']].median()
        ax.plot(g.index, g.tv_fields, color=ORANGE, ls=ls, marker='o', ms=3.5, lw=1.6)
        ax.plot(g.index, g.tv_cells, color=BLUE, ls=ls, marker='s', ms=3.5, lw=1.6)
    ax.text(3.3, 0.24, 'Random fields', color='#b4471b', fontsize=5.8); ax.text(2.2, 0.065, 'Matched random cells', color=BLUE, fontsize=5.8)
    ax.text(4.0, 0.36, 'solid: small intestine\ndashed: colon', color=INK2, fontsize=5.6)
    ax.set_xscale('log'); ax.set_xticks([1, 3, 5, 10]); ax.set_xticklabels(['1', '3', '5', '10']); ax.set_ylim(0, 0.45)
    ax.set_xlabel('Fields sampled'); ax.set_ylabel('Median TV vs whole region')
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Human intestine, 57 regions', loc='left', color=INK)
    # b design effects across organs
    ax = fig.add_subplot(gs[0, 1:]); tag(ax, 'b', -0.12)
    a1 = pd.read_csv(R + 'A1_design_effect_atlas.tsv', sep='\t'); a1 = a1[(a1.level == 'label') & (a1.p >= .01)]
    m = pd.read_csv(R + 'C3_mouse_F500_atlas.tsv', sep='\t'); m = m[m.p_true >= .01]
    h = pd.read_csv(G + 'C4_hickey_atlas.tsv', sep='\t'); h = h[h.p_true >= .01]
    groups = [('Kidney DKD\n(obs. fields)', a1[a1.source == 'GSE325587'].D), ('Kidney lupus\n(obs. fields)', a1[a1.source == 'Lupus'].D),
              ('Kidney\nGSE211785\n(obs. fields)', a1[a1.source == 'GSE211785'].D), ('Mouse\nkidney\n(whole)', m.D),
              ('Small\nintestine\n(whole)', h[h.organ == 'human small intestine'].D), ('Colon\n(whole)', h[h.organ == 'human colon'].D)]
    cols = [BLUE, BLUE, BLUE, ORANGE, AQUA, AQUA]
    rng = np.random.default_rng(1)
    for i, ((name, v), c) in enumerate(zip(groups, cols)):
        v = v.dropna()
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=4, color=c, alpha=0.45, lw=0)
        ax.plot([i - 0.25, i + 0.25], [v.median()] * 2, color=INK, lw=1.6)
        ax.text(i + 0.28, v.median(), f'{v.median():.0f}', fontsize=5.8, color=INK, va='center')
    ax.set_yscale('log'); ax.axhline(1, color=INK2, lw=0.6, ls='--')
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups], fontsize=5.6)
    ax.set_ylabel('Design effect D (cell types ≥1%)'); ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title('Design effects across organs, platforms and species (bars = medians)', loc='left', color=INK)
    # c coverage k=3 by organ
    ax = fig.add_subplot(gs[1, 0]); tag(ax, 'c', -0.32)
    it = pd.concat([pd.read_csv(G + 'C4_hickey_intervals.tsv', sep='\t'), pd.read_csv(R + 'C3_mouse_F500_intervals.tsv', sep='\t')])
    it = it[it.k == 3].groupby(['organ', 'method']).coverage.median().unstack()
    orgs = ['mouse kidney', 'human small intestine', 'human colon']; x = np.arange(3); w = 0.26
    for j, (mth, col, name) in enumerate([('Wilson', ORANGE, 'Cell Wilson'), ('ratio_t', BLUE, 'Ratio t'), ('quasibinomial', AQUA, 'Quasi-binomial')]):
        ax.bar(x + (j - 1) * w, it.loc[orgs, mth], width=w * 0.9, color=col, label=name)
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks(x); ax.set_xticklabels(['Mouse\nkidney', 'Small\nintestine', 'Colon'], fontsize=6)
    ax.set_ylim(0, 1.25); ax.set_yticks([0, .2, .4, .6, .8, 1]); ax.set_ylabel('Coverage of true proportion (3 fields)')
    ax.legend(loc='upper center', ncol=3, fontsize=5.2, handlelength=0.9, columnspacing=0.5, bbox_to_anchor=(0.5, 1.02))
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Intervals vs truth', loc='left', color=INK)
    # d planning
    ax = fig.add_subplot(gs[1, 1]); tag(ax, 'd', -0.32)
    c5 = pd.read_csv(G + 'C5_pilot_sensitivity.tsv', sep='\t'); c5 = c5[c5.pilot == 5]
    pl = pd.concat([pd.read_csv(G + 'C4_hickey_planning.tsv', sep='\t'), pd.read_csv(R + 'C3_mouse_F500_planning.tsv', sep='\t')])
    vals = pd.DataFrame({'Cell count': pl.groupby('organ').cells_assur.mean(),
                         'fovplan': c5[c5.q == 0].groupby('organ').assur.mean(),
                         'fovplan 80%': c5[c5.q == 0.8].groupby('organ').assur.mean(),
                         'fovplan 90%': c5[c5.q == 0.9].groupby('organ').assur.mean()}).loc[orgs]
    w = 0.2
    for j, (cname, col) in enumerate(zip(vals.columns, [ORANGE, BLUE, '#6da7ec', AQUA])):
        ax.bar(x + (j - 1.5) * w, vals[cname], width=w * 0.9, color=col, label=cname)
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks(x); ax.set_xticklabels(['Mouse\nkidney', 'Small\nintestine', 'Colon'], fontsize=6)
    ax.set_ylim(0, 1.25); ax.set_yticks([0, .2, .4, .6, .8, 1]); ax.set_ylabel('Share within ±25% of truth')
    ax.legend(loc='upper center', ncol=2, fontsize=5.2, handlelength=0.9, columnspacing=0.5, bbox_to_anchor=(0.5, 1.04))
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Planning vs truth', loc='left', color=INK)
    # e donor-level FPR
    ax = fig.add_subplot(gs[1, 2]); tag(ax, 'e', -0.32)
    f = pd.read_csv(G + 'C4_hickey_da_fpr.tsv', sep='\t').set_index('scenario')
    TESTS = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']
    TN = ['Cells χ²', 'Fields Wilcoxon', 'Fields quasi-bin.', 'Field model,\nclustered SE', 'Donors Wilcoxon', 'Donors quasi-bin.']
    TC = [ORANGE, BLUE, BLUE, BLUE, AQUA, AQUA]
    y = np.arange(len(TESTS))[::-1]
    for t, yy, c in zip(TESTS, y, TC):
        ax.barh(yy + 0.18, f.loc['all_fields', t], height=0.34, color=c)
        ax.barh(yy - 0.18, f.loc['10_fields', t], height=0.34, color='white', edgecolor=c, lw=0.8)
    ax.axvline(0.05, color=INK2, lw=0.6, ls='--'); ax.set_yticks(y); ax.set_yticklabels(TN, fontsize=5.6); ax.set_xlim(0, 1)
    ax.set_xlabel('False-positive rate, random donor splits'); ax.text(0.98, 0.02, 'filled: all fields\nopen: 10 fields/region', transform=ax.transAxes, ha='right', fontsize=5.4, color=INK2)
    ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('8 real donors', loc='left', color=INK)
    save(fig, 'Figure8')


if __name__ == '__main__':
    fig8(); print('done')
