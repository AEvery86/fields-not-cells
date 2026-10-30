"""G21 figures 6 (simulation with known truth) and 7 (whole-section truth, mouse kidney)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a5_figures import BLUE, ORANGE, AQUA, INK, INK2, GRID, tag, save  # shared style
YELLOW = '#eda100'
R = '../results/'
SIGCOL = {0: '#86b6ef', 0.5: '#3987e5', 1.0: '#1c5cab', 1.5: '#0d366b'}  # ordinal blue ramp (steps 250-700)
AQ_TXT = '#12855c'


def lab(ax, x, y, s, c, **kw):
    ax.text(x, y, s, color=c, fontsize=5.8, **kw)


def fig6():
    s1 = pd.read_csv(R + 'B1_sim_design_effect.tsv', sep='\t')
    s2 = pd.read_csv(R + 'B1_sim_planning.tsv.gz', sep='\t'); s2 = s2[~s2.pilot_zero]
    s3 = pd.read_csv(R + 'B1_sim_intervals.tsv', sep='\t')
    b2 = pd.read_csv(R + 'B2_sim_da_rejection.tsv', sep='\t')
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42)
    # a D analytic vs empirical
    ax = fig.add_subplot(gs[0, 0]); tag(ax, 'a', -0.3)
    for sg, c in SIGCOL.items():
        g = s1[s1.sigma == sg]
        ax.scatter(g.D_true, g.D_emp3, s=5, color=c, lw=0, alpha=0.8, label=f'σ = {sg:g}')
    ax.plot([0.5, 500], [0.5, 500], color=INK2, lw=0.6, ls='--')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.5, 500); ax.set_ylim(0.5, 500)
    ax.set_xlabel('Analytic design effect D'); ax.set_ylabel('Empirical design effect\n(3 fields vs true proportion)')
    ax.legend(loc='upper left', handletextpad=0.1, markerscale=1.6, title='Clustering', title_fontsize=6)
    ax.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('D is accurate against truth', loc='left', color=INK)
    # b planning
    ax = fig.add_subplot(gs[0, 1]); tag(ax, 'b', -0.3)
    g = s2.groupby('sigma')[['cells_assur', 'fovplan_assur', 'fovplan80_assur']].mean()
    x = np.arange(len(g))
    for c, col, mk, name, dy in [('cells_assur', ORANGE, 'o', 'Cell-count plan', -0.07), ('fovplan_assur', BLUE, 's', 'fovplan', -0.065),
                                 ('fovplan80_assur', AQUA, 'D', 'fovplan, 80% assurance', 0.04)]:
        ax.plot(x, g[c], color=col, lw=2, marker=mk, ms=4)
        pos = {'cells_assur': (1.6, 0.56), 'fovplan_assur': (1.5, 0.86), 'fovplan80_assur': (0.9, 1.03)}[c]
        lab(ax, pos[0], pos[1], name, col if col != AQUA else AQ_TXT, va='center')
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--')
    ax.set_xticks(x); ax.set_xticklabels([f'{v:g}' for v in g.index]); ax.set_xlim(-0.2, 3.2); ax.set_ylim(0, 1.1)
    ax.set_xlabel('Spatial clustering σ'); ax.set_ylabel('Share within ±25% of true proportion'); ax.set_ylim(0, 1.05)
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Planning to a precision target', loc='left', color=INK)
    # c interval coverage at k=3
    ax = fig.add_subplot(gs[0, 2]); tag(ax, 'c', -0.3)
    g = s3[s3.k == 3].groupby(['sigma', 'method']).coverage.mean().unstack()
    for m, col, mk, name, dy in [('Wilson', ORANGE, 'o', 'Cell-level Wilson', 0), ('ratio_t', BLUE, 's', 'Ratio t', 0.035),
                                 ('Fieller', YELLOW, '^', 'Fieller', -0.035), ('quasibinomial', AQUA, 'D', 'Quasi-binomial', 0.0)]:
        ax.plot(x, g[m], color=col, lw=1.8, marker=mk, ms=4)
        pos = {'Wilson': (1.75, 0.38), 'ratio_t': (1.55, 0.80), 'Fieller': (1.55, 0.74), 'quasibinomial': (0.3, 1.03)}[m]
        lab(ax, pos[0], pos[1], name, col if col not in (AQUA, YELLOW) else (AQ_TXT if col == AQUA else '#9a6a00'), va='center')
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks(x); ax.set_xticklabels([f'{v:g}' for v in g.index]); ax.set_xlim(-0.2, 3.0)
    ax.set_xlabel('Spatial clustering σ'); ax.set_ylabel('Coverage of true proportion'); ax.set_ylim(0, 1.1)
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Three-field intervals', loc='left', color=INK)
    # d false positives under true null
    ax = fig.add_subplot(gs[1, 0:2]); tag(ax, 'd', -0.1)
    TESTS = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']
    TN = ['Cells\n(χ²)', 'Fields\n(Wilcoxon)', 'Fields\n(quasi-bin.)', 'Field model,\nclustered SE', 'Individuals\n(Wilcoxon)', 'Individuals\n(quasi-bin.)']
    TC = [ORANGE, BLUE, BLUE, BLUE, AQUA, AQUA]
    nul = b2[(b2.delta == 0) & (b2.n_per_group == 6)].set_index('k')
    w = 0.26
    for j, (k, hatch) in enumerate([(3, None), (10, None), (36, None)]):
        for i, t in enumerate(TESTS):
            ax.bar(i + (j - 1) * w, nul.loc[k, t], width=w * 0.92, color=TC[i], alpha=[0.45, 0.72, 1.0][j], lw=0)
    for i in range(len(TESTS)):
        pass
    ax.axhline(0.05, color=INK2, lw=0.7, ls='--'); ax.text(5.55, 0.07, '5%', fontsize=5.8, color=INK2)
    ax.set_xticks(range(len(TESTS))); ax.set_xticklabels(TN, fontsize=5.8); ax.set_ylim(0, 1)
    ax.set_ylabel('False-positive rate (true null)'); ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.text(0.99, 0.95, 'bars left to right: 3, 10, 36 fields per individual\n6 vs 6 individuals, 300 replicates', transform=ax.transAxes, ha='right', va='top', fontsize=5.8, color=INK2)
    ax.set_title('Simulated case-control studies without a true difference', loc='left', color=INK)
    # e power
    ax = fig.add_subplot(gs[1, 2]); tag(ax, 'e', -0.3)
    pw = b2[b2.delta == 0.4]
    for k, mk, a in [(3, 'o', 0.5), (10, 's', 0.75), (36, 'D', 1.0)]:
        g = pw[pw.k == k].set_index('n_per_group')
        ax.plot(g.index, g.indiv_quasibin, color=AQUA, alpha=a, marker=mk, ms=4, lw=1.8)
        lab(ax, 12.4, g.indiv_quasibin.iloc[-1], f'{k} fields', AQ_TXT, va='center')
    ax.set_xticks([3, 6, 12]); ax.set_xlim(2.5, 15.5); ax.set_ylim(0, 1)
    ax.set_xlabel('Individuals per group'); ax.set_ylabel('Power, individual-level test')
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Power needs individuals', loc='left', color=INK)
    save(fig, 'Figure6')


def fig7():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.45))
    sa = pd.read_csv(R + 'C3_mouse_F500_sampling.tsv', sep='\t')
    ax = axs[0]; tag(ax, 'a', -0.3)
    for sec, ls in [('ShamR', '-'), ('Day14R', '--')]:
        g = sa[sa.section == sec]
        ax.plot(g.k, g.tv_fields, color=ORANGE, ls=ls, marker='o', ms=3.5, lw=1.6)
        ax.plot(g.k, g.tv_cells, color=BLUE, ls=ls, marker='s', ms=3.5, lw=1.6)
    lab(ax, 4.2, 0.20, 'Random fields', '#b4471b'); lab(ax, 4.2, 0.045, 'Matched random cells', BLUE)
    lab(ax, 6.0, 0.30, 'solid: sham\ndashed: day 14', INK2)
    ax.set_xscale('log'); ax.set_xticks([1, 3, 5, 10]); ax.set_xticklabels(['1', '3', '5', '10'])
    ax.set_xlabel('Fields sampled'); ax.set_ylabel('Median TV vs whole section'); ax.set_ylim(0, 0.36)
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Whole-section truth', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.3)
    it = pd.read_csv(R + 'C3_mouse_F500_intervals.tsv', sep='\t')
    g = it.groupby(['k', 'method']).coverage.median().unstack()
    for m, col, mk, name in [('Wilson', ORANGE, 'o', 'Cell-level Wilson'), ('ratio_t', BLUE, 's', 'Ratio t'), ('quasibinomial', AQUA, 'D', 'Quasi-binomial')]:
        ax.plot(range(3), g[m], color=col, marker=mk, ms=4, lw=1.8)
        lab(ax, 0.05, g[m].iloc[0] + (0.05 if m != 'ratio_t' else -0.08), name, col if col != AQUA else AQ_TXT)
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks(range(3)); ax.set_xticklabels(['3', '5', '10'])
    ax.set_xlabel('Fields sampled'); ax.set_ylabel('Coverage of whole-section proportion'); ax.set_ylim(0, 1.05)
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Interval coverage', loc='left', color=INK)
    ax = axs[2]; tag(ax, 'c', -0.3)
    rows = []
    for F in [500, 250]:
        pl = pd.read_csv(R + f'C3_mouse_F{F}_planning.tsv', sep='\t')
        rows.append((F, pl.cells_assur.mean(), pl.fovplan_assur.mean(), pl.fovplan80_assur.mean()))
    t = pd.DataFrame(rows, columns=['F', 'cells', 'fov', 'fov80'])
    w = 0.25
    for j, (c, col, name) in enumerate([('cells', ORANGE, 'Cell-count plan'), ('fov', BLUE, 'fovplan'), ('fov80', AQUA, 'fovplan 80%')]):
        ax.bar(np.arange(2) + (j - 1) * w, t[c], width=w * 0.9, color=col, label=name)
        for i, v in enumerate(t[c]):
            ax.text(i + (j - 1) * w, v + 0.02, f'{v:.2f}', ha='center', fontsize=5.4, color=INK2)
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks([0, 1]); ax.set_xticklabels(['500-unit fields', '250-unit fields'])
    ax.set_ylim(0, 1.28); ax.set_yticks([0, .2, .4, .6, .8, 1]); ax.set_ylabel('Share within ±25% of true proportion')
    ax.legend(loc='upper center', ncol=3, fontsize=5.4, handlelength=1, columnspacing=0.6, bbox_to_anchor=(0.5, 1.02))
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True); ax.set_title('Planning vs truth', loc='left', color=INK)
    save(fig, 'Figure7')


if __name__ == '__main__':
    fig6(); fig7(); print('done')
