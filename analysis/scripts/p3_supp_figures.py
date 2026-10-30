"""G23: supplementary figures redrawn in the unified style (S1, S3-S10). S2 needs cell coordinates and is kept."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a5_figures import BLUE, ORANGE, AQUA, INK, INK2, GRID, tag
YELLOW = '#eda100'; VIOLET = '#6250d6'; AQ_TXT = '#12855c'; GREY = '#c3c2b7'
D = '../data/'; OUT = '../../figures/supplementary/'
import os; os.makedirs(OUT, exist_ok=True)
METH = ['cell_Wilson', 'region_ratio_t', 'Fieller', 'region_percentile', 'percentile_FPC', 'random_cell_Wilson']
MNAME = {'cell_Wilson': 'Cell Wilson', 'region_ratio_t': 'Ratio t', 'Fieller': 'Fieller', 'region_percentile': 'Region bootstrap',
         'percentile_FPC': 'Bootstrap + FPC', 'random_cell_Wilson': 'Random-cell Wilson'}


def save(fig, name):
    fig.savefig(OUT + name + '.png', dpi=400, bbox_inches='tight', facecolor='white'); plt.close(fig)


def grid(ax, axis='y'):
    getattr(ax, axis + 'axis').grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)


def s1():
    p = pd.read_csv(D + 'sampling_by_sample.tsv', sep='\t'); p = p[p.k == 3].pivot(index='sample', columns='strategy', values='TV_median')
    p = p.sort_values('random_FOV'); y = np.arange(len(p))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(wspace=0.35))
    ax = axs[0]; tag(ax, 'a', -0.3)
    for i, (s, r) in enumerate(p.iterrows()):
        ax.plot([r.random_FOV, r.dispersed_FOV], [i, i], color=GREY, lw=0.8)
    ax.scatter(p.random_FOV, y, s=14, color=ORANGE, label='3 random fields', zorder=3)
    ax.scatter(p.dispersed_FOV, y, s=14, color=VIOLET, label='3 dispersed fields', zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(p.index, fontsize=5.6); ax.set_xlabel('Median TV'); grid(ax, 'x')
    ax.legend(loc='lower right'); ax.set_title('Random versus geometrically dispersed fields', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.3)
    d = (p.dispersed_FOV - p.random_FOV)
    ax.barh(y, d, color=[VIOLET if v < 0 else ORANGE for v in d], height=0.6)
    ax.axvline(0, color=INK2, lw=0.6); ax.set_yticks(y); ax.set_yticklabels([]); grid(ax, 'x')
    ax.set_xlabel('Dispersed − random median TV'); ax.set_title('Paired P = 0.68', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S1')


def s3():
    s = pd.read_csv(D + 'small_region_by_unit.tsv', sep='\t')
    s = s[(s.k == 3) & (s.method == 'cell_Wilson')].pivot(index='unit', columns='design', values='coverage')
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    for u, r in s.iterrows():
        ax.plot([0, 1], [r.observed, r.label_shuffle], color=GREY, lw=0.8)
    ax.scatter(np.zeros(len(s)), s.observed, color=ORANGE, s=14, zorder=3); ax.scatter(np.ones(len(s)), s.label_shuffle, color=AQUA, s=14, zorder=3)
    ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xticks([0, 1]); ax.set_xticklabels(['Observed labels', 'Labels redistributed']); ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel('Cell-level Wilson coverage (3 fields)'); ax.set_ylim(0, 1.02); grid(ax)
    ax.set_title('Label-redistribution control, DKD cohort', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S3')


def s4():
    m = pd.read_csv(D + 'mouse_example_sampling_summary.tsv', sep='\t')
    order = ['main', 'grid250', 'grid1000', 'shiftx', 'shifty', 'shiftxy']
    lab = ['500 (main)', '250', '1,000', 'shift x', 'shift y', 'shift xy']
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True, gridspec_kw=dict(wspace=0.08))
    for j, s in enumerate(['ShamR', 'Day14R']):
        ax = axs[j]; tag(ax, 'ab'[j], -0.12 if j else -0.2)
        g = m[m['sample'] == s].set_index('setting').loc[order]; x = np.arange(len(order))
        ax.bar(x - 0.2, g.region_TV, 0.38, color=ORANGE, label='3 regions'); ax.bar(x + 0.2, g.cell_TV, 0.38, color=BLUE, label='Matched observations')
        for xi, r in zip(x, g.ratio):
            ax.text(xi, g.region_TV.iloc[list(x).index(xi)] + 0.006, f'×{r:.1f}', ha='center', fontsize=5.6, color=INK2)
        ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=5.8); grid(ax)
        ax.set_title(f'{s}: grid width and origin', loc='left', color=INK)
        ax.set_ylim(0, 0.29)
        if j == 0: ax.set_ylabel('Median 20-class TV'); ax.legend(loc='upper left', fontsize=5.6, ncol=2)
    save(fig, 'Supplementary_Figure_S4')


def s5():
    p = pd.read_csv(D + 'sampling_by_sample.tsv', sep='\t')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(wspace=0.35))
    ax = axs[0]; tag(ax, 'a', -0.25)
    a = p[p.k == 3].pivot(index='sample', columns='strategy', values='TV_median')
    for s, r in a.iterrows():
        ax.plot([0, 1], [r.random_cells, r.random_FOV], color=GREY, lw=0.8)
    ax.scatter(np.zeros(len(a)), a.random_cells, color=BLUE, s=12, zorder=3); ax.scatter(np.ones(len(a)), a.random_FOV, color=ORANGE, s=12, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Matched random cells', '3 random fields']); ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel('Median TV'); grid(ax); ax.set_title('Paired by tissue (n = 17)', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.25)
    g = p.groupby(['k', 'strategy']).TV_median.median().unstack(); n = p[p.strategy == 'random_FOV'].groupby('k').size()
    for c, col, name in [('random_FOV', ORANGE, 'Random fields'), ('dispersed_FOV', VIOLET, 'Dispersed fields'), ('random_cells', BLUE, 'Matched random cells')]:
        ax.plot(g.index, g[c], color=col, marker='o', ms=3.5, lw=1.6, label=name)
    ax.set_xscale('log'); ax.set_xticks(g.index); ax.set_xticklabels([f'{k}\n(n={n[k]})' for k in g.index], fontsize=5.8)
    ax.set_xlabel('Fields sampled (tissues eligible)'); ax.set_ylabel('Median of tissue medians (TV)'); grid(ax)
    ax.legend(loc='upper right', fontsize=5.6); ax.set_title('Error by number of fields', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S5')


def strip_methods(ax, df, col, log=False):
    y = np.arange(len(METH))[::-1]
    for m, yy in zip(METH, y):
        v = df[df.method == m][col]
        ax.scatter(v, np.full(len(v), yy) + np.random.default_rng(yy).uniform(-0.12, 0.12, len(v)), s=6, color=BLUE, alpha=0.6, lw=0)
        ax.plot([v.median()] * 2, [yy - 0.22, yy + 0.22], color=INK, lw=1.6)
    ax.set_yticks(y); ax.set_yticklabels([MNAME[m] for m in METH]); grid(ax, 'x')
    if log: ax.set_xscale('log')


def s6():
    s = pd.read_csv(D + 'small_region_by_unit.tsv', sep='\t'); s = s[(s.k == 3) & (s.design == 'observed')]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(wspace=0.08))
    tag(axs[0], 'a', -0.38); strip_methods(axs[0], s, 'coverage'); axs[0].axvline(0.95, color=INK2, lw=0.6, ls='--'); axs[0].set_xlim(0, 1.02)
    axs[0].set_xlabel('Coverage (500 draws of 3 fields)'); axs[0].set_title('DKD cohort, 17 tissues', loc='left', color=INK)
    tag(axs[1], 'b', -0.05); strip_methods(axs[1], s, 'mean_width', log=True); axs[1].set_yticklabels([]); axs[1].set_xlabel('Mean width (log scale)')
    axs[1].set_title('Bars = medians', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S6')


def s7():
    g = pd.read_csv(D + 'GSE211785_interval_summary.tsv', sep='\t'); g = g[g.k == 3]
    t = pd.read_csv(D + 'GSE211785_sampling_summary.tsv', sep='\t')
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(width_ratios=[0.65, 1.1, 1.0], wspace=0.95))
    ax = axs[0]; tag(ax, 'a', -0.45); x = np.arange(2)
    ax.bar(x - 0.2, t.region_TV, 0.38, color=ORANGE, label='3 FOVs'); ax.bar(x + 0.2, t.cell_TV, 0.38, color=BLUE, label='Matched cells')
    for xi, r in zip(x, t.ratio): ax.text(xi - 0.2, t.region_TV.iloc[xi] + 0.004, f'×{r:.1f}', ha='center', fontsize=5.6, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(t['sample']); ax.set_ylabel('Median 13-class TV'); grid(ax); ax.legend(fontsize=5.4, loc='upper center', bbox_to_anchor=(0.5, 1.0))
    ax.set_ylim(0, 0.2); ax.set_title('Composition', loc='left', color=INK)
    for j, (col, xl) in enumerate([('coverage', 'Coverage'), ('mean_width', 'Mean width')]):
        ax = axs[j + 1]; tag(ax, 'bc'[j], -0.55 if j == 0 else -0.1)
        y = np.arange(len(METH))[::-1]
        for s_, mk, dy in [('HK2844', 'o', 0.12), ('HK3039', 's', -0.12)]:
            v = g[g['sample'] == s_].set_index('method').loc[METH, col]
            ax.scatter(v, y + dy, marker=mk, s=14, color=BLUE if s_ == 'HK2844' else AQUA, label=s_)
        ax.set_yticks(y); ax.set_yticklabels([MNAME[m] for m in METH] if j == 0 else []); grid(ax, 'x'); ax.set_xlabel(xl)
        if j == 0: ax.axvline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xlim(0, 1.02); ax.legend(fontsize=5.4, loc='lower left')
        else:
            ax.set_xscale('log'); ax.set_xticks([0.02, 0.05, 0.1, 0.2]); ax.set_xticklabels(['0.02', '0.05', '0.1', '0.2']); ax.minorticks_off()
        ax.set_title('Immune proportion, 3 FOVs', loc='left', color=INK) if j == 0 else None
    save(fig, 'Supplementary_Figure_S7')


def s8():
    ss = pd.read_csv(D + 'sampling_summary.tsv', sep='\t'); ss = ss[ss.k == 3].sort_values('ratio')
    ex = pd.read_csv(D + 'lupus_full_exact_three_field.tsv', sep='\t').set_index('tissue').loc[ss.tissue]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.0), sharey=True, gridspec_kw=dict(wspace=0.12))
    y = np.arange(len(ss))
    ax = axs[0]; tag(ax, 'a', -0.5); ax.barh(y, ss.ratio, color=[AQUA if t.startswith('SLE') else BLUE for t in ss.tissue], height=0.6)
    ax.axvline(1, color=INK2, lw=0.6, ls='--'); ax.set_yticks(y); ax.set_yticklabels(ss.tissue, fontsize=5.8); ax.set_xlabel('TV ratio, 3 FOVs / matched cells'); grid(ax, 'x')
    ax.set_title('Composition', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.08)
    ax.scatter(ex.cell_wilson_exact, y, color=ORANGE, s=14, label='Cell Wilson'); ax.scatter(ex.fieller_exact, y, color=BLUE, s=14, marker='s', label='Fieller')
    ax.axvline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_xlim(0, 1.02); ax.set_xlabel('Exact coverage (all 3-FOV subsets)'); grid(ax, 'x'); ax.legend(fontsize=5.6, loc='lower left')
    ax.set_title('Immune proportion', loc='left', color=INK)
    ax = axs[2]; tag(ax, 'c', -0.08); ax.barh(y, ex.fieller_full_exact, color=BLUE, height=0.6); ax.set_xlim(0, 1); ax.set_xlabel('Fieller sets equal to [0,1]'); grid(ax, 'x')
    ax.set_title('Uninformative sets', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S8')


def s9():
    e = pd.read_csv(D + 'external_sampling_by_entry.tsv', sep='\t')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.35))
    ax = axs[0]; tag(ax, 'a', -0.25)
    a = e[e.setting == 'primary'].pivot(index='gsm', columns='strategy', values='TV_median')
    for gsm, r in a.iterrows(): ax.plot([0, 1], [r.random_cells, r.random_regions], color=GREY, lw=0.8)
    ax.scatter(np.zeros(len(a)), a.random_cells, color=BLUE, s=12, zorder=3); ax.scatter(np.ones(len(a)), a.random_regions, color=ORANGE, s=12, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Matched random cells', '3 random 500-μm regions']); ax.set_xlim(-0.3, 1.3); ax.set_ylabel('Median 50-gene profile TV'); grid(ax)
    ax.set_title('21 nonduplicate entries', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.25)
    order = ['grid250', 'primary', 'grid1000', 'counts10']; lab = ['250 μm', '500 μm\n(main)', '1,000 μm', '500 μm,\n≥10 transcripts']
    g = e.groupby(['setting', 'strategy']).TV_median.median().unstack().loc[order]; x = np.arange(len(order))
    ax.bar(x - 0.2, g.random_regions, 0.38, color=ORANGE, label='3 regions'); ax.bar(x + 0.2, g.random_cells, 0.38, color=BLUE, label='Matched cells')
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=5.6); ax.set_ylabel('Median of entry medians'); grid(ax); ax.legend(fontsize=5.6)
    ax.set_title('Pre-specified settings', loc='left', color=INK)
    save(fig, 'Supplementary_Figure_S9')


def s10():
    m = pd.read_csv(D + 'mouse_example_interval_summary.tsv', sep='\t')
    cols = {'cell_Wilson': ORANGE, 'region_ratio_t': BLUE, 'Fieller': YELLOW, 'region_percentile': VIOLET, 'percentile_FPC': '#86b6ef', 'random_cell_Wilson': AQUA}
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True, gridspec_kw=dict(hspace=0.15, wspace=0.12))
    for j, s in enumerate(['ShamR', 'Day14R']):
        g = m[m['sample'] == s]
        for i, (col, yl) in enumerate([('coverage', 'Coverage'), ('mean_width', 'Mean width (log)')]):
            ax = axs[i, j]
            if i == 0: tag(ax, 'ab'[j], -0.15 if j == 0 else -0.05)
            for meth in METH:
                h = g[g.method == meth].sort_values('k')
                ax.plot(h.k, h[col], color=cols[meth], marker='o', ms=3, lw=1.4, label=MNAME[meth])
            if i == 0: ax.axhline(0.95, color=INK2, lw=0.6, ls='--'); ax.set_ylim(0.3, 1.02); ax.set_title(s, loc='left', color=INK)
            else: ax.set_yscale('log'); ax.set_xlabel('Regions sampled')
            if j == 0: ax.set_ylabel(yl)
            else: ax.set_yticklabels([])
            ax.set_xticks([3, 5, 10]); grid(ax)
    axs[0, 1].legend(fontsize=5.4, loc='lower right', ncol=2)
    save(fig, 'Supplementary_Figure_S10')


if __name__ == '__main__':
    for f in (s1, s3, s4, s5, s6, s7, s8, s9, s10):
        f()
    print('done')
