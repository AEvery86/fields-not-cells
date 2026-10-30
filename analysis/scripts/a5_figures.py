"""G20 main figures 1-5 (static, publication). Palette: reference categorical slots 1-3
(blue #2a78d6, orange #eb6834, aqua #1baf7a; validated all-pairs) + blue sequential ramp; grey ink.
Aqua is below 3:1 contrast, so every coloured series is also direct-labelled."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle

R = '../results/'; D = '../data/'; OUT = '../../figures/main/'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
SRC_COL = {'GSE325587': BLUE, 'GSE211785': ORANGE, 'Lupus': AQUA}
SRC_NAME = {'GSE325587': 'DKD (GSE325587)', 'GSE211785': 'GSE211785', 'Lupus': 'Lupus nephritis'}
SEQ = LinearSegmentedColormap.from_list('blue', ['#cde2fb', '#86b6ef', '#3987e5', '#256abf', '#184f95', '#0d366b'])
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
                     'xtick.labelsize': 6.3, 'ytick.labelsize': 6.3, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'legend.fontsize': 6.3, 'legend.frameon': False, 'savefig.dpi': 400,
                     'axes.spines.top': False, 'axes.spines.right': False})


def tag(ax, s, x=-0.16, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left', color=INK)


def save(fig, name):
    fig.savefig(OUT + name + '.png', bbox_inches='tight', facecolor='white')
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ---------------- Figure 1: concept + composition error -----------------
def fig1():
    fig = plt.figure(figsize=(7.2, 3.15))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.75], wspace=0.28)
    ax = fig.add_subplot(gs[0]); ax.set_anchor('N'); tag(ax, 'a', -0.05)
    rng = np.random.default_rng(7)
    # schematic tissue with spatially clustered cell types
    xs = rng.uniform(0, 10, 2600); ys = rng.uniform(0, 6, 2600)
    imm = (((xs - 2.4) ** 2 + (ys - 4.3) ** 2) < 1.6) | (((xs - 7.6) ** 2 + (ys - 1.7) ** 2) < 1.1)
    imm = imm & (rng.random(len(xs)) < 0.75) | (rng.random(len(xs)) < 0.03)
    ax.scatter(xs[~imm], ys[~imm], s=1.2, color='#c3c2b7', lw=0)
    ax.scatter(xs[imm], ys[imm], s=1.6, color=ORANGE, lw=0)
    for (x0, y0) in [(1.4, 3.3), (5.6, 3.8), (6.6, 0.7)]:
        ax.add_patch(Rectangle((x0, y0), 2.0, 2.0, fill=False, ec=BLUE, lw=1.1))
    pick = rng.choice(len(xs), 70, replace=False)
    ax.scatter(xs[pick], ys[pick], s=7, facecolor='none', edgecolor=INK, lw=0.5)
    ax.set_xlim(-0.2, 10.2); ax.set_ylim(-1.9, 6.2); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0, -0.45, 'Blue squares: 3 fields (clusters of ~1,400 cells)', color=BLUE, fontsize=6.2, va='top')
    ax.text(0, -1.05, 'Black circles: equal number of random cells', color=INK, fontsize=6.2, va='top')
    ax.text(0, -1.65, 'Orange: spatially clustered cell type', color='#b4471b', fontsize=6.2, va='top')
    ax.set_title('Fields are clusters, not independent cells', loc='left', color=INK)
    # panel b: per-tissue TV
    ax = fig.add_subplot(gs[1]); tag(ax, 'b', -0.1)
    p = pd.read_csv(D + 'sampling_by_sample.tsv', sep='\t'); p = p[p.k == 3]
    a = p.pivot(index='sample', columns='strategy', values='TV_median')
    rows = [('GSE325587', s, a.loc[s, 'random_FOV'], a.loc[s, 'random_cells']) for s in a.index]
    g = pd.read_csv(D + 'GSE211785_sampling_summary.tsv', sep='\t')
    rows += [('GSE211785', r.sample, r.region_TV, r.cell_TV) for r in g.itertuples()]
    l = pd.read_csv(D + 'sampling_summary.tsv', sep='\t'); l = l[l.k == 3]
    rows += [('Lupus', r.tissue, r.region_median_tv, r.cell_median_tv) for r in l.itertuples()]
    t = pd.DataFrame(rows, columns=['source', 'tissue', 'field', 'cell'])
    order = ['GSE325587', 'GSE211785', 'Lupus']
    t['o'] = t.source.map({s: i for i, s in enumerate(order)})
    t = t.sort_values(['o', 'field']).reset_index(drop=True)
    xpos = []; x = 0; prev = None; centers = {}
    for r in t.itertuples():
        if prev is not None and r.source != prev:
            x += 1.6
        xpos.append(x); centers.setdefault(r.source, []).append(x); prev = r.source; x += 1
    t['x'] = xpos
    for r in t.itertuples():
        ax.plot([r.x, r.x], [r.cell, r.field], color='#c3c2b7', lw=0.8, zorder=1)
    ax.scatter(t.x, t.field, s=12, color=ORANGE, zorder=3, lw=0.4, edgecolor='white', label='3 random fields')
    ax.scatter(t.x, t.cell, s=12, color=BLUE, zorder=3, lw=0.4, edgecolor='white', label='Matched random cells')
    for s, xs_ in centers.items():
        ax.text(np.mean(xs_), -0.006, SRC_NAME[s] + f'\n(n = {len(xs_)})', ha='center', va='top', fontsize=6.3, color=INK2)
    ax.set_xticks([]); ax.spines['bottom'].set_visible(False)
    ax.set_ylim(0, 0.18); ax.set_xlim(-1, t.x.max() + 1)
    ax.set_ylabel('Median composition error (TV)')
    ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend(loc='upper left', ncol=2, handletextpad=0.2, columnspacing=1.0)
    ax.set_title('Field sampling exceeded matched-cell error in all 33 tissues', loc='left', color=INK)
    save(fig, 'Figure1')


# ---------------- Figure 2: design-effect atlas -----------------
def fig2():
    a = pd.read_csv(R + 'A1_design_effect_atlas.tsv', sep='\t')
    lab = a[a.level == 'label']
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.25, 1], hspace=0.62, wspace=0.42)
    # a heatmap primary
    ax = fig.add_subplot(gs[0, :]); tag(ax, 'a', -0.03, 1.03)
    pr = lab[lab.source == 'GSE325587']
    keep = pr.groupby('label').p.median()
    keep = keep[keep >= 0.005].sort_values(ascending=False).index
    H = pr.pivot(index='label', columns='tissue', values='D').loc[keep]
    tis_order = pr.groupby('tissue').D.median().sort_values().index
    H = H[tis_order]
    im = ax.imshow(H.to_numpy(), aspect='auto', cmap=SEQ, norm=LogNorm(1, 300), interpolation='nearest')
    ax.set_yticks(range(len(H))); ax.set_yticklabels([f'{l}  ({keep_p:.1%})' for l, keep_p in zip(H.index, pr.groupby("label").p.median().loc[keep])], fontsize=5.6)
    ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns, rotation=60, ha='right', fontsize=5.8)
    for s in ax.spines.values(): s.set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label('Design effect D (log scale)', fontsize=6.3)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.8)
    ax.set_title('Design effect by cell type and tissue, GSE325587 (cell types with median proportion ≥0.5%)', loc='left', color=INK)
    # b D vs p
    ax = fig.add_subplot(gs[1, 0]); tag(ax, 'b', -0.2)
    for s in ['GSE325587', 'Lupus', 'GSE211785']:
        g = lab[lab.source == s]
        ax.scatter(g.p, g.D, s=6, color=SRC_COL[s], alpha=0.75, lw=0.25, edgecolor='white', label=SRC_NAME[s])
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-4, 1); ax.set_ylim(0.3, 500)
    ax.axhline(1, color=INK2, lw=0.6, ls='--'); ax.text(0.9, 0.62, 'D = 1: cells independent', fontsize=5.8, color=INK2, ha='right')
    ax.set_xlabel('Cell-type proportion in tissue'); ax.set_ylabel('Design effect D')
    ax.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.legend(loc='upper left', handletextpad=0.1, markerscale=1.6)
    ax.set_title('Abundant types are most penalized', loc='left', color=INK)
    # c lineage D by source
    ax = fig.add_subplot(gs[1, 1]); tag(ax, 'c', -0.42)
    ln = a[a.level == 'lineage']
    lins = ['Other epithelial', 'Immune', 'Stromal', 'Endothelial', 'Glomerular epithelial']
    off = {'GSE325587': -0.22, 'Lupus': 0.0, 'GSE211785': 0.22}
    for s in ['GSE325587', 'Lupus', 'GSE211785']:
        for i, L in enumerate(lins):
            v = ln[(ln.source == s) & (ln.label == L)].D
            y = i + off[s]
            ax.scatter(v, np.full(len(v), y) + np.random.default_rng(i).uniform(-0.05, 0.05, len(v)), s=4, color=SRC_COL[s], alpha=0.55, lw=0)
            ax.plot([v.median()] * 2, [y - 0.09, y + 0.09], color=SRC_COL[s], lw=2.0, solid_capstyle='round')
    ax.set_xscale('log'); ax.set_yticks(range(len(lins))); ax.set_yticklabels(lins); ax.invert_yaxis()
    ax.axvline(1, color=INK2, lw=0.6, ls='--'); ax.set_xlabel('Design effect D (lineage)')
    ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_xlim(0.25, 500)
    ax.set_title('Lineages (colours as in b); bars = medians', loc='left', color=INK)
    save(fig, 'Figure2')


# ---------------- Figure 3: false positives -----------------
TESTS = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']
TNAME = {'cell_chisq': 'Cells as replicates (χ²)', 'field_wilcoxon': 'Fields as replicates (Wilcoxon)',
         'field_quasibin': 'Fields as replicates (quasi-binomial)', 'field_cluster': 'Field model, individual-clustered SE',
         'indiv_wilcoxon': 'Individuals as replicates (Wilcoxon)', 'indiv_quasibin': 'Individuals as replicates (quasi-binomial)'}
TCOL = {'cell_chisq': ORANGE, 'field_wilcoxon': BLUE, 'field_quasibin': BLUE, 'field_cluster': BLUE,
        'indiv_wilcoxon': AQUA, 'indiv_quasibin': AQUA}


def fig3():
    f = pd.read_csv(R + 'A3_false_positive_rates.tsv', sep='\t')
    fam = pd.read_csv(R + 'A3_familywise_any_discovery.tsv', sep='\t')
    real = pd.read_csv(R + 'A3_real_discoveries.tsv', sep='\t')
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.35, 1, 0.8], wspace=0.12))
    ys = np.arange(len(TESTS))[::-1]
    # a FPR
    ax = axs[0]; tag(ax, 'a', -1.02)
    for s, mk, dy in [('GSE325587', 'o', 0.13), ('Lupus', 's', -0.13)]:
        for scen, fill in [('all_fields', True), ('3_fields', False)]:
            v = f[(f.source == s) & (f.scenario == scen)].iloc[0]
            for t, y in zip(TESTS, ys):
                ax.scatter(v[t], y + dy, marker=mk, s=16, color=TCOL[t] if fill else 'white', edgecolor=TCOL[t], lw=0.8, zorder=3)
    ax.axvline(0.05, color=INK2, lw=0.7, ls='--'); ax.text(0.065, ys[-1] + 0.95, 'nominal 5%', fontsize=5.8, color=INK2)
    ax.set_yticks(ys); ax.set_yticklabels([TNAME[t] for t in TESTS]); ax.set_xlim(0, 1)
    ax.set_xlabel('False-positive rate (null splits)'); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title('Per-test error rate', loc='left', color=INK)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker='o', ls='', color=INK2, mfc=INK2, ms=4, label='DKD, all fields'),
          Line2D([], [], marker='s', ls='', color=INK2, mfc=INK2, ms=4, label='Lupus, all fields'),
          Line2D([], [], marker='o', ls='', color=INK2, mfc='white', ms=4, label='3 fields per tissue')]
    ax.legend(handles=hs, loc='lower right', handletextpad=0.1, borderaxespad=0.2)
    # b family-wise
    ax = axs[1]; tag(ax, 'b', -0.08)
    for s, mk, dy in [('GSE325587', 'o', 0.13), ('Lupus', 's', -0.13)]:
        for scen, fill in [('all_fields', True), ('3_fields', False)]:
            v = fam[(fam.source == s) & (fam.scenario == scen)].iloc[0]
            for t, y in zip(TESTS, ys):
                ax.scatter(v[t], y + dy, marker=mk, s=16, color=TCOL[t] if fill else 'white', edgecolor=TCOL[t], lw=0.8, zorder=3)
    ax.set_yticks(ys); ax.set_yticklabels([]); ax.set_xlim(0, 1.02)
    ax.set_xlabel('Null splits with ≥1 BH discovery'); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title('Spurious "findings" per study', loc='left', color=INK)
    # c real contrast
    ax = axs[2]; tag(ax, 'c', -0.08)
    for s, dy in [('GSE325587', 0.17), ('Lupus', -0.17)]:
        v = real[real.source == s].iloc[0]
        for t, y in zip(TESTS, ys):
            ax.barh(y + dy, v[t + '_q'], height=0.3, color=TCOL[t] if s == 'GSE325587' else 'white', edgecolor=TCOL[t], lw=0.8)
            ax.text(v[t + '_q'] + 0.4, y + dy, str(int(v[t + '_q'])), va='center', fontsize=5.6, color=INK2)
    ax.set_yticks(ys); ax.set_yticklabels([]); ax.set_xlim(0, 28)
    ax.set_xlabel('BH-significant categories')
    ax.text(0.98, 0.02, 'filled: DKD vs control (21)\nopen: SLE vs control (25)', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.6, color=INK2)
    ax.set_title('Real disease contrasts', loc='left', color=INK)
    save(fig, 'Figure3')


# ---------------- Figure 4: interval methods -----------------
def fig4():
    t = pd.read_csv(R + 'A4_model_interval_by_tissue.tsv', sep='\t'); t = t[t.k == 3]
    meths = ['cell_Wilson', 'region_ratio_t', 'Fieller', 'quasibinomial']
    mname = {'cell_Wilson': 'Cell-level Wilson', 'region_ratio_t': 'Region-ratio t', 'Fieller': 'Fieller', 'quasibinomial': 'Quasi-binomial\nfield model'}
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw=dict(wspace=0.08))
    ys = np.arange(len(meths))[::-1]
    off = {'GSE325587': 0.2, 'Lupus': 0.0, 'GSE211785': -0.2}
    for j, (col, xl) in enumerate([('coverage', 'Coverage of nominal 95% interval'), ('width', 'Mean interval width (log scale)')]):
        ax = axs[j]; tag(ax, 'ab'[j], -0.3 if j == 0 else -0.05)
        for s in ['GSE325587', 'Lupus', 'GSE211785']:
            for m, y in zip(meths, ys):
                v = t[(t.source == s) & (t.method == m)][col]
                yy = y + off[s]
                ax.scatter(v, np.full(len(v), yy), s=6, color=SRC_COL[s], alpha=0.55, lw=0)
                ax.plot([v.median()] * 2, [yy - 0.08, yy + 0.08], color=SRC_COL[s], lw=2.2, solid_capstyle='round')
        ax.set_yticks(ys); ax.set_yticklabels([mname[m] for m in meths] if j == 0 else [])
        ax.set_xlabel(xl); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
        if j == 0:
            ax.axvline(0.95, color=INK2, lw=0.7, ls='--'); ax.set_xlim(0, 1.02)
            ax.set_title('Immune proportion, three fields per draw', loc='left', color=INK)
        else:
            ax.set_xscale('log')
            for s, yy in zip(['GSE325587', 'Lupus', 'GSE211785'], [0.97, 0.90, 0.83]):
                ax.text(0.99, yy, SRC_NAME[s], transform=ax.transAxes, ha='right', fontsize=5.8, color=SRC_COL[s] if s != 'Lupus' else '#12855c')
            ax.set_title('Bars = medians across tissues', loc='left', color=INK)
    save(fig, 'Figure4')


# ---------------- Figure 5: fovplan -----------------
def fig5():
    v = pd.read_csv(R + 'A2_fovplan_validation_summary.tsv', sep='\t'); v = v[v.source == 'All']
    req = pd.read_csv(R + 'A2_fields_required_lineage.tsv', sep='\t')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1, 1.25], wspace=0.45))
    ax = axs[0]; tag(ax, 'a', -0.2)
    ks = [3, 5, 10]; x = np.arange(3)
    for lab, col, mk, c in [('Cells independent (D = 1)', ORANGE, 'o', 'naive_mean'), ('fovplan, pilot D', BLUE, 's', 'fovplan_mean'),
                            ('fovplan, 80% assurance', AQUA, 'D', 'fovplan80_mean')]:
        ax.plot(x, v[c], color=col, lw=2, marker=mk, ms=4.5, label=lab)
        ypos = {'naive_mean': v[c].iloc[1] + 0.05, 'fovplan_mean': v[c].iloc[1] - 0.06, 'fovplan80_mean': v[c].iloc[1] + 0.045}[c]
        ax.text(1.0, ypos, lab, color=col if col != AQUA else '#12855c', fontsize=5.8, ha='center', va='center')
    ax.axhline(0.95, color=INK2, lw=0.7, ls='--'); ax.text(0.3, 0.905, 'target 95%', fontsize=5.8, color=INK2, ha='center')
    ax.set_xticks(x); ax.set_xticklabels([f'{k} fields' for k in ks]); ax.set_ylim(0, 1.05); ax.set_xlim(-0.2, 2.1)
    ax.set_ylabel('Share of held-out samples within\npredicted 95% half-width')
    ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title('Held-out validation (5 pilot fields)', loc='left', color=INK)
    ax = axs[1]; tag(ax, 'b', -0.42)
    lins = ['Immune', 'Stromal', 'Endothelial', 'Glomerular epithelial', 'Other epithelial']
    off = {'GSE325587': 0.22, 'Lupus': 0.0, 'GSE211785': -0.22}
    ys = np.arange(len(lins))[::-1]
    for s in ['GSE325587', 'Lupus', 'GSE211785']:
        for L, y in zip(lins, ys):
            r = req[(req.source == s) & (req.cell_type == L)].iloc[0]
            ax.barh(y + off[s], r.rel25, height=0.2, color=SRC_COL[s])
            ax.text(r.rel25 * 1.08, y + off[s], f'{r.rel25:g}', va='center', fontsize=5.4, color=INK2)
    ax.axvline(1, color=ORANGE, lw=0.8, ls=':')
    ax.text(1.1, ys[0] + 0.55, 'cell-count planning: 1 field', fontsize=5.6, color='#b4471b')
    ax.set_xscale('log'); ax.set_xlim(0.8, 300); ax.set_yticks(ys); ax.set_yticklabels(lins)
    ax.set_xlabel('Fields required for ±25% relative precision (median tissue)')
    ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    for s, yy in zip(['GSE325587', 'Lupus', 'GSE211785'], [0.17, 0.10, 0.03]):
        ax.text(0.99, yy, SRC_NAME[s], transform=ax.transAxes, ha='right', fontsize=5.8, color=SRC_COL[s] if s != 'Lupus' else '#12855c')
    ax.set_title('Fields needed per lineage', loc='left', color=INK)
    save(fig, 'Figure5')


if __name__ == '__main__':
    for f in (fig1, fig2, fig3, fig4, fig5):
        f()
    print('done')
