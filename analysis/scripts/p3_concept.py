"""G23: conceptual Figure 1 and graphical abstract (unified style)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from a5_figures import BLUE, ORANGE, AQUA, INK, INK2, GRID, tag, save, SRC_NAME
D = '../data/'
AQ_TXT = '#12855c'; GREY = '#c3c2b7'; LIGHT = '#f0efec'


def tissue(ax, seed=7, show_random=True, title=None):
    rng = np.random.default_rng(seed)
    xs = rng.uniform(0, 10, 2600); ys = rng.uniform(0, 6, 2600)
    clus = (((xs - 2.4) ** 2 + (ys - 4.3) ** 2) < 1.6) | (((xs - 7.6) ** 2 + (ys - 1.7) ** 2) < 1.1)
    clus = clus & (rng.random(len(xs)) < 0.75) | (rng.random(len(xs)) < 0.03)
    ax.scatter(xs[~clus], ys[~clus], s=1.2, color=GREY, lw=0)
    ax.scatter(xs[clus], ys[clus], s=1.6, color=ORANGE, lw=0)
    for (x0, y0) in [(1.4, 3.3), (5.6, 3.8), (6.6, 0.7)]:
        ax.add_patch(Rectangle((x0, y0), 2.0, 2.0, fill=False, ec=BLUE, lw=1.2))
    if show_random:
        pick = rng.choice(len(xs), 70, replace=False)
        ax.scatter(xs[pick], ys[pick], s=7, facecolor='none', edgecolor=INK, lw=0.5)
    ax.set_xlim(-0.2, 10.2); ax.set_ylim(-0.2, 6.2); ax.set_aspect('equal'); ax.axis('off')
    if title:
        ax.set_title(title, loc='left', color=INK)


def box(ax, x, y, w, h, text, fc, ec, tc=INK, fs=6.2, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.012,rounding_size=0.025', fc=fc, ec=ec, lw=0.8, transform=ax.transAxes))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs, color=tc, transform=ax.transAxes, fontweight='bold' if bold else None, linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, c=INK2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes, arrowstyle='-|>', mutation_scale=7, color=c, lw=0.9))


def fig1():
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], width_ratios=[1.05, 0.95, 1.25], hspace=0.38, wspace=0.28)
    # a schematic
    ax = fig.add_subplot(gs[0, 0]); tag(ax, 'a', -0.10, 1.02)
    tissue(ax, title='Fields are clusters')
    ax.text(0, -0.9, '■', color=BLUE, fontsize=7); ax.text(0.45, -0.9, '3 fields, ~1,600 cells each', fontsize=5.8, color=INK2)
    ax.text(0, -1.6, '○', color=INK, fontsize=7); ax.text(0.45, -1.6, 'same number of random cells', fontsize=5.8, color=INK2)
    ax.text(0, -2.3, '●', color=ORANGE, fontsize=7); ax.text(0.45, -2.3, 'spatially clustered cell type', fontsize=5.8, color=INK2)
    ax.set_ylim(-2.6, 6.2)
    # b information content
    ax = fig.add_subplot(gs[0, 1]); tag(ax, 'b', -0.2, 1.02)
    ax.barh([1], [4886], color=BLUE, height=0.5); ax.barh([0], [321], color=AQUA, height=0.5)
    ax.text(4886 * 0.98, 1, '4,886 cells\nin 3 fields', ha='right', va='center', fontsize=5.8, color='white')
    ax.text(321 + 120, 0, '≈ 321 random cells\n(effective sample size)', ha='left', va='center', fontsize=5.8, color=AQ_TXT)
    ax.set_yticks([]); ax.set_xlim(0, 5200); ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel('Cells'); ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.spines['left'].set_visible(False)
    ax.text(0.30, 0.47, 'D ≈ 1 + (m − 1)ρ\nm ≈ 1,600, ρ ≈ 0.01 → D ≈ 15', transform=ax.transAxes, fontsize=6.0, color=INK, va='center')
    ax.set_title('Same cells, less information', loc='left', color=INK)
    # c workflow
    ax = fig.add_subplot(gs[0, 2]); tag(ax, 'c', -0.10, 1.02); ax.axis('off')
    ax.set_title('Consequences and remedy', loc='left', color=INK)
    box(ax, 0.00, 0.70, 0.44, 0.22, 'Cells as replicates\n88–94% false positives\nintervals cover 25–37%', '#fde3d8', ORANGE, tc='#8a3412', fs=5.4)
    box(ax, 0.52, 0.70, 0.44, 0.22, 'Fields or individuals\nas units\n3–6% false positives', '#d9f2e8', AQUA, tc=AQ_TXT, fs=5.4)
    arrow(ax, 0.455, 0.81, 0.51, 0.81)
    box(ax, 0.00, 0.36, 0.28, 0.22, 'Pilot\n5 fields', LIGHT, INK2, fs=5.6)
    box(ax, 0.34, 0.36, 0.28, 0.22, 'Estimate D\nper cell type', '#dbe8fa', BLUE, tc='#184f95', fs=5.6)
    box(ax, 0.68, 0.36, 0.28, 0.22, 'fovplan\nfields needed', '#d9f2e8', AQUA, tc=AQ_TXT, fs=5.6, bold=True)
    arrow(ax, 0.29, 0.47, 0.335, 0.47); arrow(ax, 0.63, 0.47, 0.675, 0.47)
    box(ax, 0.00, 0.02, 0.96, 0.22, 'Validated against held-out fields, simulated truth\nand whole-tissue references in kidney and intestine', LIGHT, GREY, tc=INK2, fs=5.3)
    arrow(ax, 0.82, 0.355, 0.82, 0.245)
    # d dumbbell
    ax = fig.add_subplot(gs[1, :]); tag(ax, 'd', -0.045, 1.02)
    p = pd.read_csv(D + 'sampling_by_sample.tsv', sep='\t'); p = p[p.k == 3]
    a = p.pivot(index='sample', columns='strategy', values='TV_median')
    rows = [('GSE325587', s, a.loc[s, 'random_FOV'], a.loc[s, 'random_cells']) for s in a.index]
    g = pd.read_csv(D + 'GSE211785_sampling_summary.tsv', sep='\t')
    rows += [('GSE211785', r.sample, r.region_TV, r.cell_TV) for r in g.itertuples()]
    l = pd.read_csv(D + 'sampling_summary.tsv', sep='\t'); l = l[l.k == 3]
    rows += [('Lupus', r.tissue, r.region_median_tv, r.cell_median_tv) for r in l.itertuples()]
    t = pd.DataFrame(rows, columns=['source', 'tissue', 'field', 'cell'])
    order = ['GSE325587', 'GSE211785', 'Lupus']
    t['o'] = t.source.map({s: i for i, s in enumerate(order)}); t = t.sort_values(['o', 'field']).reset_index(drop=True)
    xpos = []; x = 0; prev = None; centers = {}
    for r in t.itertuples():
        if prev is not None and r.source != prev: x += 1.6
        xpos.append(x); centers.setdefault(r.source, []).append(x); prev = r.source; x += 1
    t['x'] = xpos
    for r in t.itertuples():
        ax.plot([r.x, r.x], [r.cell, r.field], color=GREY, lw=0.8, zorder=1)
    ax.scatter(t.x, t.field, s=12, color=ORANGE, zorder=3, lw=0.4, edgecolor='white', label='3 random fields')
    ax.scatter(t.x, t.cell, s=12, color=BLUE, zorder=3, lw=0.4, edgecolor='white', label='Matched random cells')
    for s, xs_ in centers.items():
        ax.text(np.mean(xs_), -0.007, SRC_NAME[s] + f' (n = {len(xs_)})', ha='center', va='top', fontsize=6.0, color=INK2)
    ax.set_xticks([]); ax.spines['bottom'].set_visible(False); ax.set_ylim(0, 0.16); ax.set_xlim(-1, t.x.max() + 1)
    ax.set_ylabel('Median composition\nerror (TV)'); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend(loc='upper left', ncol=2, handletextpad=0.2, columnspacing=1.0)
    ax.set_title('Field sampling exceeded matched-cell error in all 33 human kidney tissues', loc='left', color=INK)
    save(fig, 'Figure1')


def graphical_abstract():
    fig = plt.figure(figsize=(7.0, 3.2))
    ax0 = fig.add_axes([0.01, 0.12, 0.30, 0.74]); tissue(ax0, show_random=False)
    fig.text(0.02, 0.92, '1  A few fields of view', fontsize=9, fontweight='bold', color=INK)
    fig.text(0.02, 0.06, 'Cells within a field are correlated:\nfields, not cells, are the sampling units', fontsize=6.6, color=INK2)
    ax = fig.add_axes([0.38, 0.20, 0.26, 0.58])
    ax.barh([1], [4886], color=BLUE, height=0.55); ax.barh([0], [321], color=AQUA, height=0.55)
    ax.text(4800, 1, '4,886 cells', ha='right', va='center', color='white', fontsize=7)
    ax.text(420, 0, '≈ 321 effective', ha='left', va='center', color=AQ_TXT, fontsize=7)
    ax.set_xlim(0, 5200); ax.set_yticks([]); ax.spines['left'].set_visible(False); ax.tick_params(labelsize=6)
    ax.set_xlabel('Cells in 3 fields', fontsize=6.6)
    fig.text(0.37, 0.92, '2  Design effect ≈ 15–60', fontsize=9, fontweight='bold', color=INK)
    fig.text(0.37, 0.06, 'Kidney (CosMx, Xenium) and intestine (CODEX)', fontsize=6.6, color=INK2)
    fig.text(0.70, 0.92, '3  Consequences and remedy', fontsize=9, fontweight='bold', color=INK)
    ax2 = fig.add_axes([0.70, 0.10, 0.28, 0.76]); ax2.axis('off')
    box(ax2, 0.0, 0.70, 0.97, 0.24, 'Cells as replicates: 87–94% false positives\nCell-level intervals cover 20–37%', '#fde3d8', ORANGE, tc='#8a3412', fs=6.6)
    box(ax2, 0.0, 0.38, 0.97, 0.24, 'Individuals as units: 3–6% false positives\nField-model intervals cover 91–94%', '#d9f2e8', AQUA, tc=AQ_TXT, fs=6.6)
    box(ax2, 0.0, 0.04, 0.97, 0.26, 'fovplan: pilot fields → number of fields\nfor a target precision (≈95% within\ntarget with the 90% assurance plan)', '#dbe8fa', BLUE, tc='#184f95', fs=6.6, bold=True)
    fig.savefig('../../figures/main/Graphical_Abstract.png', dpi=400, facecolor='white')
    fig.savefig('../../figures/main/Graphical_Abstract.pdf', facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    fig1(); graphical_abstract(); print('done')
