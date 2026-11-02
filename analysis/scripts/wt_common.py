import os
"""Shared loaders for whole-tissue references (mouse kidney sections, human intestine regions)."""
import numpy as np, pandas as pd

PX = 0.37744  # assumed CODEX um per pixel (see Methods)


def load_cells():
    d = pd.read_pickle(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hickey_7sites.pkl')).dropna(subset=['cell_type', 'unique_region', 'x', 'y'])
    cells = [pd.DataFrame(dict(section=d.unique_region, donor=d.donor, organ=d.tissue.map({'SB': 'human small intestine', 'CL': 'human colon'}),
                               x_um=d.x * PX, y_um=d.y * PX, label=d.cell_type))]
    for s in ['ShamR', 'Day14R']:
        m = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', f'mouse_example_{s}.tsv.gz'), sep='\t')
        cells.append(pd.DataFrame(dict(section=s, donor=s, organ='mouse kidney', x_um=m.x, y_um=m.y, label=m.celltype)))
    return pd.concat(cells, ignore_index=True)


def field_tables(cells, F=500.0, min_cells=100, min_fields=20):
    """Yield (section, organ, truth Series, W DataFrame fields x labels, fx, fy) for every section."""
    for sec, g in cells.groupby('section'):
        truth = g.label.value_counts(normalize=True)
        fx = np.floor((g.x_um - g.x_um.min()) / F).astype(int)
        fy = np.floor((g.y_um - g.y_um.min()) / F).astype(int)
        key = fx * 100000 + fy
        W = pd.crosstab(key, g.label).reindex(columns=truth.index, fill_value=0)
        W = W[W.sum(axis=1) >= min_cells]
        if len(W) < min_fields:
            continue
        idx = W.index.to_numpy()
        yield sec, g.organ.iloc[0], truth, W, idx // 100000, idx % 100000
