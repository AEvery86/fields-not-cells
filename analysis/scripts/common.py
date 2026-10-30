"""Shared loaders for the G20 analyses (three annotated human CosMx sources).

Each loader returns a long table of per-field counts with columns
source, tissue, patient, group, field, label, lineage, count
and a wide helper is provided. Group labels are used only for the differential-abundance analysis.
"""
import json, os
import numpy as np, pandas as pd

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')

LUPUS_LINEAGE = {
    'Ascending.vasa.recta.endothelium': 'Endothelial', 'Peritubular.capillary.endothelium.1': 'Endothelial',
    'Peritubular.capillary.endothelium.2': 'Endothelial', 'Descending.vasa.recta.endothelium': 'Endothelial',
    'Glomerular.endothelium': 'Endothelial', 'Fibroblast': 'Stromal', 'Myofibroblast': 'Stromal',
    'Vascular.pericyte': 'Stromal', 'Mesangial.cell': 'Stromal', 'Podocyte': 'Glomerular epithelial',
    'Parietal.epithelium': 'Glomerular epithelial'}
GSE211785_LINEAGE = {'DCT': 'Other epithelial', 'Endothelium': 'Endothelial', 'Fibro': 'Stromal',
                     'IC/DCT': 'Other epithelial', 'Immune': 'Immune', 'Injured TAL': 'Other epithelial',
                     'PC': 'Other epithelial', 'PEC': 'Glomerular epithelial', 'PT': 'Other epithelial',
                     'Podo': 'Glomerular epithelial', 'TAL': 'Other epithelial', 'VSMC/Mes': 'Stromal',
                     'iPT': 'Other epithelial'}


def load_primary():
    w = pd.read_csv(f'{DATA}/fov_celltype_counts.tsv', sep='\t')
    m = pd.read_csv(f'{DATA}/sample_mapping.tsv', sep='\t').set_index('sample')
    lin = pd.read_csv(f'{DATA}/harmonized_category_mapping.tsv', sep='\t')
    lin = dict(zip(lin.original_label, lin.coarse))
    w['field'] = w.slide.astype(str) + ':' + w.fov.astype(str)
    long = w.drop(columns=['slide', 'fov']).melt(id_vars=['sample', 'field'], var_name='label', value_name='count')
    long = long.rename(columns={'sample': 'tissue'})
    long['patient'] = long.tissue
    long['group'] = np.where(m.loc[long.tissue, 'stage'].values == 'Control', 'Control', 'DKD')
    long['stage'] = m.loc[long.tissue, 'stage'].values
    long['lineage'] = long.label.map(lin)
    long['source'] = 'GSE325587'
    return long


def load_gse211785():
    w = pd.read_csv(f'{DATA}/GSE211785_fov_celltype_counts.tsv', sep='\t')
    w['field'] = w.slide.astype(str) + ':' + w.fov.astype(str)
    long = w.drop(columns=['slide', 'fov']).melt(id_vars=['sample', 'field'], var_name='label', value_name='count')
    long = long.rename(columns={'sample': 'tissue'})
    long['patient'] = long.tissue; long['group'] = 'NA'; long['stage'] = 'NA'
    long['lineage'] = long.label.map(GSE211785_LINEAGE)
    long['source'] = 'GSE211785'
    return long


def load_lupus():
    d = pd.read_csv(f'{DATA}/Danaher_fov_label_counts.tsv', sep='\t')
    d = d.rename(columns={'fovkey': 'field'})
    d['patient'] = d.tissue.str.replace(r'^(SLE8)\.\d$', r'\1', regex=True)
    d['group'] = np.where(d.tissue.str.startswith('Control'), 'Control', 'SLE')
    d['stage'] = d.group
    imm = d.is_immune.astype(str).str.upper().eq('TRUE')
    d['lineage'] = np.where(imm, 'Immune', d.label.map(LUPUS_LINEAGE).fillna('Other epithelial'))
    d['source'] = 'Lupus'
    return d[['tissue', 'field', 'label', 'count', 'patient', 'group', 'stage', 'lineage', 'source']]


def load_all():
    return pd.concat([load_primary(), load_gse211785(), load_lupus()], ignore_index=True)


def wide(long, unit='tissue'):
    """Return dict unit -> (labels, fields x labels count array, field ids)."""
    out = {}
    for u, g in long.groupby(unit):
        t = g.pivot_table(index='field', columns='label', values='count', aggfunc='sum', fill_value=0)
        out[u] = (list(t.columns), t.to_numpy(float), list(t.index))
    return out


def immune_labels(long):
    return sorted(long.loc[long.lineage == 'Immune', 'label'].unique())
