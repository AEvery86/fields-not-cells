"""Assemble the extracted Hickey et al. 2023 CODEX table (seven complete intestinal sites).
Input: data/hickey_extract/part00-06.tsv.gz (columns row, x, y, Tissue_location, tissue, donor, unique_region,
region, cell_type, tissue_unit), produced from the Dryad CSV by extract_range.py. The transverse colon is
excluded because the downloaded file was truncated within that site (see manuscript Methods)."""
import glob, os, pandas as pd
here = os.path.dirname(__file__)
fs = sorted(glob.glob(os.path.join(here, '..', 'data', 'hickey_extract', 'part0*.tsv.gz')))
d = pd.concat([pd.read_csv(f, sep='\t', dtype={'tissue_unit': str, 'region': str}) for f in fs], ignore_index=True)
d = d.dropna(subset=['unique_region', 'cell_type', 'x', 'y'])
d = d[d.Tissue_location != 'Transverse']
d.to_pickle(os.path.join(here, '..', 'data', 'hickey_7sites.pkl'))
print(len(d), d.unique_region.nunique())
