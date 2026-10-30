"""Extract selected columns from a byte range of the Hickey CSV (line-aligned).
usage: extract_range.py SRC START END OUT"""
import sys, io, pandas as pd
src, start, end, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
with open(src, 'rb') as f:
    header = f.readline()
    if start == 0:
        start = len(header)
    f.seek(start - 1)
    if f.read(1) != b'\n':
        f.readline()
    s = f.tell()
    f.seek(end - 1)
    if f.read(1) != b'\n':
        f.readline()
    e = f.tell()
    f.seek(s); data = f.read(e - s)
cols = ['Unnamed: 0', 'x', 'y', 'Tissue_location', 'tissue', 'donor', 'unique_region', 'region', 'Cell Type', 'Tissue Unit']
d = pd.read_csv(io.BytesIO(header + data), usecols=cols, dtype={'Tissue Unit': str, 'region': str})
d = d.rename(columns={'Unnamed: 0': 'row', 'Cell Type': 'cell_type', 'Tissue Unit': 'tissue_unit'})
d.to_csv(out, sep='\t', index=False, compression='gzip')
print(out, len(d))
