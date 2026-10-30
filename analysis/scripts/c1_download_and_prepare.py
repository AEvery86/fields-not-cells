"""G21-C1. Download and prepare whole-section datasets once network access to the hosts is allowed.

Run:  python c1_download_and_prepare.py hickey        (Dryad, ~3 GB CSV; human small intestine + colon, CODEX)
      python c1_download_and_prepare.py check         (only tests whether hosts are reachable)
Every download is logged with URL, size, ETag/Last-Modified and SHA256 in ../data/raw/download_log.json.
Prepared cell tables follow the c2_whole_section.py schema: section, donor, organ, x_um, y_um, label.
"""
import sys, os, json, hashlib, time
import pandas as pd

RAW = os.path.join(os.path.dirname(__file__), '..', 'data', 'raw'); os.makedirs(RAW, exist_ok=True)
LOG = os.path.join(RAW, 'download_log.json')
HOSTS = ['https://datadryad.org', 'https://ndownloader.figshare.com', 'https://cf.10xgenomics.com', 'https://zenodo.org', 'https://ftp.ncbi.nlm.nih.gov']


def log(entry):
    d = json.load(open(LOG)) if os.path.exists(LOG) else []
    d.append(entry); json.dump(d, open(LOG, 'w'), indent=1)


def fetch(url, dest):
    import requests
    h = hashlib.sha256(); n = 0
    with requests.get(url, stream=True, timeout=60, allow_redirects=True) as r:
        r.raise_for_status()
        with open(dest + '.part', 'wb') as f:
            for chunk in r.iter_content(1 << 22):
                f.write(chunk); h.update(chunk); n += len(chunk)
        meta = dict(url=url, final_url=r.url, bytes=n, sha256=h.hexdigest(), etag=r.headers.get('ETag'),
                    last_modified=r.headers.get('Last-Modified'), accessed=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    os.replace(dest + '.part', dest); log(meta); return meta


def check():
    import requests
    for u in HOSTS:
        try:
            print(u, requests.head(u, timeout=15, allow_redirects=True).status_code)
        except Exception as e:
            print(u, 'BLOCKED', type(e).__name__)


def hickey():
    """Hickey et al. 2023 Nature, CODEX human intestine (Dryad doi:10.5061/dryad.pk0p2ngrf, 2023-09-13 version)."""
    import requests
    api = 'https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.pk0p2ngrf'
    ds = requests.get(api, timeout=60).json()
    ver = requests.get('https://datadryad.org' + ds['_links']['stash:version']['href'], timeout=60).json()
    files = requests.get('https://datadryad.org' + ver['_links']['stash:files']['href'], timeout=60).json()
    out = {}
    for f in files['_embedded']['stash:files']:
        name = f['path']
        dest = os.path.join(RAW, name)
        if not os.path.exists(dest):
            print('downloading', name, f.get('size'))
            fetch('https://datadryad.org' + f['_links']['stash:download']['href'], dest)
        out[name] = dest
    csv = [p for n, p in out.items() if n.endswith('.csv') and 'donor' not in n][0]
    cols = pd.read_csv(csv, nrows=5).columns.tolist(); print('columns:', cols)

    def pick(*cands):
        for c in cands:
            for col in cols:
                if col.lower() == c.lower():
                    return col
        raise KeyError(f'none of {cands} in {cols}')
    cx, cy = pick('x', 'X', 'x_um'), pick('y', 'Y', 'y_um')
    ct = pick('Cell Type', 'cell_type', 'CellType', 'Cell_Type')
    reg = pick('unique_region', 'region')
    don = pick('donor', 'Donor')
    org = pick('Tissue_Segment', 'tissue_segment', 'Tissue Segment', 'tissue')
    parts = []
    for chunk in pd.read_csv(csv, usecols=[cx, cy, ct, reg, don, org], chunksize=500000):
        parts.append(pd.DataFrame(dict(section=chunk[reg].astype(str), donor=chunk[don].astype(str),
                                       organ=chunk[org].astype(str), x_px=chunk[cx], y_px=chunk[cy], label=chunk[ct].astype(str))))
    d = pd.concat(parts)
    # CODEX pixel size: 377.44 nm per pixel is assumed (typical 20x CODEX); VERIFY against the Dryad README before analysis
    d['x_um'] = d.x_px * 0.37744; d['y_um'] = d.y_px * 0.37744
    d.drop(columns=['x_px', 'y_px']).to_parquet(os.path.join(RAW, '..', 'hickey_cells.parquet'))
    print(d.groupby(['organ']).agg(sections=('section', 'nunique'), donors=('donor', 'nunique'), cells=('label', 'size')))


if __name__ == '__main__':
    {'check': check, 'hickey': hickey}[sys.argv[1]]()
