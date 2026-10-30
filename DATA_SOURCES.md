# Data sources

| Source | Access | Licence / terms | Derived files in `analysis/data/` |
|---|---|---|---|
| GSE325587 (CosMx, diabetic kidney disease; Meadows et al. 2026) | GEO GSE325587 | GEO terms | `fov_celltype_counts.tsv`, `sample_mapping.tsv`, `harmonized_category_mapping.tsv` |
| GSE211785 (CosMx; Abedini et al. 2024) | GEO GSE211785 | GEO terms | `GSE211785_fov_celltype_counts.tsv` |
| Pediatric lupus nephritis (CosMx; Danaher et al. 2024) | Figshare 10.6084/m9.figshare.26409307.v1 | CC BY 4.0 | `Danaher_fov_label_counts.tsv`, `Danaher_fov_mapping.tsv` |
| Mouse kidney Xenium ShamR, Day14R (Xuanyuan et al. 2025) | Figshare 10.6084/m9.figshare.26049403.v2; 26053363.v1 | CC BY 4.0 | `mouse_example_*.tsv.gz` |
| Human intestine CODEX (Hickey et al. 2023) | Dryad 10.5061/dryad.pk0p2ngrf | Public domain | `hickey_extract/part00-06.tsv.gz` (seven complete sites; transverse colon truncated in our download) |
| GSE286184 (Xenium) | GEO GSE286184 | GEO terms | summaries only |

The Hickey extract was produced with `extract_range.py` from a browser download of
`23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv` that stopped at 2,793,116,391 bytes
(SHA256 b01b5ad468da9f387766e148fc57c0f6b696d52c83d33a983a59b5dec9f2aa6b). Records are stored in contiguous blocks
by intestinal site; the first seven sites are complete. To include the transverse colon, download the complete file,
rerun `extract_range.py` over its full byte range and remove the site filter in `c0_hickey_assemble.py`.
Pixel size is assumed to be 0.37744 um per pixel; verify before use.
Original-analysis provenance (download hashes, audits) is documented in the G18 research package.
