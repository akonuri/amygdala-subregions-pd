# Zenodo deposit guide

Two records, both in the ASAP community (https://zenodo.org/communities/asaphub). Do the code record first so the data record can link to it.

## 1. Code record (automatic from GitHub)
1. Log in to Zenodo with GitHub (akonuri). Settings > GitHub > switch on `akonuri/amygdala-subregions-pd`.
2. Fill the placeholders that remain in the repo before releasing.
   - `KeyResourcesTable.csv` FreeSurfer, R and Python versions
   - `README.md` upstream steps (FreeSurfer, longCombat, subtype classification scripts)
   - `.zenodo.json` add `"orcid": "0000-..."` to each creator who has one
   - `requirements.txt` pin exact versions (`pip freeze | grep -iE "pandas|numpy|scipy|statsmodels|lifelines|scikit-learn|scikit-posthocs|matplotlib|openpyxl"`)
3. On GitHub, Releases > Draft a new release, tag `v1.0.0`, title `v1.0.0 (manuscript submission)`, publish.
4. Zenodo mints the DOI within minutes (Zenodo > Uploads). Open the record, check metadata and that it sits in the ASAP community. Note two DOIs, the version DOI and the concept DOI (all versions). Cite the version DOI in the paper.

## 2. Data record (manual upload)
1. Run the full analysis, then `bash zenodo/make_data_record.sh "$AMYG_PD_OUTPUT_DIR"`.
2. In the output folder, fill blank descriptions in `data_dictionary.csv` and the placeholders in `README.md` (code DOI from step 1).
3. Check that no file holds participant-level values. The script blocks identifier columns, but look over `figure_source_data/` and small counts.
4. Zenodo > New upload. Upload the folder contents (or one zip). Fill the form.

| Field | Value |
|---|---|
| Resource type | Dataset |
| Title | Source data for: Amygdala subregion atrophy in Parkinson's disease is associated with visual hallucinations and cognitive impairment |
| Creators | Same order and affiliations as the code record, with ORCIDs |
| Description | First paragraph of the record README |
| Licence | Creative Commons Attribution 4.0 International |
| Keywords | Parkinson's disease; amygdala; MRI; visual hallucinations; cognitive impairment; PPMI; Aligning Science Across Parkinson's |
| Version | 1.0.0 |
| Communities | Aligning Science Across Parkinson's (asaphub) |
| Funding | Award: ASAP-027449, funder Aligning Science Across Parkinson's (via the Michael J. Fox Foundation) |
| Related works | `is supplement to` code DOI · `is supplement to` preprint DOI (add when posted) |

5. Publish. Note the data DOI.

## 3. Back-fill
1. Put both DOIs in `KeyResourcesTable.csv`, `README.md` and the manuscript Data and Code Availability statements.
2. Edit the code record on Zenodo to add `is supplemented by` the data DOI.
3. Commit, then release `v1.0.1` only if code changed. Metadata-only edits can be made on Zenodo without a new release.
4. Upload `KeyResourcesTable.csv` with the manuscript and to the ASAP compliance review.
