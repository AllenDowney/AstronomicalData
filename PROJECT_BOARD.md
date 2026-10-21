# AstronomicalData Project Board

Numbered tasks for tracking work. Each task has a permanent number; add new tasks at the end. Update status as work progresses.

**Status values:** `Pending`, `In progress`, `DONE`

---

## Task 1: Refresh

**Status:** In progress (2026-10-21)

**Context:** AstronomicalData is a **self-study** workshop (GreenTeaPress Task 24): not offered live, but the notebooks are maintained and run on Colab. The guide is the "Refresh checklist" in `~/GreenTeaPress/plans/workshops.md` (self-study variant). The maintained version of this material is the Data Carpentry lesson ([site](https://datacarpentry.github.io/astronomy-python/), [source](https://github.com/datacarpentry/astronomy-python)); the README points to it. Where the lesson already fixed something that broke here, follow its fix and note differences (below). 113 stars; last content change 2021-03-26. No known recordings.

**0. Set up**
- [x] `git fetch`/pull (2 commits behind: README venues edit; a stray 2022 `Solar_Flare_Case_Study.ipynb` upload). Local state sorted with Allen (Decisions); discarded files backed up to the session scratchpad
- [x] `PROJECT_BOARD.md` and `CLAUDE.md`. Pages builds from `gh-pages` (Jupyter Book 1, last published 2021-08-30), so `.md` files on `main` are not published
- [x] Tags (pushed 2026-10-21): `carpentries-2020` at `094e93c` (2020-11-16 10:18, last commit before the Carpentries test drive, Nov 16–20, 2020; no recording); `pre-refresh-2026` at `1641df8`
- [x] Inventory (below)
- [x] Environment: mamba env `AstronomicalData` (2026-10-21): Python 3.13.16, NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1, Matplotlib 3.11.2, Astropy 8.0.2, astroquery 0.4.11, gala 1.12.0, PyTables 3.11.1

**1. Runs in a current environment**
- [x] `environment.yml` rewritten (conda-forge, `python=3.13`, pytest, nbmake, Jupyter Book 1; dropped `defaults`/`astropy` channels, `h5py`, `gsl`, `python-wget`, which no notebook imports, only `test_setup.ipynb`), `requirements.txt` (unpinned), `requirements-dev.txt`. gala 1.12.0 has a cp313 manylinux wheel, so `pip install gala` works on Colab
- [ ] First code cell: `%pip install` what Colab lacks (astroquery, gala; check pytables); download data and helper files if missing
- [ ] Fix breakage and deprecation warnings (warnings-as-errors scratch runs). First run (2026-10-21, warnings shown, live Gaia): every error in 03–07 traces to three causes: positional `to_hdf` key (TypeError in pandas 3; 03, 05, 06); GD-1 Quantity columns added to the Table before `to_pandas()` (pandas 3 fails in `describe`, `between`, `to_csv`; 03, 04, 05); and a cascade, because 03 left an empty `gd1_data.hdf` that made later `download()` calls skip the real file. Plus SyntaxWarnings for `'$\phi_1$'` labels (07), the 'partition' UserWarning in 05, and ResourceWarnings for unclosed sockets (inside astroquery; not ours). Known from the lesson: `to_hdf(..., key=...)`; raw strings for `$\phi_1$` labels (SyntaxWarning on 3.12+); add GD-1 columns to the DataFrame after `to_pandas()`, not to the Table; `seaborn-*` style names (removed in Matplotlib 3.8; the lesson hasn't fixed this)
- [ ] gala: keep it and `pip install gala` (Allen, 2026-10-21; the lesson vendored `gd1.py`/`reflex.py` instead). Check it installs on 3.13 and on Colab; drop the `astro-gala` TODO
- [ ] `np.set_printoptions(legacy='1.25')` in the cell that imports NumPy; re-run the solution notebooks
- [ ] Data URLs: `data/` on GitHub (works); MIST `read_mist_models.py` from `jieunchoi/MIST_codes` (check); `join.png` from the lesson's `gh-pages` (works 2026-10-21)
- [ ] Replace Anaconda/Binder setup instructions (README, `index.md`, `jb/README.md`) with Colab first, then pip or mamba

**2. Tests and CI**
- [ ] Load-if-exists in the query notebooks (Allen, 2026-10-21): each big query cell loads its saved file if present, otherwise queries Gaia, so CI doesn't depend on the Gaia service. Notebooks 03, 06, 07 already read saved files
- [ ] `pytest --nbmake` on the solutions; student versions too if their blanks don't break later cells
- [ ] `tests.yml` (push, PR, monthly, manual; Python 3.13; pip; timeout). Ask Allen: should the monthly run delete the saved files and query Gaia live?
- [ ] Badge in the README

**3. Colab**
- [ ] House-wording Colab cell after each title (solutions link to the student version); Colab links in the README table
- [ ] Allen runs one notebook on Colab end to end

**4. Text**
- [ ] Typos and broken links only; everything else under Content issues
- [ ] Notebook layout: one `#` title, Colab cell, copyright and license (code MIT, text CC) at the end; same in student and solution versions

**5. README in the standard shape**
- [ ] Self-study line at the top; the Data Carpentry lesson as the maintained version; recordings: none (event and tag listed); length (ask Allen); license

**6. Site and visibility**
- [ ] Rebuild the existing `gh-pages` site from `jb/` (delete `jb/_build/` first) and do the SEO checklist; Allen reviews descriptions
- [ ] Hub entry, repo description and topics

**7. Close**
- [ ] Tag `2026-refresh`; note in GreenTeaPress Task 24

### Inventory (2026-10-21)

**The workshop:** 7 notebooks. The sources are `soln/0*.ipynb`. The student versions at the root are generated by `build.sh` → `remove_soln.py`, which replaces `# Solution` cells with `# Solution goes here`, strips outputs, and drops the first raw cell. `soln/build.sh` commits the solutions and copies their output data into `data/`.

| Notebook | Topic | Gaia queries | Reads saved files |
|---|---|---|---|
| `01_query` | Gaia archive, ADQL, `load_tables`/`load_table`, `launch_job(_async)` | 4 jobs + metadata | — |
| `02_coords` | Units, SkyCoord, GD-1 frame (gala), polygon query; writes `gd1_results.fits` | 4 | — |
| `03_motion` | Transform, reflex correction, pandas, proper-motion selection; writes `gd1_data.hdf` | — | `data/gd1_results.fits` |
| `04_select` | Convex hull, polygon + proper-motion query; writes `point_series` | 1 | `data/gd1_data.hdf` |
| `05_join` | Gaia × Pan-STARRS joins; writes `candidate_df` | 6 | `gd1_data.hdf` |
| `06_photo` | Isochrone (MIST), CMD selection; writes `loop_df`, `winner_df`, `gd1_isochrone.hdf5` | — | `data/gd1_data.hdf`, `data/MIST_iso_*.cmd`, `read_mist_models.py` (jieunchoi/MIST_codes), `data/gd1_isochrone.hdf5` |
| `07_plot` | Publication-quality figure, style sheets | — | `data/gd1_data.hdf`, `az-paper-twocol.mplstyle` |

Every notebook has an `IN_COLAB` `!pip install` cell (astroquery, `astro-gala`) and a `download()` helper. Still on Gaia DR2. All three DR2 tables answer queries (curl, 2026-10-21).

**Also on `main`:** `test_setup.ipynb` (pre-workshop setup check; the 2020 event linked it), `index.md` (old copy of the README with install instructions), root `_config.yml`/`_toc.yml` (old Jupyter Book config; `jb/` is the current one), `environment.yml` (2020, `defaults` + `astropy` channels, `python>=3.6`), `images/` (HR diagram), `run_on_colab_small.png`, `az-paper-twocol.mplstyle`, `LICENSE`.

**`data/`:** the files the notebooks download (`gd1_results.fits`, `gd1_data.hdf`, `gd1_isochrone.hdf5`, `MIST_iso_5fd2532653c27.iso.cmd`) plus 2020 files no notebook reads now (`gd1_dataframe.hdf5` 20 MB, `gd1_candidates.hdf5`, `gd1_merged.hdf5`, `gd1_photo.fits`, `gd1_polygon`, `gd1_polygon.hdf5`). Kept unchanged: the lesson's instructor page downloads the MIST file from here, and old links may use the others.

**Site:** `jb/` (now tracked): `_config.yml`, `_toc.yml`, `README.md` (landing page), `prep_notebooks.py` (adds `hide-cell` to solutions), `build.sh` (copies `soln/0*.ipynb`, `jb build`, `ghp-import`). Execution is off.

**Removed from `main` 2026-10-21 (in `pre-refresh-2026` and history):** `Solar_Flare_Case_Study.ipynb` (a *Probably Overthinking It* draft uploaded here by mistake), `last_resort.ipynb` (Colab fallback; superseded), `class_notebooks/lesson_7.ipynb` (A. Bostroem's plotting functions from the 2020 event), root duplicate of `MIST_iso_5fd2532653c27.iso.cmd`.

**Discarded local state (backed up to the session scratchpad):** scratch data copies at the root and in `soln/` (2020–2021), `_build/`, `jb/_build/` and copied notebooks, `jb.zip`, `day1.ipynb`, a second `lesson_7.ipynb`, `test.*`, `introduction.md` (draft intro, 2020-11), `process.py`, `add_tags.py`, `notes.txt~`, a kernel-metadata edit to `class_notebooks/lesson_7.ipynb`, saved outputs in `test_setup.ipynb`.

**Kept and committed:** uncommitted fixes (`pytables` in the `index.md` install line and a pytables check in `test_setup.ipynb`; `soln/01_query.ipynb` ran `query` where `query2` was meant); `jb/`; `build.sh`, `remove_soln.py`, `soln/build.sh`; `figs/join.odg`/`join.png` (source of the join diagram in `05_join`).

### Data Carpentry lesson: what it changed (checked 2026-10-21, lesson at `25b78f6`, 2026-10-02)

- **Episodes:** 8, generated from these notebooks; `03_motion` split in two (2021-11). The MIST isochrone computation moved to an instructor page; learners read a prebuilt `gd1_isochrone.hdf5`. Helpers collected in `episode_functions.py`.
- **Still Gaia DR2**, same queries and counts. Nothing to follow for DR3; staying on DR2 here too.
- **gala removed** (2025-12): vendored `gd1.py` (`GD1Koposov10`) and `reflex.py` (`reflex_correct`). **Difference:** we keep gala (Allen, 2026-10-21).
- **pandas:** `to_hdf(filename, key=...)`; columns added to the DataFrame after `to_pandas()` (2026-05).
- **Matplotlib:** raw strings for mathtext labels. Both of us still use `seaborn-*` style names (gone in Matplotlib 3.8).
- **Data:** a figshare zip (`student_download.zip`) with the data files, so learners can skip queries. **Difference:** here, notebooks download from `data/` on GitHub, and the query notebooks load-if-exists.
- **Testing:** none in CI; episodes are Markdown with pasted outputs. **Difference:** we run the notebooks with nbmake.
- **Setup:** Miniforge + `environment.yml` (`python>=3.11`, conda-forge). No Colab notes.

### Decisions (Allen, 2026-10-21)

- Tags `carpentries-2020` and `pre-refresh-2026`, pushed.
- Local state and leftovers as listed above; `data/` unchanged.
- gala: `pip install gala`, not the lesson's vendored files.
- CI: load-if-exists in the query notebooks.
- Stay on Gaia DR2 (as the lesson and the paper do).
- Keep HDF (pandas + PyTables) for saved DataFrames, with the lesson's `key=` fix (Allen, 2026-10-21; Parquet would mean ~7 files per notebook, rewritten text, and divergence from the lesson).
- **Proper-motion polygon** (Allen, 2026-10-21): the archive rejects `CONTAINS(POINT(pmra, pmdec), POLYGON(...))` ("Coordinate 0 of POLYGON is out of range": vertices must be sky coordinates, and `pmra` is negative; astroquery shows only "HTTP 500: null"). The lesson hit this in 2022 and replaced the convex hull with a box chosen by eye (PR #135, `7358a73`). **Difference:** we keep the convex hull, query its bounding rectangle with `BETWEEN`, then select the stars inside the hull with `matplotlib.patches.Polygon.contains_points` (the idiom 06 uses). 04 saves the hull vertices as `pm_df`; `point_series` now holds only `point_list`. Counts: 13,052 in the rectangle, 7,350 in the hull (was 7,345; straight edges instead of great circles); 6,680 from the join query, 3,726 after the hull (was 3,725); 454 winners (unchanged).
- `data/gd1_data.hdf` regenerated in place (2026-10-21) from a full run of 03–06 with these fixes; same filename and keys plus `pm_df`. `centerline_df`/`selected_df` gained `parallax_error` and `radial_velocity` (they come from the saved `gd1_results.fits`, 140,340 rows, 8 columns, Dec 2020). `gd1_isochrone.hdf5` regenerated identical; not replaced.
- **Gaia row cap** (2026-10-21): since the archive release of 2026-10-14, anonymous queries return at most 50,000 rows (registered users: unlimited) and anonymous async jobs time out after 1 minute ([quotas](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas)). Truncation is silent. Only `02_coords` `query5` (140,339 rows) is affected; the others return under 9,000. Fix: 02 still runs the query and explains the cap (`len(results)` is 50,000), then loads the full saved `data/gd1_results.fits`. The Data Carpentry lesson (last commit 2026-10-02) has the same problem; no issue filed (Allen).

### Content issues

(For Allen to decide.)

1. Gaia DR2 vs. DR3: the lesson and these notebooks are DR2. DR3 would change column names, the Pan-STARRS cross-match table, and every count.
2. `index.md` duplicates the README with old install instructions; the root `_config.yml`/`_toc.yml` duplicate `jb/`. Delete them once the README and `jb/` are updated?
3. The paper is described as "recent" (2018).
4. `jb/README.md` and the README have drifted apart; decide which is the landing page text.
5. Gaia **1-minute anonymous time-out**: the cone-search join in 05 (the 490-row exercise) takes about a minute and failed 3 of 4 tries on 2026-10-21 (HTTP 408/500). Other queries finished in time, but all are at the mercy of server load. Options: make the slow queries cheaper, tell learners to retry, or point to the saved results.
