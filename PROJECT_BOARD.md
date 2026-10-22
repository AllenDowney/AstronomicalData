# AstronomicalData Project Board

Numbered tasks for tracking work. Each task has a permanent number; add new tasks at the end. Update status as work progresses.

**Status values:** `Pending`, `In progress`, `DONE`

---

## Task 1: Refresh

**Status:** In progress (2026-10-21)

**Context:** AstronomicalData is a **self-study** workshop (GreenTeaPress Task 24): not offered live, but the notebooks are maintained and run on Colab. The guide is the "Refresh checklist" in `~/GreenTeaPress/plans/workshops.md` (self-study variant). This material became the Data Carpentry lesson *Foundations of Astronomical Data Science* ([site](https://datacarpentry.github.io/astronomy-python/), [source](https://github.com/datacarpentry/astronomy-python)). **The Carpentries retired the lesson on 2026-10-01** (`bc3994d`; repository archived, pages stay online, "no longer actively maintained"), so these notebooks are now the maintained version (Allen, 2026-10-22); the README says so and links the lesson pages. Where the lesson already fixed something that broke here, follow its fix and note differences (below). 113 stars; last content change 2021-03-26. No known recordings.

**0. Set up**
- [x] `git fetch`/pull (2 commits behind: README venues edit; a stray 2022 `Solar_Flare_Case_Study.ipynb` upload). Local state sorted with Allen (Decisions); discarded files backed up to the session scratchpad
- [x] `PROJECT_BOARD.md` and `CLAUDE.md`. Pages builds from `gh-pages` (Jupyter Book 1, last published 2021-08-30), so `.md` files on `main` are not published
- [x] Tags (pushed 2026-10-21): `carpentries-2020` at `094e93c` (2020-11-16 10:18, last commit before the Carpentries test drive, Nov 16–20, 2020; no recording); `pre-refresh-2026` at `1641df8`
- [x] Inventory (below)
- [x] Environment: mamba env `AstronomicalData` (2026-10-21): Python 3.13.16, NumPy 2.5.4, pandas 3.0.6, SciPy 1.18.1, Matplotlib 3.11.2, Astropy 8.0.2, astroquery 0.4.11, gala 1.12.0, PyTables 3.11.1

**1. Runs in a current environment**
- [x] `environment.yml` rewritten (conda-forge, `python=3.13`, pytest, nbmake, Jupyter Book 1; dropped `defaults`/`astropy` channels, `h5py`, `gsl`, `python-wget`, which no notebook imports, only `test_setup.ipynb`), `requirements.txt` (unpinned), `requirements-dev.txt`. gala 1.12.0 has a cp313 manylinux wheel, so `pip install gala` works on Colab
- [x] First code cell: `try: import …` / `except ImportError: %pip install -q …` for astroquery and/or gala, as each notebook imports them (05 was missing gala). Colab has PyTables. Notebooks download data and helper files if missing (05 now downloads `gd1_data.hdf` too)
- [x] Breakage and deprecation warnings: all seven solutions run with Deprecation/Future/PendingDeprecation/User/Runtime/Syntax warnings as errors (2026-10-21). Fixes: `to_hdf(..., key=...)`; GD-1 columns added with `.value` (03, `make_dataframe` in 04/05; the lesson's wording in 03); centerline bounds `phi2_min = -1.0` without units, since the columns are now plain numbers (the lesson's episode 4 still has `* u.degree`, which fails; its code isn't executed); `results.to_pandas()[...]` instead of `pd.Series(masked_column)` (05, the 'partition' UserWarning); raw strings for mathtext labels (07); `seaborn-v0_8-*` style names in the 07 text. ResourceWarnings for unclosed sockets come from inside astroquery. Gaia archive changes: the row cap (02) and the proper-motion polygon (04/05) (Decisions)
- [x] gala from pip (`gala` 1.12.0, cp313 wheel); `astro-gala` and its TODO gone
- [x] `np.set_printoptions(legacy='1.25')` where NumPy is first imported (03, 04, 06, 07; 01, 02 and 05 print no NumPy scalars); solutions re-run live (2026-10-21), outputs refreshed, execution timestamps stripped; student versions regenerated (`cp soln/0*.ipynb .` + `remove_soln.py`)
- [x] Data URLs: `data/` on GitHub, `read_mist_models.py` from `jieunchoi/MIST_codes`, `join.png` from the lesson's `gh-pages`: all work (2026-10-21)
- [ ] Replace Anaconda/Binder setup instructions (README, `index.md`, `jb/README.md`) with Colab first, then pip or mamba
- [x] `test_setup.ipynb`: removed `from wget import download` (python-wget is no longer installed); the rest is content issue 6

**2. Tests and CI**
- [x] Replay saved Gaia results in CI (Allen, 2026-10-21; chosen over load-if-exists once the 1-minute time-out made every query unreliable). `tests/gaia_replay.py` patches `GaiaClass.launch_job`, `launch_job_async`, `load_tables`, `load_table`; `tests/ipython/profile_default/startup/` installs it when `IPYTHONDIR=tests/ipython`. `GAIA_REPLAY=record` makes the real calls and saves the results in `tests/gaia_cache` (34 files, 8.2 MB, mostly the 50,000-row capped query; tables as ECSV, printed text and column names as JSON); `replay` (the default) returns them and fails on a call with no saved result. Notebooks unchanged. Recorded 2026-10-21 in the same live run as the saved outputs. Re-record when a query changes: `GAIA_REPLAY=record IPYTHONDIR=$PWD/tests/ipython pytest --nbmake soln/0*.ipynb`
- [x] `pytest --nbmake soln/0*.ipynb`: 7 passed in 19 s with networking disabled (`unshare -rn`), vs. ~6 min live. Student versions don't run top to bottom: their blanks include the cells that save `centerline_df` etc., so later cells and notebooks fail. Not tested
- [x] `.github/workflows/tests.yml`: push to `main`, PR, monthly, manual; Python 3.13; pip; 30-minute timeout; replay, except the monthly run and a manual run with `live` checked, which query the archive (so archive changes show up; expect occasional failures from time-outs)
- [x] Badge in the README (the README itself is rewritten in step 5)

**3. Colab**
- [x] House-wording Colab cell after each title; solutions also link to the version with blanks (except 06, which has no exercises); `remove_soln.py` rewrites the cell in the student versions to link to themselves. Colab links in the README table (2026-10-22)
- [ ] Allen runs one notebook on Colab end to end

**4. Text**
- [x] Typos (codespell on the cell sources, plus a search for doubled words): "unnused", "coordindates", "the figure we've made have been", "remove `TOP 10` run", "capabilities that (almost) a superset", "identity the stars", "Now can transform", "except that have fewer", "in in Lesson", "the data and and analysis"; "a recent paper" → "a 2018 paper" (Allen). Links: http → https for SIMBAD; MIST → mist.science; NumPy indexing, GeeksforGeeks join and Matplotlib `plot` links moved. Not changed: `ivoa.net` (no response, 2026-10-21), `spiff.rit.edu` (no https); bot-blocked (403/405: Stack Overflow, ADS, STScI) and briefly down (503: Gaia GDR2 docs, IERS) links left as they are. Dead links with no replacement: content issues 8, 9
- [x] Notebook layout (2026-10-22): Carpentries front-matter raw cells deleted (Allen); the numbered title alone in the first cell (Allen: keep the numbers), then the Colab cell, then the intro; copyright 2020–2026 and license (code MIT, text [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), Allen: to match the Carpentries lesson) in a markdown cell at the end; empty trailing cells removed. Same in student versions (regenerated). `LICENSE` year 2020–2026
- [x] 01: a "Note about limits" after the first query: the anonymous limits (1 minute, 50,000 rows), "run the cell again" on a time-out, and how to register and `Gaia.login()` (Allen)
- [x] `test_setup.ipynb` rewritten as a check for local installs (Allen): import checks only, a pointer to the README, copyright cell; the 2020 event material (Jupyter intro, Slack, survey, old pandas check) dropped. Tested in CI
- [x] Deleted `index.md` and the root `_config.yml`/`_toc.yml` (Allen)

**5. README in the standard shape**
- [x] README rewritten (2026-10-22): badge; author; self-study line; abstract and credits (kept); the lesson's retirement and that these notebooks are maintained; level (intermediate, Allen); length (about 20 hours, five 4-hour sessions, Allen); outline; notebook table with Colab links (exercises and solutions; 06 solutions only); link to the site; Gaia archive limits; recordings: none, with the 2020 event and the `carpentries-2020` tag; running instructions (Colab, then pip or mamba, then `test_setup.ipynb`); next steps (hub; no related book); license (MIT, CC BY 4.0)

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
- **Data:** a figshare zip (`student_download.zip`) with the data files, so learners can skip queries. **Difference:** here, notebooks download from `data/` on GitHub, and CI replays saved Gaia results.
- **Testing:** none in CI; episodes are Markdown with pasted outputs. **Difference:** we run the notebooks with nbmake.
- **Setup:** Miniforge + `environment.yml` (`python>=3.11`, conda-forge). No Colab notes.

### Decisions (Allen, 2026-10-21)

- Tags `carpentries-2020` and `pre-refresh-2026`, pushed.
- Local state and leftovers as listed above; `data/` unchanged.
- gala: `pip install gala`, not the lesson's vendored files.
- CI: replay saved Gaia results (Allen, 2026-10-21; first load-if-exists, changed after the 1-minute time-out showed that every query is unreliable).
- Stay on Gaia DR2 (as the lesson and the paper do).
- Keep HDF (pandas + PyTables) for saved DataFrames, with the lesson's `key=` fix (Allen, 2026-10-21; Parquet would mean ~7 files per notebook, rewritten text, and divergence from the lesson).
- **Proper-motion polygon** (Allen, 2026-10-21): the archive rejects `CONTAINS(POINT(pmra, pmdec), POLYGON(...))` ("Coordinate 0 of POLYGON is out of range": vertices must be sky coordinates, and `pmra` is negative; astroquery shows only "HTTP 500: null"). The lesson hit this in 2022 and replaced the convex hull with a box chosen by eye (PR #135, `7358a73`). **Difference:** we keep the convex hull, query its bounding rectangle with `BETWEEN`, then select the stars inside the hull with `matplotlib.patches.Polygon.contains_points` (the idiom 06 uses). 04 saves the hull vertices as `pm_df`; `point_series` now holds only `point_list`. Counts: 13,052 in the rectangle, 7,350 in the hull (was 7,345; straight edges instead of great circles); 6,680 from the join query, 3,726 after the hull (was 3,725); 454 winners (unchanged).
- `data/gd1_data.hdf` regenerated in place (2026-10-21) from a full run of 03–06 with these fixes; same filename and keys plus `pm_df`. `centerline_df`/`selected_df` gained `parallax_error` and `radial_velocity` (they come from the saved `gd1_results.fits`, 140,340 rows, 8 columns, Dec 2020). `gd1_isochrone.hdf5` regenerated identical; not replaced.
- **Gaia row cap** (2026-10-21): since the archive release of 2026-10-14, anonymous queries return at most 50,000 rows (registered users: unlimited) and anonymous async jobs time out after 1 minute ([quotas](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas)). Truncation is silent. Only `02_coords` `query5` (140,339 rows) is affected; the others return under 9,000. Fix: 02 still runs the query and explains the cap (`len(results)` is 50,000), then loads the full saved `data/gd1_results.fits`. The Data Carpentry lesson (last commit 2026-10-02) has the same problem; no issue filed (Allen).

### Content issues

(For Allen to decide.)

Resolved 2026-10-22 (Allen): Gaia DR2 kept (DR3 would change column names, the Pan-STARRS cross-match table, and every count; the lesson is DR2 too); `index.md` and the root Jupyter Book config deleted; "a recent paper" → "a 2018 paper"; HDF kept; `test_setup.ipynb` rewritten; the time-out note with registering added to 01.

Open:

1. `jb/README.md` (the site's landing page) has drifted from the README; step 6 should rebuild it from the new README.
2. The cone-search join in 05 (the 490-row exercise) takes about a minute and often hits the anonymous time-out. Could it be made cheaper (e.g. a smaller cone), or is the retry note enough?
3. The [gaia.ac.uk ADQL cookbook](https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook) (01) is gone (404); no obvious replacement.
4. The ESA "archive tips" link about negative parallax (03) now redirects to a general [FAQ page](https://www.cosmos.esa.int/web/gaia/faqs); the anchor is lost.
5. Several links point to old versioned Matplotlib docs (3.1.0, 3.1.1, 3.3.1); they still work.
6. GreenTeaPress `plans/workshops.md` still says the lesson is the maintained version; update it there (and the hub entry) to say it was retired.
