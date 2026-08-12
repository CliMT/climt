# climt Website (Quarto Migration + RT Walkthrough) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace climt's Sphinx/ReadTheDocs docs with a single Quarto GitHub Pages site that unifies the existing user/API docs, an 8-chapter radiative-transfer walkthrough, three pedagogical notebooks, the picket-fence Experiments post, and an auto-generated API reference — deployed via GitHub Actions.

**Architecture:** The Quarto project roots at `docs/`. It builds on the minimal scaffolding sub-project B already shipped (`docs/_quarto.yml`, `references.bib`, `scripts/build_experiments.py`, the `make experiments` target, `docs/experiments/`). Existing `.rst` is pandoc-converted to `.qmd`; autodoc is replaced by `quartodoc`; every walkthrough figure is a regenerated artifact wired through B's `sources.yml` + `build_experiments.py` pipeline (not a one-off script); code excerpts are pulled live from climt source so they cannot drift; citations resolve against the single root `references.bib`. The 8 RT chapters' prose is adapted from the committed phase-4 plan (`docs/superpowers/plans/2026-04-20-picket-fence-radiation-phase4.md`, Part F, Tasks 17–25) — that document is the prose source of record; this plan specifies the Quarto transformation of each.

**Tech Stack:** Quarto 1.9.38 (installed), `quartodoc` (to be added), Python in the `climt` conda env, `pandoc` (bundled with Quarto), `linepyline` (for some LBL figures; degrades gracefully when absent), GitHub Actions + `peaceiris/actions-gh-pages`.

**Conventions for every task:**
- Run Python via the `climt` env: `conda run --no-capture-output -n climt python …` (works under `sh -c`; if your sandbox blocks it, use `/Users/joymonteiro/miniconda3/envs/climt/bin/python …`).
- Render with `quarto`; point it at the env interpreter for notebook/execute cells: `QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render docs/`.
- Commit at the end of each task. Do not touch unrelated working-tree changes.

**Branch dependency:** This work depends on sub-project B's scaffolding, which currently lives on `feature/picket-fence-radiation` (PR #207). Branch this work off `feature/picket-fence-radiation` (or off `develop` only after #207 merges). Do **not** start from a clean `develop` that lacks `docs/_quarto.yml` / `scripts/build_experiments.py`.

---

## File Structure

```
docs/
├── _quarto.yml                      # MODIFY: minimal seed -> full site config (nav, sections, theme, quartodoc)
├── index.qmd                        # CREATE: landing page (from index.rst)
├── get-started/                     # CREATE: install, quickstart, first model (from installation/quickstart .rst)
│   ├── installation.qmd
│   └── quickstart.qmd
├── user-guide/                      # CREATE: components, init, config, etc. (from existing .rst)
│   ├── introduction.qmd ... (one per migrated .rst)
│   └── plug-and-play.qmd            # from PLUG_AND_PLAY_ARCHITECTURE.md
├── radiative-transfer/              # CREATE: the 8-chapter walkthrough
│   ├── index.qmd
│   ├── 01-why-nongrey.qmd ... 08-multiplanet.qmd
│   ├── performance.qmd
│   ├── table-generation.qmd         # from docs/radiative_transfer/table_generation.rst
│   ├── _figures.py                  # CREATE: walkthrough figure functions (ported from phase4 _generate_figures.py)
│   ├── sources.yml                  # CREATE: artifact manifest for the chapter figures
│   └── _artifacts/                  # generated figures (committed)
├── experiments/                     # EXISTS (sub-project B): index.qmd + the picket-fence post
├── api/                             # CREATE: quartodoc-generated reference (_metadata + generated .qmd)
│   └── _metadata.yml
├── notebooks/ (authored under examples/, embedded into chapters)
└── _site/                           # gitignored build output (B added the ignore)

examples/
├── k_distribution_demo.ipynb        # CREATE (phase4 Task 13)
├── spectral_radiation_anatomy.ipynb # CREATE (phase4 Task 14)
└── picket_fence_vs_rrtmg.ipynb      # CREATE (phase4 Task 15)

references.bib                       # MODIFY: append every newly cited work (each with url=)
.github/workflows/docs.yml           # CREATE: build + deploy
scripts/build_experiments.py         # EXISTS (sub-project B) — reused as-is
```

Deletions (Task 5): `docs/conf.py`, `docs/Makefile`, `docs/make.bat`, `docs/modules.rst`, `docs/climt.rst`, `readthedocs.yml`, the 6 `docs/*.pdf`, and every `docs/*.rst` once its content is migrated.

---

## PART A — Quarto project foundation

### Task 1: Add quartodoc and expand the site config

**Files:**
- Modify: `docs/_quarto.yml`
- Modify: `setup.py` (docs extra) or `requirements_dev.txt` (whichever the repo uses for dev deps)

- [ ] **Step 1: Install quartodoc into the env**

Run: `conda run --no-capture-output -n climt pip install quartodoc`
Then verify: `conda run --no-capture-output -n climt python -c "import quartodoc; print(quartodoc.__version__)"`
Expected: a version prints.

- [ ] **Step 2: Record the dev dependency**

Add `quartodoc` to the project's dev/docs dependency list. Find it first: `grep -rn "sphinx" setup.py requirements*.txt docs/requirements*.txt 2>/dev/null`. Replace the Sphinx docs dependencies (`sphinx`, `sphinx_rtd_theme`, etc.) with `quartodoc` in that same list. If a `docs` extras_require group exists in `setup.py`, edit it there.

- [ ] **Step 3: Expand `docs/_quarto.yml` to the full site config**

Replace the minimal B seed with (preserving the `project`/`format` blocks B established):

```yaml
project:
  type: website
  output-dir: _site
  render:
    - "*.qmd"
    - "experiments/**/*.qmd"
    - "radiative-transfer/**/*.qmd"
    - "!_artifacts/"

website:
  title: "climt"
  description: "A climate modelling toolkit"
  repo-url: https://github.com/CliMT/climt
  navbar:
    left:
      - text: "Get Started"
        href: get-started/installation.qmd
      - text: "User Guide"
        href: user-guide/introduction.qmd
      - text: "Radiative Transfer"
        href: radiative-transfer/index.qmd
      - text: "Experiments"
        href: experiments/index.qmd
      - text: "API"
        href: api/index.qmd
    right:
      - icon: github
        href: https://github.com/CliMT/climt
  sidebar:
    - title: "User Guide"
      contents: user-guide/*
    - title: "Radiative Transfer"
      contents:
        - radiative-transfer/index.qmd
        - radiative-transfer/01-why-nongrey.qmd
        - radiative-transfer/02-line-by-line.qmd
        - radiative-transfer/03-k-distribution.qmd
        - radiative-transfer/04-correlated-k.qmd
        - radiative-transfer/05-gas-overlap.qmd
        - radiative-transfer/06-picket-fence.qmd
        - radiative-transfer/07-two-stream.qmd
        - radiative-transfer/08-multiplanet.qmd
        - radiative-transfer/table-generation.qmd
        - radiative-transfer/performance.qmd

format:
  html:
    theme: cosmo
    toc: true
    code-fold: false
    code-tools: true

bibliography: references.bib

metadata-files:
  - api/_metadata.yml
```

- [ ] **Step 4: Verify the config parses**

Run: `conda run --no-capture-output -n climt python -c "import yaml; yaml.safe_load(open('docs/_quarto.yml')); print('ok')"`
Expected: `ok`. (Full `quarto render` is deferred until content exists; it will fail on missing hrefs now, which is expected.)

- [ ] **Step 5: Commit**

```bash
git add docs/_quarto.yml setup.py requirements_dev.txt 2>/dev/null
git commit -m "build(docs): add quartodoc dep + full Quarto site config"
```

---

## PART B — Migrate existing content (RST → QMD)

### Task 2: Pandoc-convert the user-guide and get-started pages

**Files:**
- Create: `docs/get-started/installation.qmd`, `docs/get-started/quickstart.qmd`
- Create: `docs/user-guide/*.qmd` (introduction, interaction, realistic, component_types, configuration, components, component_manual, initialisation, utilities, naming, memory_management, emanuel_convection, rrtmg_clouds, contributing, authors, history)
- Create: `docs/user-guide/plug-and-play.qmd` (from `docs/PLUG_AND_PLAY_ARCHITECTURE.md` if present)

- [ ] **Step 1: Batch-convert with pandoc (Quarto bundles it)**

For each source `.rst`, convert to `.qmd`. Example for one file (repeat per file, routing to the right subdir):

```bash
cd /Users/joymonteiro/github/climt/docs
quarto pandoc installation.rst -f rst -t markdown -o get-started/installation.qmd
quarto pandoc quickstart.rst    -f rst -t markdown -o get-started/quickstart.qmd
for f in introduction interaction realistic component_types configuration components \
         component_manual initialisation utilities naming memory_management \
         emanuel_convection rrtmg_clouds contributing authors history; do
  quarto pandoc "$f.rst" -f rst -t markdown -o "user-guide/$f.qmd"
done
```

- [ ] **Step 2: Add front-matter + tidy directives**

For each converted `.qmd`, prepend a YAML title block and fix directives pandoc didn't translate:
- Add `---\ntitle: "<Human Title>"\n---` at the top (derive the title from the original `.rst` heading).
- Replace any `.. autoclass::` / `.. automodule::` blocks with a note + cross-link to the API page (these become quartodoc, Task 6) — do NOT hand-write API docs.
- Replace `:ref:` / `:doc:` cross-references with Quarto relative links (`[text](other-page.qmd)`).
- Replace `.. code-block:: python` with fenced ```` ```python ```` blocks (pandoc usually does this; verify).

- [ ] **Step 3: Verify each page renders in isolation**

Run (one example): `QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render docs/get-started/installation.qmd`
Expected: HTML produced, no fatal directive errors. Fix any leftover `.. raw::` or unknown directives.

- [ ] **Step 4: Commit**

```bash
git add docs/get-started docs/user-guide
git commit -m "docs: migrate user-guide and get-started pages rst->qmd"
```

### Task 3: Migrate the landing page and the two existing notebooks

**Files:**
- Create: `docs/index.qmd` (from `docs/index.rst`, minus the Sphinx toctree — nav now lives in `_quarto.yml`)
- Create: `docs/user-guide/second-best.qmd`, `docs/user-guide/implicit-scheme.qmd` (embedding `docs/Description_of_SecondBEST.ipynb` and `docs/Implicit_Scheme_Derivation.ipynb`)

- [ ] **Step 1: Landing page**

`quarto pandoc docs/index.rst -f rst -t markdown -o docs/index.qmd`, then delete the `.. toctree::` block and add a `---\ntitle: "climt"\n---` front-matter with a short intro and links to the four sections.

- [ ] **Step 2: Embed the two docs notebooks as pages**

Create `docs/user-guide/second-best.qmd`:

```markdown
---
title: "Description of SecondBEST"
---

{{< embed ../Description_of_SecondBEST.ipynb echo=true >}}
```

Create `docs/user-guide/implicit-scheme.qmd` analogously embedding `../Implicit_Scheme_Derivation.ipynb`. (The notebooks must execute cleanly first — if either errors on render, run `conda run -n climt jupyter nbconvert --to notebook --execute --inplace docs/<nb>.ipynb` once and commit the executed version.)

- [ ] **Step 3: Add these two pages to the User Guide sidebar** in `docs/_quarto.yml` (`contents: user-guide/*` already globs them — verify they appear).

- [ ] **Step 4: Commit**

```bash
git add docs/index.qmd docs/user-guide/second-best.qmd docs/user-guide/implicit-scheme.qmd docs/Description_of_SecondBEST.ipynb docs/Implicit_Scheme_Derivation.ipynb
git commit -m "docs: migrate landing page + embed existing notebooks as pages"
```

### Task 4: Migrate citations and delete PDFs

**Files:**
- Modify: `references.bib`
- Delete: `docs/aa22342-13.pdf`, `docs/J Adv Model Earth Syst - 2019 - Pincus …​.pdf`, `docs/Louis.79.pdf`, `docs/stab1851.pdf`, `docs/aa23127-13.pdf`, `docs/Emanuel and Živković-Rothman - 1999 …​.pdf`

- [ ] **Step 1: Identify what cites the PDFs**

Run: `grep -rniE "pincus|louis|emanuel|zivkovic|209458|stab1851|aa22342|aa23127" docs/*.rst docs/user-guide/*.qmd 2>/dev/null`
For each PDF that is actually referenced in prose, add a BibTeX entry to `references.bib` with a `url=` to the publisher/DOI/arXiv page (no PDF in repo). Known mappings:
- `Pincus 2019` → J. Adv. Model. Earth Syst., doi 10.1029/2019MS001621.
- `Louis 1979` → Boundary-Layer Meteorol., doi 10.1007/BF00117978.
- `Emanuel & Živković-Rothman 1999` → J. Atmos. Sci., doi 10.1175/1520-0469(1999)056<1766:DAEOAC>2.0.CO;2.
- `aa22342-13` / `aa23127-13` → the two A&A picket-fence papers (Parmentier & Guillot 2014 is already in `references.bib` as `parmentier2014`; add the companion `parmentier2015` A&A 574, A35, doi 10.1051/0004-6361/201323147 if cited).
- `stab1851` → Lee et al. (2021) MNRAS HD 209458b paper (doi 10.1093/mnras/stab1851).

Add each as `@article{...}` with `url=` mirroring the existing entries' style.

- [ ] **Step 2: Replace any in-text PDF references with `@key` citations** in the migrated `.qmd` pages.

- [ ] **Step 3: Delete the PDFs**

```bash
git rm "docs/aa22342-13.pdf" "docs/Louis.79.pdf" "docs/stab1851.pdf" "docs/aa23127-13.pdf" \
       "docs/J Adv Model Earth Syst - 2019 - Pincus - Balancing Accuracy Efficiency and Flexibility in Radiation Calculations for.pdf" \
       "docs/Emanuel and Živković-Rothman - 1999 - Development and Evaluation of a Convection Scheme for Use in Climate Models.pdf"
```

- [ ] **Step 4: Commit**

```bash
git add references.bib docs
git commit -m "docs: migrate PDF citations to references.bib (external links), delete PDFs"
```

### Task 5: Remove Sphinx machinery

**Files:**
- Delete: `docs/conf.py`, `docs/Makefile`, `docs/make.bat`, `docs/modules.rst`, `docs/climt.rst`, `readthedocs.yml`, and every remaining migrated `docs/*.rst`

- [ ] **Step 1: Confirm every `.rst` has a `.qmd` counterpart**

Run: `for f in docs/*.rst; do base=$(basename "$f" .rst); find docs -name "$base.qmd" -o -name "${base//_/-}.qmd" | grep -q . && echo "OK $base" || echo "MISSING $base"; done`
Resolve any `MISSING` before deleting (some, like `modules.rst`/`climt.rst`, are autodoc and are intentionally replaced by quartodoc — mark those handled).

- [ ] **Step 2: Delete Sphinx files and the migrated `.rst`**

```bash
cd /Users/joymonteiro/github/climt
git rm docs/conf.py docs/Makefile docs/make.bat docs/modules.rst docs/climt.rst readthedocs.yml
git rm docs/*.rst
git rm -r docs/radiative_transfer 2>/dev/null  # old Sphinx figures dir; superseded by docs/radiative-transfer
```

(Keep `docs/radiative_transfer/table_generation.rst` content — migrate it to `docs/radiative-transfer/table-generation.qmd` in Task 14 before deleting; if not yet migrated, hold this specific delete until Task 14.)

- [ ] **Step 3: Verify the project still renders what exists**

Run: `QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render docs/ 2>&1 | tail -20`
Expected: renders the migrated pages + experiments; broken links only to not-yet-created radiative-transfer/api pages (acceptable at this stage — note them).

- [ ] **Step 4: Commit**

```bash
git add -A docs readthedocs.yml
git commit -m "docs: remove Sphinx machinery (conf.py, Makefiles, autodoc rst, readthedocs.yml)"
```

---

## PART C — API reference (quartodoc)

### Task 6: Generate the API reference with quartodoc

**Files:**
- Create: `docs/api/_metadata.yml`
- Create: `docs/api/index.qmd` (overview; quartodoc fills the rest)

- [ ] **Step 1: Write the quartodoc config**

Create `docs/api/_metadata.yml`:

```yaml
quartodoc:
  package: climt
  dir: api
  title: "API Reference"
  style: pkgdown
  sections:
    - title: "Radiation"
      desc: "Longwave and shortwave radiation components."
      contents:
        - RRTMGLongwave
        - RRTMGShortwave
        - GrayLongwaveRadiation
        - Frierson06LongwaveOpticalDepth
    - title: "Picket-fence radiation"
      desc: "Non-grey correlated-k / Parmentier scheme."
      contents:
        - name: PicketFenceLongwave
          package: climt._components.picket_fence
        - name: PicketFenceShortwave
          package: climt._components.picket_fence
    - title: "Convection & surface"
      contents:
        - EmanuelConvection
        - DryConvectiveAdjustment
        - SlabSurface
        - IceSheet
    - title: "Dynamics & utilities"
      contents:
        - get_default_state
        - get_grid
```

(Confirm each symbol resolves: `conda run -n climt python -c "import climt; [getattr(climt, n) for n in ['RRTMGLongwave','EmanuelConvection','SlabSurface','get_default_state']]"`. Drop or rename any that error; add others present in `climt.__all__`.)

- [ ] **Step 2: Write the API landing page** `docs/api/index.qmd`:

```markdown
---
title: "API Reference"
---

Auto-generated reference for climt's public components and helpers. Use the
sidebar to browse by category, or jump to a specific component below.
```

- [ ] **Step 3: Build the reference**

Run: `cd docs && conda run --no-capture-output -n climt quartodoc build`
Expected: generates `docs/api/*.qmd` stubs for each listed symbol. Then `quartodoc interlinks` if the version supports it.

- [ ] **Step 4: Render to confirm the API section appears**

Run: `QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render docs/ 2>&1 | tail -10`
Expected: API pages render; navbar "API" link resolves.

- [ ] **Step 5: Commit**

```bash
git add docs/api
git commit -m "docs(api): quartodoc-generated API reference replacing autodoc"
```

---

## PART D — Walkthrough figure pipeline (reuse B's build_experiments)

### Task 7: Port the figure functions and wire them through build_experiments

**Files:**
- Create: `docs/radiative-transfer/_figures.py`
- Create: `docs/radiative-transfer/sources.yml`
- Create: `docs/radiative-transfer/_artifacts/` (regenerated, committed)

- [ ] **Step 1: Port the figure functions**

Create `docs/radiative-transfer/_figures.py` containing the figure functions from phase4 `_generate_figures.py` (lines 1593–1715 of `docs/superpowers/plans/2026-04-20-picket-fence-radiation-phase4.md`): `fig_01_mean_of_exp`, `fig_02_lbl_spectrum`, `fig_03_k_distribution_construction`, `fig_04_correlation_across_T`, `fig_06_picket_fence_opacity`, `fig_07_two_stream_phases`. Adapt the module so each function takes an explicit `--figure NAME --out PATH` via an `argparse` CLI (mirroring `scripts/experiments/make_subproject_B_figures.py`), writing one PNG per call, so it can be driven by `sources.yml`. Keep the `linepyline`-dependent figures (`fig_02`, `fig_04`) guarded with `try/import linepyline` → print-and-skip when absent.

- [ ] **Step 2: Write the chapter `sources.yml`**

Create `docs/radiative-transfer/sources.yml`, one artifact per figure, e.g.:

```yaml
artifacts:
  - out: _artifacts/01_mean_of_exp.png
    cmd: conda run --no-capture-output -n climt python docs/radiative-transfer/_figures.py --figure 01_mean_of_exp --out docs/radiative-transfer/_artifacts/01_mean_of_exp.png
    deps:
      - docs/radiative-transfer/_figures.py
  - out: _artifacts/03_k_distribution_construction.png
    cmd: conda run --no-capture-output -n climt python docs/radiative-transfer/_figures.py --figure 03_k_distribution_construction --out docs/radiative-transfer/_artifacts/03_k_distribution_construction.png
    deps:
      - docs/radiative-transfer/_figures.py
  - out: _artifacts/06_picket_fence_opacity.png
    cmd: conda run --no-capture-output -n climt python docs/radiative-transfer/_figures.py --figure 06_picket_fence_opacity --out docs/radiative-transfer/_artifacts/06_picket_fence_opacity.png
    deps:
      - docs/radiative-transfer/_figures.py
      - climt/_components/picket_fence/optics/parmentier.py
  - out: _artifacts/07_two_stream_phases.png
    cmd: conda run --no-capture-output -n climt python docs/radiative-transfer/_figures.py --figure 07_two_stream_phases --out docs/radiative-transfer/_artifacts/07_two_stream_phases.png
    deps:
      - docs/radiative-transfer/_figures.py
      - climt/_components/picket_fence/sw/kernels.py
```

Add the `02_lbl_H2O_1000_1200` and `04_correlation_across_T` entries too; mark in a `# linepyline` comment that those regenerate only where `linepyline` is installed (the cmd still succeeds — the function prints-and-skips, but the figure won't be produced; for committed artifacts, generate them once on a machine with `linepyline` and commit the PNGs).

- [ ] **Step 3: Regenerate the figures**

Run: `make experiments` (B's driver walks `docs/**/sources.yml`, so the new manifest is picked up automatically).
Expected: PNGs appear under `docs/radiative-transfer/_artifacts/`. Confirm: `ls docs/radiative-transfer/_artifacts/`.

- [ ] **Step 4: Verify idempotence**

Run: `make experiments-check`
Expected: exit 0.

- [ ] **Step 5: Commit**

```bash
git add docs/radiative-transfer/_figures.py docs/radiative-transfer/sources.yml docs/radiative-transfer/_artifacts/
git commit -m "docs(rt): walkthrough figure pipeline wired through build_experiments"
```

---

## PART E — Radiative-transfer walkthrough chapters

Each chapter task creates one `.qmd`. The **prose is adapted from the committed phase-4 plan** (`docs/superpowers/plans/2026-04-20-picket-fence-radiation-phase4.md`, Part F) — that is the source draft; this plan specifies the Quarto transformation. Transformation rules applied to every chapter:
- `.rst` math `.. math::` → Quarto `$$ … $$`; inline `:math:` → `$…$`.
- `.. figure:: figures/NAME.png` → `![caption](_artifacts/NAME.png)`.
- `.. literalinclude:: … :pyobject: FUNC` → an executed code cell that prints the live source, so excerpts cannot drift:
  ````markdown
  ```{python}
  #| echo: false
  import inspect
  from climt._components.picket_fence.optics.correlated_k import load_k_table
  print(inspect.getsource(load_k_table))
  ```
  ````
- "Further reading" bullets → `@key` citations resolved against `references.bib` (add any missing entry with a `url=`).
- "Try it yourself" pointer → a `.callout-tip` linking to the relevant `examples/*.ipynb` (Part F) or test.
- Graduate-student physics derivations stay in the body; pitfalls go in `.callout-warning`.

### Task 8: Walkthrough index + chapter stubs that render

**Files:**
- Create: `docs/radiative-transfer/index.qmd`
- Create: `docs/radiative-transfer/{01-why-nongrey,02-line-by-line,03-k-distribution,04-correlated-k,05-gas-overlap,06-picket-fence,07-two-stream,08-multiplanet,table-generation,performance}.qmd` as minimal stubs (front-matter + one-line summary), so the site renders and the sidebar is complete from the start (spec success criterion: all chapters render with at least placeholders).

- [ ] **Step 1: Index page** — create `docs/radiative-transfer/index.qmd` adapting the phase4 Task 16 index prose (lines 1530–1546) into Quarto, with a chapter list.

- [ ] **Step 2: Stub each chapter** — for each filename above, write:

```markdown
---
title: "<Chapter Title>"
---

::: {.callout-note}
This chapter is being written. See the [walkthrough index](index.qmd).
:::
```

- [ ] **Step 3: Render the whole site** — `QUARTO_PYTHON=… quarto render docs/`; expected: clean, every sidebar link resolves.

- [ ] **Step 4: Commit** — `git add docs/radiative-transfer && git commit -m "docs(rt): walkthrough index + rendering chapter stubs"`.

### Task 9: Chapter 1 — "Why non-grey?"

**Files:** Create `docs/radiative-transfer/01-why-nongrey.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 17 prose (lines 1749–1837). Sections: "The mean of an exponential is not the exponential of the mean" (with `$$\langle e^{-\sigma L}\rangle$$` math), figure `![](_artifacts/01_mean_of_exp.png)`, "Non-grey phenomena you get for free" (stratospheric cooling, CO₂ forcing, solar stratospheric heating), a `.callout-tip` "Try it yourself" → `examples/spectral_radiation_anatomy.ipynb`, citations `@goody1989`, `@pierrehumbert2010` (add to `references.bib` with url=).
- [ ] **Step 2: Render** — `QUARTO_PYTHON=… quarto render docs/radiative-transfer/01-why-nongrey.qmd`; expected: figure + citations resolve.
- [ ] **Step 3: Commit** — `git commit -m "docs(rt): chapter 1 — why non-grey"`.

### Task 10: Chapter 2 — "Line-by-line physics"

**Files:** Create `docs/radiative-transfer/02-line-by-line.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 18 (lines 1862–1944). Sections: Voigt line shape (`$$\phi(\nu)$$` convolution), "What linepyline does" (a fenced `python` usage block — NOT executed, since `linepyline` may be absent), figure `![](_artifacts/02_lbl_H2O_1000_1200.png)`, "Why not just ship LBL?" (the 22 GB argument), `.callout-tip` → `examples/k_distribution_demo.ipynb`, citations `@rothman2013` (+ Goody & Yung).
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 2 — line-by-line physics"`.

### Task 11: Chapter 3 — "The k-distribution"

**Files:** Create `docs/radiative-transfer/03-k-distribution.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 19 (lines 1965–2054). Sections: "Band-averaging re-ordered" (the `$$\langle T(L)\rangle=\int_0^1 e^{-k(g)L}dg$$` derivation), "Constructing k(g) in practice" (sort → quadrature), figure `![](_artifacts/03_k_distribution_construction.png)`, "The code in climt" — live excerpt of `load_k_table` via the `inspect.getsource` cell pattern, `.callout-tip` → `examples/k_distribution_demo.ipynb` cells 5–7, citations `@lacis1991`, `@fu1992`.
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 3 — the k-distribution"`.

### Task 12: Chapter 4 — "Correlated-k"

**Files:** Create `docs/radiative-transfer/04-correlated-k.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 20 (lines 2075–2135). Sections: the correlated-k assumption (rank-order preserved with T,p), figure `![](_artifacts/04_correlation_across_T.png)`, "When it holds / breaks down", "What it gains", "Accuracy" (8 g-pts 1–3% vs 2 g-pts 5–10%), citations `@lacis1991`, `@mlawer1997` (already present).
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 4 — correlated-k"`.

### Task 13: Chapter 5 — "Gas overlap: additive vs ESFT"

**Files:** Create `docs/radiative-transfer/05-gas-overlap.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 21 (lines 2157–2219). Sections: additive overlap (`$$k_\text{total}(g)=\sum_i k_i(g) q_i$$`), ESFT (`$$k_{ij}, w_{ij}$$` outer product, `$G^2$` cost), live excerpt of `_esft_combine` via `inspect.getsource`, worked example (4 g-pts → 16; additive over-absorbs 20–30%), citations `@mlawer1997`, `@mitsel1995`. **Verify** `_esft_combine` exists: `conda run -n climt python -c "from climt._components.picket_fence.optics.correlated_k import _esft_combine"`; if the name differs, use the actual ESFT-combination function and update the excerpt import.
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 5 — gas overlap"`.

### Task 14: Chapter 6 — "The picket-fence model" + table-generation page

**Files:** Create `docs/radiative-transfer/06-picket-fence.qmd`, `docs/radiative-transfer/table-generation.qmd`

- [ ] **Step 1: Author chapter 6** — adapt phase4 Task 22 (lines 2239–2340; read this range for the full prose). Sections: the Parmentier picket-fence formulation, figure `![](_artifacts/06_picket_fence_opacity.png)`, live excerpt of the relevant `parmentier.py` function, `.callout-tip` → `examples/spectral_radiation_anatomy.ipynb`, citations `@parmentier2014` (and `parmentier2015` if cited).
- [ ] **Step 2: Migrate table-generation** — `quarto pandoc docs/radiative_transfer/table_generation.rst -f rst -t markdown -o docs/radiative-transfer/table-generation.qmd`, add front-matter, then `git rm -r docs/radiative_transfer` (the old Sphinx dir, now fully migrated).
- [ ] **Step 3: Render. Step 4: Commit** — `git commit -m "docs(rt): chapter 6 — picket-fence model + table-generation page"`.

### Task 15: Chapter 7 — "The two-stream solver"

**Files:** Create `docs/radiative-transfer/07-two-stream.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 23 (lines 2348–2455; read for full prose). Sections: the Meador & Weaver two-stream equations, figure `![](_artifacts/07_two_stream_phases.png)`, live excerpt of `_sw_dif_and_source` (or the LW solver kernel) via `inspect.getsource`, citations `@meador1980`. **Verify** the kernel import path `from climt._components.picket_fence.sw.kernels import _sw_dif_and_source`; correct it if the symbol moved.
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 7 — two-stream solver"`.

### Task 16: Chapter 8 — "Switching planets"

**Files:** Create `docs/radiative-transfer/08-multiplanet.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 24 (lines 2463–2580; read for full prose + figures). Sections: how the same scheme retargets to Mars/Venus/Titan/hot-Jupiters via swapped tables and Parmentier coefficients; include any multiplanet comparison figure the phase4 task defines (add it to `_figures.py` + `sources.yml` if it introduces a new figure; regenerate via `make experiments`), citations as listed in phase4 Task 24.
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): chapter 8 — switching planets"`.

### Task 17: Performance appendix

**Files:** Create `docs/radiative-transfer/performance.qmd`

- [ ] **Step 1: Author** — adapt phase4 Task 25 (lines 2581–2697; read for full prose). Reuse the throughput artifact B already produced (`docs/experiments/2026-06-05-picket-fence-co2-bands/_artifacts/throughput.png`) or add a dedicated perf figure to `_figures.py`/`sources.yml`. Cover the njit performance arc (faithful scheme faster than RRTMG).
- [ ] **Step 2: Render. Step 3: Commit** — `git commit -m "docs(rt): performance appendix"`.

---

## PART F — Pedagogical notebooks

These three notebooks are referenced by the chapters' "Try it yourself" callouts and embedded at chapter ends. Each: author with `nbformat`/Jupyter, execute in the `climt` env, commit pre-executed (same pattern proven in sub-project B).

### Task 18: `examples/k_distribution_demo.ipynb`

**Files:** Create `examples/k_distribution_demo.ipynb`

- [ ] **Step 1: Author** the 8 cells from phase4 Task 13 (lines 1441–1450): LBL step (`linepyline`), mean-of-exp vs exp-of-mean, sort+CDF, 2-gpoint quadrature, transmission reconstruction (assert max rel error < 2%), link to `earth_low_res_lw.nc` `k_coefficients`. **Guard** the `linepyline` cells with `try/except` + a pre-baked fallback array so the notebook executes even where `linepyline` is absent (mirror B's `RRTMG_AVAILABLE` pattern).
- [ ] **Step 2: Execute** — `cd examples && conda run -n climt jupyter nbconvert --to notebook --execute --inplace k_distribution_demo.ipynb`; expected exit 0.
- [ ] **Step 3: Embed** — add `{{< embed ../../examples/k_distribution_demo.ipynb echo=true >}}` at the end of `03-k-distribution.qmd`.
- [ ] **Step 4: Commit** — `git commit -m "docs(rt): k_distribution_demo notebook + embed"`.

### Task 19: `examples/spectral_radiation_anatomy.ipynb`

**Files:** Create `examples/spectral_radiation_anatomy.ipynb`

- [ ] **Step 1: Author** the cells from phase4 Task 14 (lines 1475–1484): load a standard Earth profile, run `PicketFenceLongwave(table="earth_low_res_lw")`, plot per-band τ(p), transmittance, net flux, heating-rate contributions, CO₂-doubling Δτ/Δheating/ΔOLR. Use the component's per-band diagnostics if exposed; otherwise compute from the public outputs. **Verify** the per-band diagnostic names against the shipped component before relying on them: `conda run -n climt python -c "from climt._components.picket_fence import PicketFenceLongwave; c=PicketFenceLongwave(); print([k for k in c.diagnostic_properties])"`.
- [ ] **Step 2: Execute. Step 3: Embed** at the end of `06-picket-fence.qmd`. **Step 4: Commit** — `git commit -m "docs(rt): spectral_radiation_anatomy notebook + embed"`.

### Task 20: `examples/picket_fence_vs_rrtmg.ipynb`

**Files:** Create `examples/picket_fence_vs_rrtmg.ipynb`

- [ ] **Step 1: Author** the three sections from phase4 Task 15 (lines 1500): Parmentier-mode hot-Jupiter T-p vs reference; correlated-k Earth clear-sky broadband flux vs RRTMG; discussion of agreement/divergence. Guard the RRTMG section with `RRTMG_AVAILABLE` (Fortran) so it degrades under Pyodide.
- [ ] **Step 2: Execute. Step 3: Embed** at the end of `08-multiplanet.qmd` (or chapter 4). **Step 4: Commit** — `git commit -m "docs(rt): picket_fence_vs_rrtmg notebook + embed"`.

---

## PART G — Deployment

### Task 21: GitHub Actions build + deploy workflow

**Files:** Create `.github/workflows/docs.yml`

- [ ] **Step 1: Write the workflow** (mirrors the spec's Deployment section):

```yaml
name: docs
on:
  push:
    branches: [develop, main]
  pull_request:
    branches: [develop, main]
  workflow_dispatch:

jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - uses: quarto-dev/quarto-actions/setup@v2
      - name: Install climt + docs deps
        run: |
          pip install quartodoc pyyaml jupyter matplotlib numpy
          pip install -e . --no-deps
      - name: Verify experiment artifacts are not stale
        run: python scripts/build_experiments.py --check
      - name: Build API reference
        run: cd docs && quartodoc build
      - name: Render site
        run: quarto render docs/
      - name: Upload site artifact
        uses: actions/upload-artifact@v4
        with:
          name: site
          path: docs/_site/
      - name: Deploy to gh-pages
        if: github.event_name == 'push'
        uses: peaceiris/actions-gh-pages@v4
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: docs/_site
          destination_dir: ${{ github.ref_name == 'main' && '.' || 'dev' }}
```

- [ ] **Step 2: Lint the workflow YAML** — `conda run -n climt python -c "import yaml; yaml.safe_load(open('.github/workflows/docs.yml')); print('ok')"`.

- [ ] **Step 3: Note the manual step** — GitHub Pages must be enabled once in repo settings (Source = `gh-pages` branch). This is a human action; record it in the PR description (cannot be automated here).

- [ ] **Step 4: Commit** — `git add .github/workflows/docs.yml && git commit -m "ci(docs): Quarto build + gh-pages deploy workflow"`.

### Task 22: Full-site verification

- [ ] **Step 1: Clean render from scratch**

Run: `QUARTO_PYTHON=/Users/joymonteiro/miniconda3/envs/climt/bin/python quarto render docs/ 2>&1 | tee /tmp/quarto_render.log | tail -20`
Expected: "Output created: …"; **no** unresolved `@key` (`grep -i "could not" /tmp/quarto_render.log` returns nothing), every navbar/sidebar link resolves.

- [ ] **Step 2: Artifact + test gates**

Run: `make experiments-check && conda run -n climt python -m pytest tests/test_tropopause.py tests/test_build_experiments.py -q`
Expected: experiments-check exit 0; tests pass.

- [ ] **Step 3: Link/citation audit**

Run: `grep -rn "???" docs/_site 2>/dev/null | head` (Quarto renders unresolved cross-refs as `???`). Expected: no hits in the chapter/API pages.

- [ ] **Step 4: Finish** — invoke `superpowers:finishing-a-development-branch` (this is a large change; a PR off `develop` is the natural choice, with the GH-Pages-enablement manual step called out in the PR body).

---

## Self-Review

**Spec coverage** (against `2026-05-18-climt-website-design.md`):
- Site structure (Get Started / User Guide / RT Walkthrough / Experiments / API / Contributing / References) → Tasks 1–6, 8–17; Experiments already shipped by B; References = `references.bib` (Tasks 4, 9–16).
- Quarto generator + GH Pages (dev/main) → Tasks 1, 21.
- RT walkthrough Ch 1–8 + perf + table-generation → Tasks 8–17.
- Three pedagogical notebooks, embedded in motivating chapters → Tasks 18–20.
- Artifact regeneration (`sources.yml` + `build_experiments.py` + `make experiments`) → reuses B's pipeline; extended in Task 7.
- quartodoc API reference replacing autodoc → Task 6.
- PDFs deleted, citations external-link only → Task 4.
- Sphinx machinery deleted → Task 5.
- CI fails on stale artifacts (`build_experiments.py --check`) → Task 21 Step + Task 22.
- Pyodide live cells → correctly OUT of scope (deferred by the spec; notebooks authored with degradation guards so a later follow-up can wire them).

**Decomposition note:** Per the user's choice, foundation and all 8 RT chapters are in this one plan. The deferred live-cells follow-up remains a separate future spec/plan (already noted in the website spec).

**Placeholder scan:** Chapter tasks point to the committed phase-4 Part F as the prose source of record (a stable reference document, not a TBD) and specify the exact Quarto transformation, figure, excerpt mechanism, and citations for each — no "fill in later" steps. Infra tasks carry concrete commands/code.

**Consistency:** figure filenames in `_figures.py`/`sources.yml` (Task 7) match the `![](_artifacts/NAME.png)` embeds in the chapter tasks; the `inspect.getsource` excerpt pattern and `{{< embed >}}` notebook pattern are used identically across chapters; `references.bib` `@key`s introduced in chapter tasks are added in the same task that cites them.
