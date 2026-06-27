# Contributing to PrecisionNutrition

Welcome! This repository holds the raw data, analysis scripts, and manuscript
files for precision-nutrition studies from the [Bridges Lab](http://bridgeslab.sph.umich.edu)
and our collaborators. Whether you're a new lab member finding your way around,
a long-time collaborator adding an analysis, or someone from outside the lab who
spotted a problem, thank you for being here — this guide will get you oriented.

If anything below is unclear or out of date, that itself is worth an issue or a
pull request. Improving these guidelines counts as a contribution.

---

## How this repository is organized

Each top-level folder is a self-contained line of work, and **every folder has a
`README.md` that orients you** before you open a single script. Start there.

| Folder | What lives here |
|---|---|
| `Manuscripts` | Text and figure files for manuscripts. Each subfolder's README points to the underlying data and scripts. |
| `Mouse Genetics` | Mouse datasets with their analysis and figure-generation scripts. |
| `Meta Analysis` | Meta-analyses, generally of human data. |

The general pattern throughout: **raw data is provided in raw form, and an R
script (as an `.Rmd`/`.qmd`) analyzes it, rendering to `.html` and `.md` outputs
alongside the figures.** When you add new work, follow this pattern rather than
inventing a new one — predictability is what makes the repo navigable for the
next person.

---

## The workflow at a glance

We do active work on **branches** and keep `main` in a state we'd be comfortable
sharing. The short version:

1. **Branch** off `main` for your analysis or change.
2. **Work** on that branch — commit as you go.
3. **Open a pull request** when it's ready for a second set of eyes.
4. **Review and merge** into `main`.
5. **Tag and release** when the work reaches a milestone (typically preprint or
   publication).

The sections below expand on each step.

### 1. Work on a branch

Create a branch named for what you're doing, e.g. `calcium-mr-leave-one-out` or
`fix-locf-cap-bug`:

```bash
git checkout main
git pull
git checkout -b short-descriptive-name
```

Keep a branch alive until the work is ready to share — for a new analysis, that
usually means **ready for preprint submission**. This keeps `main` clean and
citable while exploratory work is in flight.

### 2. Commit as you go

- Commit in logical chunks rather than one giant commit at the end.
- Write commit messages that say *what changed and why*, e.g.
  `Add leave-one-out sensitivity for calcium MR instrument` rather than
  `update`.
- Don't edit raw data in place to "fix" an analysis — see
  [Data discipline](#data-discipline) below.

### 3. Open a pull request

When the work is ready, open a PR into `main`. A good PR description says what
the analysis does, what changed, and anything a reviewer should look at closely.
Even within the lab, **one other person reviewing before merge** catches a
surprising number of issues — please don't merge your own PR straight to `main`
without a look unless it's trivial (a typo, a README clarification).

Use [issues](https://github.com/BridgesLab/PrecisionNutrition/issues) to track
TODOs, bugs, and open questions, and link PRs to the issues they close.

### 4. Tag and release at milestones

This repo pins **dataset tags to publications** — the tag records the state of
the data and analysis at the moment a manuscript was submitted or published. When
work reaches that milestone, the release ritual is:

1. Merge the PR into `main`.
2. Tag the release (e.g. `1.1.0`) — see [Versioning](#versioning).
3. Cut a [Zenodo](https://zenodo.org) release so the snapshot gets a DOI.
4. Add a row to the publication table in the top-level `README.md` linking the
   publication, the dataset, and the tag.
5. Update `CITATION.cff` so citation metadata stays in sync.

Doing these together is what keeps the tag, the DOI, the README table, and the
citation file from drifting apart.

---

## Reproducibility

The whole point of this repo is that someone can re-run the analysis and get the
same result. To keep that true:

- **Render from a fresh clone before you submit.** If it only renders on your
  machine, it isn't reproducible yet.
- **Use relative paths.** No absolute paths like `/Users/yourname/...` — they
  break for everyone else.
- **Record your environment.** Include `sessionInfo()` at the end of each
  rendered document (already standard in our QMDs); for heavier dependency
  needs, consider [`renv`](https://rstudio.github.io/renv/) to pin package
  versions.
- **Set seeds** for anything stochastic, so resampling, permutation, and
  simulation results are reproducible.
- **One analysis = one `.Rmd`/`.qmd`** that renders to `.html` + `.md`. We
  commit the rendered outputs alongside the source on purpose, so results are
  visible without re-running — please keep rendered outputs current when you
  change the source.

---

## Data discipline

- **Raw data is read-only.** Treat anything in a raw-data location as immutable.
  If a value is wrong, document and correct it in the analysis script (or a
  clearly named cleaning step), not by hand-editing the raw file.
- **Derived data is regenerated, not hand-patched.** Cleaned/derived files
  should be produced by a script that anyone can re-run, following the lab
  convention of cleaning scripts that write to CSV.
- **Be explicit about provenance.** Each folder README should say which raw
  files feed which scripts and produce which outputs.

---

## Naming and hygiene

- **File and folder names:** prefer lowercase with hyphens or underscores and no
  spaces for anything new. Spaces in names (e.g. `Meta Analysis`) occasionally
  trip up scripts and URLs; we won't rename existing folders gratuitously, but
  new additions should avoid them.
- **Keep OS and editor cruft out of the repo.** Files like `.DS_Store`,
  `.Rhistory`, and `.Rproj.user/` shouldn't be committed — they belong in
  `.gitignore`. If you spot one already tracked, removing it is a welcome PR.
- **Large or sensitive data:** don't commit anything that can't be shared under
  the repository's license, and check with the lab before adding large binary
  data.

---

## Versioning

Tags follow `MAJOR.MINOR.PATCH` (e.g. `1.0.0`):

- **MAJOR** — a new study/dataset or a breaking reorganization.
- **MINOR** — a new analysis or substantive addition to existing work.
- **PATCH** — corrections and small fixes that don't change conclusions.

Each tag should correspond to a Zenodo release and, where applicable, a row in
the README publication table.

---

## Pre-submission checklist

Before opening a PR that's meant to accompany a preprint or publication:

- [ ] Renders cleanly from a fresh clone
- [ ] All paths are relative
- [ ] Rendered `.html`/`.md` outputs and figures are up to date
- [ ] `sessionInfo()` (or `renv` lockfile) is current
- [ ] The relevant folder README is updated
- [ ] Seeds set for any stochastic steps
- [ ] No raw data edited by hand
- [ ] No stray OS/editor files added

---

## License

Data in this repository is released under the
[Open Data Commons Attribution License](http://opendatacommons.org/licenses/by/1.0).
By contributing, you agree that your contributions will be shared under the same
terms. See [`LICENSE.txt`](LICENSE.txt).

---

## Contributors

We credit contributions using the [CRediT taxonomy](https://credit.niso.org/)
(Contributor Roles Taxonomy). This makes the repository a running record of who
did what, and maps directly onto the author-contributions statement for the
resulting manuscripts. Please keep this section consistent with `CITATION.cff`.

When you make a meaningful contribution, add yourself (or ask a maintainer to add
you) with your ORCID and your role(s). Common CRediT roles:
*Conceptualization, Data curation, Formal analysis, Funding acquisition,
Investigation, Methodology, Project administration, Resources, Software,
Supervision, Validation, Visualization, Writing – original draft, Writing –
review & editing.*

| Contributor | ORCID | Role(s) |
|---|---|---|
| Dave Bridges | [0000-0002-5334-972X](https://orcid.org/0000-0002-5334-972X) | Conceptualization, Supervision, Funding acquisition, Methodology, Software, Writing – review & editing |
| *Your name here* | `0000-0000-0000-0000` | *e.g. Formal analysis, Software, Visualization* |

> Note: confirm the ORCID above is correct before relying on it — it's filled in
> as a placeholder from the maintainer's name and should be verified.

---

## Questions

Open an [issue](https://github.com/BridgesLab/PrecisionNutrition/issues) or reach
out to the lab. No contribution is too small — fixing a broken path or clarifying
a README genuinely helps.
