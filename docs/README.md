# Applying the climt design system to the Quarto docs site

The live docs are migrating to **Quarto with the `cosmo` theme** (`docs/_quarto.yml`).
Quarto themes are just **Bootstrap + SASS**, so adopting this design system means
layering one SCSS file on top of `cosmo` that maps Bootstrap's SASS variables onto
the climt tokens. No fork of the theme needed.

## What's here
- `climt.scss` — the drop-in theme layer. A `scss:defaults` block sets brand color,
  fonts, surfaces, radii and the four Quarto callout colors from the climt tokens;
  a `scss:rules` block adds the warm-dark code blocks, the ochre sidebar rail and
  rounded cards.
- `index.qmd` — the **marketing landing page** as a Quarto custom-layout page.
  Self-contained: design-system tokens are inlined into a scoped `<style>` block, no
  React/Babel, no dependency on the design-system bundle. Drop it in and it builds.
- `_preview.html` / `_preview_shot.png` — a standalone render of `index.qmd` for review
  here (not for the repo).

## Adding the landing page (`index.qmd`)

1. **Copy** `index.qmd` to the docs root (e.g. `docs/index.qmd`) — it becomes the site home.
2. **Copy the logo** to `docs/assets/climt-logo.jpg` (the page references `assets/climt-logo.jpg`).
3. Fix up the nav/footer links to your real page slugs (`quickstart.html`, `user-guide.html`, …).
4. If you use a Quarto **`website`** project with a top navbar, either set
   `page-navbar: false` in the front matter to let this page own its chrome, or delete the
   page's own `<nav class="ch-nav">` block and keep the site navbar. Both work.

The page already pulls IBM Plex Sans + JetBrains Mono from Google Fonts via `include-in-header`;
self-host them the same way as `climt.scss` if the site must work offline.

## How to wire it up (2 edits in the repo)

1. **Copy** `climt.scss` into the docs folder, e.g. `docs/climt.scss`.
2. **Edit `docs/_quarto.yml`** — point the HTML format at it:

   ```yaml
   format:
     html:
       theme:
         - cosmo
         - climt.scss      # <-- layer climt on top of cosmo
       toc: true
       highlight-style: a11y-dark   # warm-dark code; pairs with climt.scss
   ```

That's the whole adoption. Callouts (`::: {.callout-note}` etc.) automatically pick up
the climt note/tip/warning/important colors; headings, links, code and surfaces follow.

## About the branch / PR

I can read and import from the repo, but I **can't push commits or open a PR from here** —
GitHub access in this workspace is read-only. So the flow is:

1. I generate / iterate on `climt.scss` (and any extra CSS) here until you're happy.
2. You (or I hand you the exact files + the `_quarto.yml` diff) commit them on a branch —
   e.g. `docs/climt-theme` — and open the PR against `develop`.

If you'd like, I can also produce: a ready-to-paste **PR description**, a logo asset sized
for the navbar/favicon, and a matching `theme-dark.scss` for a dark-mode toggle. Just say so.

## Notes & caveats
- `climt.scss` is a faithful first pass from the tokens; exact callout/border/code shades
  are easy to tune once you see it rendered by Quarto.
- The DS code blocks are dark; I set `highlight-style: a11y-dark` as the closest Quarto
  built-in. A bespoke climt Pygments/highlight theme can be generated if you want the
  ochre/sage/ocean syntax accents exactly.
- Fonts load from Google Fonts via `@import` in the SCSS. If the site must work offline,
  self-host IBM Plex Sans + JetBrains Mono and swap the `@import` for `@font-face`.
