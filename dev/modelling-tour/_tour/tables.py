"""Locate the spectrum table, in the browser and natively.

The high-resolution table is too large to ship inside the wheel, so it lives
beside these pages as a static site asset. Pages 1-3 declare it in their
``pyodide: resources:`` front matter exactly the way they declare
``_tour/*.py``, so quarto-live stages it into the Pyodide filesystem at document
setup, at the same relative path it has on disk. That is why the lookup below is
a plain ``os.path.isfile`` in both environments -- the browser case needs no
network code of its own, and no CORS negotiation, because the docs site is
same-origin with the page. (A GitHub release asset would not work: it sends no
CORS headers.)

Pages call :func:`spectrum_table` and pass the result straight to
``CorkLongwaveRadiation(table=...)``. If the asset cannot be found the shipped
14-band table is returned instead, so a page degrades to a coarser spectrum
rather than failing -- which is also what happens on the two pages that never
ask for it.
"""
import os

FALLBACK = "earth_low_res_lw"
ASSET = "earth_spectrum_lw.npz"


def _candidates(base_url):
    """Paths to try, in order, for the staged asset."""
    yield os.path.join(base_url, ASSET)
    # `_tour/` is one level below the page, so the asset sits next to our own
    # parent directory. This is what makes the function work from a native
    # render or a test run, where the working directory is not the page's.
    here = os.path.dirname(os.path.abspath(__file__))
    yield os.path.join(here, os.pardir, base_url, ASSET)


def spectrum_table(prefer_hires=True, base_url="_data"):
    """Return a table name or path for ``CorkLongwaveRadiation``.

    Args:
        prefer_hires: if False, always return the shipped 14-band table.
        base_url: directory (or URL prefix) holding the asset, relative to the
            page.

    Returns:
        A table name or filesystem path, always usable as ``table=``.
    """
    if not prefer_hires:
        return FALLBACK

    for path in _candidates(base_url):
        if os.path.isfile(path):
            return path

    # Last resort: the page did not declare the asset as a resource, but we are
    # in the browser and it is same-origin, so fetch it once into the Pyodide
    # filesystem. Synchronous by necessity -- quarto-live does not await a
    # cell's async work before starting the next one.
    try:
        from pyodide.http import open_url

        data = open_url(f"{base_url}/{ASSET}")
        with open(ASSET, "wb") as handle:
            handle.write(data.getvalue().encode("latin-1"))
        return ASSET
    except Exception:
        return FALLBACK
