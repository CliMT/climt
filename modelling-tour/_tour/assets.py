"""Locate a Modelling Tour site asset, in the browser and natively.

Some of what these pages need is too large or too configuration-specific to
ship inside the climt wheel: the 56-band spectrum table (~5 MB) and the two
radiative-convective equilibrium states. Those live beside the pages as static
site assets under ``_data/``.

``docs/_quarto.yml`` lists ``modelling-tour/_data/*.npz`` under
``project: resources:`` so Quarto publishes them at all (it skips
underscore-prefixed paths otherwise), and each page lists the ones it uses
under ``pyodide: resources:`` so quarto-live stages them into the Pyodide
filesystem at document setup, at the same relative path they have on disk.

That is why the lookup below is a plain ``os.path.isfile`` in both
environments: the browser case needs no network code of its own and no CORS
negotiation, because the docs site is same-origin with the page. (A GitHub
release asset would not work -- it sends no CORS headers.)

This module is the single place that knows any of the above.
``tables.py`` and ``states.py`` are its callers.
"""
import os

DEFAULT_BASE_URL = "_data"


def _candidates(name, base_url):
    """Paths to try, in order, for a staged asset."""
    yield os.path.join(base_url, name)
    # `_tour/` is one level below the page, so the asset sits next to our own
    # parent directory. This is what makes resolution work from a native
    # render or a test run, where the working directory is not the page's.
    here = os.path.dirname(os.path.abspath(__file__))
    yield os.path.join(here, os.pardir, base_url, name)


def resolve(name, base_url=DEFAULT_BASE_URL):
    """Return a filesystem path to ``name``, or None if it is not staged.

    Args:
        name: the asset's file name, e.g. ``"rce_dry_equilibrium.npz"``.
        base_url: directory holding it, relative to the page.

    Returns:
        A path that ``open()`` will accept, or ``None``.
    """
    for path in _candidates(name, base_url):
        if os.path.isfile(path):
            return path
    return None


def fetch_same_origin(name, base_url=DEFAULT_BASE_URL):
    """Browser-only last resort: fetch a same-origin asset into the FS.

    Used when a page forgot to declare the asset under ``pyodide: resources:``.
    Synchronous by necessity -- quarto-live does not await a cell's async work
    before starting the next one, so an ``await``-based fetch would race.

    Returns the local file name on success, or ``None`` outside Pyodide or if
    the fetch fails for any reason.
    """
    try:
        from pyodide.http import open_url

        data = open_url("{}/{}".format(base_url, name))
        with open(name, "wb") as handle:
            handle.write(data.getvalue().encode("latin-1"))
        return name
    except Exception:
        return None


def locate(name, base_url=DEFAULT_BASE_URL):
    """``resolve``, then the browser fetch. ``None`` if neither works."""
    return resolve(name, base_url) or fetch_same_origin(name, base_url)
