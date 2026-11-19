"""Choose the longwave absorption table a Modelling Tour page runs on.

The high-resolution 56-band table is a site asset rather than wheel data --
``_tour/assets.py`` explains the mechanism and owns the path resolution. This
module is the policy on top of it: which table a page gets, and what happens
when the asset is absent.

Pages call :func:`spectrum_table` and pass the result straight to
``CorkLongwaveRadiation(table=...)``. If the asset cannot be found, the shipped
14-band table is returned instead, so a page degrades to a coarser spectrum
rather than failing -- which is also what happens on the pages that never ask
for it.
"""
import assets

FALLBACK = "earth_low_res_lw"
ASSET = "earth_spectrum_lw.npz"


def spectrum_table(prefer_hires=True, base_url=assets.DEFAULT_BASE_URL):
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
    return assets.locate(ASSET, base_url) or FALLBACK
