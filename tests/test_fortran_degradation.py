import climt


def test_has_fortran_extensions_is_bool():
    assert isinstance(climt.has_fortran_extensions(), bool)


def test_no_import_warning(recwarn):
    import importlib
    importlib.reload(climt)
    assert not [w for w in recwarn if "compiled" in str(w.message).lower()]
