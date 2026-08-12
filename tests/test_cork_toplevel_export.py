import climt


def test_cork_is_top_level():
    assert hasattr(climt, "CorkLongwaveRadiation")
    assert hasattr(climt, "CorkShortwaveRadiation")
