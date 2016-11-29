from .exceptions import SharedKeyException


def add_dicts_inplace(dict1, dict2):
    """
    Takes two dictionaries. For any keys in both dictionary, it adds the value
    in dict2 to the value in dict1. This is done in-place if the values are
    array-like, to avoid data copying. None is returned.
    """
    for key in dict2.keys():
        if key not in dict1:
            dict1[key] = dict2[key]
        else:
            try:
                # works for array-like objects, in-place
                dict1[key][:] += dict2[key][:]
            except TypeError:
                dict1[key] += dict2[key]
    return  # not returning anything emphasizes that this is in-place


def ensure_no_shared_keys(dict1, dict2):
    """
    Raises SharedKeyException if there exists a key present in both
    dictionaries.
    """
    if len(set(dict1.keys()).intersection(dict2.keys())) > 0:
        raise SharedKeyException()
