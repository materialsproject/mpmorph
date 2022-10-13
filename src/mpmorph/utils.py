import collections
from datetime import datetime


def recursive_update(orig_dict, new_dict):
    for key, val in new_dict.items():
        if isinstance(val, collections.Mapping):
            tmp = recursive_update(orig_dict.get(key, {}), val)
            orig_dict[key] = tmp
        elif isinstance(val, list):
            orig_dict[key] = orig_dict.get(key, []) + val
        else:
            orig_dict[key] = new_dict[key]
    return orig_dict


def datetime_str() -> str:
    """
    Get a string representation of the current time. Borrowed from atomate2.
    """
    return str(datetime.utcnow())
