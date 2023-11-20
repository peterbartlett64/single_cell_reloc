def remove_crap(str_in):
    if "test" in str_in:
        return None
    elif '.sync' in str_in:
        return None
    else:
        return (str_in)

