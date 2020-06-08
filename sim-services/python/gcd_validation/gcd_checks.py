'''
A collection of higher level GCD validation tests.
'''

def all_standard_inice_strings_exist(dom_geo_map):
    '''
    Returns True if at least one DOM from the standard IC86
    configuration is found.  Returns False if whole strings
    are missing from the geometry.
    '''
    all_strings_exist = True
    for string in range(1,87):
        found = False
        for omkey, i3omgeo in dom_geo_map:
            if omkey.string == string :
                found = True
        if not found:
            print("string %d does not exist in the geometry" % string)
            all_strings_exist = False
    return all_strings_exist

