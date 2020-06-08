"""
copyright  (c) 2005 the icecube collaboration
@version: $Revision: $
@date: $Date: $
@author: Juan Carlos Diaz Velez <juancarlos@icecube.wisc.edu>
@author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
"""

class Option:
    def __init__(self, name=None,
                 short_name=None, long_name=None,
                 type=str, action='store', default=None,
                 help='', const=None,
                 group=None, choices=None):
        self.name=name
        self.short_name=short_name
        self.long_name=long_name
        self.action=action
        self.default=default
        self.help=help
        self.const=const
        self.group=group
        self.choices=choices
        self.type=type
        if action in ['store_const', 'store_true', 'store_false', 'count', 'callback', 'append_const']:
            self.type=None
        if self.choices:
            self.type=None
        if not self.short_name and not self.long_name and not self.name:
            raise Exception('name, long_name and short_name are false')
        if not self.name:
            self.name = self.long_name[2:].replace('-','_') if self.long_name else self.short_name[1:].replace('-','_')
class OptionGroup:
    def __init__(self, name, title, help=''):
        self.name = name
        self.title = title
        self.help = help

class OptionSuite:
    '''
    FIXME : Add docs.
    '''
    def __init__(self, options):
        self.options = options

    def configure(self, options):
        '''
        This method will be called after "parse_args" and is meant
        to be used to set defaults that depend on other user defined
        parameters that could potentially belong to a different policy.
        The result of the parse_args, 'options' is what's passed here.
        '''
        self.options = options

    def append_to_steering(self):
        '''
        Derived modules should return a string that will get appended
        directly to the INPUTS steering file.
        '''
        return ''

    def validate(self):
        '''
        Check to see that the input is valid and self-consistent.
        '''
        return True


