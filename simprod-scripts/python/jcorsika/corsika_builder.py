#!/bin/env python

import os, sys
import corsika_options

class CorsikaBuilder(object):
    def __init__(self, version):
        from corsika_options.coconut_options import CoconutOptionCollection
        self.options = CoconutOptionCollection(version)

    def coconut(self, options, log_file, suffix=None, cwd=None):
        """
        Spawn a process to compile CORSIKA.

        This function calls coconut and passes the options required through coconut's stdin.
        If a suffix is provided, the executable name is renamed by replacing everything
        after the first occurrence of '_'.

        The options parameter must be a dictionary like this:
            {'float': '2',
             'geometry': '1',
             'he_model': '6',
             'le_model': '2',
             'other': ['9', '6', 'i', 's']}
        Note that each key corresponds to a section in coconut and the values are the ones that correspond to each option.
        """
        import subprocess
        import re, os, shutil
        options['other'] = ' '.join(sorted(options['other']))
        option_string = "{float}\n{he_model}\n{le_model}\n{geometry}\n{other}\n\n\n\n"
        coconut_in = option_string.format(**options)
        if cwd is None:
            cwd = os.getcwd()
        if suffix is None:
            suffix = '_'.join(coconut_in.strip().split())

        if self.options.option_from_name['coast'][1] in options['other']:
            if not 'COAST_USER_LIB' in os.environ:
                raise Exception('Trying to compile with COAST, but COAST_USER_LIB is not in the environment')
            if not os.path.exists(os.path.join(os.environ['COAST_USER_LIB'], 'libCOAST.so')):
                raise Exception('Trying to compile with COAST, but libCOAST.so is not in COAST_USER_LIB')
        if options['le_model'] == self.options.option_from_name['fluka'][1] and (not 'FLUPRO' in os.environ or not os.path.exists(os.path.join(os.environ['FLUPRO'], 'libflukahp.a'))):
            raise Exception('Trying to build with Fluka but either FLUPRO is not defined or the Fluka library is not there')

        #subprocess.call(['env'], stdout=log_file, stderr=log_file)
        p = subprocess.Popen(['./coconut'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=log_file, cwd=cwd)
        coconut_out = p.communicate(input=coconut_in)[0]
        log_file.write(coconut_out)
        m = re.search('"(\S+)" successfully installed', coconut_out)
        if not m:
            print("could not figure out the name of binary file")
            sys.exit(1)
        binary_file = m.group(1)
        if suffix:
            # executable needs to be renamed
            binary_file_2 = m.group(1).split('_')[0] + '_' + suffix
            log_file.write('{cwd}/run/{binary_file} -> {cwd}/run/{binary_file_2}'.format(**{'cwd':cwd, 'binary_file_2':binary_file_2, 'binary_file':binary_file}))
            shutil.move(os.path.join(cwd, 'run', binary_file), os.path.join(cwd, 'run', binary_file_2))
            binary_file = binary_file_2
        subprocess.call(['./coconut', '-d'], stdout=log_file, stderr=log_file, cwd=cwd)
        return binary_file

    def build_fluka(self, fluka_dir, fluka_tarball, fortran='gfortran', force=False, log_file=sys.stdout, **kwargs):
        """
        Create fluka directory, cd into it and build fluka. If force==True, previous build is removed.

        This is not a general purpose script!
        """
        import subprocess
        import shutil
        old_path = os.getcwd()
        if not force and os.path.exists(fluka_dir) and  os.path.exists(os.path.join(fluka_dir, 'libflukahp.a')):
            return
        os.environ['FLUFOR']=fortran
        if os.path.exists(fluka_dir):
            shutil.rmtree(fluka_dir)
        os.makedirs(fluka_dir)
        os.chdir(fluka_dir)
        subprocess.call(['tar', 'xzf', fluka_tarball], stdout=log_file, stderr=log_file)
        subprocess.call(['make'], stdout=log_file, stderr=log_file)
        del os.environ['FLUFOR']
        os.chdir(old_path)

    def build(self, options, suffix=None, cwd=None, fluka_dir=None, fluka_tarball=None, log_file=sys.stdout):
        """
        convenience function that just compiles Fluka and then calls coconut
        """
        if fluka_dir:
            if not fluka_tarball:
                raise Exception('If you want Fluka, you have to specify the --fluka-tarball option.')
            if not os.path.exists(fluka_tarball):
                raise Exception('The fluka tarball is not found.')
            os.environ['F77']='gfortran'
            print("Building Fluka")
            self.build_fluka(fluka_dir, fluka_tarball, log_file=log_file)
            os.environ['FLUPRO'] = fluka_dir
        print("Running coconut")
        self.coconut(options, log_file, cwd=cwd, suffix=suffix)



if __name__ == "__main__":
    try:
        from corsika_options import get_version
        version = get_version.get_version()
    except Exception as e:
        print("Exception: " + str(e))
        sys.exit(1)

    builder = CorsikaBuilder(version)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', required=True)
    parser.add_argument('--directory', default=os.environ['PWD'])
    parser.add_argument('--fluka-tarball')

    required = ('float', 'he_model', 'le_model', 'geometry')
    action = {'float':'store_const', 'he_model':'store_const', 'le_model':'store_const', 'geometry':'store_const', 'other':'append_const'}
    option_groups = {}
    for g in builder.options.option_groups:
        option_groups[g] = parser.add_argument_group(g)
        if g in required:
            option_groups[g] = option_groups[g].add_mutually_exclusive_group(required=True)
    for n,o in builder.options.options.iteritems():
        option_groups[o.group].add_argument('--'+o.name, dest=o.group, action=action[o.group], const=o.const)

    args = parser.parse_args()
    if not os.path.isabs(args.directory):
        args.directory = os.path.join(os.environ['PWD'], args.directory)
    if args.other is None:
        args.other = []

    log_file = open("corsika_build.log", "w")
    fluka_dir = os.path.join(args.directory, 'fluka') if args.fluka_tarball else None
    builder.build(vars(args), cwd=args.directory, fluka_tarball=args.fluka_tarball, fluka_dir=fluka_dir, log_file=log_file)
    log_file.close()
