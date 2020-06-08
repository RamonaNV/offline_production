#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/snobo/py3-v4.0.1/

import os
import warnings

# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu
# bsmithers@mail.icecube.wisc.edu

# this script can be used to create a series of dags you can submit
# it produces individual ones to be run with individual scripts

# TODO:
#       add argparse like interface for better UX
#       allow the output files to be saved somewhere that's not my personal directory 

# to someone who might want to use this:
#   + you may want to change which gcd file is specified in the 'write_dag_chain' function to one that's not in my /data/user folder
#   + you will need to change the 'workpath' variable so somewhere you have write permission
#   + depending on what kind of generation you want to run, you might have to change the hardcoded arguments in the 'write_dag_chain' function
# Not sure where these arguments come from? You may want to check the scripts that are called like
#   /path/to/script/1-process_gen.py --help

def make_condor_submit(execPath, name, gpus, memory, cpus=1):
    """
    writes the job.submit file to go with the dag. 
    execPath is which executable we will call, this writes the submit file with name "$name.submit" 
    will call for $gpus GPUs, $memory megabytes of memory, and $cpus cpus 
    """
    
    
    # just open a file, write in the file, and close the file...
    condor_file = open("{}.submit".format(name), 'w' )
    condor_file.write("executable = {}\n".format(execPath))
    condor_file.write("arguments  = $(arg)\n" )
    condor_file.write("\n")
    
    condor_file.write("should_transfer_files = true\n")
    condor_file.write("when_to_transfer_output = ON_EXIT\n")
    condor_file.write("getenv = False \n")
    condor_file.write("output = /scratch/bsmithers/{}.$(cluster).out\n".format(name))
    condor_file.write("error  = /scratch/bsmithers/{}.$(cluster).err\n".format(name))
    condor_file.write("log    = /scratch/bsmithers/{}.$(cluster).log\n".format(name))
    condor_file.write("\n")
    condor_file.write("request_memory = {}\n".format(memory))
    condor_file.write("request_gpus   = {}\n".format(gpus))
    condor_file.write("request_cpus   = {}\n".format(cpus))
    condor_file.write("queue 1\n")
    condor_file.close()




def write_dag(execPath,gen_args, outname, nJobs, memory = 1000, gpus = 0, seedroot = 5000, inName=None):
    """
    this script writes the dag. It's the one you'll call
    1.  execPath - where is the executable file?
    2.  gen_args - generic arguments that all of dag jobs will have (things like energies and spectal indices)
    3.  outname  - what pattern should be used to name the files. 
                        --->  $outname00000.i3.zst
    4.  nJobs    - how many jobs to run
    5.  meory    - amount of memory to request in the condor job
    6.  gpus     - number of GPUs to request
    7.  seedroot - some integer. Each job will be given a seed equal to '$seedroot + job number'
    8.  inName   - the name of the file to use as input. It is assumed these will have a similar naming scheme to the one used here. May be the previous dag's 'outname' argument

    OUTPUT: creates two text files. A $name.submit script and a $name_dag.dag script

    Note: you may want to use 'write_dag_chain', it adds the parent-child relations so you only have to submit one dag. 
    """
    # first try casting some arguments as numbers. Get all pissy if they aren't the right type
    if not (isinstance(memory, int) or isinstance(memory, float)):
        raise Exception("Expected type {} for memory. Got {}.",format(int, type(memory)))
    # why would you do this to me
    if not (isinstance(gpus, int) or isinstance(gpus, float)):
        raise Exception("Expected type {} for gpus. Got {}.".format(int, type(gpus)))
    # just in case the user ties something mean
    if not isinstance(execPath, str):
        raise Exception("execPath is type {}. Why would you do this to me".format(type(execPath)))
    # if the specified executable doesn't exist, something is very wrong. Tell the user.
    if not os.path.isfile(execPath):
        raise Exception("Executable doesn't exist at {}".format(execPath))
    # the specified executable is not an executable, the job will fail unless this is changed. Warn the user. 
    if not os.access(execPath, os.X_OK):
        warnings.warn("WARNING! CANNOT EXECUTE {}".format(execPath))

    job_name = outname.split("/")[-1]
    
    dag_file = open("{}_dag.dag".format(job_name),'w')
    
    suffix = ".i3.zst"
    print("Using file suffix: {}".format(suffix))
    workpath = "/data/user/bsmithers/runs/snobo_large/"

    
    for job in range(nJobs):
        # JOB thing00002 thing.submit
        dag_file.write("JOB  {}{}  {}.submit\n".format(job_name, get_5d_string(job), job_name ))
    for job in range(nJobs):
        # -o /path/to/output/thingout_00002.suffix 
        outfile = "-o {}{}_{}{}".format(workpath,outname,get_5d_string(job),suffix) 
        
        if inName is not None:
            # -i /path/to/output/thingin_000002.suffix
            infile = "-i {}{}_{}{}".format(workpath, inName, get_5d_string(job), suffix)
        else:
            infile = ""
        
        args =  "-s {} {} {} {}".format( seedroot+job, gen_args, outfile, infile ) 
        
        #VARS thing00002 -arg1 A -arg2 B -arg3 C [...]
        dag_file.write('VARS  {}{}  arg="{}"\n'.format(job_name, get_5d_string(job), args)) 

    dag_file.close()

    make_condor_submit(execPath, job_name , gpus, memory)


def write_job_script(outname, nJobs,execPath):
    """
    DEPRECATED
    """
    raise Exception("deprecated")
    
    # keeping this here because it's horrible and I love it 
    job_sub = open("job_submit.sh", 'w')
    job_sub.write("export execPath={}\n".format(execPath))
    job_sub.write("for (( FILE=1; FILE<={}; FILE++ ))\n".format(nJobs))
    job_sub.write("do\n")
    job_sub.write('    export word=""\n')
    job_sub.write("    for (( constructor=1; constructor<=5; constructor++ ))\n")
    job_sub.write("    do\n")
    #that's a bash implementation of the magic python script
    job_sub.write('        export word="$word$(($(($FILE/$((10**$((5-$constructor)))) ))%10))"\n')
    job_sub.write("    done\n")
    job_sub.write("    condor_submit {}-$word.submit\n".format(outname))
    job_sub.write("done\n")
    job_sub.close()
#    os.chmod("job_submit.sh",777)


# yes there are better ways of doing this.
# no, I don't care
# 
# to get the i'th digit, you shift the decimal point over a few orders of magnitude,
#   round the number (lobbing off the stuff to the right of the dot),
#   and find the modulus of 10, therefore lobbing off everything before the ones place
def get_5d_string( number ):
    """
    takes an integer < 100000, returns a string of the number that's the 5 characters long.
    If the number is < 10000, it puts zeros at the front to make the string 5 chars long
    """
    thing = ""
    for i in range(5):
        thing+= str( int(number/( 10**(4-i))) % 10 )
    return(thing)

def write_dag_chain(execPath, dataPath, base_name, nJobs, nEvents):
    """
    writes a dag file that will generate a sample all the way from 
    generation -> particle propagation -> photon propagation -> detector sim -> (more eventually)

    Arguments:
    1. folder to executables
    2. folder to save data
    3. base name for output data
    4. number of jobs
    5. number of events per job (default 2500)

    Several arguments are hard-coded in. These are all explicitly written for SnowSuite scripts, though it might be easy to modify to suit other needs  
    """
    
    # check to make sure these are the right data types and the folders actually exist
    if not isinstance(nJobs, int):
        raise Exception("nJobs should be type {}, received type {}".format(int, type(nJobs)))
    if not isinstance(nEvents, int):
        raise Exception("nEvents should be type {}, received type {}".format(int, type(nEvents)))
    if not os.path.exists(dataPath):
        raise Exception("Folder {} does not exist".format(dataPath))
    if not os.path.exists(execPath):
        raise Exception("Executable folder {} does not exist.".format(execPath))
    if not isinstance(base_name, str):
        raise Exception("base_name should be type {}, received {}".format(str, type(base_name)))
    
    # make sure my gcd file is still there
    gcd_file    = "/data/user/bsmithers/calibration/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_wNoiseUpdate.i3.gz"  
    if not os.path.isfile(gcd_file):
        raise Exception("Could not find gcd file at {}".format(gcd_file))

    seedroot    = 5000 # doesn't really matter what this is, but it gives control over the simulation
    
    #arguments are provided as a list ["generic arguments", "input file"]
    # the input thingy may be provided as none, in which case there are no inputs
    generator_name          = "gen"
    generator               = "{}/1-process-Gen.py".format(execPath)
    
    # TODO  may need to be frequently changed! 
    #
    #
    # 
    generation_arguments    = ["-i 1 --Emin 1000 --Emax 10000000 -t All --Zmin 0 --Zmax 180 -c CC -n {}".format(nEvents), None]
    #
    #
    # TODO
    

    propagator_name         ="prop"
    propagator              = "{}/2-Propagate.py".format(execPath)
    propagation_arguments   = ["-g {}".format(gcd_file), generator_name]

    events_per_model = 100
    snobo_name              = "snobo"
    snobo_prop       = "{}/3-SNOBO_Snowstorm.py".format(execPath)
    snobo_arguments         = ["-g {} --events-per-model {}".format(gcd_file, events_per_model), propagator_name]
    
    det_name = "det"
    detector = "{}/4-process-Det.py".format(execPath)
    det_arguments           = ["-g {} --dom_oversize=5.0 ".format(gcd_file), det_name ]

    # default I3 file suffix! 
    suffix = ".i3.zst"
    print("Using file suffix: {}".format(suffix))
    
    
    # dictionary for the names of each level of the MC generation and their arguments 
    level_arguments = {base_name+generator_name : generation_arguments, base_name+propagator_name : propagation_arguments, base_name+snobo_name : snobo_arguments, base_name+det_name: det_arguments }
    
    # make sure each of the necessary executables exist and are executable. If they are not executable just issue a warning! 
    for executable in [generator, propagator, snobo_prop ,detector]:
        if not os.path.isfile(executable):
            raise Exception("Executable doesn't exist at {}".format(executable))
        # the specified executable is not an executable, the job will fail unless this is changed. Warn the user. 
        if not os.access(executable, os.X_OK):
            warnings.warn("WARNING! CANNOT EXECUTE {}".format(executable))
    
    
    dag_file = open("chain_dag.dag",'w')
    
    
    
    # use that dictionary, get the keys (which will be the data file names too)
    file_names = [base_name+generator_name, base_name+propagator_name, base_name+snobo_name, base_name+det_name ]
    for outname in file_names:
        for job in range(nJobs):
            # JOB thing00002 thing.submit
            dag_file.write("JOB  {}{}  {}.submit\n".format(outname, get_5d_string(job), outname ))
    
    # write the arguments from each dictionary entry into the dag
    for outname in file_names:
        for job in range(nJobs):
            # -o /path/to/output/thingout_00002.suffix 
            outfile = "-o {}{}_{}{}".format(dataPath,outname,get_5d_string(job),suffix) 
            
            if level_arguments[outname][1] is not None: 
                # gen is the only level that doesn't take an input
                # -i /path/to/output/thingin_000002.suffix
                infile = "-i {}{}_{}{}".format(dataPath, base_name+level_arguments[outname][1], get_5d_string(job), suffix)
            else:
                infile = ""
            
            args =  "-s {} {} {} {}".format( seedroot+job, level_arguments[outname][0] , outfile, infile ) 
            
            #VARS thing00002 -arg1 A -arg2 B -arg3 C [...]
            dag_file.write('VARS  {}{}  arg="{}"\n'.format(outname, get_5d_string(job), args)) 
    
    #write the family dag tree
    for job in range(nJobs):
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[0], get_5d_string(job), file_names[1],get_5d_string(job)))
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[1], get_5d_string(job), file_names[2],get_5d_string(job))) 
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[2], get_5d_string(job), file_names[3],get_5d_string(job)))
    dag_file.close()
    
    # make the .submit files
    make_condor_submit(generator,  file_names[0], gpus=0, memory=1000)
    make_condor_submit(propagator, file_names[1], gpus=1, memory=1000)
    make_condor_submit(snobo_prop, file_names[2], gpus=1, memory=2000, cpus=2)
    make_condor_submit(detector,   file_names[3], gpus=0, memory=2000)


nJobs = 200

root_dir = "/data/user/bsmithers/snobo/simprod-scripts/resources/scripts/SnowSuite/"
workpath = "/data/user/bsmithers/runs/snobo_large/CC_runs"

# exec dir | data output | name | jobs | events per job
write_dag_chain(root_dir, workpath, "CC", nJobs, 2500)

# write_dag(execPath,gen_args, outname, nJobs, memory = 1000, gpus = 0, seedroot = 5000, inName=None)

gcd_file    = "/data/user/bsmithers/calibration/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_wNoiseUpdate.i3.gz"  
#snobo_arguments         = ["-g {} --events-per-model {}".format(gcd_file, 100), "prop"]
#write_dag( root_dir+"3-SNOBO_Snowstorm.py",  snobo_arguments[0], "CC_hobo", nJobs, 8000, gpus=1, seedroot=5000, inName="CC_prop")

# execPath | gen_args | outname  | nJobs | memory | gpus | seedroot | inName
#det_arguments           = ["-g {} --dom_oversize=5.".format(gcd_file)]
#write_dag( root_dir+"4-process-Det.py", det_arguments[0], "/NC_runs/detector_level/NCdet", nJobs, memory=2000, seedroot=5000, inName="/NC_runs/photon_level/NCsnobo")
