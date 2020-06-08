How to run a test of production scripts
=======================================

Set up environment
~~~~~~~~~~~~~~~~~~

First, build iceprod and load it's env-shell.sh.
Then build a metaproject or download a tarball and load it's env-shell.sh.

Get a dataset config
~~~~~~~~~~~~~~~~~~~~

.. highlight:: bash

Use iceprodsh and ipstranslate to make python files that you can use::

    $ iceprodsh
    iceprodsh> download 10036
    iceprodsh> save script.xml
    iceprodsh> exit
    adios.
    $ iptranslate -l python -o script.py -v 0 script.xml
    script.xml ----> script_steering.py
    script.xml ----> script0.py
    script.xml ----> script1.py
    script.xml ----> script2.py
    script.xml ----> script3.py
    script.xml ----> script4.py
    script.xml ----> script5.py
    script.xml ----> script6.py
    script.xml ----> script7.py
    script.xml ----> script8.py
    
.. highlight:: python

Edit script_steering.py if necessary to change some variables::

    # Default values
    cmd_opts = {
       'seed':0, 'procnum':0, 'nproc':1, 'tray':0, 'iter':0,
       'dataset': 0,
       'fetch': 'http://convey.icecube.wisc.edu/data/sim/sim-new/downloads',
       # a better fetch alternative: 'file:/data/sim/sim-new/downloads'
       }

.. highlight:: bash

Set some environment variables::
    
    $ export PHOTON_TABLES=/data/sim/sim-new/PhotonTablesProduction
    
    Others that can be set:
    
        JAVA_HOME
        SCRATCH
        GLOBUS_LOCATION
        X509_USER_PROXY

Run a script
~~~~~~~~~~~~

Use those default values above as options to the script::
    
    $ python script0.py --seed=0 --procnum=0 --nproc=1 --tray=0 --iter=0 --dataset=10036
    
.. highlight:: python

Note that older datasets and occasionally parts of datasets attempt to load module scripts from the internet.  A typical example of this is::

    # Configure PreTray modules 
    pre.AddModule("i3.IceTray","generate_BG")(
        ("mjd",int(exparser.parse("$steering(mjd_11)"))),
        ("showers",int(exparser.parse("$eval($steering(CORSIKA::showers) * 2)"))),
        ("outputfile",str(exparser.parse("$steering(current_file_BG)"))),
        ("summaryfile",str(exparser.parse("$steering(summaryfile_BG)"))),
        ("IPModuleURL",str(exparser.parse("$steering(SCRIPTS::repository)/simulation/generators.py"))),
        ("gcdfile",str(exparser.parse("$steering(gcdfile_11)"))),
        ("IPModuleClass","generators.CorsikaBG"),
        ("IPModuleDependencies",[str(exparser.parse("$steering(SCRIPTS::repository)/simulation/dcorsika.py")),
                            str(exparser.parse("$steering(SCRIPTS::repository)/simulation/polygonato.py"))]),
    )
    pre.SetParser("generate_BG",exparser)

Here it loads from the SCRIPTS::repository (webdata's copy of svn most likely).

Newer datasets (icesim4 especially) have module scripts baked in::

    # Configure PreTray modules 
    pre.AddModule("icecube.simprod.generators.CorsikaUCR","generate_BG")(
        ("RunCorsika",True),
        ("mjd",int(exparser.parse("$steering(mjd_11)"))),
        ("seed",int(exparser.parse("$steering(seed)"))),
        ("procnum",int(exparser.parse("$args(procnum)"))),
        ("nproc",int(exparser.parse("$args(nproc)"))),
        ("nshowers",int(exparser.parse("$eval($steering(CORSIKA::showers) * 2)"))),
        ("outputfile",str(exparser.parse("$steering(current_file_BG)"))),
        ("corsikaVersion","v6960-5comp"),
        ("summaryfile",str(exparser.parse("$steering(summaryfile_BG)"))),
        ("gcdfile",str(exparser.parse("$steering(gcdfile_11)"))),
    )
    pre.SetParser("generate_BG",exparser)
    
Here it loads from the metaproject copy of simprod-scripts.  This is nice for testing simprod-scripts itself.

More Help
~~~~~~~~~

Feel free to email simprod@icecube.wisc.edu for more information and assistance.  We're here to help you.
