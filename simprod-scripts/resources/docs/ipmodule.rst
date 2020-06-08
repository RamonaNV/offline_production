.. highlight:: python

IceProd Module Base
===================

Basic IceProd Modules
---------------------

IceProd Modules are basic wrappers around tray segments that provide an interface for IceProd.
The interface of an IceProd module or IPModule consist of

* An __init__ method where one should initialize and define paramters
* A method to define input parameters to be passed by IceProd (Similar to I3Module)
* An Execute method that takes as input a Python dictionary for storing performance statistics.

An IceProd can basically do anything you want but often the purpose is to run IceTray code. The typical 
structure of the Execute method is the following:

* Instantiate I3Tray
* Add common services (RNG, SummaryService)
* Add I3Reader module (unless this is a generator)
* Add/configure tray segments
* Add I3Writer module
* Terminate chain (Trash,Execute,Finish)

.. code-block:: python

   from icecube.simprod import ipmodule
   class MyModule(ipmodule.IPModule):
         def __init__(self):
             # initialize parent class
             ipmodule.IPModule.__init__(self)

             # define configurable parameters
             self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
             self.AddParameter('outputfile','Output filename','')

         def Execute(self,stats):
             # Base class defines an "Excute" option which tells 
             # the module wether to run or not
             if not ipmodule.IPModule.Execute(self,stats): return 0

             # Load libraries
             from I3Tray import *
             from icecube.simprod.segments import PropagateMuons,GenerateNeutrinos

             # Instantiate a tray
             tray = I3Tray()

             # Add services
             randomService = phys_services.I3SPRNGRandomService(
                               self.rngseed, 
                               self.rngnumberofstreams, 
                               self.rngstream)

             tray.context['I3RandomService'] = randomService

             # Call tray segments
             tray.AddSegment(GenerateNeutrinos, 
                             'generator',
                             RandomService=randomService,
                             NumEvents=self.nevents,
                             SimMode=self.simmode)

             # Add IO modules
             tray.AddModule("I3Writer","writer", filename=self.outputfile)

             # Run IceTray
             tray.Execute()
             



Intatiating the IPModule
~~~~~~~~~~~~~~~~~~~~~~~~
Executing the IPModule then can be done the following way::

   #!/bin/env python
   from icecube.simprod.modules import MyModule
   if __name__ == '__main__':
      stats = {}
      mymod = MyModule()
      mymod.SetParameter("gcdfile","GCD.i3.gz")
      mymod.SetParameter("outputfile","test.i3.bz2")
      mymod.Execute(stats)
      print stats
      
ParsingModule: IceProd Module + OptionParser
--------------------------------------------
The class :py:class:`.ipmodule.ParsingModule` inherits from both :py:class:`~.ipmodule.IPModule` and :py:class:`~.optparse.OptionParser`. 
This provides an interface 
that can be run by humans (through command-line options) as well as by IceProd (through the IPModule interface).
In order to implement a ParsingModule, you need to replace the parent base class as well as the invocation of __init__
and Execute::

   class MyOptParsingModule(ipmodule.ParsingModule):

         def __init__(self):
             ipmodule.ParsingModule.__init__(self)
             ...

         def Execute(self,stats):
             if not ipmodule.ParsingModule.Execute(self,stats): return 0


And then call :py:meth:`ExecuteOpts` instead of Execute in your instantiated class::

   #!/bin/env python
   from icecube.simprod.modules import MyModule
   if __name__ == '__main__':
      stats = {}
      myopmod = MyModule()
      myopmod.ExecuteOpts(stats)
      print stats


Since :py:class:`.ipmodule.ParsingModule` inherits from both :py:class:`~.ipmodule.IPModule` and :py:class:`~.optparse.OptionParser`. You can call your script with a 
:py:class:`ParsingModule` object on the command line::

   Usage: nugen.py [options]
   Options:
      -h, --help show this help message and exit
      --execute boolean condition to execute
      --gcdfile=GCDFILE GeoCalibDetStatus filename
      --outputfile=OUTPUTFILE
      Output filename
      --summaryfile=SUMMARYFILE JSON Summary filename
      --mjd=MJD MJD for the GCD file
      --RunId=RUNID Configure run ID
      --RNGSeed=RNGSEED RNG seed
      --RNGStream=RNGSTREAM

