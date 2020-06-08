Detectors
=========
The detectors module contains a base class *ICBase* for detectors IC86, IC79, IC59.
Common parts factored out to reduce code redundancy.
Overwrite addparameter_X and segment_X in derived classes as needed
(X may be input or output, and must be main).

Like most IceProd Modules, it is basically a wrapper to :doc:`segements.DetectorSegment` that 

* Instatiantes and configures a Random Number Generator (RNG) service.
* Adds and configures I3Reader and I3Writer modules.
* Adds and configures Sanity checkers.
* Configures frame keys to omit.

The module contains the following subclasses for ICBase:

* ICBase (ipmodule.ParsingModule)
* IC86 (ICBase)
* IC79 (ICBase)
* IC59 (IC79)
* IceTop (ipmodule.ParsingModule)


