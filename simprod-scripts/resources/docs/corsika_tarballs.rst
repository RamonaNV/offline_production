CORSIKA Tarballs
================

Getting the source code
~~~~~~~~~~~~~~~~~~~~~~~
The CORSIKA code IceSim supports can currently (9/15/15) be obtained from:
http://code.icecube.wisc.edu/tools/distfiles/corsika/

Untar and unzip:

.. tar -zxvf [corsika-binary.tar.gz]

change directory to the newly create corsika installation dir:

.. cd [corsika]

The binary is created after corsika options have been configured.  To configure run the "coconut" executible in expert mode as follows:

.. ./coconut --expert

The expert flag allows you to select the date and time routine manually.  The default option is for the "coconut" script to detect it automatically and set it accordingly.  The problem is that the automatic detection does not always work.  When you get to the following screen:

.. Which routine for date and time ?
    1 - automatic detection by configure
      (only use other choices if this one fails) [DEFAULT]
    2 - new date_and_time routine
    3 - old date routine
    4 - timerc routine
    5 - date and time for IBM risc
    6 - old date routine for pgf77
..
    r - restart (reset all options to cached values)
    x - exit make
..
    (only one choice possible): 4
    SELECTED         : TIMERC 

select option 4 - timerc routine.  This date/time routine has been tested to work on cobalt.

If the date/time is routine is not the right one, at the very end, when a final binary is being compiled, CORSIKA will return an error such as:

.. corsika-corsikacompilefile.o: In function aamain:
.. corsika.F:1657: undefined reference to cpu_time


Naming the binary
~~~~~~~~~~~~~~~~~



Generating the config file
~~~~~~~~~~~~~~~~~~~~~~~~~~

