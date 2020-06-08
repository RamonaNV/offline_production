
.. _ppc-main:

Photon Propagation Code (ppc)
=============================

ppc is software program that can be used as either a stand-alone command-line executable, or within IceCube software chain as a module to generate and propagate photons from either in-situ light sources (flasher LEDs, standard candle lasers) or photons emitted or caused by high-energy particles or processes (Cherenkov photons and delayed luminescence).

Files
-----

- private/ppc

  - i3ppc.cxx

    ppc module wrapper, compiles into libppc.so

  - llh/

    llh/DirectFit code (stand-alone only)

  - ocl/

    the core of ppc, compiles into libxppc.so

  - ppc.h

    header file for code in ocl/, used by both i3ppc.cxx and llh/.

- resources/ice

  one of the ice models, fully configured for use by module, stand-alone program, or llh/DirectFit.


Configuration
-------------

ppc is conifigured via environmental variables, configuration files in the ice directory, and module parameters (or command-line parameters for stand-alone executable).

Environmental variables
+++++++++++++++++++++++

these are set as usual within a shell (with an "export" as necessary). Within a python script these are set with os.putenv(). If the ppc module is loaded with a call to load() these must be set prior to loading the module. If using an import statement instead the environmental variables can be set at any point prior to the Execute() statement.

- system-wide

  these are settings used by the GPU driver that are not directly accessed by ppc

  - **DISPLAY**

    *use: DISPLAY=:0*

  - **COMPUTE**

    *use: COMPUTE=:0*

  - **GPU_DEVICE_ORDINAL**

    *example: GPU_DEVICE_ORDINAL=0,2,3*

    sets which GPUs to use (and makes those not set invisible to the program). In the example above, assuming a system with 4 GPUS 0, 1, 2, and 3, sets to use 3 GPUs 0, 2, and 3, and makes GPU #1 invisible to the program. This re-numbers GPUs: in the example above the GPU numbers that the program sees will be 0, 1, and 2.

  - **GSL_RNG_SEED**

    *example: GSL_RNG_SEED=$RANDOM*

    sets the gsl random number generator seed, used by llh/DirectFit

- ppc

  - **PPCTABLESDIR**

    *example: PPCTABLESDIR=ice/*

    sets directory where ice and most other configuration files are located. By default uses current directory, ""

  - **WFLA**

    *example: WFLA=337*

    overrides the wavelength sampling curve conatined in wv.dat with a single wavelength in nm

  - **FLDR**

    *example: FLDR=-1*

    set FLDR=x+(n-1)*360, where 0<=x<360 and n>0 to simulate n LEDs in a symmetrical n-fold pattern, with first LED centered in the direction x. Negative or unset FLDR simulate a symmetric in azimuth pattern of light. Ignored in particle simulation.

  - **OFLA**

    *use: OFLA=0*

    setting this disables the default mode where photons that come back to the flashing DOM are omitted. 

  - **FWID**

    *example: FWID=9.7*

    sets the width (in degrees) of the 2d gaussian (von-Mieses-Fisher distribution) that determines the light emission profile of flasher LEDs. Set to -1 to simulate isotropcal emission profile

  - **NPHO/NPHO_X**

    *example: NPHO_2=512*

    sets the average number of photons to process in a single thread on a GPU. If underscore syntax is used, the number that follows the underscore sets the GPU to which to apply this setting. Setting this to 0 takes that GPU out of use. Default is 1024.

  - **XMLT/XMLT_X**

    *example: XMLT_1=4*

    oversubscribes to resouces advertised by the driver by the factor specified. Default is 1 for NVidia and 8 for AMD cards.

  - **OCPU/OGPU/OACC**

    *use: OGPU=1*

    use: only CPUs, only GPUs, or only Accelarator cards when specified. Can be combined to specify multiple devices. If not set the program will use every device available to it. In a system with GPUs it is recommended to use OGPU=1. Otherwise, if the driver advertises both GPUs and CPU, the load will be spread equally between both. GPUs, usually being faster, will complete their load quickly and then wait for the CPU device to complete, thus leading to idling GPUs and slower overall execution.

- llh/DirectFit

  - **CYLR**

    *use: CYLR=0/1*

    for cable simulation: (1) simulate straight cylindrical cable, which is faster or (0) curved gaussian-like shape of cable that curves around the DOM and asymptotically approaches the DOM axis above/below the DOM, which is slower (and is the default)

  - **ANGR**

    *example: ANGR="0 0 1"*

    sets nx ny nz components of the cascade/particle direction. At the same time the angular width of the proposal distribution is set to 0, so the direction is held fixed during iterations. This is overriden if the input file "ini" exists and is successfully read at initialization

  - **FSEP**

    *example: FSEP=1*

    llh only accepts DOMs that are more than FSEP DOMs away from flasher. Default is 1 (so if DOM 4 is flashing, DOMs 3 and 5 are not used)

  - **SREP**

    *example: SREP=10*

    simulate event this many times at each step of the calculation. Default is 1

  - **DREP**

    *example: DREP=250*

    the data file contained averages for this many events. Default is 1. Numbers above 1 are usually used only if there were multiple in-situ light source events taken with the same configuration (e.g. 250 flasher events)

  - **LOOP**

    *example: LOOP=1000*

    number of llh steps in a sub-chain. Different search methods might use this number differently. E.g., localized random search has this many evaluations. However, it is repeated 10 times in method 11 (usually used for cascade reconstruction). Default is 1000

  - **FAIL**

    *use: FAIL=0/1*

    set to cause program to fail on some warnings. Default is 0

  - **FAST**

    *use: FAST=0/1*

    1: only use time-integrated charges during simultaneous t0 (start time) and energy unfolding steps. 0: use time-binned charges in parts of the calculation when optimizing t0 and unfolding energy/flasher brightness. Default is 0

  - **MLPD**

    *use: MLPD=0/1*

    short for millipede. Enables/disables pattern unfolding: loss profile along the track (0 to reconstruct as a cascade, 1 to reconstruct as a track), or azimuthal flasher light emission profile (1 to enable unfolding into 2 up/down components and 72 azimuthal components spaced out 5 degrees apart; or 0 to use emission profile determined by the FLDR setting)

  - **FLSH**

    *example: FLSH=63,20*

    invokes flasher mode. Sets the flasher position to the value of the parameter

  - **FDUR**

    *example: FDUR=70*

    width of flasher emission pulse in ns assuming rectangular profile. Default is 70 ns

  - **QSAT**

    *example: QSAT=500*

    maximum integrated charge per DOM to accept that DOM into the calculation. Default is 500

  - **CNUM**

    *example: CNUM=40*

    number of cos(arrival angle wrt. PMT axis) bins. Default is 1

  - **LSIG**

    *example: LSIG=0.05*

    value of the sigma/model error to be used in the likelihood evaluation. 0 reverts to likelihood containing only Poisson terms. Default it 0.1 (i.e., 10%)

  - **FLOR**

    *example: FLOR=1*

    1: tilt the flasherboard in a direction consistent with the DOM tilt (only when MLPD=1). Default is 0

- inv (the code used to fit RDEs and to unfold angular sensitivity curve)

  - **SREP**
  - **DREP**

    these have the same meaning as when used with llh and described in the previous section

  - **IGEO**

    *example: IGEO=ice/geo-f2k*

    sets the geometry file

  - **IEFF**

    *example: IEFF=ice/eff-f2k*

    eff-f2k file used in the simulation

  - **IORI**

    *example: IORI=ice/eff-f2k.ori*

    sets the file specifying nominal RDE values (for use with XMAX and XSIG parameters described below)

  - **IANG**

    *example: IANG=ice/as.dat*

    sets the angular sensitivity file used in the simulation

  - **XINI**

    *example: XINI=xini*

    sets the file containing an appoximation to the unfolding result (usually a result from the previous interation)

  - **XMAX**

    *example: XMAX=1.5*

    sets the hard limits on RDEs around values contained in the IORI file. In the example above the limits are [1/1.5; 1.5] for a DOM with a nominal RDE value of 1. 0 disables the hard limits. Default is 0

  - **XSIG**

    *example: XSIG=0.01*

    adds regularization around the RDE values specified in the IORI file with width specified. Default is 0.1

Configuration files
+++++++++++++++++++

- ice (set by PPCTABLESDIR)

  - **as.dat**

    Definition of the angular sensitivity of a DOM. There are two possible variations in the format:

      - 1st number is greater than 0:

        the rest of the file contains coefficients of the polinomial expansion of the angular sensitivity curve vs. cos(photon arrival angle wrt. PMT axis). The first number is the maximum reached by this curve. Numbers lower than 1 accelerate calculation (since fewer photons need to be simulated)

      - 1st number is 0:

        This defines the "surface sensitivity" option. Second number in the file defines sensitive area: cos(angle to photon hit point on surface from center wrt. PMT axis) must be greater than this number to accept the photon

  - **cfg.txt**

    main configuration file. See example below for explanation.

    ::

      # ppc configuration file: follow strict order below
      5     # over-R: DOM radius "oversize" scaling factor
      1.0   # overall DOM efficiency correction
      0.35  # 0=HG; 1=SAM
      0.9   # g=<cos(theta)>

      130   # direction of major anisotropy axis
      -0.106 # magnitude of major anisotropy coefficient k1
      0.053  # magnitude of minor anisotropy coefficient k2

      0.5   # hole ice radius in units of [DOM radius]
      0.5   # hole ice effective scattering length [m]
      100   # hole ice absorption length [m]
      0.35  # hole ice 0=HG; 1=SAM
      0.9   # hole ice g=<cos(theta)>

    Last block or two last blocks can be optionally omitted (disabling anisotropy and hole ice parts of the calculation)

  - **cx.dat**

    DOM tilt map, each line contains: String#, OM#, nx, ny, nz, uncertainty (degrees). nx, ny, nz are components of the tilt vector that is defined as opposite of PMT axis direction

  - **dx.dat**

    Cable position map, each line contains: String#, OM#, azimuth direction to cable (degrees), uncertainty (degrees).

  - **eff-f2k**

    RDE (relative DOM efficiency) map, each line contains: String#, OM#, RDE, Type. If no entry RDE=1, Type=0 are assumed. DOMs that use corrected  wavelength acceptance from file wv.rde (for high-QE DOMs) have Type=1. It is possible to specify high-QE DOMs with Type=0 and simply a higher RDE value (nominally 1.35), of with an RDE value near 1 and Type=1. The acceptance correction curve parametrized in wv.rde file nears a value of 1.35 for wavelengths near 400 nm. RDE values taken from the GCD frame are matched with Type=0. If a corrected wavelength dependence is desired, GCD values need to be overridden by having this file (and wv.rde) present in the ice configuration directory

  - **geo-f2k**

    Geometry map, each line contains: DOM ID, Mainboard ID, x, y, z, String#, OM#. This file is necessary for running ppc from command-line. When present and running as an icetray module, will override the values from GCD

  - **hvs-f2k**

    High-voltage map, each line contains: String#, OM#, high voltage. Used only to specify that the DOM is on when HV>0. This file overrides the map of "ON" DOMs from GCD when present in the ice directory.

  - **icemodel.bbl**

    parametrization of air bubble contribution to scattering. Has 3 values: b, d1, d2. The parametrized contribution is b*(d1-d)*(d2-d) for d that specifies a shallower depth than both d1 and d2. The contribution is 0 otherwise (i.e. for deeper locations)

  - **icemodel.dat**

    main ice properties table: depth of the center of the layer, be(400), adust(400), delta tau (as defined in section 4 of the SPICE paper). All layers must be of equal width, and there must be at least 2 layers defined in the file. If the file icemodel.par contains 6 parameters, then the absorption coefficient is calculated as adust(400)=(D*[3rd element in a line]+E)*400-kappa.

    this file may contain 2 additional optional columns, containing the anisotropy coefficients k1 and k2. Ice layers defined with lines containing k1 and k2 will use these anisotropy coefficients instead of those specified in file cfg.txt

  - **icemodel.par**

    file with 4 parameters of the icemodel: alpha, kappa, A, B (as defined in section 4 of the SPICE paper). Each parameter is followed by its measurement uncertainty, which is ignored by the program. The older models (older than SPICE Lea or WHAM) have 6 parameters: alpha, kappa, A, B, D, E.

  - **rnd.txt**

    table of random number multipliers for the multiply-with-carry random number generator used by the parallelized kernel. Can have one or more elements per line, but only the first number is used (this is to make is copmatible with older formats of this file)

  - **tilt.dat**

    Describes ice layer tilt, each line contains: center depth of layer, and several depth corrections for locations specified in file tilt.par

  - **tilt.par**

    Containes a map of tabulated ice layer tilt locations, each line has: string number, and a relative distance along the gradient tilt direction (225 degrees SW)

  - **wv.dat**

    parametrization of wavelength-tabulated DOM acceptance (calculated from qe_dom2007a table of efficiency.h file of photonics), convolved with input spectrum. Each line contains: normalized integrated acceptance, and wavelength in nm.

  - **wv.rde**

    parametrization of the correction to the wavelength acceptance curve to be used for high-QE DOMs. Each line has: wavelength in nm, and correction factor (ratio of high-QE to nominal)

- llh/DirectFit additional configuration/input files, to be placed in the "current" directory

  - **as**

    this has the same format as as.dat in the ice directory. llh needs to be able to apply the angular sensitivity within its code when fitting for the DOM tilt or cable position. When using file "as" make sure to apply a uniform/flat sensitivity curve in file as.dat (e.g., by having it contain 2 numbers: 0.68 and 0.68)

  - **zs**

    contains the grid of search directions, used when fitting the DOM tilt (which is performed if this file is found). Each line contains: a unique identifying number, nx, ny, nz. This file can be generated with program "ico" in llh subdirectory

  - **cx**

    this has the same format as cx.dat in the ice directory. Make sure that only one of "cx", "cx.dat" is available at run time. If fitting for DOM tilt iteratively with llh/DirectFit, make sure that only "cx" is available.

  - **cs**

    contains the set of azimuthal positions of cable to test used when fitting for the cable position (which is performed if this file is found). Each line contains: a unique identifying number, and azimuth angle in degrees.

  - **dx**

    this has the same format as dx.dat in the ice directory. Make sure that only one of "dx", "dx.dat" is available at run time. If fitting for cable position iteratively with llh/DirectFit, make sure that only "dx" is available.

  - **bad**

    contains String#, OM# of DOMs that are to be considered bad in the fit. If this file is found the DOMs in it are excluded from the fit, and the no-hit contribution to the log likelihood sum is taken into account

  - **ert**

    contains String#, OM#, ti, tf that define the "DOM errata" list containing time intervals of bad data, which are not to be used. May define more than one interval [ti; tf) for each DOM

  - **dat**

    main data file. Each line contains: String#, OM#, time in ns, and average charge in p.e.s. The data is internally rebinned in 25 ns bins before applying the bayesian deblocking method to merge bins. If an event spans over more than 5000 ns then to avoid resizing the fixed 200 bin internal buffers the bin size is increased. It is recommended to trim events to keep then at 5000 ns or less in length but throwing away late pulses and coinsident events before or after the main event. Coincident events should be cleaned away anyway with, e.g., topological trigger. Longer events such as muons crossing the entire detector should of course not be shortened just to fit into 5000 ns, but only to remove afterpulses and coincident events.

  - **ini**

    Contains cascade/flasher parameters (to be used as best fit, or as initial approximation, or to facilitate iterations passing the solution between separate runs of llh). It may contain one or more lines, ordered as listed below:

    1) x, y, z (meters), zenith, azimuth (degrees), energy (GeV)/flasher brightness (bunches), time (ns), scattering and sbsorption scaling coefficients (last two unsupported, set both to 1.0). This line needs 5 or more elements to be accepted (some values at the end may be omittes, like the scaling coefficients)
    2) sequence numbers representing the unfolded pattern. The number of elements must match that expected by llh (usually defined by the geometry of the event) exactly, and the elements should sum up to 1. This line may be left empty if it is not needed but the following lines are.
    3) estimates of the proposal distribution parameters: rr (correlation between position and direction), dr (spacial width in m), da (angular width in degrees), di (intended for use with scattering and absorption scaling coeffifients, so should be left as 0), and optionally, lx (threshold value of likelihood, used in the ABC method for calculating uncertainties)

Module parameters
+++++++++++++++++

- ppc

  - **gpu**

    GPU to use. Default is -1, which uses all devices available to ppc

  - **fla**

    Flasher DOM (string, OM)

  - **nph**

    number of photons to simulate. A positive value enables in-situ light source simulation. Default value is 0

  - **wid**

    Flasher pulse width in ns. Default is 0

  - **MCTree**

    MCTree containing particles that ppc needs to process

  - **cyl**

    use cylinder (1) or strict +300 m (0) detector volume. Default is 1. For comparisons with clsim use 0

  - **keep**

    Keep events that don't produce hits (1). Otherwise (0, default) remove such events

  - **verbose**

    print some informational messages

  - **photons**

    save photons that cause hits into the MCTree if 1. Default is 0

  - **tau_dnde_vec**

    vector of pairs of luminescence decay time and dnde, tau in ns, dNdE in gamma/eV. Used in monopole simulation

  - **infoName**

    Name of the ppc info dictionary. Will not be created if set to empty string. Used in monopole simulation

Command-line parameters
+++++++++++++++++++++++

- ppc

  - no parameters

    Prints a summary of available tables, and an error if something is missing. If all necessary tables are found, also prints a summary of the available GPU devices within your system. These are numbered starting with 0 and must be specified with a [gpu] parameter in the examples below.

  - one parameter "-"
  - optionally "-" [x] [y]

    Print out the table of ice parameters (IceCube coordinate z of the center of the ice layer, absorption coefficient, and effective scattering coefficient) for wavelength w in [nm] (if set with WFLA=[w]) at the IceCube coordinates x and y in [m] (or 0, 0 if not specified). The parameters are computed using formulae of section 4 of the SPICE paper.

  - one integer parameter [gpu]

    Process particle simulation in f2k format from stdin. The muons must have been processed by mmc with the "-recc" option, which prints out all muon segments individually as "amu" particles. Here is an example of f2k input to ppc:

    ::

      #!/bin/awk -f

      BEGIN {
        print "V 2000.1.2"
        print "TBEGIN ? ? ?"

        srand(1);
        for(i=0; i<100; i++){
          x=(2*rand()-1)*500 # meters
          y=(2*rand()-1)*500 # meters
          z=(2*rand()-1)*500 # meters
          zenith=rand()*180  # degrees
          azimuth=rand()*360 # degrees
          l=500              # length, m
          energy=1.e5        # GeV
          t=0                # ns

          print "EM 1 1 1970 0 0 0"
          print "TR 1 0 e    ", x, y, z, zenith, azimuth, 0, energy, t
          print "TR 1 0 amu  ", x, y, z, zenith, azimuth, l, energy, t
          print "TR 1 0 hadr ", x, y, z, zenith, azimuth, 0, energy, t
          print "EE"
        }
        print "TEND ? ? ?"
        print "END"
     }


  - 4 parameters: [str] [dom] [num] [gpu]

    Simulate [num] photons emitted by a flasher or a standard candle at the position [str],[om]. Please note the following rules:

    - positive [str] simulates horizontal flashers, negative [-str] simulates tilted flashers,
    - str=0 and om=1,2 simulates standard candles 1 and 2,
    - you must set WFLA=337 before simulating the standard candles,
    - if the wv.dat file contains the flasher wavelength profile, WFLA=405 should be omitted,
    - if [num] is specified as x*y, x photons are simulated y times (with y trailing empty lines).

- llh/DirectFit

  - 1 parameter [method]

    - 0, 1 (same as 0): calculate llh. -1 additionally prints out the simulated hit data for the best solution
    - 10: applies localized random seach after calculating initial guess
    - 16: ABC (Approximate Bayesian Calculation) to estimate uncertainties



Description of output
---------------------

- ppc

  command-line ppc reports hits with "HIT" lines:

  HIT String# OM# time(ns) wavelength(nm) p_theta p_phi d_theta d_phi

  p_theta and p_phi specify direction of the photon at the point where it impacts the DOM

  d_theta and d_phi specify direction from the DOM center to the point of photon impact



- llh

  \*,?,|,-,+,0,1,... llh x y z zenith azimuth n t s a

  first element is a special character or an integer indicating the progress of the program as it goes through iterations. The second element, llh, is the saturated log likelihood (a measure of the goodness of fit), that indicates how well the simulation matches data (lower values are better). Coordinates x,y,z are in meters, zenith and azimuth in degrees. The next element, n, is either reconstructed particle deposited energy in GeV or light source brightness in photon bunches. This number maybe a constant factor off (effects like SPE mean, cable shadow, etc.) One way to figure out this factor is to reconstruct a few simulated events of known energy (i.e., calibrate the output of llh/DirectFit with the proper IceCube simulation). Next is t, the t0 (in ns) of the event. Finally, s and a are the scaling ice corrections. These are currently not used and are left at 1 each.

  If negative method is used as run-time parameter, the best match between data and simulation will be printed out in the following format:

  String# OM# bin_size charge_data(p.e.) charge_simulation(p.e.)

  internally llh/DirectFit bins the data nominally in 25 ns bins. It might be necessary to increase the bin size to a larger number if the overal length of the event is larger than 5 us (200 of the 25 ns bins). The bin size is printed out on stderr at initialization as "lbin". Then llh applies the Bayesian blocks procedure to merge some of the data bins. The number of initial bins contained in such a merged bin is indicated in the line above as bin_size. For each String#,OM# the bins are printed out in their time order, so it should be possible to infer the time structure (waveform) of the detected charge, albeit only with precision limited by the variable-size bins. This can be used to plot data and best simultion for visual inspection of the quality of fit.


Random notes on code structure
------------------------------

the point of next scatter is calculated by solving the following equation:

Exp( - integral_0^(distance) (scattering coefficient at x) dx ) = (random number)

since we have ice organized in layers of constant optical properties the integral reduces to a sum, and calculating the distance to next scatter is as simple as solving a linear equation with a couple boundary checks.


the point of absorption is calculated by solving the following equation:

Exp( - integral_0^(distance) (absorption coefficient at x) dx ) = (random number)

since we have ice organized in layers of constant optical properties the integral reduces to a sum, and calculating the distance to next scatter is as simple as solving a linear equation with a couple boundary checks. One complication compared to scattering is that the sum is done over multiple segments because of intermediate scatterings. So the code keeps subtracting the integral evaluated between successive scatters from the -log(random number) until it drops below zero. When that happens the particle does not make it to the next scatter point and the point of absorption is calculated instead.

on any segment the direction is fixed and absorption coefficient is modified according to the anisotropy model. It should be easy to do the same to the scattering coefficient, if necessary.


Specific code details requested in tickets
------------------------------------------

Cable shadow
++++++++++++

Cable shadow is implemented in ppc as a folloowing approximation: the photon landing coordinates (on a DOM) and final direction are used to "backtrack" the photon to sek whether it could have intersected with cable positioned next to the DOM before landing on the surface of the DOM. The cable shape is implemented as a vertical cylinder with a radius of 23 mm that is touching the surface of the DOM (for the nominally oriented DOM, at the equator). An implementation of the cable as a curve touching the DOM surface but asymptotically aligning with the DOM axis above and below the DOM, also exists within the llh/DirectFit code. This code can be switched on by setting CYLR=0. The location of the cable wrt. the DOM is specified via angle to the cable in the configuration file dx.dat.

The cable shadow code is implemented outside of the propagate kernel (which runs on the GPU), thus executing on the CPU side, during the post-processing of the hits (the code is in file f2k.cxx).

Hole Ice
++++++++

The hole ice is implemented to describe the following physical construction precisely: A vertical infinite cylinder column of constant ice properties, going through the center of the string (for each string). The ice properties are specified in the cfg.txt file in the optional block as described in the configuration files section above. The configuration describes the radius of the cylindrical column, scattering and absorption coefficients, and shape parameters of the scattering function: f_SL and g. If all of the DOMs of the string have exactly the same x and y coordinates, the hole ice column is simulated as concentric with the DOMs. In order to simulate the situation where the hole ice is to a side of the DOM, the DOM coordinates need to be adjusted. Keep in mind that this will in turn modify the average x and y of the string (i.e. of the center of the string), so the coordinates of the rest of the DOMs need to be adjusted in the opposite direction by a small amount (1/60th for a nominal IceCube string). Since this implementation is currently just a placeholder, still awaiting detailed calibration of the hole ice properties, a more verstile configuration has not yet been implemented. It may turn out that the future configurations fully implementing hole ice can be fully specified with the existing scheme, or that modification may be required.

DOM tilt
++++++++

DOM tilt is implemented by assuming a tilt in the DOM axis (i.e. deviation from the vertical) during application of the angular DOM sensitivity. This is done for either the "effective" angular DOM sensitivity that only depends on the direction of the photon (and then the angle on which the angular sensitivity depends is calculated wrt. the DOM tilt axis rather than the vertical), or the new DARD-style DOM sensitivity which accepts photons a certain distance down from the DOM equator (i.e., effectively simulating a sensitive survace of the PMT). The tilt directions of DOMs are specifies in the configuration file cx.dat.

The DOM tilt code is implemented outside of the propagate kernel (which runs on the GPU), thus executing on the CPU side, during the post-processing of the hits (the code is in file f2k.cxx).

DOM Oversize
++++++++++++

In order to accelerate calculation, the DOMs can be (and normally are) oversized. The latest ice models were fitted with the oversize factor of 5, and much of the nominal IceCube simulation also uses an oversize factor of 5. Simplifying the actually ipmlemented geometry (explained in the following paragraph) a little, this means the DOM radius is increased by the factor of 5, increasing the DOM cross-sectional area presented to the oncoming photons by a square of the oversize factor, i.e., by a factor of 25. This, in turn, means that a factor 25 times fewer photons need to be generated and propagated, thus accelerating the simulation by a factor of 25.

A number of variations of the oversizing geometries and simulation strategies were studied some time after the feature was introduced, and the best option was settled on during the development of the SPICE Mie ice model (and used mostly unchanged since then). The main issue with scaling the DOM as a perfect sphere was that the DOM was occupying much more space, 125 time more (a cube of the oversize factor). This space was then unavailable to the scattering and absorption processes, significantly changing the timing distributions (as well as time-integrated charge). Thus, another approach, the "pancake" scheme, was implemented instead. In the pancake construction the DOM dimensions are scaled only in the directions perpendicular to the traveling photon, while maintaining the nominal dimension in the direction along the photon travel. This maintains the factor 25 scaling in the DOM cross section presented to the photon, while also reducing the volume occupied by the oversized DOM from a factor 125 to only 25 times the nominal volume. As the photon scatters and changes direction, so does the pancake rotate so that the area presented to the photon is always 25 times the nominal. Such changes in simulated DOM geometries, as well as the larger dimensions of the DOMs (compared to nominal) do lead to some deviations in the timing distributions, and for oversize factor of 16 (used in the development of ice models up to and including SPICE Lea) lead to about 1-2 ns distortion of the timing distribution at 17 m, and 3-9 ns distortion at 125 m (as measured in the positions of the leading edge and peak in the clear deep ice). SPICE 3.x ice models used the entire volume of IceCube detector in the ice model fits, including the more dense part DeepCore. It was thought that the oversize factor 16 was too big to adequately approximate the physics of the denser parts of the array in cleaner ice, so a smaller factor of 5 (also matching the factor used in most of the simluation produced at the time) was settled on and used in SPICE 3.x models.

Additional feature in the strategy of the oversize implementation is to continue propagating the photon after it hits the DOM. This is disabled when the oversize factor is set to 1, and the photon is stopped and disappears once it hits the surface of the DOM. Such a strategy is thought to be more correct, as in the equivalent nominal-size DOM treatement, while propgating a bunch of 25 photons that would otherwise hit theoversized DOM, only one of them will hit the nominal-size DOM, while 24 will continue unimpeded. This strategy was also included into the timing distribution distortion numbers stated in the previous paragraph.



References:
-----------

SPICE models
++++++++++++

`Measurement of South Pole ice transparency with the IceCube LED calibration system (SPICE Paper): arXiv:1301.5361 <http://icecube.wisc.edu/~dima/work/WISC/ppc/spice/new/paper/a.pdf>`_

`Evidence of optical anisotropy of the South Pole ice (Ice Anisotropy Paper): arXiv:1309.7010 (pp. 17-20) <http://icecube.wisc.edu/~dima/work/WISC/papers/2013_ICRC/ice/icrc2013-0580.pdf>`_

llh/DirectFit
+++++++++++++

`Event reconstruction in IceCube based on direct event re-simulation (DirectFit paper): arXiv:1309.7010 (pp. 21-24) <http://icecube.wisc.edu/~dima/work/WISC/papers/2013_ICRC/dir/icrc2013-0581.pdf>`_

`Likelihood description for comparing data with simulation of limited statistics (Likelihood Paper): arXiv:1304.0735 <http://icecube.wisc.edu/~dima/work/WISC/papers/2013/llh/a.pdf>`_

`Likelihood description for comparing data to simulation of limited statistics (LLH ICRC Paper) <http://icecube.wisc.edu/~dima/work/WISC/papers/2013_ICRC/llh/icrc2013-0582.pdf>`_

PPC
+++

`Photon tracking with GPUs in IceCube, Nuclear Inst. and Methods in Physics Research, A, Volume 725, pp. 141-143. <http://icecube.wisc.edu/~dima/work/BKP/DCS/VLVNT11/paper/ppc.pdf>`_

`Photon Propagation with GPUs in IceCube, Proceedins to GPUHEP2014, DESY-PROC-2014-05, pp. 217-220 <http://icecube.wisc.edu/~dima/work/WISC/new/2014/gpu2014/chirkin_dmitry.pdf>`_

Miscellaneous
+++++++++++++

`Older ppc information pages on my homepage <http://icecube.wisc.edu/~dima/work/WISC/ppc/>`_ and `readme file <http://icecube.wisc.edu/~dima/work/WISC/ppc/readme.html>`_

`AMANDA file format definition (used as input for command-line ppc in particle mode) <https://www-zeuthen.desy.de/~steffenp/f2000/>`_

`Muon Monte Carlo (MMC), a java program that can be used to process muons for use with command-line ppc <http://icecube.wisc.edu/~dima/work/MUONPR/>`_
