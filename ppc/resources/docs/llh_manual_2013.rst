*this is a llh/DirectFit manual written by Jacob Morrison circa 2013*

**Setup Process**

This is the process to setup the llh software. It shows you what directories and tools you need to run the llh software on your gpu. It also allows you to check and make sure you are properly setup.

1. To be able to run the llh software, you need the /gpu, /ice, /ocl, /llh, and /resources directories

2. After finding these directories in the svn directory of IceCube, place them where you will be working with them.

3. Double check src file in /gpu and /ocl directories - /gpu if you are using CUDA or /ocl if you are using opencl

4. For /llh/src, make sure I3\_PORTS points towards your local copy of IceCube tools – tools should have lib/gsl

5. Set environment variables with **. src** in /gpu (or /ocl) and /llh

6. Compile /gpu with **make go** and **make gpu** (or if using /ocl – **make oo** and **make ocl**)

7. Compile /llh with **make gpu** (or if using /ocl – **make ocl**)

8. Run **./ppc** in /gpu to make sure everything is running properly

Requires – [str] [om] [num] [device] (flasher)

[str] – string

[om] – DOM number

[num] – number of individual photons

[device] – GPU that ./ppc runs on

**./ppc [str] [om] [num] [device] (flasher)**

9. To run **./llh** on a specific GPU – edit src file in /gpu (or /ocl) directory and change CUDA\_VISIBLE\_DEVICES=”gpu,you,want”. Be sure to **. src** again after changing the src file

**Typical Command Sequence**

This is the typical sequence of commands I entered in my terminal to run the llh software. It also contains explanations as to what information is outputted with each command.

1. **PPCTABLESDIR=../resources/ice/ice model [FLSH=str,dom] [MLPD=0,1] [SREP] [DREP] ./llh [method=1]**

2. **PPCTABLESDIR=../resources/ice/ice model [FLSH=str,dom] [MLPD=0,1] [SREP] [DREP] ./llh [method=2]**

1. For the first command, you don’t want ./llh to read in the ini file, so change the name of the ini file. This first command runs llh for a specific flashing string and DOM with specifications of the ice model being used, whether the millipede function is used, and the number of simulated events compared to the events in the input file. The input file is a text file that is soft linked to the dat file in the /llh directory. The text file includes the string and DOM number of the receiving string and DOM, as well as the time, in nanoseconds, for the hit on the receiving DOM. Finally, it also includes the charge per hit in units of photoelectrons.

    This process gives both light yield (e.) and llh

2. echo command fills the ini file with x,y,z,th,ph,e,t,s,a (x,y,z,th,ph,t are all set to 0 and s,a are set to 1, e is equal to the value found in 1.), it then runs llh with the same string and DOM

This process gives only llh and reads in the ini file

Run for each of the ice models and each DOM file

Each DOM file must be soft-linked to dat file in /llh directory

**rm dat ; ln –s /where/data/file/lives dat**

**Ice Model Directories**

These are the four ice models that are used in the llh software. A short description of each ice model and where they all lie is provided in this section.

There are four important ice models (directories in /resources directory)

WHAM – ice model derived from dim (SPE) data

Spice 1 – sp1 – derived from bright (MPE) and old scattering function

Spice Mie – mie – MPE data with new scattering function

Spice Lea – lea – new scattering function, MPE data, anisotropy

**Files in Ice Model Directories**

These are the files that should be in each of the four ice model directories. A short description of each file is provided, as well as some helpful information.

as.dat – angular sensitivity – how the DOM is sensitive to the angle which the photon hits the DOM at (as.holeice=as.h2-50cm)

cfg.txt – configuration of the ice, unique to each ice model

icemodel.dat – scattering as a function of depth, unique to each ice model

icemodel.par – don’t all link to the same file

rnd.txt – random number files used with simulation

tilt.dat – ice layers aren’t flat, need to change in /gpu/llh.cxx from *#define T2LT* to *#define TILT*

tilt.par – see tilt.dat

wv.dat – response to the wavelength (nm)

eff-f2k, geo-f2k, hvs-f2k files and wv.rde must be in or copied into the four ice model directories (copied from /ice directory)

wv.rde file must be soft-linked to /resources/All/wv.rde once in ice model directory

(ln –s ../All/wv.rde wv.rde)

**Required Command Arguments**

This is the command and the various arguments that are needed for llh to run properly.

The command – **PPCTABLESDIR=../resources/ice/icemodel [FLDR] [FLSH] [MLPD] [SREP] [DREP] ./ llh [method]**

PPCTABLESDIR=../resources/ice/icemodel – This points to where the ice properties are stored that **./llh** uses when running.

[FLDR] = **FLDR=(angle in degrees)** – number of LEDs simulated and the azimuthal direction, left out of the command if using MLPD option or if there is no directional dependence

[FLSH] = **FLSH=str,dom** – the string and DOM that is being looked at

[MLPD] = **MLPD=0/1** – 0=off, 1=on - enables track reconstruction with detailed loss pattern (along the track) unfolding when FLSH is unset, or azimuthal unfolding for the flasher LEDs, when FLSH is set. In the latter case the FLDR will be ignored

[SREP] = **SREP=(number of repetitions)** – how many repetitions of simulation you want to compare

[DREP] = **DREP=(number of repetitions)** – how many events are in the input file

[method] = **(method number)** – argument after **./llh**

**0/1** – 2-step likelihood evaluation – reconstructs energy (light yield) and time, as well as calculating the likelihood

**2** – calculates the likelihood, tends to be used when running **./llh** reading in the ini file

**3** – reconstructs the energy and time

**-1** – same as **1** plus makes and output file

**-2** – same as **2** plus makes and output file

    **-3** – same as **3** plus makes and output file

**The ini file**

The llh software reads in this file if there is an ini file living in the /llh directory. If you do not want llh to read in the ini file, all you need to do is change the name of the file before running **./llh** and switch the new file back to ini if you want to run it with llh reading in the file. These are the variables that are listed in the ini file and what they stand for. The section also includes the variables printed out by **./llh** once it is done running.

llh – reduced log likelihood – the probability that the seen outcome or event was “caused” by the hypothesized process, only printed out by **./llh**, not a part of the ini file

x,y,z – position variables

th, ph – θ,φ – directional variables – θ=zenith angle, φ=azimuthal angle

e – light yield/energy – given in bunches of 2.5x10\ :sup:`9` photons

t – time

s,a – scale factors – set to one for our case

**Example**

Here is an example of how the llh software can be used the direction of a simulated flasher.

In this example, I ran a simulated flasher on DOM 45 of string 80 with 250 events and an azimuthal angle of ninety-three (93) degrees. To determine the angle using the llh software, I ran **PPCTABLESDIR=../ice FLSH=80,45 MLPD=0 SREP=1 DREP=250 ./llh 1**. I ran this without llh reading in the ini file to find both the light yield and the likelihood for this text file. Once I had a determined light yield or 0.385 for this file, I ran the same command, except with method 2, where only the likelihood is found. I had **./llh** read in the ini file with variables **0 0 0 0 0 0.385 0 1 1 (x y z th ph e t s a)**. I also added the FLDR argument to the command line and ran **./llh** for FLDR’s ranging from 0 to 360 in increments of 15, recording the likelihood for each increment with the set light yield. After collecting the likelihoods, I plotted a graph of likelihood versus FLDR and found the best likelihood to approximated 90 degrees.

**Troubleshooting**

Here is a list of common problems and solutions to those problems.

Be sure to source (**. src**) after each login in both /gpu and /llh.

Problem: *CUDA Error: setting the device when a process is active is not allowed*.

Description: The problem is that you are trying to run multiple GPU’s on a single CPU thread. To be able to do this you need a driver and runtime of at least 4.0 on your GPU machine.

Solution: Go to CUDA\_VISIBLE\_DEVICES in the /gpu/src file and switch to only a single GPU, instead of multiple GPUs.
