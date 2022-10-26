# Notes on use of Garfieldpp

## Installation
Installation was relatively simple, with a few hiccups here or there installing the dependencies such as..
* ROOT
* GEANT4
* GSL
* CMake

The following PDF by Irina Kempf was useful in understanding how to compile and build Garfieldpp and all of the dependencies.
https://garfieldpp.web.cern.ch/getting-started/Garfield_Installation_Ubuntu_by_Irina_Kempf_20210428.pdf

\
Garfieldpp gitlab site: https://gitlab.cern.ch/garfield/garfieldpp \
Garfield Documentation: https://garfieldpp.web.cern.ch/garfieldpp/  

Useful things to note were making sure to **source** root and garfieldpp when starting a new terminal session so the following lines of code needed to be added to the *./bashrc* file:
```
source /$path-to-directory/root/bin/thisroot.sh
source /$path-to-directory/garfieldpp/install/share/Garfield/setupGarfield.sh
```
## Calculating Weighting Potential

I'm attempting to calculate the weighting potential and simulate the retrieval of a signal from a semiconductor solid with anode/cathode electrodes.

The following example documentation is a good start in learning how to use Garfieldpp to simulate the signal in a silicon sensor: https://garfieldpp.web.cern.ch/examples/silicon/
