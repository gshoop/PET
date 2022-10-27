# Notes on use of Garfieldpp

## __TO-DO LIST__

- [ ] Dig into the `Medium` class and understand how it derives the fundamental material properties such as electron/hole velocities/mobilities, density of states, etc..


## Installation
Installation was relatively simple, with a few hiccups here or there installing the dependencies such as..
* [ROOT](https://root.cern.ch/)
* [GEANT4](https://geant4.web.cern.ch/)
* GSL
* CMake

The following PDF by Irina Kempf was useful in understanding how to compile and build Garfieldpp and all of the dependencies.
https://garfieldpp.web.cern.ch/getting-started/Garfield_Installation_Ubuntu_by_Irina_Kempf_20210428.pdf


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

## Transport properties

Getting the transport properties for electrons and holes is simple. We need a medium object which is defined in Garfield. In the example provided we use a Silicon object derived from the class `MediumSilicon`. 

## Notes and Understanding of the example


### **Future Considerations**

1. If Garfield is to be used we need to create our own `MediumCZT` class. In the example above we are using Si however there *is* a `MediumCdTe` class that could serve as a good starting point.
    * It seems that a `MediumCZT` class can be inherited and created on top of the `Medium` class provided by Garfield.