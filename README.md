# SPCI-Reconstruction

![Version](https://img.shields.io/github/tag/dkornek/SPCI-Reconstruction.svg?style=flat-square) ![License](https://img.shields.io/github/license/dkornek/SPCI-Reconstruction.svg?style=flat-square)

Tool for backprojection of activity distributions measured with *Single Plane Compton Cameras* (SPCC) in high energy physics. Typical applications are molecular imaging or range verification in hadron therapy.

## Image reconstruction
SPCC consist of two distinct detector materials arranged in a checkerboard pattern in a plane. Energy spectra due to incoherent scattering can be obtained for any detector pair. These spectra are used to backproject the probable activity distribution. The algorithms used are:
* [maximum likelihood expectation maximization](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm)
* [origin ensemble](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2590772/) -- *yet to be implemented*

---

## Getting Started
Use these instructions to get a copy of this repo on your local machine for research or development.

### Prerequisites
You need:
* [**ROOT**](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html) -- data analysis framework
* [**boost**](https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html) -- C++ libraries


### Installing
To clone this repo, use: `$ git clone https://github.com/dkornek/SPCI-Reconstruction.git`

To compile your own binary file, type:
```
$ cd SPCI-Reconstruction/
$ mkdir build
$ cd build/
$ cmake ../
$ make
```

To execute the binary file, type:
```
$ ./SPCI-Reconstruction
```

---

## Usage
* You need a .root file containing the energy spectra for each considered detector for each voxel in an ascending order. This file will be used to generate the system matrix.
* You need a .root file containing the measured spectra for each considered detector. This file will be used to generate the 3D activity distribution with means of the system matrix.
> Both the choice of considered detectors as well as the binning pattern of the spectra must be consistent. The binning pattern can be changed with the [RebinningMacro](macros/RebinningMacro.cpp).

---

## Authors
* **Dominik Kornek** -- [TU Dresden](https://tu-dresden.de/)

---

## License
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.
