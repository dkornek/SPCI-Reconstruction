# SPCI-Reconstruction

![Version](https://img.shields.io/github/tag/dkornek/SPCI-Reconstruction.svg?style=flat-square) ![License](https://img.shields.io/github/license/dkornek/SPCI-Reconstruction.svg?style=flat-square)

Tool for image reconstruction of emission densities measured with *Single Plane Compton Cameras* (SPCC). Typical applications are nuclear medicine imaging or range verification in hadron therapy.

## Image reconstruction
SPCC combine scatter and absorption detectors into one plane. By means of electronic collimation, the emission density can be backprojected. The algorithms used are:
* [maximum-likelihood expectation-maximization](https://www.ncbi.nlm.nih.gov/pubmed/18238264)
* [origin ensemble](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2590772/)

---

## Getting Started
Use these instructions to get a copy of this repo on your local machine for research or development.

### Prerequisites
You need:
* [**ROOT**](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html)
* [**boost**](https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html)


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
* You need a .root file containing the measured spectra for each considered detector. This file will be used to backproject the emission density.
> Both the choice of considered detectors as well as the binning pattern of the spectra must be consistent. The binning pattern can be changed with the [RebinningMacro](macros/RebinningMacro.cpp).

---

## Authors
* **Dominik Kornek** -- [TU Dresden](https://tu-dresden.de/)

---

## License
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.
