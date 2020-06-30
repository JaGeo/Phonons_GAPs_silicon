# Repository to reproduce "Combining phonon accuracy with high transferability in Gaussian approximation potential models"

The following Python packages have been used to produce the data: 
The first five can be installed via conda/pip, the other one has to be installed manually.

pymatgen 2019.12.22

Atomic  Simulation  Environment  (ASE) 3.19.0

phonopy 2.4.2

phono3py 1.18.2

quippy (including  the  GAP  code,  development  version  of  10  Jan 2020, see https://github.com/libAtoms/QUIP for help with the installation)

You have to download all potentials from zenodo.org (https://doi.org/10.5281/zenodo.3924402), unzip the files and adapt the filepaths in the two scripts of the folder accordingly. Please make sure that they are correct. Otherwise, the code won't work.

To test the the core class of the code, you have to add the FakeCalculator.py to core from the following repo: https://github.com/JaGeo/FakeCalculator/ (https://zenodo.org/record/3888134)

This code is licensed under a BSD 3-Clause License. 
