# ElectronicSurfaces

Calculates potential energy surfaces for any tri-atomic molecule given a grid a points, can also work for sections of surfaces for larger molecules (nosymm).
Can be run using no symmetry or with symmetry on - maximum symmetry is Cs. 
Currently works for Molpro. 

Ver2.0 - Splits state averaging of spin multiplicities.

### To run:

`./production_surfaces.py -inp input.json -g geom.xyz &`

For help: `./production_surfaces.py --help`

OR `python3` without shebangs.

Ensure inputs are specified in Angstroms - programme will convert to bohr.

If using symmetry, ensure the unique symmetry axis chosen is compatable with the plane you set the molecule to in xyz format. 

You should have extensively tested the levels of theory you wish to use for fitting surfaces prior to running.

### Features:

Single point energies - CASSCF, CASPT2, MRCI.

NACMEs (by finite differences) - CASSCF/ MRCI.

Gradients - CASSCF (Analytical). Note analytical gradients require the use of spin-adapted CSFs and thus for gradients state-averaing is split over each spin multiplicity.

Numerical MRCI gradients to be added in if required.

Spin-orbit Couplings - CASSCF/ MRCI.

CASPT2 couplings and gradients can be calculated by a transformation of CASSCF level NACMEs/ Gradients using the PT2 mixing matrix.

The progamme will save all output data to .npy files. The index for each grid point corrosponds to the index of that geometry in the initial list fed in.

If a calculation fails, elements of the output arrays for that grid point are set to NaN.

### TO DO:

- [x] Routine for setting up geometries at each point.

- [x] Routines for settuping up HF/CASSCF/CASPT2/MRCI blocks.

- [x] Routine for calculting single point energies (CASSCF/CASPT2/MRCI). 

- [x] Routine for calculating SOC (CASSCF/MRCI).

- [x] Routine for displacing molecule along x, y, and z coordinate by distance dr (depending on unique axis chosen).

- [x] Routine for calculating all possible combinations of coupled states.

- [x] Routine for incrementing record files.

- [x] Routine for calculating transition dipole matrices.

- [x] Routine for calculating NACMEs (CASSCF/ MRCI).

- [x] Routine for calculating analytical gradients (CASSCF).

- [ ] Routine for calculating numerical gradients (MRCI).

- [x] Logic throughout whole programme to determine symmetry vs no symmetry calculations.

- [x] Routines to store output data in arrays - energies, gradients, NACMEs.

- [x] Parallel execution - Must force sharing of array memory for global data storage.

- [ ] Routine to extract SOC.

- [ ] Internal coordinate routine for spherical coordinates. Option to use this or feed a predefined list of arrays to the programme.

- [ ] Routine to extract CASPT2 mixing matrix.

- [ ] Add functionality to run on a HPC queuing system.

- [ ] Extend to Molcas - Framework already in place, just rewrite molpro.py in molcas format in molcas.py. 
