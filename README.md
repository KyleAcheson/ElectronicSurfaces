# ElectronicSurfaces

Calculates potential energy surfaces for any tri-atomic molecule given a grid a points, can also work for sections of surfaces for larger molecules (nosymm).
Can be run using no symmetry or with symmetry on - maximum symmetry is Cs. 
Currently works for Molpro. 

Runs on a local linux system (tested on Ubuntu) or any HPC system running either Sun Grid Engine or PBS Pro queueing systems.

Requires Python 3.5.2+ and Numpy.

### To run:

`python production_surfaces.py -inp input.json -g [geometry_file] &`

`[geometry_file]` = .npy or .mat file containing array of geometries over configuration space (in Angstrom)! 

Running on a HPC queueing system: `./production_surfaces.py -inp input.json -g geom.xyz -q [template filename]` 

For help: `python production_surfaces.py --help`

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

##### File System Structure:

Parent Directory (where script is executed)
    |
    |-> GP_XXXXX/
        |
        |-> SPE_GP_XXXXX.input
            SPE_GP_XXXXX.out
            SPE_GP_XXXXX.molden
            GRADS/
            |
            |-> GRAD_GP_XXXXX.input
                GRAD_GP_XXXXX.out
            NACMES/
            |
            |-> NACME_GP_XXXXX.input
                NACME_GP_XXXXX.out
            

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

- [x] Logic throughout whole programme to determine symmetry vs no symmetry calculations.

- [x] Routines to store output data in arrays - energies, gradients, NACMEs.

- [x] Parallel execution - Must force sharing of array memory for global data storage.

- [x] Routine to extract SOC.

- [x] Caclulations in CSF basis.

- [x] Add functionality to run on a HPC queuing system.

- [ ] Routine for calculating numerical gradients (MRCI).

- [ ] Routine to extract CASPT2 mixing matrix.

- [ ] Extend to Molcas - Framework already in place, just rewrite molpro.py in molcas format in molcas.py. 
