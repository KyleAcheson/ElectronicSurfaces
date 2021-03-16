#!/usr/bin/python3
import os
import re
import sys
import multiprocessing
import numpy as np
import argparse
import textwrap
import json
import scipy.io as sio
import setup
import utilities as util
import molpro
import molcas

ANG2AU = 1.88973

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-inp", dest="inputfile", required=True, type=str,
                    help=textwrap.dedent('''\
                    Input options in Json format\n
                    Required Input Specifications:\n
                    "code": molpro/molcas
                    "hpc" : HPC system - 'local', 'Sun Grid Engine' or 'PBS'
                    "mem": file memory
                    "symm": Symmetry generator (X, Y, Z) or  'nosymm' - max Cs symmetry
                    "basis": basis set in chosen programme format
                    "soc": yes/no
                    "nacme": yes/no
                    "nacme_level": casscf/mrci
                    "dr": displacement value (ANGSTROM) if nacmes are numerical
                    "paxis": principle axis
                    "grad": yes/no - gradient calculation
                    "spe": casscf/mrci/caspt2 - single point energy calculation
                    "atoms": list of atoms (list of str)
                    "nelec": number electrons (int)
                    "singlets": yes/no
                    "triplets": yes/no
                    "occ": nlists of occupied orbitals - [[Irrep1, Irrep2], ... ]
                    "closed": nlists of closed orbitals - [[Irrep1, Irrep2], ... ]
                    "states": nists of no. states - [[Irrep1(S=0), Irrep2(S=0), Irrep1(S=2), Irrep2(S=2)], ... ]\n
                    NB:
                    Can provide n lists (n = 1, 2, 3...) of orbitals and states
                    to project through active spaces at CASSCF level.
                    If a MR method is chosen, the MR calculation will corrospond
                    to the last list of active spaces and states.\n
                    If running with  symmetry off, set the orbitals and states
                    that corrospond to Irrep2 to 0.\n
                    If running singlet or triplet states only, set the multiplicity
                    you wish to ommit to 0.\n'''))
parser.add_argument("-g", dest="geom_inputfile", required=True, type=str,
                    help=textwrap.dedent('''\
                    Nd Array of Geometry in xyz coordinates (ANGSTROM) over all space\n
                    Must be of type .mat or .npy\n
                    Must have dimensions (Natom, 3, Ngrid DoF1, Ngrid DoF2, ..., Ngrid DoFN)\n
                    Degrees of freedom will be flattened from left to right to give a list of geometries\n'''))
parser.add_argument("-q", dest="hpc_script", required=False, type=str,
                    help="Submission script template for running on a HPC queueing system.")

args = parser.parse_args()
inp = args.inputfile
input_geom = args.geom_inputfile
submit_script = args.hpc_script

# Load inputs
f = open(inp, 'r')
inputs = json.load(f)
f.close()

file_type = input_geom.split('.')[-1]
if file_type == 'npy':
    referance_geom_raw = np.load(input_geom)
elif file_type == 'mat':
    referance_mat_raw = sio.loadmat(input_geom)
    keys = [k for k in referance_mat_raw if '__' not in k]
    if len(keys) > 1:
        raise TypeError('.Mat file must have only one variable')
    else:
        key_val = keys[0]
    referance_geom_raw = referance_mat_raw[key_val]
else:
    raise TypeError('Input goemetry file must be .npy or .mat')


try:
    f = open(submit_script, 'r')
    submit_script_lines = f.readlines()
    f.close()
    submission_script = ''.join(submit_script_lines)
except TypeError:
    submission_script = None


##################
#  Sanity Checks #
##################


required_inputs = ["code", "hpc", "mem", "symm", "basis", "soc", "dr", "paxis",
                  "nacme", "nacme_level", "grad", "spe", "nelec", "occ",
                  "atoms", "closed", "states", "singlets", "triplets"]

code_except = ['molpro', 'molcas']
yn_except = ['yes', 'no']  # For all y/n inputs
symm_except = ['x', 'y', 'z', 'nosymm']
nacme_except = ['casscf', 'mrci']
paxis_except = ['x', 'y', 'z']
except_spe = ['casscf', 'mrci', 'caspt2']
except_hpc = ['local', 'Sun Grid Engine', 'PBS']

InputError = ("ERROR: Missing required input fields from input.json.\n"
              "RUN: %s --help for a list of fields." % sys.argv[0])


class InputGeomError(Exception):
    def __init__(self, geom, message='''Input geometry array has incorrect dimensions
                                       or does not match number of user specified atoms'''):
        self.geom = geom
        self.message = message
        super().__init__(self.message)


# Make sure all keys are correct case

inputs = {k.lower(): v for k, v in inputs.items()}

# Check all required input params provided

for k in required_inputs:
    if k not in inputs.keys():
        raise ValueError(InputError)

# Check input params have correct data type

if not isinstance(inputs['code'], str):
    raise ValueError(" '%s' must be type str" % inputs['code'])
if not isinstance(inputs['hpc'], str):
    raise ValueError(" '%s' must be type str" % inputs['hpc'])
if not isinstance(inputs['mem'], str):
    raise ValueError(" '%s' must be type str" % inputs['mem'])
if not isinstance(inputs['symm'], str):
    raise ValueError(" '%s' must be type str" % inputs['symm'])
if not isinstance(inputs['singlets'], str):
    raise ValueError(" '%s' must be type str" % inputs['singlets'])
if not isinstance(inputs['triplets'], str):
    raise ValueError(" '%s' must be type str" % inputs['triplets'])
if not isinstance(inputs['basis'], str):
    raise ValueError(" '%s' must be type str" % inputs['basis'])
if not isinstance(inputs['soc'], str):
    raise ValueError(" '%s' must be type str" % inputs['soc'])
if not isinstance(inputs['nacme'], str):
    raise ValueError(" '%s' must be type str" % inputs['nacme'])
if not isinstance(inputs['nacme_level'], str):
    raise ValueError(" '%s' must be type str" % inputs['nacme_level'])
if not isinstance(inputs['grad'], str):
    raise ValueError(" '%s' must be type str" % inputs['grad'])
if not isinstance(inputs['spe'], str):
    raise ValueError(" '%s' must be type str" % inputs['spe'])
if not isinstance(inputs['paxis'], str):
    raise ValueError(" '%s' must be type str" % inputs['paxis'])
if not isinstance(inputs['atoms'], list):
    raise ValueError(" '%s' must be type list of str" % inputs['atoms'])
if not isinstance(inputs['nelec'], int):
    raise ValueError(" '%s' must be type int" % inputs['nelec'])
if not isinstance(inputs['dr'], float):
    raise ValueError(" '%s' must be type float" % inputs['dr'])
if not isinstance(inputs['occ'], list):
    raise ValueError(" '%s' must be type list of lists of int" % inputs['occ'])
if not isinstance(inputs['closed'], list):
    raise ValueError(" '%s' must be type list of lists of int" % inputs['closed'])
if not isinstance(inputs['states'], list):
    raise ValueError(" '%s' must be type list of lists of int" % inputs['states'])

# Ensure values of type str are lowercase

for k, v in inputs.items():
    if isinstance(v, str):
        inputs[k] = v.lower()
inputs['dr'] = inputs['dr']*ANG2AU

# Check inputs have acceptable values

if inputs['code'] not in code_except:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['code'], code_except))
if inputs['hpc'] not in except_hpc:
    raise ValueError(" choose from %s." % except_hpc)
if inputs['symm'] not in symm_except:
    raise ValueError(" choose from %s." % (symm_except))
if inputs['soc'] not in yn_except:
    raise ValueError(" choose from %s." % (yn_except))
if inputs['nacme'] not in yn_except:
    raise ValueError(" choose from %s." % (yn_except))
if inputs['grad'] not in yn_except:
    raise ValueError(" choose from %s." % (yn_except))
if inputs['singlets'] not in yn_except:
    raise ValueError(" choose from %s." % (yn_except))
if inputs['triplets'] not in yn_except:
    raise ValueError(" choose from %s." % (yn_except))
if inputs['nacme_level'] not in nacme_except:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['nacme_level'], nacme_except))
if inputs['paxis'] not in paxis_except:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['paxis'], paxis_except))
if inputs['spe'] not in except_spe:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['spe'], except_spe))


# Check input state and orbitals are correct

if inputs['singlets'] == 'yes':
    for nstates in inputs['states']:
        if inputs['symm'] != 'nosymm':
            if not isinstance(nstates[0], int) or not isinstance(nstates[1], int):
                raise ValueError('Please set number of singlet states')
        elif inputs['symm'] == 'nosymm':
            if not isinstance(nstates[0], int) or nstates[1] != 0:
                raise ValueError('Please set number of singlet states for irrep 1 to an integer and irrep 2 to 0.')
elif inputs['singlets'] == 'no':
    for nstates in inputs['states']:
        if nstates[0] != 0 and nstates[1] != 0:
            raise ValueError('Singlets input set to no, please set singlet states for both irreps to 0.')

if inputs['triplets'] == 'yes':
    for nstates in inputs['states']:
        if inputs['symm'] != 'nosymm':
            if not isinstance(nstates[2], int) or not isinstance(nstates[3], int):
                raise ValueError('Please set number of triplet states')
        elif inputs['symm'] == 'nosymm':
            if not isinstance(nstates[2], int) or nstates[3] != 0:
                raise ValueError('Please set number of triplet states for irrep 1 to an integer and irrep 2 to 0.')
elif inputs['triplets'] == 'no':
    for nstates in inputs['states']:
        if nstates[2] != 0 and nstates[3] != 0:
            raise ValueError('Triplets input set to no, please set triplet states for both irreps to 0.')

for occ_orb in inputs['occ']:
    if inputs['symm'] != 'nosymm':
        if not isinstance(occ_orb[0], int) or not isinstance(occ_orb[1], int):
            raise ValueError('Ensure occupied orbitals are integers for both irreps.')
    elif inputs['symm'] == 'nosymm':
        if not isinstance(occ_orb[0], int) or occ_orb[1] != 0:
            raise ValueError('Ensure occupied orbitals of irrep 1 is an integer and irrep 2 is 0.')

for closed_orb in inputs['closed']:
    if inputs['symm'] == 'yes':
        if not isinstance(closed_orb[0], int) or not isinstance(closed_orb[1], int):
            raise ValueError('Ensure closed orbitals are integers for both irreps.')
    elif inputs['symm'] == 'nosymm':
        if not isinstance(closed_orb[0], int) or closed_orb[1] != 0:
            raise ValueError('Ensure closed orbitals of irrep 1 is an integer and irrep 2 is 0.')


# Ensure input geometry and HPC queue submission script is correct

if referance_geom_raw.shape[0:2] != (len(inputs['atoms']), 3):
    raise InputGeomError(referance_geom_raw)
if submit_script is None and inputs['hpc'] != 'local':
    raise ValueError(" If not running via. a queueing system set hpc input to local")
elif submit_script is not None and inputs['hpc'] == 'none':
    raise ValueError(" If running via. a queueing system set hpc input to 'local', 'Sun Grid Engine' or 'PBS'")

if submission_script is not None and not re.search(r'template', submission_script):
    raise ValueError(" Place the keyword 'template' where your input file goes in the HPC submission script")


# While waiting to add new functionality

if inputs['spe'] == 'mrci' and inputs['grad'] == 'yes':
    raise ValueError(" Gradients for MRCI are yet to be added")
if inputs['code'] == 'molcas':
    raise ValueError(" Molcas is not yet available")


#########################
# Quantum Code KEYWORDS #
#########################

molpro_keys = {'cas_prog': '1PROGRAM * MULTI',
              'nacme_regex': '!Total NACME.*$',
              'termination_code': 'Variable memory released',
              'submit_command': 'molpro',
              'output_ext': 'out'}

molcas_keys = {}  # TO DO


if inputs['code'] == 'molpro':  # Set input dependent keys
    if inputs['spe'] == 'casscf':
        molpro_keys['energy_regex'] = "MCSCF STATE [0-9]{1,2}\.\d Energy.*$"
        molpro_keys['grad_regex'] = 'RSPT2 GRADIENT FOR STATE'
        molpro_keys['numerical'] = False
    elif inputs['spe'] == 'mrci':
        molpro_keys['energy_regex'] = "MRCI STATE [0-9]{1,2}\.\d Energy.*$"
        molpro_keys['grad_regex'] = 'Numerical gradient'
        molpro_keys['numerical'] = True
    elif inputs['spe'] == 'caspt2':
        molpro_keys['energy_regex'] = "RS2 STATE [0-9]{1,2}\.\d Energy.*$"
        molpro_keys['grad_regex'] = 'RSPT2 GRADIENT FOR STATE'
        molpro_keys['numerical'] = False


elif inputs['code'] == 'molcas':
    if inputs['spe'] == 'casscf':
        pass
    elif inputs['spe'] == 'mrci':
        pass
    elif inputs['spe'] == 'caspt2':
        pass


############################
# Coordinate Reading       #
############################


#def coordinateGenerator(referance_geom_raw):
#    list_geom = []
#    for i in np.arange(1.40, 1.85, 0.05):
#        new_geom = referance_geom_raw.copy()
#        new_geom[1, 2] = i
#        new_geom = new_geom*ANG2AU
#        list_geom.append(new_geom)
#    return list_geom


def coordinate_reader(raw_geom):
    geom_array = raw_geom*ANG2AU
    list_geom = []
    dims = geom_array.shape[2:]
    ngrid = int(np.product(dims))
    geom_flat = geom_array.reshape(len(inputs['atoms']), 3, ngrid)
    for i in range(ngrid):
        list_geom.append(geom_flat[:, :, i])
    return list_geom

##################################
# MAIN - Run surface calculation #
##################################


def run(grid_point):
    """Takes each geometry and an index for each grid point and runs the series
       of calculations specified by the user for either quantum chemistry code
       in parallel. Data for each grid point is sorted in global arrays according
       to the index."""
    index, geom = grid_point
    parent = os.getcwd()
    checkfile = "%s/Check_GP%s.chk" % (parent, index)
    setup.logging(checkfile, message='Calculation started!')
    logfile = "%s/Err_GP%s.log" % (parent, index)

    if inputs['code'] == 'molpro':

        [spe_wdir, spe_inputfile] = molpro.setup_spe(inputs, geom, pwd, index)
        [NormTerm, spe_output] = util.run_calculation(inputs['hpc'], molpro_keys, pwd, spe_wdir, spe_inputfile, submission_script, index)
        if NormTerm:
            data.extract_energies(spe_wdir+spe_output, inputs['spe'], molpro_keys['energy_regex'],
                                  molpro_keys['cas_prog'], index)
            if inputs['soc'] == 'yes':
                data.extract_soc(spe_wdir+spe_output, index)

            if inputs['nacme'] == 'yes':
                [nacme_wdir, nacme_inputfile, displaced_axes] = molpro.main_nacme(inputs, geom, spe_wdir, index)
                [NacmeNormTerm, nacme_output] = util.run_calculation(inputs['hpc'], molpro_keys, pwd, nacme_wdir, nacme_inputfile, submission_script, index)
                if NacmeNormTerm:
                    data.extract_nacme(nacme_wdir+nacme_output, molpro_keys['nacme_regex'], index, displaced_axes)
                else:
                    setup.logging(logfile, message='NACME FAILED')
                    data.nacmes[:, :, :, index] = 'NaN'

            if inputs['grad'] == 'yes':
                [grad_wdir, grad_inputfile] = molpro.setup_gradient(inputs, geom, spe_wdir, index)
                [GradNormTerm, grad_output] = util.run_calculation(inputs['hpc'], molpro_keys, pwd, grad_wdir, grad_inputfile, submission_script, index)
                if GradNormTerm:
                    data.extract_grad_molpro(grad_wdir+grad_output, molpro_keys['grad_regex'], molpro_keys['numerical'], index)
                else:
                    setup.logging(logfile, message='GRAD FAILED')
                    data.grads[:, :, :, index] = 'NaN'

        else:
            data.energies[index, :] = 'NaN'
            setup.logging(logfile, message='SPE FAILED')
            if inputs['nacme'] == 'yes':
                data.nacmes[:, :, :, index] = 'NaN'
            if inputs['grad'] == 'yes':
                data.grads[:, :, :, index] = 'NaN'
            if inputs['soc'] == 'yes':
                data.socs[:, :, :, index] = 'NaN'


    elif inputs['code'] == 'molcas':   # TO DO
        pass


if __name__ == "__main__":
    pwd = os.getcwd()  # parent dir for all calculations
    list_geom = coordinate_reader(referance_geom_raw)
    if inputs['code'] == 'molpro':
        couplings = molpro.state_couplings(inputs['states'][-1])  # All possible state couplings
    elif inputs['code'] == 'molcas':
        pass
    pmanager = setup.ProccessManager()  # Force global mem sharing for ouput data
    pmanager.start()
    data = setup.DataStore(len(inputs['atoms']), inputs['states'][-1], couplings, len(list_geom), pmanager)
    cpus = int((multiprocessing.cpu_count())/2)
    pool = multiprocessing.Pool(processes=cpus)
    pool.map_async(run, [(k, v) for k, v in enumerate(list_geom)])
    pool.close()
    pool.join()
    np.save('ENERGIES.npy', data.energies)  # Extract data
    if inputs['nacme'] == 'yes':
        np.save('NACMES.npy', data.nacmes)
    if inputs['grad'] == 'yes':
        np.save('GRADS.npy', data.grads)
    if inputs['soc'] == 'yes':
        np.save('SOCS.npy', data.socs)
