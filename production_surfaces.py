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

ang2au = 1.88973

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
parser.add_argument("-g", dest="inputGeomFile", required=True, type=str,
                    help=textwrap.dedent('''\
                    Nd Array of Geometry in xyz coordinates (ANGSTROM) over all space\n
                    Must be of type .mat or .npy\n
                    Must have dimensions (Natom, 3, Ngrid DoF1, Ngrid DoF2, ..., Ngrid DoFN)\n
                    Degrees of freedom will be flattened from left to right to give a list of geometries\n'''))
parser.add_argument("-q", dest="scriptHPC", required=False, type=str,
                    help="Submission script template for running on a HPC queueing system.")

args = parser.parse_args()
inp = args.inputfile
inputGeom = args.inputGeomFile
submitScript = args.scriptHPC

# Load inputs
f = open(inp, 'r')
inputs = json.load(f)
f.close()

fileType = inputGeom.split('.')[-1]
if fileType == 'npy':
    refGeom = np.load(inputGeom)
elif fileType == 'mat':
    refMat = sio.loadmat(inputGeom)
    keys = [k for k in refMat if '__' not in k]
    if len(keys) > 1:
        raise TypeError('.Mat file must have only one variable')
    else:
        keyval = keys[0]
    refGeom = refMat[keyval]
else:
    raise TypeError('Input goemetry file must be .npy or .mat')


try:
    f = open(submitScript, 'r')
    submissionScriptLines = f.readlines()
    f.close()
    submissionScript = ''.join(submissionScriptLines)
except TypeError:
    submissionScript = None


##################
#  Sanity Checks #
##################


requiredInputs = ["code", "hpc", "mem", "symm", "basis", "soc", "dr", "paxis",
                  "nacme", "nacme_level", "grad", "spe", "nelec", "occ",
                  "atoms", "closed", "states", "singlets", "triplets"]

codeExcept = ['molpro', 'molcas']
ynExcept = ['yes', 'no']  # For all y/n inputs
symmExcept = ['x', 'y', 'z', 'nosymm']
nacmeExcept = ['casscf', 'mrci']
paxisExcept = ['x', 'y', 'z']
exceptSPE = ['casscf', 'mrci', 'caspt2']
exceptHPC = ['local', 'Sun Grid Engine', 'PBS']

inputError = ("ERROR: Missing required input fields from input.json.\n"
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

for k in requiredInputs:
    if k not in inputs.keys():
        raise ValueError(inputError)

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
inputs['dr'] = inputs['dr']*ang2au

# Check inputs have acceptable values

if inputs['code'] not in codeExcept:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['code'], codeExcept))
if inputs['hpc'] not in exceptHPC:
    raise ValueError(" choose from %s." % exceptHPC)
if inputs['symm'] not in symmExcept:
    raise ValueError(" choose from %s." % (symmExcept))
if inputs['soc'] not in ynExcept:
    raise ValueError(" choose from %s." % (ynExcept))
if inputs['nacme'] not in ynExcept:
    raise ValueError(" choose from %s." % (ynExcept))
if inputs['grad'] not in ynExcept:
    raise ValueError(" choose from %s." % (ynExcept))
if inputs['singlets'] not in ynExcept:
    raise ValueError(" choose from %s." % (ynExcept))
if inputs['triplets'] not in ynExcept:
    raise ValueError(" choose from %s." % (ynExcept))
if inputs['nacme_level'] not in nacmeExcept:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['nacme_level'], nacmeExcept))
if inputs['paxis'] not in paxisExcept:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['paxis'], paxisExcept))
if inputs['spe'] not in exceptSPE:
    raise ValueError(" '%s' not available, choose from %s." % (inputs['spe'], exceptSPE))


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

for occOrb in inputs['occ']:
    if inputs['symm'] != 'nosymm':
        if not isinstance(occOrb[0], int) or not isinstance(occOrb[1], int):
            raise ValueError('Ensure occupied orbitals are integers for both irreps.')
    elif inputs['symm'] == 'nosymm':
        if not isinstance(occOrb[0], int) or occOrb[1] != 0:
            raise ValueError('Ensure occupied orbitals of irrep 1 is an integer and irrep 2 is 0.')

for closedOrb in inputs['closed']:
    if inputs['symm'] == 'yes':
        if not isinstance(closedOrb[0], int) or not isinstance(closedOrb[1], int):
            raise ValueError('Ensure closed orbitals are integers for both irreps.')
    elif inputs['symm'] == 'nosymm':
        if not isinstance(closedOrb[0], int) or closedOrb[1] != 0:
            raise ValueError('Ensure closed orbitals of irrep 1 is an integer and irrep 2 is 0.')


# Ensure input geometry and HPC queue submission script is correct

if refGeom.shape[0:2] != (len(inputs['atoms']), 3):
    raise InputGeomError(refGeom)
if submitScript is None and inputs['hpc'] != 'local':
    raise ValueError(" If not running via. a queueing system set hpc input to local")
elif submitScript is not None and inputs['hpc'] == 'none':
    raise ValueError(" If running via. a queueing system set hpc input to 'local', 'Sun Grid Engine' or 'PBS'")

if submissionScript is not None and not re.search(r'template', submissionScript):
    raise ValueError(" Place the keyword 'template' where your input file goes in the HPC submission script")


# While waiting to add new functionality

if inputs['spe'] == 'mrci' and inputs['grad'] == 'yes':
    raise ValueError(" Gradients for MRCI are yet to be added")
if inputs['code'] == 'molcas':
    raise ValueError(" Molcas is not yet available")


#########################
# Quantum Code KEYWORDS #
#########################

molproKeys = {'cas_prog': '1PROGRAM * MULTI',
              'nacme_regex': '!Total NACME.*$',
              'termination_code': 'Variable memory released',
              'submit_command': 'molpro',
              'output_ext': 'out'}

molcasKeys = {}  # TO DO


if inputs['code'] == 'molpro':  # Set input dependent keys
    if inputs['spe'] == 'casscf':
        molproKeys['energy_regex'] = "MCSCF STATE [0-9]{1,2}\.\d Energy.*$"
        molproKeys['grad_regex'] = 'RSPT2 GRADIENT FOR STATE'
        molproKeys['numerical'] = False
    elif inputs['spe'] == 'mrci':
        molproKeys['energy_regex'] = "MRCI STATE [0-9]{1,2}\.\d Energy.*$"
        molproKeys['grad_regex'] = 'Numerical gradient'
        molproKeys['numerical'] = True
    elif inputs['spe'] == 'caspt2':
        molproKeys['energy_regex'] = "RS2 STATE [0-9]{1,2}\.\d Energy.*$"
        molproKeys['grad_regex'] = 'RSPT2 GRADIENT FOR STATE'
        molproKeys['numerical'] = False


elif inputs['code'] == 'molcas':
    if inputs['spe'] == 'casscf':
        pass
    elif inputs['spe'] == 'mrci':
        pass
    elif inputs['spe'] == 'caspt2':
        pass


############################
# Coordinate System - Temp #
############################


#def coordinateGenerator(refGeom):
#    listGeom = []
#    for i in np.arange(1.40, 1.85, 0.05):
#        new_geom = refGeom.copy()
#        new_geom[1, 2] = i
#        new_geom = new_geom*ang2au
#        listGeom.append(new_geom)
#    return listGeom


def coordinateReader(GeomIn):
    GeomArr = GeomIn*ang2au
    listGeom = []
    dims = GeomArr.shape[2:]
    ngrid = int(np.product(dims))
    GeomFlat = GeomArr.reshape(len(inputs['atoms']), 3, ngrid)
    for i in range(ngrid):
        listGeom.append(GeomFlat[:, :, i])
    return listGeom

##################################
# MAIN - Run surface calculation #
##################################


def run(gridPoint):
    """Takes each geometry and an index for each grid point and runs the series
       of calculations specified by the user for either quantum chemistry code
       in parallel. Data for each grid point is sorted in global arrays according
       to the index."""
    index, geom = gridPoint

    if inputs['code'] == 'molpro':

        [workdirSPE, inputSPE] = molpro.setupSPE(inputs, geom, pwd, index)
        [normalTermination, outputSPE] = util.runCalculation(inputs['hpc'], molproKeys, pwd, workdirSPE, inputSPE, submissionScript, index)
        if normalTermination:
            data.energiesExtract(workdirSPE+outputSPE, inputs['spe'], molproKeys['energy_regex'],
                                  molproKeys['cas_prog'], index)
            if inputs['soc'] == 'yes':
                data.extractSOC(workdirSPE+outputSPE, index)

            if inputs['nacme'] == 'yes':
                [nacmeWorkdir, nacmeInput, daxes] = molpro.nacmeSetup(inputs, geom, workdirSPE, index)
                [nacmeNormalTermination, nacmeOutput] = util.runCalculation(inputs['hpc'], molproKeys, pwd, nacmeWorkdir, nacmeInput, submissionScript, index)
                if nacmeNormalTermination:
                    data.nacmeExtract(nacmeWorkdir+nacmeOutput, molproKeys['nacme_regex'], index, daxes)
                else:
                    data.nacmes[:, :, :, index] = 'NaN'

            if inputs['grad'] == 'yes':
                [gradWorkdir, gradInput] = molpro.gradientSetup(inputs, geom, workdirSPE, index)
                [gradNormalTermination, gradOutput] = util.runCalculation(inputs['hpc'], molproKeys, pwd, gradWorkdir, gradInput, submissionScript, index)
                if gradNormalTermination:
                    data.gradExtractMolpro(gradWorkdir+gradOutput, molproKeys['grad_regex'], molproKeys['numerical'], index)
                else:
                    data.grads[:, :, :, index] = 'NaN'

        else:
            data.energies[index, :] = 'NaN'
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
    listGeom = coordinateReader(refGeom)
    if inputs['code'] == 'molpro':  # All possible state couplings
        couplings = molpro.stateCouplings(inputs['states'][-1])
    elif inputs['code'] == 'molcas':
        pass
    pmanager = setup.ProccessManager()  # Force global mem sharing for ouput data
    pmanager.start()
    data = setup.datastore(refGeom, inputs['states'][-1], couplings, len(listGeom), pmanager)
    cpus = int((multiprocessing.cpu_count())/2)
    pool = multiprocessing.Pool(processes=cpus)
    pool.map_async(run, [(k, v) for k, v in enumerate(listGeom)])
    pool.close()
    pool.join()
    np.save('ENERGIES.npy', data.energies)  # Extract data
    if inputs['nacme'] == 'yes':
        np.save('NACMES.npy', data.nacmes)
    if inputs['grad'] == 'yes':
        np.save('GRADS.npy', data.grads)
    if inputs['soc'] == 'yes':
        np.save('SOCS.npy', data.socs)
