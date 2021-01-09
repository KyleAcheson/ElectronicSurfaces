#!/usr/bin/python3
import os
import sys
import multiprocessing
import numpy as np
import scipy.io
import argparse
import textwrap
import json
import copy
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
                    "mem": file memory
                    "symm": Symmetry generator (X, Y, Z) or  'nosymm' - max Cs symmetry
                    "basis": basis set in chosen programme format
                    "soc": yes/no
                    "nacme": yes/no
                    "nacme_level": casscf/mrci
                    "dr": displacement value (angstrom) if nacmes are numerical
                    "paxis": principle axis
                    "grad": yes/no - gradient calculation
                    "spe": casscf/mrci/caspt2 - single point energy calculation
                    "atoms": list of atoms (list of str)
                    "nelec": number electrions (int)
                    "occ": nlists of occupied orbitals - [[Irrep1, Irrep2], ... ]
                    "closed": nlists of closed orbitals - [[Irrep1, Irrep2], ... ]
                    "states": nists of no. states - [[Irrep1, Irrep2], ... ]\n
                    Optional Input Specifications:\n
                    "pspace": molpro pspace thresh\n
                    NB:
                    Can provide n lists (n = 1, 2, 3...) of orbitals and states
                    to project through active spaces at CASSCF level.
                    If a MR method is chosen, the MR calculation will corrospond
                    to the last list of active space and states.\n
                    If running with  symmetry off, set the orbitals and states
                    that corrospond to Irrep2 to 0.\n'''))
parser.add_argument("-g", dest="inputGeomFile", required=True, type=str,
                    help="Geometry in xyz coordinates at equalibrium.")
args = parser.parse_args()
inp = args.inputfile
inputGeom = args.inputGeomFile

##################
#  Sanity Checks #
##################

requiredInputs = ["code", "mem", "symm", "basis", "soc", "dr", "paxis",
                   "nacme", "nacme_level", "grad", "spe", "nelec", "occ",
                   "atoms", "closed", "states"]

inputError = ("ERROR: Missing required input fields from input.json.\n"
               "RUN: %s --help for a list of fields." % sys.argv[0])

# Load input params in json format into dict
f = open(inp, 'r')
inputs = json.load(f)
f.close()
# Make sure all keys are correct case
inputs = {k.lower(): v for k, v in inputs.items()}

# Check all required input params provided
for k in requiredInputs:
    if k not in inputs.keys():
        sys.exit(inputError)

# Check input params have correct data type
if not isinstance(inputs['code'], str):
    sys.exit("Error: '%s' must be type str" % inputs['code'])
if not isinstance(inputs['mem'], str):
    sys.exit("Error: '%s' must be type str" % inputs['mem'])
if not isinstance(inputs['symm'], str):
    sys.exit("Error: '%s' must be type str" % inputs['symm'])
if not isinstance(inputs['singlets'], str):
    sys.exit("Error: '%s' must be type str" % inputs['singlets'])
if not isinstance(inputs['triplets'], str):
    sys.exit("Error: '%s' must be type str" % inputs['triplets'])
if not isinstance(inputs['basis'], str):
    sys.exit("Error: '%s' must be type str" % inputs['basis'])
if not isinstance(inputs['soc'], str):
    sys.exit("Error: '%s' must be type str" % inputs['soc'])
if not isinstance(inputs['nacme'], str):
    sys.exit("Error: '%s' must be type str" % inputs['nacme'])
if not isinstance(inputs['nacme_level'], str):
    sys.exit("Error: '%s' must be type str" % inputs['nacme_level'])
if not isinstance(inputs['grad'], str):
    sys.exit("Error: '%s' must be type str" % inputs['grad'])
if not isinstance(inputs['spe'], str):
    sys.exit("Error: '%s' must be type str" % inputs['spe'])
if not isinstance(inputs['paxis'], str):
    sys.exit("Error: '%s' must be type str" % inputs['paxis'])
if not isinstance(inputs['atoms'], list):
    sys.exit("Error: '%s' must be type list of str" % inputs['atoms'])
if not isinstance(inputs['nelec'], int):
    sys.exit("Error: '%s' must be type int" % inputs['nelec'])
if not isinstance(inputs['dr'], float):
    sys.exit("Error: '%s' must be type float" % inputs['dr'])
if not isinstance(inputs['occ'], list):
    sys.exit("Error: '%s' must be type list of lists of int" % inputs['occ'])
if not isinstance(inputs['closed'], list):
    sys.exit("Error: '%s' must be type list of lists of int" % inputs['closed'])
if not isinstance(inputs['states'], list):
    sys.exit("Error: '%s' must be type list of lists of int" % inputs['states'])

# Ensure values of type str are lowercase
for k, v in inputs.items():
    if isinstance(v, str):
        inputs[k] = v.lower()
inputs['dr'] = inputs['dr']*ang2au

# Check inputs have acceptable values
codeExcept = ['molpro', 'molcas']
ynExcept = ['yes', 'no']
symmExcept = ['x', 'y', 'z', 'nosymm']
nacmeExcept = ['casscf', 'mrci']
paxisExcept = ['x', 'y', 'z']
exceptSPE = ['casscf', 'mrci', 'caspt2']

if inputs['code'] not in codeExcept:
    sys.exit("Error: '%s' not available, choose from %s." % (inputs['code'], codeExcept))
if inputs['symm'] not in symmExcept:
    sys.exit("Error: choose from %s." % (symmExcept))
if inputs['soc'] not in ynExcept:
    sys.exit("Error: choose from %s." % (ynExcept))
if inputs['nacme'] not in ynExcept:
    sys.exit("Error: choose from %s." % (ynExcept))
if inputs['grad'] not in ynExcept:
    sys.exit("Error: choose from %s." % (ynExcept))
if inputs['singlets'] not in ynExcept:
    sys.exit("Error: choose from %s." % (ynExcept))
if inputs['triplets'] not in ynExcept:
    sys.exit("Error: choose from %s." % (ynExcept))
if inputs['nacme_level'] not in nacmeExcept:
    sys.exit("Error: '%s' not available, choose from %s." % (inputs['nacme_level'], nacmeExcept))
if inputs['paxis'] not in paxisExcept:
    sys.exit("Error: '%s' not available, choose from %s." % (inputs['paxis'], paxisExcept))
if inputs['spe'] not in exceptSPE:
    sys.exit("Error: '%s' not available, choose from %s." % (inputs['spe'], exceptSPE))

# While waiting to add new functionality
if inputs['spe'] == 'mrci' and inputs['grad'] == 'yes':
    sys.exit("Error: Gradients for MRCI are yet to be added")
if inputs['code'] == 'molcas':
    sys.exit("Error: Molcas is not yet available")


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


def coordinateGenerator(refGeom):
    listGeom = []
    for i in np.arange(1.60, 2.1, 0.1):
        new_geom = refGeom.copy()
        new_geom[1, 2] = i
        new_geom = new_geom*ang2au
        listGeom.append(new_geom)
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
        [normalTermination, outputSPE] = util.runCalculation(molproKeys, pwd, workdirSPE, inputSPE, index)
        if normalTermination:
            data.energiesExtract(workdirSPE+outputSPE, inputs['spe'], molproKeys['energy_regex'],
                                  molproKeys['cas_prog'], index)

            if inputs['nacme'] == 'yes':
                [nacmeWorkdir, nacmeInput, daxes] = molpro.nacmeSetup(inputs, geom, workdirSPE, index)
                [nacmeNormalTermination, nacmeOutput] = util.runCalculation(molproKeys, pwd, nacmeWorkdir, nacmeInput, index)
                if nacmeNormalTermination:
                    data.nacmeExtract(nacmeWorkdir+nacmeOutput, molproKeys['nacme_regex'], index, daxes)
                else:
                    data.nacmes[:, :, :, index] = 'NaN'

            if inputs['grad'] == 'yes':
                [gradWorkdir, gradInput] = molpro.gradientSetup(inputs, geom, workdirSPE, index)
                [gradNormalTermination, gradOutput] = util.runCalculation(molproKeys, pwd, gradWorkdir, gradInput, index)
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


    elif inputs['code'] == 'molcas':   # TO DO
        pass


if __name__ == "__main__":
    try:
        pwd = os.getcwd()  # parent dir for all calculations
        f = open(inputGeom, 'r')
        refGeom = np.genfromtxt(f, usecols=(1, 2, 3), encoding=None)
        f.close()
        listGeom = coordinateGenerator(refGeom)  # Generate coordinates

        if inputs['code'] == 'molpro':  # All possible state couplings
            couplings = molpro.stateCouplings(inputs['states'][-1])
        elif inputs['code'] == 'molcas':
            pass

        pmanager = setup.ProccessManager()  # Force global mem sharing for ouput data
        pmanager.start()
        data = setup.datastore(refGeom, inputs['states'][-1], couplings, len(listGeom), pmanager)
        runPool = multiprocessing.Pool()
        runPool.map(run, [(k, v) for k, v in enumerate(listGeom)])

        np.save('ENERGIES.npy', data.energies)  # Extract data
        if inputs['nacme'] == 'yes':
            np.save('NACMES.npy', data.nacmes)
        if inputs['grad'] == 'yes':
            np.save('GRADS.npy', data.grads)

    except StopIteration:
        print('Please ensure goemetry file is correct.')
