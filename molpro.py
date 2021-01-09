import re
import os
import numpy as np

"""This module contains all routines for setting up Molpro inputs.
    Includes - Geometry, basis, symmetry etc. initialisation,
    SPE (CASSCF/MRCI/CASPT2), NACME(CASSCF/MRCI), SOC(CASSCF/MRCI),
    and GRADIENTS (CASSCF/MRCI/CASPT2). For CASPT2 couplings/ gradients
    the CASSCF couplings/ gradients must be transformed using the CASPT2
    mixing matrix."""

#########################
# INPUT FILE GENERATORS #
#########################


def setupSPE(inputs: dict, geom, pwd: str, gp: int) -> str:
    """Sets up a single point energy calculation given inputs from input.json.
        Options - CASSCF/ MRCI/ CASPT2. Will project through active spaces provided
        sequentually. Multi-Referance methods are run using the last casscf projection.
        SOC matrix can be computed at the CASSCF or MRCI level. """
    calculationSPE = []
    workdir = '%s/GP%s/' % (pwd, gp)
    os.mkdir(workdir)
    inputfile = 'SPE_GP%s.input' % gp
    #  Initialise SPE calculation for given geometry.
    initCalc = "***PES Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    geometryBlock = geometrySetup(inputs, geom)
    geometryBlock = '\n'.join(geometryBlock)
    calculationSPE.extend([initCalc, geometryBlock])

    #  Set up chosen electronic structure theory
    casscf = casscfSetup(inputs, nacme=False, splitSpin=False)
    casscf = '\n'.join(casscf)
    if inputs['spe'] == 'casscf':
        calculationSPE.append(casscf)
    elif inputs['spe'] == 'mrci':
        mrci = mrciSetup(inputs, 6000.2, 8000.2, noexc=False, tdm=False)
        mrci = '\n'.join(mrci)
        calculationSPE.extend([casscf, mrci])
    elif inputs['spe'] == 'caspt2':
        caspt2 = caspt2Setup(inputs, grad=False)
        caspt2 = '\n'.join(caspt2)
        calculationSPE.extend([casscf, caspt2])

    #  Set up SOC calculation - if EST is not MRCI - pass CAS wf through CI
    #  programme using noexc flag
    if inputs['soc'] == 'yes':
        if inputs['spe'] != 'mrci':
            ciwf = mrciSetup(inputs, 6000.2, 8000.2, noexc=True, tdm=False)
            ciwf = '\n'.join(ciwf)
            records = re.findall(r"[0-9]{4}\.2", ciwf)
            soc = "{ci;hlsmat,amfi,%s\nprint,hls=1};" % (','.join(records))
            calculationSPE.extend([ciwf, soc])
        else:
            records = re.findall(r"[0-9]{4}\.2", mrci)
            soc = "{ci;hlsmat,amfi,%s\nprint,hls=1};" % (','.join(records))
            calculationSPE.append(soc)

    calculationSPE = '\n'.join(calculationSPE)
    f = open(workdir+inputfile, 'w+')
    f.write(calculationSPE)
    f.close()
    return workdir, inputfile


def nacmeSetup(inputs: dict, refGeom, workdir: str, gp: int) -> (str, list):
    """Sets up NACME calculation at CASSCF or MRCI level. Done using
       the molpro 3-point DDR procedure, one calculation at a referance geometry
       and two (positive & negative) at displaced geometries - done for each atom
       along each x, y and z axis."""
    nacmeInput = []
    states = inputs['states'][-1]  # Only need final CAS projection for NACME.
    nacmeWorkdir = '%s/NACMEs/' % workdir
    inputfile = "NACME_GP%s.input" % (gp)
    os.mkdir(nacmeWorkdir)  # Setup dir for nacmes as subdir of grid point dir
    if inputs['nacme_level'] == 'casscf':  # If not MRCI, pass CAS wf to CI programme using noexc.
        noexc = True
    elif inputs['namcme_level'] == 'mrci':
        noexc = False
    [nacmeRef, wfsRef, tdmsRef] = nacmeReferanceSetup(inputs, refGeom, noexc)

    # Initialise displaced geoms to be used in following calculation
    casscf = casscfSetup(inputs, nacme=True, splitSpin=False)
    casscfIndexs = [ind for ind, line in enumerate(casscf) if 'casscf' in line]
    casscf = casscf[casscfIndexs[-1]:]  # For displaced calc - only need final CAS projection
    casscf = '\n'.join(casscf)
    [nd, dx, dy, dz] = displace(inputs['symm'], len(inputs['atoms']), inputs['dr'], refGeom)
    [orbitals, dwfs, dtdms] = recordIncrementer(2140.2, wfsRef, tdmsRef, nd, np.count_nonzero(states))

    # Generate displaced calculations in each x, y, z coord
    if dx is not None:
        dx = nacmeDisplacedSetup(inputs, dx, wfsRef, tdmsRef, dwfs, dtdms, orbitals, casscf, noexc, states, nd)
    if dy is not None:
        dy = nacmeDisplacedSetup(inputs, dy, wfsRef, tdmsRef, dwfs, dtdms, orbitals, casscf, noexc, states, nd)
    if dz is not None:
        dz = nacmeDisplacedSetup(inputs, dz, wfsRef, tdmsRef, dwfs, dtdms, orbitals, casscf, noexc, states, nd)
    daxis = [True if dnacme is not None else False for dnacme in (dx, dy, dz)]
    dnacmes = ['\n'.join(dnacme) for dnacme in (dx, dy, dz) if dnacme is not None]
    nacmeInput.extend(['\n'.join(nacmeRef), '\n'.join(dnacmes)])

    f = open(nacmeWorkdir+inputfile, 'w+')
    f.write('\n'.join(nacmeInput))
    f.close()
    return nacmeWorkdir, inputfile, daxis


def nacmeReferanceSetup(inputs: dict, refGeom, noexc: bool) -> (list, list, list):
    """Procedure for setting up referance calculation to be used in
       NACME calculation. The referance calculation must project through all
       user specified active spaces."""
    #  Initialise reference calculation where no atoms are displaced - project
    #  through all specified CAS calculations
    nacmeRef = []
    initCalc = "***NACME Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    geometryBlock = geometrySetup(inputs, refGeom)
    casscf = casscfSetup(inputs, nacme=False, splitSpin=False)
    nacmeRef.extend([initCalc, '\n'.join(geometryBlock), '\n'.join(casscf), '\n'])
    # Set up following CI calculation for reference - do tdm calculation
    ciwf = mrciSetup(inputs, 6000.2, 8000.2, noexc, tdm=True)
    ciwf = '\n'.join(ciwf)
    nacmeRef.append(ciwf)
    wfsRef = re.findall(r"6[0-9]{3}\.2", ciwf)   # Get list of reference wf
    tdmsRef = re.findall(r"8[0-9]{3}\.2", ciwf)  # and tdm record files
    return nacmeRef, wfsRef, tdmsRef


def nacmeDisplacedSetup(inputs, geom, wfsRef, tdmsRef, dwfs, dtdms, orbitals, casscf, noexc, states, nd):
    """Sets up calculations at displaced coordinates for NACME calculation.
       Here only the final active space in the preceeding referance calculation
       is required as the orbitals/CI coeffs from reference calculation are
       read into each displaced calcualtion as a starting guess"""
    displacedCalculation = []
    orbitalsGrouped = splitList(orbitals, 2)  # Groups records for +ve and -ve
    tdmsGrouped = splitList(dtdms, 2)      # displacements - easier indexing
    count = 0

    for i in range(nd):
        dwf = dwfs[i]
        dtdm = dtdms[i]
        dorbital = orbitals[i]
        geomDisplaced = geometrySetup(inputs, geom[i, :, :])
        dcasscf = re.sub(r'orbital,2140.2', 'orbital,%s' % dorbital, casscf)  # displaced casscf
        displacedCalculation.extend(['\n', '\n'.join(geomDisplaced), dcasscf, '\n'])
        dciwf = mrciSetup(inputs, dwf[0], dtdm[0], noexc, tdm=False)  # displaced CI
        displacedCalculation.extend(['\n'.join(dciwf), '\n'])
        citdm = citdmSetup(wfsRef, dwf, dtdm, np.count_nonzero(states))  # Overlap of displaced w/ reference
        displacedCalculation.extend(['\n'.join(citdm), '\n'])
        couplings = stateCouplings(states)  # Get all possible coupling combinations

        if i % 2 != 0:
            ddr = ddrSetup(inputs['dr'], 2140.2, tdmsRef, orbitalsGrouped[count], tdmsGrouped[count], states, couplings)  # NACME for each combination
            displacedCalculation.append('\n'.join(ddr))
            count += 1

    return displacedCalculation


def gradientSetup(inputs: dict, geom, workdir: str, gp: int):
    """Sets up gradient calculation, gradients are calculated for all states.
       If a CASSCF or CASPT2 SPE energy calculation preceeds then analytical
       gradients are calculated using the rs2 programme. For MRCI numerical
       gradientds are used."""
    gradientCalculation = []
    states = inputs['states'][-1]
    gradWorkdir = "%s/GRAD/" % workdir
    inputfile = "GRAD_GP%s.input" % gp
    os.mkdir(gradWorkdir)
    initCalc = "***GRAD Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    geometryBlock = geometrySetup(inputs, geom)
    gradientCalculation.extend([initCalc, '\n'.join(geometryBlock)])
    if inputs['spe'] != 'mrci':
        if inputs['singlets'] and inputs['triplets'] == 'yes':
            nSinglets = sum(states[:2])
            nTriplets = sum(states[2:])
            casscf = casscfSetup(inputs, nacme=False, splitSpin=True)
            casscfIndex = [ind for ind, line in enumerate(casscf) if 'casscf' in line]
            casscfInit = casscf[:casscfIndex[-2]]
            casscfSinglet = casscf[casscfIndex[-2]:casscfIndex[-1]]
            casscfTriplet = casscf[casscfIndex[-1]:]
            gradients = caspt2Setup(inputs, grad=True)
            singletGradients = gradients[:nSinglets]
            tripletGradients = gradients[nSinglets:]
            gradientCalculation.extend(['\n'.join(casscfInit), '\n'.join(casscfSinglet),
                                        '\n'.join(singletGradients), '\n'.join(casscfTriplet),
                                        '\n'.join(tripletGradients)])
        else:
            casscf = casscfSetup(inputs, nacme=False, splitSpin=False)
            gradients = caspt2Setup(inputs, grad=True)
            gradientCalculation.extend(['\n'.join(casscf), '\n'.join(gradients)])

    elif inputs['spe'] == 'mrci':
        pass  # TO DO - Numerical grads

    f = open(gradWorkdir+inputfile, 'w+')
    f.write('\n'.join(gradientCalculation))
    f.close()
    return gradWorkdir, inputfile

###########################
# NACME SETUP UTILITIES #
###########################


def stateCouplings(states: list) -> list:
    """Determines the possible combinations of coupled states for a given
       symmetry/ multiplicity. Returns a list of lists, where each sublist
       corrosponds to a different symm/spin. Each sublist contains N tuples
       of coupled state combinations.."""
    couplingCombinations = []
    for i, nstate in enumerate(states):
        temp = []
        if i % 2 == 0:  # Even elements corrspond to A' symm
            for j in range(1, nstate+1):
                for k in range(1, j+1):  # Get lower triangle of matrix only
                    if j != k:
                        temp.append((j+0.1, k+0.1))
            couplingCombinations.append(temp)

        else:           # Odd elements corrospond to A" symm
            for j in range(1, nstate+1):
                for k in range(1, j+1):
                    if j != k:
                        temp.append((j+0.2, k+0.2))
            couplingCombinations.append(temp)
    return couplingCombinations


def recordIncrementer(orbital: float, wfsRef: list, tdmsRef: list, nd: int, nblocks: int) -> (list, list, list):
    """Given a list of referance records for each symm/ multiplicity - increment
       the record number for n displacements. Returns lists of tuples, each tuple
       corrosponds to records for one displacement - elements of tuples are
       the symmetry/ multiplicity. """
    dwfs, dtdms, orbitals = [], [], []
    wf = float(wfsRef[-1])
    tdm = float(tdmsRef[-1])
    orbital = float(orbital)

    for disp in range(nd):      # Loop over num. displacements
        orbital = orbital + 1   # For each displacement increment orbital rec by one
        orbitals.append(orbital)
        wfTemp, tdmTemp = [], []
        for irrep in range(nblocks):  # Loop over symm/ spin blocks
            wf = wf + 1
            tdm = tdm + 1             # Increment each wf/ tdm rec by one for each symm/ spin
            wfTemp.append(wf)
            tdmTemp.append(tdm)
        dwfs.append(tuple(wfTemp))   # symm/ spin blocks for each displacement grouped into tuples
        dtdms.append(tuple(tdmTemp))
    return orbitals, dwfs, dtdms


def displace(uaxis: str, natom: int, dr: float, refGeom) -> (int, list, list, list):
    """Generates geometries displaced by an infinitesimal ammount dr along
       the x, y, and z axis in cartesian coordinates. Returns three np.arrays
       with dimensions [num. displacements, 3, natom] for each axis
       and the number of displacements. """
    numDisplacements = 6  # num axis (3) * 2
    displacedGeoms = np.zeros((numDisplacements*natom,  3, natom))
    c = 0

    for i in range(3):  # Loop over axes
        for j in range(natom):  # Loop over atoms
            positiveDisplacement = refGeom.copy()  # positive and negative
            negativeDisplacement = refGeom.copy()  # displacement for each element
            positiveDisplacement[j, i] = positiveDisplacement[j, i] + dr
            negativeDisplacement[j, i] = negativeDisplacement[j, i] - dr
            displacedGeoms[c, :, :] = positiveDisplacement
            displacedGeoms[c+1, :, :] = negativeDisplacement
            c += 2

    dx, dy, dz = np.split(displacedGeoms, 3)
    if uaxis == 'x':  # Do not want to displace along unique symmetry axis
        dx = None
    elif uaxis == 'y':
        dy = None
    elif uaxis == 'z':
        dz = None
    return numDisplacements, dx, dy, dz


def citdmSetup(wfsRef: list, wfs: list, tdms: list, nblocks: int) -> list:
    """ This function sets up a transition dipole moment calculation using
        the MRCI programme in Molpro. """
    citdmBlock = []
    for i in range(nblocks):  # Loop over symm/ spin blocks
        citdm = "{ci;trans,%s,%s;dm,%s}" % (wfsRef[i], wfs[i], tdms[i])
        citdmBlock.append(citdm)
    return citdmBlock


def ddrSetup(dr: float, orbitalRef: float, tdmsRef: list, orbitalDisp: list, tdmsDisp: list, states: list, couplingCombinations: list) -> list:
    """Sets up nacme computation using finite differences. Loops over all
       states of each symm/spin and generates a ddr calculation for each
       possible combination of states in each symm/spin. """
    ddrBlock = []
    count = 0
    for i in range(len(states)):  # Loop over symm/ spin blocks
        if states[i] != 0:
            for coupledPair in couplingCombinations[i]:
                ddr = "{ddr,2*%s\norbital,%s,%s,%s;\ndensity,%s,%s,%s;\nstates,%s,%s}" % (dr,
                        orbitalRef, orbitalDisp[0], orbitalDisp[1], tdmsRef[count], tdmsDisp[0][count], tdmsDisp[1][count], coupledPair[0], coupledPair[1])
                ddrBlock.append(ddr)
            count += 1
    return ddrBlock


##############################
# Ab-Initio Methods Setup    #
##############################


def geometrySetup(inputs: dict, geom) -> list:
    """ Sets up the geometry block for a calculation given inputs.
        Block includes basis, symmetry generators and xyz coordinates"""
    geometryBlock = []
    basis = "basis=%s" % inputs['basis']

    if inputs['symm'] != 'nosymm':
        symmetry = "symmetry,%s" % inputs['symm']
    else:
        symmetry = "symmetry,nosymm"

    geom = formatGeometry(geom, inputs['atoms'])  # Format array into coords w/ atomic symbols in column 0
    geometry = "geomtyp=xyz\nbohr\ngeometry={\n%s\nSurf Calc\n%s\n}" % (len(inputs['atoms']), geom)
    geometryBlock.extend((basis, symmetry, geometry))
    return geometryBlock


def casscfSetup(inputs: dict, nacme: bool, splitSpin: bool) -> list:
    """Sets up N single point calculations at the CASSCF level given N active spaces
    and states. Works for any combination of symmetry and multiplicities."""
    casscfBlock = []
    if not nacme:
        startwf = ''
    if nacme:
        startwf = '\nstart,2140.2'
    casscfBlock.append('{hf}')

    if inputs['symm'] != 'nosymm':
        for i in range(len(inputs['states'])):  # Loop over projections
            occ = inputs['occ'][i]
            closed = inputs['closed'][i]        # Init active spaces and num. states
            states = inputs['states'][i]        # for each projection
            casscf = "{casscf\nfrozen,0\nclosed,%s,%s\nocc,%s,%s%s\norbital,2140.2" % (closed[0], closed[1], occ[0], occ[1], startwf)
            casscfBlock.append(casscf)
            if inputs['singlets'] == 'yes' and inputs['triplets'] == 'yes':
                wf1 = "wf,%s,1,0\nstate,%s" % (inputs['nelec'], states[0])  # Init wf depending on states
                wf2 = "wf,%s,2,0\nstate,%s" % (inputs['nelec'], states[1])
                wf3 = "wf,%s,1,2\nstate,%s" % (inputs['nelec'], states[2])
                wf4 = "wf,%s,2,2\nstate,%s\n};" % (inputs['nelec'], states[3])
                if splitSpin:  # For analytical gradient calculations
                    casscfBlock.extend([wf1, wf2+'};', casscf, wf3, wf4])  # can not SA over diff spin states here
                else:
                    casscfBlock.extend([wf1, wf2, wf3, wf4])  # Do SA over different spins
            elif inputs['singlets'] == 'yes' and inputs['triplets'] == 'no':
                wf1 = "wf,%s,1,0\nstate,%s" % (inputs['nelec'], states[0])
                wf2 = "wf,%s,2,0\nstate,%s\n};" % (inputs['nelec'], states[1])
                wf3 = ""
                wf4 = ""
                casscfBlock.extend([wf1, wf2, wf3, wf4])
            elif inputs['singlets'] == 'no' and inputs['triplets'] == 'yes':
                wf1 = ""
                wf2 = ""
                wf3 = "wf,%s,1,2\nstate,%s" % (inputs['nelec'], states[2])
                wf4 = "wf,%s,2,2\nstate,%s\n};" % (inputs['nelec'], states[3])
                casscfBlock.extend([wf1, wf2, wf3, wf4])

    elif inputs['symm'] == 'nosymm':
        for i in range(len(states)):
            occ = inputs['occ'][i]
            closed = inputs['closed'][i]
            states = inputs['states'][i]
            casscf = "{casscf\nfrozen,0\nclosed,%s\nocc,%s%s\norbital,2140.2" % (closed[0], occ[0], startwf)
            casscfBlock.append(casscf)
            if inputs['singlets'] == 'yes' and inputs['triplets'] == 'yes':
                wf1 = "wf,%s,1,0\nstate,%s" % (inputs['nelec'], states[0])  # Only one irrep for each spin
                wf2 = "wf,%s,1,2\nstate,%s};" % (inputs['nelec'], states[2])  # if no symm is chosen
                if splitSpin:  # do not SA over different spin multiplicities
                    casscfBlock.extend([wf1+'};', casscf, wf2])
                else:
                    casscfBlock.extend([wf1, wf2])
            elif inputs['singlets'] == 'yes' and inputs['triplets'] == 'no':
                wf1 = "wf,%s,1,0\nstate,%s};" % (inputs['nelec'], states[0])
                wf2 = ""
                casscfBlock.extend([wf1, wf2])
            elif inputs['singlets'] == 'no' and inputs['triplets'] == 'yes':
                wf1 = ""
                wf2 = "wf,%s,1,2\nstate,%s};" % (inputs['nelec'], states[2])
                casscfBlock.extend([wf1, wf2])

    return casscfBlock


def mrciSetup(inputs: dict, orbitalRecord: float, tdmRecord: float, noexc: bool, tdm: bool) -> list:
    """Sets up an MRCI calculation (based on final CASSCF projection)
       for each symmetry/ multiplicity provided. No excitations (noexc) input
       bool allows to setup actual MRCI calculation or just to pass a CASSCF wf
       through the CI programme for use in other programmes eg. SOC, NACME.
       Transition dipole moment (tdm) bool allows subsequent calculation and
       storage of the tdm for each state."""
    mrciBlock = []
    states = inputs['states'][-1]  # Only need final CAS projections states
    wfs = ["wf,%s,1,0" % inputs['nelec'], "wf,%s,2,0" % inputs['nelec'],
           "wf,%s,1,2" % inputs['nelec'], "wf,%s,2,2" % inputs['nelec']]  # wfs corrosponding to list of states

    for i in range(len(states)):
        if states[i] != 0:  # Loop over symm/ spin blocks
            if not noexc:  # Include excitations - run MRCI
                if not tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s};" % (wfs[i], orbitalRecord, states[i])
                if tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s\ndm,%s};" % (wfs[i], orbitalRecord, states[i], tdmRecord)

            else:          # Exclude excitations - pass CAS wf through CI programme
                if not tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s;NOEXC};" % (wfs[i], orbitalRecord, states[i])
                if tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s;NOEXC\ndm,%s};" % (wfs[i], orbitalRecord, states[i], tdmRecord)

            orbitalRecord += 1  # Increment records
            tdmRecord += 1
            mrciBlock.append(mrci)

    return mrciBlock


def caspt2Setup(inputs: dict, grad: bool) -> list:
    """Sets up N XMS-CASPT2 calculations for each symmetry/ multiplicity
       with referance to the final CASSCF projection. The level shift is
       specified by user in input.json. If grad == True then an analytical
       gradient calculation is set up instead."""
    caspt2Block = []
    states = inputs['states'][-1]  # Final CAS projection
    wfs = ["wf,%s,1,0" % inputs['nelec'], "wf,%s,2,0" % inputs['nelec'],  # Corrosponding wfs
           "wf,%s,1,2" % inputs['nelec'], "wf,%s,2,2" % inputs['nelec']]

    for i in range(len(states)):  # Loop over symm/ spin
        if states[i] != 0:
            if not grad:
                caspt2 = "{rs2,shift=%s,MIX=%s,XMS=1\n%s\nstate,%s};" % (inputs['level_shift'], states[i], wfs[i], states[i])
                caspt2Block.append(caspt2)
            else:
                for root in range(states[i]):
                    caspt2 = "{rs2;NOEXC\n%s\nOPTION,NSTATI=%s;state,1,%s};\nforces" % (wfs[i], states[i], root+1)
                    caspt2Block.append(caspt2)

    return caspt2Block


def formatGeometry(xyzGeom, atoms: list) -> str:
    """Turns coorindates in a np.array into a string with correct format. """
    geomList = [",\t".join(row) for row in xyzGeom.astype(str)]
    geomList = [atoms[i] + ",\t" + row for i, row in enumerate(geomList)]
    geom = '\n'.join(geomList)
    return geom


def splitList(alist, chunks):
    """Splits a list into a list of lists with equal no. elements (chunks)"""
    newList = []
    for i in range(0, len(alist), chunks):
        newList.append(alist[i:i+chunks])

    return newList
