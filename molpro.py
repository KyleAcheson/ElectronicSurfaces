import re
import os
import numpy as np
import setup

"""This module contains all routines for setting up Molpro inputs.
    Includes - Geometry, basis, symmetry etc. initialisation,
    SPE (CASSCF/MRCI/CASPT2), NACME(CASSCF/MRCI), SOC(CASSCF/MRCI),
    and GRADIENTS (CASSCF/MRCI/CASPT2). For CASPT2 couplings/ gradients
    the CASSCF couplings/ gradients must be transformed using the CASPT2
    mixing matrix."""


class FileTreeExistsError(Exception):
    def __init__(self, directory, message='Delete, move or spawn new parent directory to run script from.'):
        self.wdir = directory
        self.message = message
        super().__init__(self.message)

#########################
# INPUT FILE GENERATORS #
#########################


def setup_spe(inputs: dict, geom, pwd: str, grid_index: int) -> str:
    """Sets up a single point energy calculation given inputs from input.json.
        Options - CASSCF/ MRCI/ CASPT2. Will project through active spaces provided
        sequentually. Multi-Referance methods are run using the last casscf projection.
        SOC matrix can be computed at the CASSCF or MRCI level.

        Inputs:
        inputs - dict of user specified inputs from input.json.
        geom - geometry for each calculation (np.array).
        pwd - parrent working directory for all calculations.
        grid_index - index for grid point according to position in list_geom.

        Outputs:
        workdir - new working directory for this grid points SPE calculation.
        inputfile - SPE inputfile name."""

    try:
        workdir = '%s/GP%s/' % (pwd, grid_index)
        os.mkdir(workdir)
    except FileExistsError:
        raise FileTreeExistsError(workdir)
    except:
        #setup.Logging(logfile, message="Can not create SPE Directory.")
        raise OSError('Can not create SPE directory')

    spe_calculation_input = []
    states = inputs['states']
    singlets = inputs['singlets']
    triplets = inputs['triplets']
    inputfile = 'SPE_GP%s.input' % grid_index
    molden = inputfile.split('.')[0]+'.molden'
    #  Initialise SPE calculation for given geometry.
    init_calculation = "***PES Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    geometry_block = geometry_setup(inputs, geom)
    geometry_block = '\n'.join(geometry_block)
    spe_calculation_input.extend([init_calculation, geometry_block, 'hf'])

    #  Set up chosen electronic structure theory
    casscf = setup_casscf(inputs, states, singlets, triplets, nacme=False)
    if inputs['spe'] == 'casscf':
        spe_calculation_input.append('\n'.join(casscf))
    elif inputs['spe'] == 'mrci':
        mrci = setup_mrci(inputs, states[-1], 6000.2, 8000.2, noexc=False, tdm=False)
        spe_calculation_input.extend(['\n'.join(casscf), '\n'.join(mrci)])
    elif inputs['spe'] == 'caspt2':
        caspt2 = setup_caspt2(inputs, states[-1], grad=False)
        spe_calculation_input.extend(['\n'.join(casscf), '\n'.join(caspt2)])

    #  Set up SOC calculation - if EST is not MRCI - pass CAS wf through CI
    #  programme using noexc flag
    if inputs['soc'] == 'yes':
        if inputs['spe'] != 'mrci':
            ciwf = setup_mrci(inputs, states[-1], 6000.2, 8000.2, noexc=True, tdm=False)
            ciwf = '\n'.join(ciwf)
            records = re.findall(r"[0-9]{4}\.2", ciwf)
            soc = "{ci;hlsmat,amfi,%s\nprint,hls=1};" % (','.join(records))
            spe_calculation_input.extend([ciwf, soc])
        else:
            records = re.findall(r"[0-9]{4}\.2", mrci)
            soc = "{ci;hlsmat,amfi,%s\nprint,hls=1};" % (','.join(records))
            spe_calculation_input.append(soc)

    molden_str = "put,molden,%s;" % molden
    spe_calculation_input.append(molden_str)
    spe_calculation_input = '\n'.join(spe_calculation_input)
    f = open(workdir+inputfile, 'w+')
    f.write(spe_calculation_input)
    f.close()
    return workdir, inputfile


def main_nacme(inputs: dict, referance_geom, workdir: str, grid_index: int):
    """Master function for setting up NACME calculation. Collects input for
       each spin and writes to file. Splitting the treatment of spin into
       two seperate calculations if singlets and triplets are both requested.

       Inputs:
       inputs - dict of input params from json file.
       referance_geom - geometry used to generate displacements for numerical calculation.
       workdir - parent directory for that grid point.
       grid_index - grid point number.

       Outputs:
       nacme_dir - directory path for NACME - subdir of workdir (str).
       inputfile - inputfile name for that grid point (str).
       displaced_axes_list - record of which axes have been displaced (list of bools). """

    try:
        nacme_dir = '%s/NACMEs/' % workdir
        os.mkdir(nacme_dir)  # Setup dir for nacmes as subdir of grid point dir
    except FileExistsError:
        raise FileTreeExistsError(nacme_dir)
    except:
        #setup.Logging(logfile, message="Can not create NACME Directory.")
        raise OSError('Can not create NACME directory')

    nacme_calculation_input = []
    singlets = inputs['singlets']
    triplets = inputs['triplets']
    states = inputs['states']
    inputfile = "NACME_GP%s.input" % (grid_index)
    if inputs['nacme_level'] == 'casscf':  # If not MRCI, pass CAS wf to CI programme using noexc.
        noexc = True
    elif inputs['namcme_level'] == 'mrci':
        noexc = False

    init_calculation = "***NACME Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    nacme_calculation_input.append(init_calculation)

    if inputs['singlets'] == 'yes' and inputs['triplets'] == 'yes':
        singlet_states = format_states(states, singlet=True, triplet=False)
        singlets, triplets = 'yes', 'no'
        [singlet_calculation_input, displaced_axes_list] = setup_nacme(inputs, singlet_states, singlets, triplets, referance_geom, noexc)
        nacme_calculation_input.append(singlet_calculation_input)

        triplet_states = format_states(states, singlet=False, triplet=True)
        singlets, triplets = 'no', 'yes'
        [triplet_calculation_input, displaced_axes_list] = setup_nacme(inputs, triplet_states, singlets, triplets, referance_geom, noexc)
        triplet_calculation_inputB = re.sub(r'hf', '', triplet_calculation_input)
        lines = triplet_calculation_inputB.split('\n')
        nacme_calculation_input.append(triplet_calculation_inputB)

    else:
        [nacme, displaced_axes_list] = setup_nacme(inputs, states, singlets, triplets, referance_geom, noexc)
        nacme_calculation_input.append(nacme)

    f = open(nacme_dir+inputfile, 'w+')
    f.write('\n'.join(nacme_calculation_input))
    f.close()
    return nacme_dir, inputfile, displaced_axes_list



def setup_nacme(inputs: dict, states: list, singlets: str, triplets: str, referance_geom, noexc: bool):
    """ Generates NACME input for the spin requested by user. If multiple
        active spaces are requested, only the final active space of the
        projection series is used for displaced calculations. Each atom
        in the molecule is displaced by +ve and -ve value of user input.

        Inputs (additional):
        states - list of states for one multiplicity only.
        singlets/ triplets - specifies which spin is being setup.
        noexc - flag for NOEXC option of MRCI programme - for passing CASSCF wf through CI programme.

        Outputs:
        nacme_calculation_input - complete input of nacme calculation for one spin (str).
        displaced_axes_list - record of displaced axes."""

    nacme_calculation_input = []
    [referance_nacme, referance_records] = setup_referance_nacme(inputs, states, referance_geom, singlets, triplets,  noexc) # at non-displaced geom
    casscf = setup_casscf(inputs, states, singlets, triplets, nacme=True)
    casscf_block_idxs = [i for i, line in enumerate(casscf) if 'casscf' in line]
    casscf = '\n'.join(casscf[casscf_block_idxs[-1]:]) # Only need final projection for displaced calculation
    [nd, dx, dy, dz] = displace(inputs['symm'], len(inputs['atoms']), inputs['dr'], referance_geom)
    states = states[-1] # Only need final set after ref setup
    incremented_records = record_incrementer(2140.2, referance_records, nd, np.count_nonzero(states))
    if dx is not None:
        dx = setup_displaced_nacme(inputs, dx, referance_records, incremented_records, casscf,
                                 noexc, states, nd)
    if dy is not None:
        dy = setup_displaced_nacme(inputs, dy, referance_records, incremented_records, casscf,
                                 noexc, states, nd)
    if dz is not None:
        dz = setup_displaced_nacme(inputs, dz, referance_records, incremented_records, casscf,
                                 noexc, states, nd)
    displaced_axes_list = [True if dnacme is not None else False for dnacme in (dx, dy, dz)]

    dnacmes = ['\n'.join(dnacme) for dnacme in (dx, dy, dz) if dnacme is not None]
    nacme_calculation_input.extend(['\n'.join(referance_nacme), '\n'.join(dnacmes)])
    return '\n'.join(nacme_calculation_input), displaced_axes_list


def setup_referance_nacme(inputs: dict, states: list, referance_geom, singlets: str, triplets: str, noexc: bool) -> (list, tuple):
    # Need to modify to split into spin multis
    """Procedure for setting up referance calculation to be used in
       NACME calculation. The referance calculation must project through all
       user specified active spaces.

       Outputs:
       referance_nacme - complete input for referance calculation prior to displacing molecule.
       (ref_wf, ref_tdm) - tuples of lists of record file numbers for wf and tdm at referance geometry.
       """

    #  Initialise reference calculation where no atoms are displaced - project
    #  through all specified CAS calculations
    referance_nacme = []
    geometry_block = geometry_setup(inputs, referance_geom)
    casscf = setup_casscf(inputs, states,  singlets, triplets, nacme=False)
    referance_nacme.extend(['\n'.join(geometry_block), 'hf', '\n'.join(casscf), '\n'])
    # Set up following CI calculation for reference - do tdm calculation
    states = states[-1]
    ciwf = setup_mrci(inputs, states,  6000.2, 8000.2, noexc, tdm=True)  # ADD new state functionality
    ciwf = '\n'.join(ciwf)
    referance_nacme.append(ciwf)
    ref_wf = re.findall(r"6[0-9]{3}\.2", ciwf)   # Get list of reference wf
    ref_tdm = re.findall(r"8[0-9]{3}\.2", ciwf)  # and tdm record files
    return referance_nacme, (ref_wf, ref_tdm)


def setup_displaced_nacme(inputs: dict, geom, referance_records: tuple, incremented_records: tuple,
                        casscf: str, noexc: bool, states: list, nd: int):
    """Sets up calculations at displaced coordinates for NACME calculation.
       Here only the final active space in the preceeding referance calculation
       is required as the orbitals/CI coeffs from reference calculation are
       read into each displaced calcualtion as a starting guess.

       Inputs (additional):
       geom - displaced geometry (array).
       referance_records - contains referance wf and tdm record file numbers.
       incremented_records - contains orbital, wf and tdm record file numbers
                            for n displaced geometries.
       casscf - block of casscf input for the spin requested.
       nd - number of displacements for each axis i.e. natom*2.

       Outputs:
       displaced_calculation_input - input for every displaced geometry."""

    displaced_calculation_input = []
    ref_wf, ref_tdm = referance_records[0], referance_records[1]
    orbitals, displaced_wf, displaced_tdm = incremented_records[0], incremented_records[1], incremented_records[2]
    grouped_orbitals = split_list(orbitals, 2)  # Groups records for +ve and -ve
    grouped_tdms = split_list(displaced_tdm, 2)      # displacements - easier indexing
    count = 0
    for i in range(nd):
        dwf = displaced_wf[i]
        dtdm = displaced_tdm[i]
        dorbital = orbitals[i]
        displaced_geom = geometry_setup(inputs, geom[i, :, :])
        dcasscf = re.sub(r'orbital,2140.2', 'orbital,%s\ndiab,2140.2' % dorbital, casscf)  # displaced casscf
        displaced_calculation_input.extend(['\n', '\n'.join(displaced_geom), dcasscf, '\n'])
        dciwf = setup_mrci(inputs, states, dwf[0], dtdm[0], noexc, tdm=False)  # displaced CI
        displaced_calculation_input.extend(['\n'.join(dciwf), '\n'])
        citdm = setup_ci_tdm(ref_wf, dwf, dtdm, np.count_nonzero(states))  # Overlap of displaced w/ reference
        displaced_calculation_input.extend(['\n'.join(citdm), '\n'])
        couplings = state_couplings(states)  # Get all possible coupling combinations

        if i % 2 != 0:
            ddr = setup_ddr(inputs['dr'], 2140.2, ref_tdm, grouped_orbitals[count], grouped_tdms[count], states, couplings)  # NACME for each combination
            displaced_calculation_input.append('\n'.join(ddr))
            count += 1

    return displaced_calculation_input



def setup_gradient(inputs: dict, geom, workdir: str, grid_index: int):
    """Sets up gradient calculation, gradients are calculated for all states.
       If a CASSCF or CASPT2 SPE energy calculation preceeds then analytical
       gradients are calculated using the rs2 programme. For MRCI numerical
       gradientds are used."""

    try:
        grad_dir = "%s/GRAD/" % workdir
        os.mkdir(grad_dir)
    except FileExistsError:
        raise FileTreeExistsError(grad_dir)
    except:
        #setup.Logging(logfile, message="Can not create GRAD Directory.")
        raise OSError('Can not create GRAD directory')

    gradient_calculation_input = []
    singlets = inputs['singlets']
    triplets = inputs['triplets']
    states = inputs['states']
    inputfile = "GRAD_GP%s.input" % grid_index
    init_calculation = "***GRAD Calculation\nmemory,%s\nprint,orbitals,civectors;\n" % (inputs['mem'])
    geometry_block = geometry_setup(inputs, geom)
    gradient_calculation_input.extend([init_calculation, '\n'.join(geometry_block), 'hf'])
    if inputs['spe'] != 'mrci':
        if inputs['singlets'] and inputs['triplets'] == 'yes':
            singlet_states = format_states(states, singlet=True, triplet=False)
            singlets, triplets = 'yes', 'no'
            singlet_casscf = setup_casscf(inputs, singlet_states, singlets, triplets, nacme=False)
            singlet_grads = setup_caspt2(inputs, singlet_states[-1], grad=True) # edit
            triplet_states = format_states(states, singlet=False, triplet=True)
            singlets, triplets = 'no', 'yes'
            triplet_casscf = setup_casscf(inputs, triplet_states, singlets, triplets, nacme=False)
            triplet_grads = setup_caspt2(inputs, triplet_states[-1], grad=True) # edit too
            gradient_calculation_input.extend(['\n'.join(singlet_casscf), '\n'.join(singlet_grads), '\n'.join(triplet_casscf), '\n'.join(triplet_grads)])
        else:
            casscf = setup_casscf(inputs, states, singlets, triplets, nacme=False)
            grads = setup_caspt2(inputs, states[-1], grad=True)
            gradient_calculation_input.extend(['\n'.join(casscf), '\n'.join(grads)])

    elif inputs['spe'] == 'mrci':
        pass  # TO DO - Numerical grads

    f = open(grad_dir+inputfile, 'w+')
    f.write('\n'.join(gradient_calculation_input))
    f.close()
    return grad_dir, inputfile


###########################
# NACME SETUP UTILITIES #
###########################


def record_incrementer(orbital: float, referance_records: tuple, nd: int, nblocks: int) -> tuple:
    """Given a list of referance records for each symm/ multiplicity - increment
       the record number for n displacements. Returns lists of tuples, each tuple
       corrosponds to records for one displacement - elements of tuples are
       the symmetry/ multiplicity.

       Inputs (additional):
       orbital - Referance orbital record to start incrementing from.
       referance_records - tuple of lists containing ref wf and tdm record file numbers.
       nblocks - number symmetry/ spin blocks in each displaced calculation.

       Outputs:
       (orbitals, displaced_wf, displaced_tdm) - tuple of lists of record file numbers for displaced calculation."""

    displaced_wf, displaced_tdm, orbitals = [], [], []
    ref_wf, ref_tdm = referance_records[0], referance_records[1]
    wf = float(ref_wf[-1])
    tdm = float(ref_tdm[-1])
    orbital = float(orbital)

    for disp in range(nd):      # Loop over num. displacements
        orbital = orbital + 1   # For each displacement increment orbital rec by one
        orbitals.append(orbital)
        temp_wf, temp_tdm = [], []
        for irrep in range(nblocks):  # Loop over symm/ spin blocks
            wf = wf + 1
            tdm = tdm + 1             # Increment each wf/ tdm rec by one for each symm/ spin
            temp_wf.append(wf)
            temp_tdm.append(tdm)
        displaced_wf.append(tuple(temp_wf))   # symm/ spin blocks for each displacement grouped into tuples
        displaced_tdm.append(tuple(temp_tdm))
    return (orbitals, displaced_wf, displaced_tdm)


def state_couplings(states: list) -> list:
    """Determines the possible combinations of coupled states for a given
       symmetry/ multiplicity.

       Outputs:
       coupling_combinations - list of lists, each sublist corrosponds to
                              a different symmetry/spin. Each of which contains
                              N tuples of coupling combinations."""

    coupling_combinations = []
    for i, nstate in enumerate(states):
        temp = []
        if i % 2 == 0:  # Even elements corrspond to A' symm
            for j in range(1, nstate+1):
                for k in range(1, j+1):  # Get lower triangle of matrix only
                    if j != k:
                        temp.append((j+0.1, k+0.1))
            coupling_combinations.append(temp)

        else:           # Odd elements corrospond to A" symm
            for j in range(1, nstate+1):
                for k in range(1, j+1):
                    if j != k:
                        temp.append((j+0.2, k+0.2))
            coupling_combinations.append(temp)
    return coupling_combinations


def displace(uaxis: str, natom: int, dr: float, referance_geom) -> (int, list, list, list):
    """Generates geometries displaced by an infinitesimal ammount dr along
       the x, y, and z axis in cartesian coordinates. Returns three np.arrays
       with dimensions [num. displacements, 3, natom] for each axis
       and the number of displacements.

       Inputs (additional):
       uaxis - unique axis under symmetry generator specified by user.
       natom - number of atoms.
       dr - amount to displace geometry by along each axis in +ve and -ve direction.

       Outputs:
       displacements_per_atom - num. displacements along each axis.
       dx (or dy or dz) - list of arrays containing displaced geometries along
                          each axis, one will be None if symmetry is on.
       """

    displacements_per_atom = 6  # num axis (3) * 2
    displaced_geoms = np.zeros((displacements_per_atom*natom,  3, natom))
    c = 0

    for i in range(3):  # Loop over axes
        for j in range(natom):  # Loop over atoms
            pos_displacement = referance_geom.copy()  # positive and negative
            neg_displacement = referance_geom.copy()  # displacement for each element
            pos_displacement[j, i] = pos_displacement[j, i] + dr
            neg_displacement[j, i] = neg_displacement[j, i] - dr
            displaced_geoms[c, :, :] = pos_displacement
            displaced_geoms[c+1, :, :] = neg_displacement
            c += 2

    dx, dy, dz = np.split(displaced_geoms, 3)
    if uaxis == 'x':  # Do not want to displace along unique symmetry axis
        dx = None
    elif uaxis == 'y':
        dy = None
    elif uaxis == 'z':
        dz = None
    return displacements_per_atom, dx, dy, dz


def setup_ci_tdm(ref_wf: list, wfs: list, tdms: list, nblocks: int) -> list:
    """ This function sets up a transition dipole moment calculation using
        the MRCI programme in Molpro.

        Inputs:
        ref_wf - record file numbers for wf at reference geometry for each spin/symm..
        wfs - record file numbers for wf at displaced geometry for each spin/symm.
        tdms - record file numbers for tdm at displaced geometry for each spin/symm.
        nblocks - num. symmetry/ spin blocks in each calculation.

        Outputs:
        ci_tdm_calculation - input for tdm calculation for every spin/symmetry using
        MRCI programme."""

    ci_tdm_calculation = []
    for i in range(nblocks):  # Loop over symm/ spin blocks
        citdm = "{ci;trans,%s,%s;dm,%s}" % (ref_wf[i], wfs[i], tdms[i])
        ci_tdm_calculation.append(citdm)
    return ci_tdm_calculation


def setup_ddr(dr: float, ref_orbs: float, ref_tdm: list, disp_orbs: list, disp_tdm: list, states: list, coupling_combinations: list) -> list:
    """Sets up nacme computation using finite differences. Loops over all
       states of each symm/spin and generates a ddr calculation for each
       possible combination of states in each symm/spin.

       Inputs:
       dr - amount geometry is displaced by in +ve and -ve direction.
       ref_orbs - orbital record file number at referance geometry.
       ref_tdm - tdm record file number at referance geometry for each symm/ spin.
       disp_orbs - orbital record file number at +ve/-ve displacement for every
                     symm/ spin.
       disp_tdm - tdm record file number at each displacement for each symm/ spin.
       states - list of states for each symm/ spin.
       coupling_combinations - every possible state coupling.

       Output:
       ddr_calculation - input for 3-point DDR calculation for every state combination."""

    ddr_calculation = []
    count = 0
    for i in range(len(states)):  # Loop over symm/ spin blocks
        if states[i] != 0:
            for coupled_pair in coupling_combinations[i]:
                ddr = "{ddr,2*%s\norbital,%s,%s,%s;\ndensity,%s,%s,%s;\nstates,%s,%s}" % (dr,
                        ref_orbs, disp_orbs[0], disp_orbs[1], ref_tdm[count], disp_tdm[0][count], disp_tdm[1][count], coupled_pair[0], coupled_pair[1])
                ddr_calculation.append(ddr)
            count += 1
    return ddr_calculation


##############################
# Ab-Initio Methods Setup    #
##############################


def geometry_setup(inputs: dict, geom) -> list:
    """ Sets up the geometry block for a calculation given inputs.
        Block includes basis, symmetry generators and xyz coordinates"""

    geometry_block = []
    basis = "basis=%s" % inputs['basis']

    if inputs['symm'] != 'nosymmetry':
        symmetry = "symmetry,%s" % inputs['symm']
    else:
        symmetry = "symmetry,nosymmetry"

    geom = format_geometry(geom, inputs['atoms'])  # Format array into coords w/ atomic symbols in column 0
    geometry = "geomtyp=xyz\nbohr\ngeometry={\n%s\nSurf Calc\n%s\n}" % (len(inputs['atoms']), geom)
    geometry_block.extend((basis, symmetry, geometry))
    return geometry_block


def setup_casscf(inputs: dict, total_states: list, singlets: str, triplets: str, nacme: bool) -> list:
    """Sets up N single point calculations at the CASSCF level given N active spaces
       and states. Works for any combination of symmetry and multiplicities.

       Inputs:
       inputs - user defined inputs from json file.
       total_states - List of states for each spin/symm.
       singlets/ triplets - yes/no str to decide which wfs to use.
       nacme - bool to add reading of orbitals from ref calculation.

       Outputs:
       casscf_calculation_input - block of casscf input for each symm/spin."""

    casscf_calculation_input = []
    if not nacme:
        init_singlet_wf, init_triplet_wf = '', ''
    if nacme:
        init_singlet_wf, init_triplet_wf = '\nstart,2140.2', '\nstart,2141.2'

    if inputs['symm'] != 'nosymmetry':
        for i in range(len(total_states)):  # Loop over projections
            occ = inputs['occ'][i]
            closed = inputs['closed'][i]    # Init active spaces and num. states
            states = total_states[i]        # for each projection
            singlet_casscf = "{casscf\nfrozen,0\nclosed,%s,%s\nocc,%s,%s%s\norbital,2140.2\nconfig,csf" % (
                        closed[0], closed[1], occ[0], occ[1], init_singlet_wf)
            triplet_casscf = "{casscf\nfrozen,0\nclosed,%s,%s\nocc,%s,%s%s\norbital,2141.2\nconfig,csf" % (
                        closed[0], closed[1], occ[0], occ[1], init_triplet_wf)
            if singlets == 'yes' and triplets == 'yes':
                wf1 = "wf,%s,1,0\nstate,%s" % (inputs['nelec'], states[0])  # Init wf depending on states
                wf2 = "wf,%s,2,0\nstate,%s\n};" % (inputs['nelec'], states[1])
                wf3 = "wf,%s,1,2\nstate,%s" % (inputs['nelec'], states[2])
                wf4 = "wf,%s,2,2\nstate,%s\n};" % (inputs['nelec'], states[3])
                casscf_calculation_input.extend([singlet_casscf, wf1, wf2, triplet_casscf, wf3, wf4])  # Do SA over different spins
            elif singlets == 'yes' and triplets == 'no':
                wf1 = "wf,%s,1,0\nstate,%s" % (inputs['nelec'], states[0])
                wf2 = "wf,%s,2,0\nstate,%s\n};" % (inputs['nelec'], states[1])
                casscf_calculation_input.extend([singlet_casscf, wf1, wf2])
            elif singlets == 'no' and triplets == 'yes':
                wf3 = "wf,%s,1,2\nstate,%s" % (inputs['nelec'], states[2])
                wf4 = "wf,%s,2,2\nstate,%s\n};" % (inputs['nelec'], states[3])
                casscf_calculation_input.extend([triplet_casscf, wf3, wf4])

    elif inputs['symm'] == 'nosymmetry':
        for i in range(len(total_states)):
            occ = inputs['occ'][i]
            closed = inputs['closed'][i]
            states = total_states[i]
            singlet_casscf = "{casscf\nfrozen,0\nclosed,%s\nocc,%s%s\norbital,2140.2\nconfig,csf" % (
                      closed[0], occ[0], init_singlet_wf)
            triplet_casscf = "{casscf\nfrozen,0\nclosed,%s\nocc,%s%s\norbital,2141.2\nconfig,csf" % (
                      closed[0], occ[0], init_triplet_wf)
            if singlets == 'yes' and triplets == 'yes':
                wf1 = "wf,%s,1,0\nstate,%s};" % (inputs['nelec'], states[0])  # Only one irrep for each spin
                wf2 = "wf,%s,1,2\nstate,%s};" % (inputs['nelec'], states[2])  # if no symm is chosen
                casscf_calculation_input.extend([singlet_casscf, wf1, triplet_casscf, wf2])
            elif singlets == 'yes' and triplets == 'no':
                wf1 = "wf,%s,1,0\nstate,%s};" % (inputs['nelec'], states[0])
                casscf_calculation_input.extend([singlet_casscf, wf1])
            elif singlets == 'no' and triplets == 'yes':
                wf2 = "wf,%s,1,2\nstate,%s};" % (inputs['nelec'], states[2])
                casscf_calculation_input.extend([triplet_casscf, wf2])

    return casscf_calculation_input


def setup_mrci(inputs: dict, states: list,  orbital_record: float, tdm_record: float, noexc: bool, tdm: bool) -> list:
    """Sets up an MRCI calculation (based on final CASSCF projection)
       for each symmetry/ multiplicity provided.

       Inputs:
       orbital_record - which record file to save mrci wf to.
       tdm - bool that specifies a tdm calculation using ci programme.
       tdm_record - in event of tdm calculation - which record to save tdm to.
       noexc - bool that specifies to just pass a CASSCF wf through CI programme.

       Outputs:
       mrci_calculation_input - input for mrci calculation block."""

    mrci_calculation_input = []
    wfs = ["wf,%s,1,0" % inputs['nelec'], "wf,%s,2,0" % inputs['nelec'],
           "wf,%s,1,2" % inputs['nelec'], "wf,%s,2,2" % inputs['nelec']]  # wfs corrosponding to list of states

    for i in range(len(states)):
        if states[i] != 0:  # Loop over symm/ spin blocks
            if not noexc:  # Include excitations - run MRCI
                if not tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s};" % (wfs[i], orbital_record, states[i])
                if tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s\ndm,%s};" % (wfs[i], orbital_record, states[i], tdm_record)

            else:          # Exclude excitations - pass CAS wf through CI programme
                if not tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s;NOEXC};" % (wfs[i], orbital_record, states[i])
                if tdm:
                    mrci = "{ci;maxit,100,200\n%s;save,%s\nstate,%s;NOEXC\ndm,%s};" % (wfs[i], orbital_record, states[i], tdm_record)

            orbital_record += 1  # Increment records
            tdm_record += 1
            mrci_calculation_input.append(mrci)
    return mrci_calculation_input


def setup_caspt2(inputs: dict, states: list, grad: bool) -> list:
    """Sets up N XMS-CASPT2 calculations for each symmetry/ multiplicity
       with referance to the final CASSCF projection. The level shift is
       specified by user in input.json. If grad == True then an analytical
       gradient calculation is set up instead."""

    caspt2_calculation_input = []
    wfs = ["wf,%s,1,0" % inputs['nelec'], "wf,%s,2,0" % inputs['nelec'],  # Corrosponding wfs
           "wf,%s,1,2" % inputs['nelec'], "wf,%s,2,2" % inputs['nelec']]

    for i in range(len(states)):  # Loop over symm/ spin
        if states[i] != 0:
            if not grad:
                caspt2 = "{rs2,shift=%s,MIX=%s,XMS=1\n%s\nstate,%s};" % (inputs['level_shift'], states[i], wfs[i], states[i])
                caspt2_calculation_input.append(caspt2)
            else:
                for root in range(states[i]):
                    caspt2 = "{rs2;NOEXC\n%s\nOPTION,NSTATI=%s;state,1,%s};\nforces" % (wfs[i], states[i], root+1)
                    caspt2_calculation_input.append(caspt2)
    return caspt2_calculation_input


def format_geometry(xyz_geom, atoms: list) -> str:
    """Turns coorindates in a np.array into a string with correct format. """
    geom_list = [",\t".join(row) for row in xyz_geom.astype(str)]
    geom_list = [atoms[i] + ",\t" + row for i, row in enumerate(geom_list)]
    geom = '\n'.join(geom_list)
    return geom


def split_list(alist, chunks):
    """Splits a list into a list of lists with equal no. elements (chunks)"""
    new_list = []
    for i in range(0, len(alist), chunks):
        new_list.append(alist[i:i+chunks])
    return new_list


def format_states(states: list, singlet: bool, triplet: bool):
    """Edits list of states for both multiplicties to set one
       multiplicity to zeros - for splitting SA of each spin."""
    new_states = []
    if singlet:
        for block in states:
            new_block = block[:2]
            new_block.extend([0, 0])
            new_states.append(new_block)
    elif triplet:
        for block in states:
            new_block = [0, 0]
            new_block.extend(block[2:])
            new_states.append(new_block)
    return new_states

