import os
import numpy as np
import re
import multiprocessing.managers
from io import StringIO

""" This module contains all functions/ methods required to set up workspaces
    for either user specified quantum code, and to intialise matrices for
    storing data, as well as extracting and sorting data from
    each type of calculation"""

ANG2AU = 1.88973
flatten = lambda nestedList: [item for sublist in nestedList for item in sublist]


class ProccessManager(multiprocessing.managers.BaseManager):
    ''' Custom subclass of BaseManager to share objects between processes. '''
    pass


ProccessManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)


class DataStore:
    """Global datastoring class for Energies, NACMES, SOCs, Gradients and
       CASPT2 Mixing Matrix. This class is initialised prior to the parallel run.
       Memory sharing between processes is forced for required
       instances using ProcessManager."""

    def __init__(self, natom, states: list, couplings: list, ngrid: int, pmanager):
        self.ngrid = ngrid
        self.natom = natom
        self.nsinglets = sum(states[0:2])
        self.ntriplets = sum(states[2:4])
        self.nmultiplets = self.nsinglets + (3*self.ntriplets)
        self.nstates = sum(states)
        self.energies = pmanager.np_zeros((ngrid, self.nstates))
        self.couplings = flatten(couplings)
        self.nacmes = pmanager.np_zeros((self.natom, 3, len(self.couplings), self.ngrid))
        self.grads = pmanager.np_zeros((self.natom, 3, self.nstates, self.ngrid))
        self.socs = pmanager.np_zeros((self.nmultiplets, self.nmultiplets, 2, self.ngrid))

    def extract_energies(self, inputs: dict,  outfile: str, progregex: str, multi: str, grid_index: int):
        """Extracts energies from a single point calculation at each grid point.
        Takes an outfile path, user specified programme and corrosponding regex,
        if multi, will only extract energies of last projection. Sorts in a global
        array according to the grid point index."""
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        if inputs['spe'] == 'casscf':
            if inputs['singlets'] == 'yes' and inputs['triplets'] == 'yes':
                calculation_indexes = [idx for idx, line in enumerate(lines) if multi in line]
                lines = lines[calculation_indexes[-2]:]
            else:
                calculation_indexes = [idx for idx, line in enumerate(lines) if multi in line]
                lines = lines[calculation_indexs[-1]:]  # Only final CASSCF projection
            energies = regex(lines, progregex)
            self.energies[grid_index, :] = energies
        else:
            energies = regex(lines, progregex)
            self.energies[grid_index, :] = energies

    def extract_nacme(self, outfile: str, nacme_regex: str, grid_index: int, daxes: list):
        """Extracts NACMEs for each grid point, accorindg to each atom/ axis
           and stores in a global array"""
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        ndisplaced_axes = sum(daxes)
        displaced_axes_index = [ind for ind, axis in enumerate(daxes) if axis]
        nacme_array_raw = regex(lines, nacme_regex)  # Array with N combination couplings, N axis displacments and N atoms
        nacme_axial_sort = np.split(nacme_array_raw, ndisplaced_axes)  # Split into 1-3 arrays depending on which axis displaced
        for count, nacme_axis in enumerate(nacme_axial_sort):
            nacme_atomic_sort = np.split(nacme_axis, self.natom)  # Split into array for each atom
            ind = displaced_axes_index[count]
            for j in range(self.natom):
                self.nacmes[j, ind, :, grid_index] = nacme_atomic_sort[j]

    def extract_grad_molpro(self, outfile: str, grad_regex: str, numerical: bool, grid_index: int):
        """Extracts gradients for each atom along each axis as an array.
           Works for both numerical or analytical gradients in molpro. """
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        if numerical:
            grad_line_index = [i+5 for i, line in enumerate(lines) if grad_regex in line]
        else:
            grad_line_index = [i+4 for i, line in enumerate(lines) if grad_regex in line]
        for x, i in enumerate(grad_line_index):
            grad_matrix_list = lines[i:i+self.natom]
            grad_matrix_str = ''.join(grad_matrix_list)
            grad_array = np.genfromtxt(StringIO(grad_matrix_str), usecols=(1, 2, 3), encoding=None)
            self.grads[:, :, x, grid_index] = grad_array

    def extract_soc(self, outfile: str, grid_index: int):
        temp_soc = np.zeros((self.nmultiplets, self.nmultiplets, 2))
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        idx = [i for i, line in enumerate(lines) if 'Nr  State  S   SZ' in line]
        state_block = [int(lines[i].split()[-1]) for i in idx]
        last_index = idx[-1]+((self.nmultiplets*2)+(self.nmultiplets+2))
        idx.append(last_index)
        state_block.insert(0, 0)
        for j in range(len(idx[0:-1])):
            matrix_block_raw = lines[idx[j]+1:idx[j+1]]
            matrix_block_raw = [line.strip() for line in matrix_block_raw]
            matrix_block = list(filter(None, matrix_block_raw))
            real_valued_raw = [line for i, line in enumerate(matrix_block) if i % 2 == 0]
            real_valued = [' '.join(line.split()[4:]) for line in real_valued_raw]
            imag_valued = [line for i, line in enumerate(matrix_block) if i % 2 != 0]
            real_soc_temp = np.genfromtxt(StringIO('\n'.join(real_valued)), encoding=None)
            imag_soc_temp = np.genfromtxt(StringIO('\n'.join(imag_valued)), encoding=None)
            temp_soc[:, state_block[j]:state_block[j+1], 0] = real_soc_temp
            temp_soc[:, state_block[j]:state_block[j+1], 1] = imag_soc_temp

        self.socs[:, :, :, grid_index] = temp_soc

    def extract_mixing_matrix():
        pass


def regex(lines: list, progregex: str) -> list:
    """Searches through list of lines for a match to a specific method,
       returning the energies/ NACMEs for each local grid point."""
    values = []
    for line in lines:
        line = line.strip()
        match = re.findall(progregex, line)
        if match:
            match = match[-1]
            value = re.split(r'\s', str(match))[-1]
            values.append(float(value))

    return np.array(values)


def split_list(alist, chunks):
    """Splits a list into a list of lists with equal no. elements (chunks)"""
    new_list = []
    for i in range(0, len(alist), chunks):
        new_list.append(alist[i:i+chunks])
    return new_list


def logging(logfile, message):
    with open(logfile, 'w+') as f:
        f.write(message)


def extract_energies_test(inputs: dict,  outfile: str, progregex: str, multi: str):
    """Extracts energies from a single point calculation at each grid point.
    Takes an outfile path, user specified programme and corrosponding regex,
    if multi, will only extract energies of last projection. Sorts in a global
    array according to the grid point index."""
    f = open(outfile, 'r')
    lines = f.readlines()
    f.close()
    Energies = np.zeros((1, 19))
    if inputs['spe'] == 'casscf':
        if inputs['singlets'] == 'yes' and inputs['triplets'] == 'yes':
            calculation_indexes = [idx for idx, line in enumerate(lines) if multi in line]
            lines = lines[calculation_indexes[-2]:]
        else:
            calculation_indexes = [idx for idx, line in enumerate(lines) if multi in line]
            lines = lines[calculation_indexes[-1]:]  # Only final CASSCF projection
        energies = regex(lines, progregex)
        Energies[0, :] = energies
    else:
        energies = regex(lines, progregex)
        Energies[0, :] = energies
    return Energies


def extract_soc_test(outfile: str, grid_index: int):
    SOC = np.zeros((39, 39, 2))
    temp_soc = np.zeros((39, 39, 2))
    f = open(outfile, 'r')
    lines = f.readlines()
    f.close()
    idx = [i for i, line in enumerate(lines) if 'Nr  State  S   SZ' in line]
    state_block = [int(lines[i].split()[-1]) for i in idx]
    last_index = idx[-1]+((39*2)+(39+2))
    idx.append(last_index)
    state_block.insert(0, 0)
    for j in range(len(idx[0:-1])):
        matrix_block_raw = lines[idx[j]+1:idx[j+1]]
        matrix_block_raw = [line.strip() for line in matrix_block_raw]
        matrix_block = list(filter(None, matrix_block_raw))
        real_valued_raw = [line for i, line in enumerate(matrix_block) if i % 2 == 0]
        real_valued = [' '.join(line.split()[4:]) for line in real_valued_raw]
        imag_valued = [line for i, line in enumerate(matrix_block) if i % 2 != 0]
        real_soc_temp = np.genfromtxt(StringIO('\n'.join(real_valued)), encoding=None)
        imag_soc_temp = np.genfromtxt(StringIO('\n'.join(imag_valued)), encoding=None)
        temp_soc[:, state_block[j]:state_block[j+1], 0] = real_soc_temp
        temp_soc[:, state_block[j]:state_block[j+1], 1] = imag_soc_temp
    SOC[:, :, :] = temp_soc
    return SOC
