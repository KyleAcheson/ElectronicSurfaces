import os
import numpy as np
import re
import multiprocessing.managers
from io import StringIO

""" This module contains all functions/ methods required to set up workspaces
    for either user specified quantum code, and to intialise matrices for
    storing data, as well as extracting and sorting data from
    each type of calculation"""

ang2au = 1.88973
flatten = lambda nestedList: [item for sublist in nestedList for item in sublist]


class ProccessManager(multiprocessing.managers.BaseManager):
    ''' Custom subclass of BaseManager to share objects between processes. '''
    pass


ProccessManager.register('np_zeros', np.zeros, multiprocessing.managers.ArrayProxy)


class datastore:
    """Global datastoring class for Energies, NACMES, SOCs, Gradients and
       CASPT2 Mixing Matrix. This class is initialised prior to the parallel run.
       Memory sharing between processes is forced for required
       instances using ProcessManager."""

    def __init__(self, refGeom, states: list, couplings: list, ngrid: int, pmanager):
        self.refGeom = refGeom
        self.ngrid = ngrid
        self.Natm = len(self.refGeom)
        self.nsinglets = sum(states[0:2])
        self.ntriplets = sum(states[2:4])
        self.nmultiplets = self.nsinglets + (3*self.ntriplets)
        self.nstates = sum(states)
        self.energies = pmanager.np_zeros((ngrid, self.nstates))
        self.couplings = flatten(couplings)
        self.nacmes = pmanager.np_zeros((self.Natm, 3, len(self.couplings), self.ngrid))
        self.grads = pmanager.np_zeros((self.Natm, 3, self.nstates, self.ngrid))
        self.socs = pmanager.np_zeros((self.nmultiplets, self.nmultiplets, 2, self.ngrid))

    def energiesExtract(self, outfile: str, programme: str, progregex: str, multi: str, gp: int):
        """Extracts energies from a single point calculation at each grid point.
        Takes an outfile path, user specified programme and corrosponding regex,
        if multi, will only extract energies of last projection. Sorts in a global
        array according to the grid point index."""
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        if programme == 'casscf':
            calculationIndex = []
            for idx, line in enumerate(lines):
                if multi in line:
                    calculationIndex.append(idx)
            lines = lines[calculationIndex[-1]:]  # Only final CASSCF projection
            energies = regex(lines, progregex)
            self.energies[gp, :] = energies
        else:
            energies = regex(lines, progregex)
            self.energies[gp, :] = energies

    def nacmeExtract(self, outfile: str, nacme_regex: str, gp: int, daxes: list):
        """Extracts NACMEs for each grid point, accorindg to each atom/ axis
           and stores in a global array"""
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        numAxesDisplaced = sum(daxes)
        displacedAxesIndex = [ind for ind, axis in enumerate(daxes) if axis]
        nacmeArrayRaw = regex(lines, nacme_regex)  # Array with N combination couplings, N axis displacments and N atoms
        nacmeAxialSort = np.split(nacmeArrayRaw, numAxesDisplaced)  # Split into 1-3 arrays depending on which axis displaced
        for count, nacmeAxis in enumerate(nacmeAxialSort):
            nacmeAtomicSort = np.split(nacmeAxis, self.Natm)  # Split into array for each atom
            ind = displacedAxesIndex[count]
            for j in range(self.Natm):
                self.nacmes[j, ind, :, gp] = nacmeAtomicSort[j]

    def gradExtractMolpro(self, outfile: str, gradRegex: str, numerical: bool, gp: int):
        """Extracts gradients for each atom along each axis as an array.
           Works for both numerical or analytical gradients in molpro. """
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        if numerical:
            gradLineIndex = [i+5 for i, line in enumerate(lines) if gradRegex in line]
        else:
            gradLineIndex = [i+4 for i, line in enumerate(lines) if gradRegex in line]
        for x, i in enumerate(gradLineIndex):
            gradMatrixList = lines[i:i+self.Natm]
            gradMatrixStr = ''.join(gradMatrixList)
            gradArray = np.genfromtxt(StringIO(gradMatrixStr), usecols=(1, 2, 3), encoding=None)
            self.grads[:, :, x, gp] = gradArray

    def extractSOC(self, outfile: str, gp: int):
        TempSOC = np.zeros((self.nmultiplets, self.nmultiplets, 2))
        f = open(outfile, 'r')
        lines = f.readlines()
        f.close()
        idx = [i for i, line in enumerate(lines) if 'Nr  State  S   SZ' in line]
        StateBlock = [int(lines[i].split()[-1]) for i in idx]
        LastIdx = idx[-1]+((self.nmultiplets*2)+(self.nmultiplets+2))
        idx.append(LastIdx)
        StateBlock.insert(0, 0)
        for j in range(len(idx[0:-1])):
            MatrixBlockRaw = lines[idx[j]+1:idx[j+1]]
            MatrixBlockRaw = [line.strip() for line in MatrixBlockRaw]
            MatrixBlock = list(filter(None, MatrixBlockRaw))
            RealValuedRaw = [line for i, line in enumerate(MatrixBlock) if i % 2 == 0]
            RealValued = [' '.join(line.split()[4:]) for line in RealValuedRaw]
            ImagValued = [line for i, line in enumerate(MatrixBlock) if i % 2 != 0]
            SOCTempReal = np.genfromtxt(StringIO('\n'.join(RealValued)), encoding=None)
            SOCTempImag = np.genfromtxt(StringIO('\n'.join(ImagValued)), encoding=None)
            TempSOC[:, StateBlock[j]:StateBlock[j+1], 0] = SOCTempReal
            TempSOC[:, StateBlock[j]:StateBlock[j+1], 1] = SOCTempImag

        self.socs[:, :, :, gp] = TempSOC

    def caspt2MatrixExtract():
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
