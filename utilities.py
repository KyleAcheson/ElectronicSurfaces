import subprocess
import os
import time

"""This module contains routines for both quantum codes which submit the jobs,
    check their status and also any additional utilities for formatting text. """


def getPID(inputfile: str) -> int:
    ''' This function returns the pid of the job associated
        with an input file. '''
    a = subprocess.Popen(['ps', '-eo', 'pid,ppid,command'],
                         stdout=subprocess.PIPE)  # List all running procs
    b = subprocess.Popen(['grep', '%s' % (inputfile)], stdin=a.stdout,
                         stdout=subprocess.PIPE)  # Search for input job PID
    output, error = b.communicate()
    output = output.decode("utf-8").split('\n')
    pid = ''
    pid = int(pid.join(list(output[0])[0:6]))  # just pid number
    return pid


def runCalculation(codeKeys: dict, pwd: str, workdir: str, inputfile: str, index: str) -> (bool, str):
    '''This function submits a single point energy and SO calculation at
        each point in the grid. Starting in the parent directory it opens the
        sub-directory at the grid point and submits the job, waiting until
        the job PID no longer exists. '''
    index = "GP_"+str(index)
    outputfile = inputfile.split('.')[0]+'.out'
    os.chdir('%s' % (workdir))
    subprocess.Popen(['%s %s &' % (codeKeys['submit_command'], inputfile)],
                     shell=True, preexec_fn=os.setsid)  # Submit job

    pid = getPID(inputfile)
    time.sleep(1)
    while os.path.exists('/proc/%s' % (pid)) == 1:  # Wait while job running
        time.sleep(1)
    lastline = str(subprocess.check_output(['tail', '-1', outputfile]))
    os.chdir(pwd)

    if codeKeys['termination_code'] in lastline:  # Check termination
        normalTermination = True
    else:
        normalTermination = False
    return normalTermination, outputfile


def insertGeom(start, refGeomBlock, geom, codeKeys):
    '''This function takes a geometry block with some place holder geometry
        and inserts the new reference geometry for use in a nacme calculation
        in its place. Returns a list of strings that can be used as the
        input geometry block for the next calculation. '''
    newGeomBlock = []
    ind = [i for i, j in enumerate(refGeomBlock) if j.strip() == str(start)]
    newGeomBlock = refGeomBlock[0:ind[0]+2]
    geomInput = ["\t".join(row) for row in geom.astype(str)]
    for atom in geomInput:
        newGeomBlock.append(atom)
    newGeomBlock.append(codeKeys['block_end'])
    newGeomBlock = '\n'+'\n'.join(newGeomBlock)+'\n'
    return newGeomBlock


def insertBlock(start, end, block, insert):
    '''This function takes a geometry block and inserts it within
        the template input for a tdm calculation at the reference geometry -
        the first step in a nacme calculation.'''
    outputText = []
    copying = True
    for line in block:
        line = line.strip()
        if copying:
            if line.startswith('%s' % (start)):
                outputText.append(insert)
                copying = False
            else:
                outputText.append(line)
        elif line.startswith('%s' % (end)):
            copying = True
    outputText = '\n'.join(outputText)+'\n'
    return outputText
