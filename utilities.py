import subprocess
import os
import time
import re

''' This module contains routines for both quantum codes which submit the jobs,
    check their status and wait for termination. Compatible with any local linux
    system, and Sun Grid Engine or PBS Pro HPC systems. '''


def runCalculation(system: str, codeKeys: dict, pwd: str, workdir: str, inputfile: str, submitscript: str, index: int):
    '''Submits calculations to either a local linux OS, or a HPC running eithers Sun Grid Engine or PBS Pro.
       Waits for calculation to finish running and then returns True if it terminated with no errors. '''
    outputfile = inputfile.split('.')[0]+'.out'
    os.chdir('%s' % (workdir))

    if system == 'local':  #  Run calculation and wait for termination in each case
        runLocal(codeKeys, inputfile)
    elif system == 'sun grid engine':
        gridEngineScript = setup_submit_script(submitscript, inputfile, index)
        qsub(gridEngineScript)
    elif system == 'pbs':
        pass

    terminationCode = calculationTermination(codeKeys, outputfile)  #  check normal termination
    os.chdir(pwd)
    return terminationCode, outputfile


def calculationTermination(codeKeys, outputfile):
    ''' Check final line of output for termination code. '''
    lastline = str(subprocess.check_output(['tail', '-1', outputfile]))
    if codeKeys['termination_code'] in lastline:
        normalTermination = True
    else:
        normalTermination = False
    return normalTermination

############################
# Setup for local linux OS #
############################


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


def runLocal(codeKeys: dict, inputfile: str):
    '''This function submits a single point energy and SO calculation at
        each point in the grid. Starting in the parent directory it opens the
        sub-directory at the grid point and submits the job, waiting until
        the job PID no longer exists. '''

    subprocess.Popen(['%s %s &' % (codeKeys['submit_command'], inputfile)],
                     shell=True, preexec_fn=os.setsid)  # Submit job
    pid = getPID(inputfile)
    time.sleep(1)
    while os.path.exists('/proc/%s' % (pid)) == 1:  # Wait while job running
        time.sleep(3)

#############################
# Setup for Sun Grid Engine #
#############################


def subprocess_cmd(command, return_stdout):
    ''' Returns output for a shell command. '''
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    proc_stdout = process.communicate()[0].strip()
    if return_stdout:
        return(proc_stdout)


def setup_submit_script(submitScript: str, inputfile: str, index: int):
    ''' Set up input for Grid Engine submission script.
        Job name follows the form [calculation type][grid point index]. '''
    subjobFile = "submit_GP%s.sh" % index
    inputid = inputfile.split('.')[0]
    submitScriptGP = re.sub(r'template', inputid, submitScript)
    f = open(subjobFile, 'w+')
    f.write(submitScriptGP)
    f.close()
    return subjobFile


def submit_qsub_job(submitScript, params='-j y'):
    ''' Submits job to Grid Engine. Returns output. '''
    qsub_command = 'qsub %s %s' % (params, submitScript)
    proc_stdout = subprocess_cmd(command=qsub_command, return_stdout=True)
    return(proc_stdout)


def get_qsub_job_ID_name(proc_stdout):
    ''' Returns a tuple of job id and name.
        proc_stdout = 'Your job <id> (<name>) has been submitted. '''
    proc_stdout_list = proc_stdout.split()
    job_id = proc_stdout_list[2]
    job_name = proc_stdout_list[3]
    job_name = re.sub(r'^\("', '', str(job_name))
    job_name = re.sub(r'"\)$', '', str(job_name))
    return((job_id, job_name))


def check_qsub_job_status(job_id, desired_status="r"):
    ''' Use qstat to check run status of a job. Return True if status
        matches desired states for job, else False. '''
    job_id_pattern = r"^.*%s.*\s%s\s.*$" % (job_id, desired_status)
    qstat_stdout = subprocess_cmd('qstat', return_stdout=True)
    job_match = re.findall(str(job_id_pattern), str(qstat_stdout), re.MULTILINE)
    job_status = bool(job_match)
    if job_status is True:
        status = True
        return(status)
    elif job_status is False:
        return(job_status)


def wait_qsub_job_start(job_id, job_name, return_True=False):
    ''' Wait for job to start, print when job has started running. '''
    while check_qsub_job_status(job_id=job_id, desired_status="r") is not True:
        time.sleep(3)
    if check_qsub_job_status(job_id=job_id, desired_status="r") is True:
        print('job %s has started' % job_name)
        if return_True is True:
            return(True)


def check_qsub_job_absent(job_id):
    ''' Check if job is absent from qstat output. When finished running,
        the job id will be absent from the list. '''
    qstat_stdout = subprocess_cmd('qstat', return_stdout=True)
    job_id_pattern = r"^.*{0}.*$".format(job_id)
    job_match = re.findall(str(job_id_pattern), str(qstat_stdout), re.MULTILINE)
    job_status = bool(job_match)
    if job_status is True:
        return(False)
    elif job_status is False:
        return(True)


def check_grid_termination(job_id, job_name):
    ''' While a job is not absent from qstat output, wait.'''
    while not check_qsub_job_absent(job_id):
        time.sleep(3)
    if check_qsub_job_absent(job_id):
        print('job %s has finished' % job_name)


def qsub(submitScript):
    ''' Submit a job using a submission script, wait for job to start
        and for terminaton. '''
    proc_stdout = submit_qsub_job(submitScript)
    job_id, job_name = get_qsub_job_ID_name(proc_stdout)
    wait_qsub_job_start(job_id, job_name)
    check_grid_termination(job_id, job_name)

#####################
# Setup for PBS Pro #
#####################


def pbs():
    pass
