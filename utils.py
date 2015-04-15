import os
import datetime
import logging
import subprocess
import select

logger = logging.getLogger(__name__)


def touch(fname, times=None):
    """
Similar to nix command, touch, to create an empty file, but also added some
meta data to the touched file
"""
    with open(fname, 'a') as opf:
        opf.write('created: {0}\n'.format(datetime.datetime.now()))
        opf.write('location of code execution: {0}\n'.format(
            os.path.abspath('.')))
        os.utime(fname, times)


def execute(cmd, msg_id='', flag=None, debug=False):
    """
    This execute logs all stdout and stderr, which could look funny, especially
    when it comes to tools like aspc and wget
    
    :param flag: a flage file to indicate whether the process succeed or not.
    """
    logger.info(cmd)
    if debug:
        return
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True,
                                executable="/bin/bash")

        # stdout, stderr = [], []
        while True:
            reads = [proc.stdout.fileno(), proc.stderr.fileno()]
            ret = select.select(reads, [], [])

            for fd in ret[0]:
                if fd == proc.stdout.fileno():
                    read = proc.stdout.readline()
                    logger.info('stdout: ' + read.strip())
                    # stdout.append(read)
                if fd == proc.stderr.fileno():
                    read = proc.stderr.readline()
                    logger.info('stderr: ' + read.strip())
                    # stderr.append(read)

            if proc.poll() != None:
                break

        returncode = proc.returncode

        if returncode != 0:
            logger.error(
                '{0}, started, but then failed with returncode: {1}\n\t'
                'CMD "{2}"'.format(msg_id, returncode, cmd))
        else:
            logger.info('{0}, execution succeeded with returncode: {1}\t\n'
                        'CMD "{2}"'.format(msg_id, returncode, cmd))
            if flag is not None:
                touch(flag)
        return returncode
    except OSError as err:
        logger.exception(
            '{0}, failed to start, raising OSError {1}.\n\t'
            'CMD: "{2}"'.format(msg_id, err, cmd))
