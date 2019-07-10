import subprocess
import sys
from cancersig.utils import logger

def exec_sh(cmd, silent=False):
    logger.dbg("executing: " + repr(cmd))
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         )
    stdout_data, stderr_data = p.communicate()
    return_code = p.returncode
    if not silent:
        sys.stdout.write(stdout_data.decode("utf-8"))
    if return_code:
        logger.throw("Error found during execute command '%s' with error code: %d, %s" % (cmd, return_code, stderr_data.decode("utf-8")))
        raise
    sys.stderr.write(stderr_data.decode("utf-8"))
    return p, stdout_data


