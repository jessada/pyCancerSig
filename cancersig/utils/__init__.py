import subprocess
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
        sys.stdout.write(stdout_data)
        if return_code:
            mylogger.throw("Error found during execute command '%s' with error code: %d, %s" % (cmd, return_code, stderr_data))
        sys.stderr.write(stderr_data)
    elif stderr_data:
        sys.stderr.write(stderr_data)
    return p, stdout_data


