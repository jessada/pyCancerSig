from cancersig.utils import logger
import pkg_resources
import sys
import getpass
import socket
from collections import OrderedDict
from collections import defaultdict

param_display_fmt = "  {name:<50}{value}"

def center_txt(txt, width):
    txt = " " + txt + " "
    return txt.center(width, "*")

def new_section_txt(txt):
    logger.info("")
    logger.info(center_txt(txt, 140))

def disp_header(header_txt, new_line=True):
    if new_line:
        logger.info("")
    logger.info(header_txt)

def disp_subheader(subheader_txt):
    logger.info("  " + subheader_txt)

def disp_param(param_name, param_value):
    logger.info(param_display_fmt.format(name=param_name+": ",
                                           value=param_value))

def debug_param(param_name, param_value):
    logger.debug(param_display_fmt.format(name=param_name+": ",
                                            value=param_value))

def disp_subparam(subparam_name, subparam_value):
    disp_param("  "+subparam_name, subparam_value)

def show_config(app_description,
                required_params,
                optional_params,
                third_party_software_version=None,
                ):
    disp_header("Description")
    logger.info("  " + app_description)
    disp_header("Version and environment configuration")
    disp_param("pyCancerSig version", pkg_resources.get_distribution("pyCancerSig").version)
    disp_param("executed command", " ".join(sys.argv))
    disp_param("hostname", socket.gethostname())
    disp_param("user", getpass.getuser())
    disp_params_set("Third party software version", third_party_software_version)
    if (required_params is not None) and (len(required_params) > 0):
        disp_params_set("Required parameters", required_params) 
    if (optional_params is not None) and (len(optional_params) > 0):
        disp_params_set("Optional parameters", optional_params) 

def disp_params_set(params_name,
                    params,
                    indent="",
                    ):
    if params is None:
        return
    if len(indent) == 0:
        disp_header(params_name)
    else:
        disp_header(params_name, new_line=False)
    if type(params) is list:
        for entry_idx in range(len(params)):
            entry = params[entry_idx]
            if ((type(entry) is list) or
                (type(entry) is dict) or
                (type(entry) is OrderedDict)
                ):
                disp_params_set(indent+"  #"+str(entry_idx+1),
                                entry,
                                indent=indent+"  ")
            else:
                disp_param(indent+"#"+str(entry_idx+1), entry)
    else:
        # assuming params is dict
        for key in params:
            val = params[key]
            if ((type(val) is list) or
                (type(val) is dict) or
                (type(val) is OrderedDict) or
                (type(val) is defaultdict)
                ):
                disp_params_set("  "+indent+key, val, indent=indent+"  ")
            else:
                disp_param(indent+key, params[key])

