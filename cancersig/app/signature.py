import sys
from collections import OrderedDict
from cancersig.utils import logger
from cancersig.utils import disp
from cancersig.signature.decipher import CancerSigNMFController

APP_SIGNATURE_DESCRIPTION = "To perform unsupervised learning and decipher underlying cancer signature processes"

def app_signature_decipher(*args, **kwargs):
    logger.getLogger(__name__)
    app_name = "cancersig signature decipher"

    mutation_profiles = kwargs['mutation_profiles']
    output_prefix = kwargs['output_prefix']
    min_signatures = kwargs['min_signatures']
    max_signatures = kwargs['max_signatures']

    disp.new_section_txt("S T A R T < " + app_name + " >")
    required_params = OrderedDict()
    required_params['input mutation profiles (-i/--mutation_profiles)'] = mutation_profiles
    required_params['output file prefix (-o/--output_prefix)'] = output_prefix
    optional_params = OrderedDict()
    optional_params['minimum number of signatures (--min_signatures)'] = min_signatures
    optional_params['maximum number of signatures (--max_signatures)'] = max_signatures
    disp.show_config(app_description=APP_SIGNATURE_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + app_name + " >")
    cancersig_nmf_controller = CancerSigNMFController()
    cancersig_nmf_controller.decipher(mutation_profiles,
                                      output_prefix,
                                      min_signatures,
                                      max_signatures,
                                      )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + app_name + " >")
