import sys
from collections import OrderedDict
from cancersig.utils import logger
from cancersig.utils import disp
from cancersig.profile.msi import MSIProfiler
from cancersig.profile.sv import SVProfiler

APP_PROFILE_MSI_DESCRIPTION = "To extract MSI features from msisensor output"
APP_PROFILE_SV_DESCRIPTION = "To extract structural variantion features from FindSV output"

def app_profile_snv(*args, **kwargs):
    logger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name[4:]

    disp.new_section_txt("S T A R T <" + func_name + ">")
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

def app_profile_sv(*args, **kwargs):
    logger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name[4:]

    input_vcf_file = kwargs['input_vcf_file']
    output_file = kwargs['output_file']
    sample_id = kwargs['sample_id']

    disp.new_section_txt("S T A R T < " + func_name + " >")
    required_params = OrderedDict()
    required_params['input vcf file name (-i/--input_vcf_file)'] = input_vcf_file
    required_params['output file name (-o/--output_file)'] = output_file
    optional_params = OrderedDict()
    optional_params['sample id (--sample_id)'] = sample_id
    disp.show_config(app_description=APP_PROFILE_SV_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + func_name + " >")
    sv_profiler = SVProfiler()
    sv_profiler.extract_features(input_vcf_file,
                                 output_file,
                                 sample_id,
                                 )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + func_name + " >")

def app_profile_msi(*args, **kwargs):
    logger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name[4:]

    raw_msisensor_out = kwargs['raw_msisensor_out']
    raw_msisensor_out_somatic = kwargs['raw_msisensor_out_somatic']
    sample_id = kwargs['sample_id']
    output_file = kwargs['output_file']

    disp.new_section_txt("S T A R T < " + func_name + " >")
    required_params = OrderedDict()
    required_params['msisensor output (--raw_msisensor_out)'] = raw_msisensor_out
    required_params['msisensor somatic (--raw_msisensor_out_somatic)'] = raw_msisensor_out_somatic
    required_params['sample id (--sample_id)'] = sample_id
    required_params['output file name (-o/--output_file)'] = output_file
    optional_params = OrderedDict()
    disp.show_config(app_description=APP_PROFILE_MSI_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + func_name + " >")
    msi_profiler = MSIProfiler()
    msi_profiler.extract_features(raw_msisensor_out,
                                  raw_msisensor_out_somatic,
                                  sample_id,
                                  output_file,
                                  )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + func_name + " >")
