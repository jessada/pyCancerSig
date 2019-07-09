import sys
from collections import OrderedDict
from cancersig.utils import logger
from cancersig.utils import disp
from cancersig.profile.snv import SNVProfiler
from cancersig.profile.sv import SVProfiler
from cancersig.profile.msi import MSIProfiler
from cancersig.profile.merge import ProfileMerger

APP_PROFILE_SNV_DESCRIPTION = "To extract SNV features from somatic SNV call in vcf format"
APP_PROFILE_SV_DESCRIPTION = "To extract structural variantion features from FindSV output"
APP_PROFILE_MSI_DESCRIPTION = "To extract MSI features from msisensor output"
APP_PROFILE_MERGE_DESCRIPTION = "To merge all profiles from all samples into one single profile"

def app_profile_snv(*args, **kwargs):
    logger.getLogger(__name__)
    app_name = "cancersig profile snv"

    input_vcf_file = kwargs['input_vcf_file']
    ref_genome_file = kwargs['ref_genome_file']
    output_file = kwargs['output_file']
    gt_format = kwargs['gt_format']
    sample_id = kwargs['sample_id']

    disp.new_section_txt("S T A R T < " + app_name + " >")
    required_params = OrderedDict()
    required_params['input vcf file name (-i/--input_vcf_file)'] = input_vcf_file
    required_params['reference genome (-r/--reference)'] = ref_genome_file
    required_params['output file name (-o/--output_file)'] = output_file
    optional_params = OrderedDict()
    optional_params['sample id (--sample_id)'] = sample_id
    optional_params['genotype format (-g/--gt_format)'] = gt_format
    disp.show_config(app_description=APP_PROFILE_SNV_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + app_name + " >")
    snv_profiler = SNVProfiler()
    snv_profiler.profile(input_vcf_file,
                         ref_genome_file,
                         output_file,
                         gt_format,
                         sample_id,
                         )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + app_name + " >")

def app_profile_sv(*args, **kwargs):
    logger.getLogger(__name__)
    app_name = "cancersig profile sv"

    input_vcf_file = kwargs['input_vcf_file']
    output_file = kwargs['output_file']
    sample_id = kwargs['sample_id']

    disp.new_section_txt("S T A R T < " + app_name + " >")
    required_params = OrderedDict()
    required_params['input vcf file name (-i/--input_vcf_file)'] = input_vcf_file
    required_params['output file name (-o/--output_file)'] = output_file
    optional_params = OrderedDict()
    optional_params['sample id (--sample_id)'] = sample_id
    disp.show_config(app_description=APP_PROFILE_SV_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + app_name + " >")
    sv_profiler = SVProfiler()
    sv_profiler.profile(input_vcf_file,
                        output_file,
                        sample_id,
                        )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + app_name + " >")

def app_profile_msi(*args, **kwargs):
    logger.getLogger(__name__)
    app_name = "cancersig profile msi"

    raw_msisensor_out = kwargs['raw_msisensor_out']
    raw_msisensor_out_somatic = kwargs['raw_msisensor_out_somatic']
    sample_id = kwargs['sample_id']
    output_file = kwargs['output_file']

    disp.new_section_txt("S T A R T < " + app_name + " >")
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
    disp.new_section_txt("E X E C U T E < " + app_name + " >")
    msi_profiler = MSIProfiler()
    msi_profiler.profile(raw_msisensor_out,
                         raw_msisensor_out_somatic,
                         sample_id,
                         output_file,
                         )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + app_name + " >")

def app_profile_merge(*args, **kwargs):
    logger.getLogger(__name__)
    app_name = "cancersig profile merge"

    input_dirs = kwargs['input_dirs'].split(",")
    output_file = kwargs['output_file']
    profile_types = kwargs['profile_types']

    disp.new_section_txt("S T A R T < " + app_name + " >")
    required_params = OrderedDict()
    required_params['input directories (-i/--input_dirs)'] = input_dirs
    required_params['output file name (-o/--output_file)'] = output_file
    optional_params = OrderedDict()
    optional_params['profile types to be merged (--profile_types)'] = profile_types
    disp.show_config(app_description=APP_PROFILE_MERGE_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T E < " + app_name + " >")
    profile_merger = ProfileMerger()
    profile_merger.merge(input_dirs,
                         output_file,
                         profile_types,
                         )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + app_name + " >")
