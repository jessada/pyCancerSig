import sys
from collections import OrderedDict
from cancersig.utils import logger
from cancersig.utils import disp
from cancersig.profile.msi import MSIProfiler

APP_PROFILE_MSI_DESCRIPTION = "To extract MSI features from msisensor output"

def app_profile_snv(*args, **kwargs):
    logger.getLogger(__name__)
    func_name = sys._getframe().f_code.co_name[4:]

    disp.new_section_txt("S T A R T <" + func_name + ">")
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H <" + func_name + ">")

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
#    required_params['project output directory (-O)'] = kwargs['project_out_dir']
#    input_binary = kwargs['input_binary']
#    input_file_prefix = kwargs['input_file_prefix']
#    if kwargs['input_binary']:
#        required_params['input file prefix (--bfile)'] = input_file_prefix
#    else:
#        required_params['input file prefix (--file)'] = input_file_prefix
    optional_params = OrderedDict()
#    if kwargs['input_dna_regions'] is not None:
#        optional_params['input dna regions (-R)'] = kwargs['input_dna_regions'].split(",")
#    if kwargs['filter_criteria'] is not None:
#        optional_params['filter criteria (--filter_criteria)'] = kwargs['filter_criteria'].split(",")
#    required_params['phenotype file (--pheno)'] = kwargs['phenotype_file']
#    if kwargs['project_code'] is not None:
#        optional_params['project code (-p)'] = kwargs['project_code']
#    if kwargs['sample_info'] is not None:
#        optional_params['sample information (-p)'] = kwargs['sample_info'].split(",")
#    optional_params['job allocation time (--job_alloc_time)'] = kwargs['job_alloc_time']
#    optional_params['cut-off p-value (--cutoff_pvalue)'] = kwargs['cutoff_pvalue']
#    optional_params['cut-off odds ratio (--cutoff_ors)'] = kwargs['cutoff_ors']
#    optional_params['families haplotyep file prefix (--fam_hap)'] = kwargs['fam_hap_prefix']
#    optional_params['haplotype windows (--hap_window)'] = kwargs['hap_window_sizes'].split(",")
#    optional_params['output jobs setup file (-o)'] = kwargs['out_jobs_setup_file']
    disp.show_config(app_description=APP_PROFILE_MSI_DESCRIPTION,
                     required_params=required_params,
                     optional_params=optional_params,
                     )
    disp.new_section_txt("E X E C U T I N G < " + func_name + " >")
#    create_jobs_setup_file(project_name=kwargs['project_name'],
#                           project_out_dir=kwargs['project_out_dir'],
#                           input_file_prefix=kwargs['input_file_prefix'],
#                           phenotype_file=kwargs['phenotype_file'],
#                           input_binary=kwargs['input_binary'],
#                           input_dna_regions=kwargs['input_dna_regions'],
#                           cutoff_pvalue=kwargs['cutoff_pvalue'],
#                           cutoff_ors=kwargs['cutoff_ors'],
#                           hap_window_sizes=kwargs['hap_window_sizes'],
#                           project_code=kwargs['project_code'],
#                           sample_info=kwargs['sample_info'],
#                           job_alloc_time=kwargs['job_alloc_time'],
#                           filter_criteria=kwargs['filter_criteria'],
#                           fam_hap_prefix=kwargs['fam_hap_prefix'],
#                           out_jobs_setup_file=kwargs['out_jobs_setup_file'],
#                           )
    msi_profiler = MSIProfiler()
    msi_profiler.extract_features(raw_msisensor_out,
                                  raw_msisensor_out_somatic,
                                  sample_id,
                                  output_file,
                                  )
    logger.getLogger(__name__)
    disp.new_section_txt("F I N I S H < " + func_name + " >")
