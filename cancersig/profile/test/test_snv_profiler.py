import unittest
import filecmp
from os.path import join as join_path
from cancersig.config import ENABLE_CPU_INTENSIVE_UNITTEST
from cancersig.template import Tester
from cancersig.profile.snv import SNVProfiler

class TestSNVProfiler(Tester):

    def __init__(self, methodName):
        super(TestSNVProfiler, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        self.__snv_profiler = SNVProfiler()
        self.__snv_profiler.debug_mode = True

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_profile_1(self):
        """ test extracting snv features on standard argument (sample id from vcf header) """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_output_feature_file")
        exp_output_event_file = join_path(self.data_dir,
                                          "exp_output_event_file")
        self.assertTrue(filecmp.cmp(output_event_file, exp_output_event_file), 'Malfunction in SNVProfiler.profile()')
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_profile_2(self):
        """ test extracting snv features with sample id as an argument """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    sample_id="example_sample_id_from_argument",
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                             "exp_output_feature_file")
        exp_output_event_file = join_path(self.data_dir,
                                          "exp_output_event_file")
        self.assertTrue(filecmp.cmp(output_event_file, exp_output_event_file), 'Malfunction in SNVProfiler.profile()')
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_profile_3(self):
        """ test simple vcf """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                             "exp_output_feature_file")
        exp_output_event_file = join_path(self.data_dir,
                                          "exp_output_event_file")
        self.assertTrue(filecmp.cmp(output_event_file, exp_output_event_file), 'Malfunction in SNVProfiler.profile()')
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    @unittest.skip("No real life data to be tested for SNV called by MuTect2")
    def test_profile_4(self):
        """ test extracting snv features with custom GT field """

        self.init_test(self.current_func_name)
#        input_vcf_file = join_path(self.data_dir,
#                                   "input.vcf")
#        ref_genome_file = join_path(self.data_dir,
#                                    "ref_genome.fa")
#        output_file = join_path(self.working_dir,
#                                self.current_func_name+".txt")
#        self.__snv_profiler.profile(input_vcf_file,
#                                             ref_genome_file,
#                                             output_file,
#                                             )
#        exp_output_file = join_path(self.data_dir,
#                                    "exp_output_file")
#        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in SNVProfiler.profile()')
