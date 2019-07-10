import unittest
import filecmp
from os.path import join as join_path
from cancersig.template import Tester
from cancersig.profile.sv import SVProfiler

class TestSVProfiler(Tester):

    def __init__(self, methodName):
        super(TestSVProfiler, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        self.__sv_profiler = SVProfiler()
        self.__sv_profiler.debug_mode = True

    def test_profile_1(self):
        """ test extracting sv features on standard argument (sample id from vcf header) """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf")
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__sv_profiler.profile(input_vcf_file,
                                   output_file,
                                   )
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in SVProfiler.profile()')

    def test_profile_2(self):
        """ test extracting sv features with sample id as an argument """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf")
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__sv_profiler.profile(input_vcf_file,
                                   output_file,
                                   sample_id="example_sample_from_argument",
                                   )
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in SVProfiler.profile()')
