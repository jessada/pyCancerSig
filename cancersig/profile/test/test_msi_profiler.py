import unittest
import filecmp
from os.path import join as join_path
from cancersig.template import Tester
from cancersig.utils.logger import info
from cancersig.profile.msi import MSI_FEATURES_TEMPLATE
from cancersig.profile.msi import MSIProfiler

class TestMSIProfiler(Tester):

    def __init__(self, methodName):
        super(TestMSIProfiler, self).__init__(methodName=methodName,
                                              test_module_name=__name__,
                                              )

    def setUp(self):
        self.__msi_profiler = MSIProfiler()
        self.__msi_profiler.debug_mode = True

    def test_profile_1(self):
        """ test extracting msi features for an MSI positive sample """

        self.init_test(self.current_func_name)
        raw_msisensor_out = join_path(self.data_dir,
                                      "raw_msisensor_out")
        raw_msisensor_out_somatic = join_path(self.data_dir,
                                              "raw_msisensor_out_somatic")
        sample_id = "example_sample"
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__msi_profiler.profile(raw_msisensor_out,
                                    raw_msisensor_out_somatic,
                                    sample_id,
                                    output_file,
                                    )
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in MSIProfiler.profile()')

    def test_profile_2(self):
        """ test extracting msi features for an MSI negative sample """

        self.init_test(self.current_func_name)
        raw_msisensor_out = join_path(self.data_dir,
                                      "raw_msisensor_out")
        raw_msisensor_out_somatic = join_path(self.data_dir,
                                              "raw_msisensor_out_somatic")
        sample_id = "example_sample"
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__msi_profiler.profile(raw_msisensor_out,
                                    raw_msisensor_out_somatic,
                                    sample_id,
                                    output_file,
                                    )
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in MSIProfiler.profile()')

