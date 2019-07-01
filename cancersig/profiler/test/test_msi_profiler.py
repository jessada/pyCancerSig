import unittest
from os.path import join as join_path
from cancersig.template import Tester
from cancersig.utils.logger import info
from cancersig.profiler.msi import MSI_FEATURES_TEMPLATE

class TestMSIProfiler(Tester):

    def __init__(self, methodName):
        super(TestMSIProfiler, self).__init__(methodName=methodName,
                                         test_module_name=__name__,
                                         )

    def setUp(self):
        pass


    def test_extract_features_1(self):
        """ test extracting msi features in basic scenario """

        self.init_test(self.current_func_name)
        sample_info = join_path(self.data_dir,
                                "sample.info")
        info()
        info(MSI_FEATURES_TEMPLATE)
