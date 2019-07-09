import unittest
import gc
import inspect
import datetime
import os
import shutil
from os.path import join as join_path
from os.path import dirname
from cancersig.utils import logger

class pyCancerSigBase(object):
    """ pyCancerSig base class """

    def __init__(self, debug_mode=False, verbose=True, *args, **kwargs):
        self.__time_stamp = datetime.datetime.now()
        self.pkg_root_dir = dirname(dirname(__file__))
        self.__debug_mode = debug_mode
        self.__verbose = verbose

    @property
    def current_func_name(self):
        return inspect.stack()[1][3]

    def remove_dir(self, dir_name):
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name, ignore_errors=True)

    def create_dir(self, dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    def dbg(self, dbg_msg=""):
        if self.__debug_mode:
            frm = inspect.stack()[1]
            mod = inspect.getmodule(frm[0])
            logger.getLogger(mod.__name__)
            logger.dbg(dbg_msg)

    def info(self, info_msg=""):
        if self.__verbose:
            frm = inspect.stack()[1]
            mod = inspect.getmodule(frm[0])
            logger.getLogger(mod.__name__)
            logger.info(info_msg)

    def warning(self, warning_msg=""):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        logger.getLogger(mod.__name__)
        logger.warning(warning_msg)

    def throw(self, err_msg=""):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        logger.getLogger(mod.__name__)
        logger.throw(err_msg)


class Tester(unittest.TestCase, pyCancerSigBase):
    """ general pyCancerSig template for testing """

    individual_debug = False

    def __init__(self, test_module_name, *args, **kwargs):
        self.test_module_name = test_module_name
        super(Tester, self).__init__(*args, **kwargs)
        pyCancerSigBase.__init__(self)

    def init_test(self,
                  test_function,
                  ):
        self.test_function = test_function
        self.__set_dir()
        self.empty_working_dir()

    def __set_dir(self):
        working_subdir = "/".join(self.test_module_name.split('.')[:-2])
        working_subdir = join_path(working_subdir,
                                   self.test_module_name.split('.')[-1][5:])
        self.working_dir = os.getenv("PYTHON_TEST_DIR",
                                     join_path(self.pkg_root_dir,
                                               "tmp"))
        self.working_dir = join_path(self.working_dir,
                                     working_subdir)
        self.working_dir = join_path(join_path(self.working_dir,
                                               self.__class__.__name__[4:]),
                                     self.test_function)
        data_subdir = "/".join(self.test_module_name.split('.')[:-1])
        data_subdir = join_path(data_subdir,
                                "data")
        data_subdir = join_path(data_subdir,
                                self.test_module_name.split('.')[-1][5:])
        self.data_dir = join_path(self.pkg_root_dir,
                                  data_subdir)
        self.data_dir = join_path(join_path(self.data_dir,
                                            self.__class__.__name__[4:]),
                                  self.test_function)

    def empty_working_dir(self):
        if not self.individual_debug:
            self.remove_dir(self.working_dir)
        self.create_dir(self.working_dir)

