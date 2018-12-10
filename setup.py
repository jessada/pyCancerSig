from setuptools import setup

setup(
    name='pyCancerSig',
    version='0.0.1',
    author='Jessada Thutkawkorapin',
    author_email='jessada.thutkawkorapin@gmail.com',
    packages=['cancersig',
              'cancersig.features',
              ],
    package=['pyCancerSig'],
    scripts=['bin/extract_msi_features',
             ],
    url='http://pypi.python.org/pypi/pyCancerSig/',
    license='LICENSE.txt',
    description='Python packages for deciphering cancer signature processes',
    long_description=open('README.md').read(),
)
