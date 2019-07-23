from setuptools import setup

setup(
    name='pyCancerSig',
    version='0.0.1',
    author='Jessada Thutkawkorapin',
    author_email='jessada.thutkawkorapin@gmail.com',
    packages=['cancersig',
              'cancersig.profile',
              'cancersig.app',
              'cancersig.utils',
              'cancersig.signature',
              ],
    scripts=['bin/cancersig',
             ],
    package=['pyCancerSig'],
    url='http://pypi.python.org/pypi/pyCancerSig/',
    license='LICENSE.txt',
    description='Python packages for deciphering cancer signature processes',
    long_description=open('README.md').read(),
    install_requires=[
        "Biopython >= 1.72",
        "pandas >= 0.23.4",
        "matplotlib >= 3.0.2",
        "scikit-learn >= 0.20.1",
        "seaborn >= 0.9.0",
        ],
)
