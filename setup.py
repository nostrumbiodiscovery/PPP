import PPP
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from distutils.extension import Extension

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
     long_description = f.read()

setup(
    name="PPPele",
    version=PPP.__version__,
    description='Protein PELE preparation',
    long_description=long_description,
    url="https://github.com/NostrumBioDiscovery/PPP",
    author='Nostrum Biodiscovery',
    author_email='it@nostrumbiodiscovery.com',
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={"PPP/Data/Templates/OPLS2005/HeteroAtoms/": ['*'],
                  "PPP/Data/Templates/OPLS2005/Protein": ['*'],
                  "PPP/Data/RotamerLibs/": ['*']},
    include_package_data=True,
    install_requires=["biopython", "prody", "pytest", "scipy"],
    cmdclass=cmdclass,
    ext_modules=ext_modules  # accepts a glob pattern
)
