from setuptools import setup, find_packages


setup(
    name        = 'MultilayerOptics',
    version     = '0.1.0',
    author      = 'Kanglin Xiong',
    packages    = find_packages(),

    package_data = {"MultilayerOptics" : ["*/*.txt"]},

    url         = 'http://pypi.python.org/pypi/MultilayerOptics/',
    license     = 'LICENSE.txt',
    description = 'Calculate optical properties of layered structures.',
    long_description = open('README.txt').read()
)
