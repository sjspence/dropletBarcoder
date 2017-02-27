from distutils.core import setup

setup(
    name='epicBarcoder',
    version='0.1dev',
    packages=['epicBarcoder', ],
    license='MIT',
    install_requires=['biopython',
                      'matplotlib',
                      'pandas'],
    long_description=open('README.md').read(),
)
