from distutils.core import setup

setup(
    name='dropletBarcoder',
    version='0.1dev',
    packages=['dropletBarcoder', ],
    license='MIT',
    install_requires=['biopython',
                      'matplotlib',
                      'pandas'],
    long_description=open('README.md').read(),
)
