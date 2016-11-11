from distutils.core import setup

setup(
    name='spenceOTU',
    version='0.1dev',
    packages=['spenceOTU', ],
    license='MIT',
    install_requires=['biopython',
                      'matplotlib',
                      'pandas'],
    long_description=open('README.md').read(),
)
