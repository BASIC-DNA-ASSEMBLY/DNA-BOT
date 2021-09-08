from setuptools import setup, find_packages

setup(
    name='dnabot',
    version='v1.0.0',
    description='DNA assembly using BASIC on OpenTrons (DNA-BOT)',
    license='MIT',
    author='Matthew C Haines',
    author_email='matthew.haines10@imperial.ac.uk',
    packages=find_packages(),
    include_package_data=True,
    keywords=['dnabot'],
    url='https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT',
    classifiers=[
        'Topic :: Scientific/Engineering',
    ]
)