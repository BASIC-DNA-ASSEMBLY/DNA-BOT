from setuptools import setup, find_packages

setup(
    name='dnabot',
    version='0.0.0',
    description='DNA assembly using BASIC on OpenTrons (DNA-BOT)',
    license='MIT',
    author='Matthew C Haines',
    author_email='thomas.duigou@inrae.fr',
    packages=find_packages(),
    include_package_data=True,
    keywords=['dnabot'],
    url='https://github.com/brsynth/DNA-BOT',
    classifiers=[
        'Topic :: Scientific/Engineering',
    ]
)