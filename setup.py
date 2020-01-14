import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dnabot",
    version="1.0.1",
    author="Matthew Haines",
    author_email="hainesm6@gmail.com",
    description="Code for DNA assembly using BASIC on OpenTrons (DNA-BOT)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BASIC-DNA-ASSEMBLY/DNA-BOT",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
