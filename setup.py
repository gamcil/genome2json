import setuptools


with open("README.md", "r") as fp:
    long_description = fp.read()


setuptools.setup(
    name="genome2json",
    version="0.0.1",
    author="Cameron Gilchrist",
    description="Parse genomes in GenBank/GFF3 formats to JSON",
    long_description=long_description,
    url="https://www.github.com/gamcil/genome2json",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">3.6",
)
