import setuptools
import re

from pathlib import Path


with open("README.md", "r") as fp:
    long_description = fp.read()


def get_version():
    """Get version number from __init__.py"""
    version_file = Path(__file__).resolve().parent / "g2j" / "__init__.py"
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read_text(), re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


setuptools.setup(
    name="genome2json",
    version="0.0.1",
    author="Cameron Gilchrist",
    description="Parse genomes in GenBank/GFF3 formats to JSON",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.github.com/gamcil/genome2json",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">3.6",
    entry_points={"console_scripts": ["g2j = g2j.main:run"]},
)
