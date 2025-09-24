"""A setuptools based setup module."""

# Always prefer setuptools over distutils
# To use a consistent encoding
from codecs import open
from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

with open("requirements.txt", "r") as fh:
    requirements = fh.read().splitlines()

with open("requirements_test.txt", "r") as fh:
    test_requirements = fh.read().splitlines()

setup(
    name="pygdsm",
    version="1.6.4",
    description="Python Global Sky Model of diffuse Galactic radio emission",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # The project's main homepage.
    url="https://github.com/telegraphic/pygdsm",
    # Author details
    author="Danny C. Price",
    author_email="daniel.price@skao.int",
    # Choose your license
    license="MIT",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    # What does your project relate to?
    keywords="radio astronomy sky model galactic diffuse emission",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=["contrib", "docs", "tests*"]),
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=requirements,
    tests_require=test_requirements,
)
