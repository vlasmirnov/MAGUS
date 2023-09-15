import setuptools
from setuptools import find_packages
import os

# load version number, follows https://stackoverflow.com/a/16084844
exec(open('./version.py').read())

# https://stackoverflow.com/a/36693250/13241395
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="magus-msa",
    version=__version__,
    author="vlasmirnov",
    description="Multiple Sequence Alignment using Graph Clustering",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vlasmirnov/MAGUS",
    project_urls={
        "Bug Tracker": "https://github.com/vlasmirnov/MAGUS/issues",
    },
    entry_points={
        'console_scripts': [
            'magus=magus:main',
        ],
    },
    py_modules=['magus', 'magus_configuration', 'version'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    package_dir={"": "."},
    package_data={'': package_files('magus_tools')},
    packages=find_packages(),
    install_requires=["dendropy>=4.5.2"],
    python_requires=">=3.6",
)
