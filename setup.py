import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="magus-msa",
    version="0.1.0a0",
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    package_dir={"": "."},
    packages=[],
    install_requires=["dendropy>=4.5.2"],
    python_requires=">=3.6",
)