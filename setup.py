from setuptools import setup, find_namespace_packages
import sys

if (sys.version_info.major, sys.version_info.minor) < (3, 9) or (sys.version_info.major, sys.version_info.minor) > (3, 11):
    print("This Python version is not supported:")
    print(sys.version)
    print("pVACtools supports Python versions 3.9, 3.10, and 3.11")
    sys.exit(1)

def readme():
    with open("README.md") as f:
        return f.read()

setup(
    name="pvactools",
    version="7.0.1",
    entry_points={
        "console_scripts":[
            "pvactools = pvactools.tools.main:main",
            "pvacseq = pvactools.tools.pvacseq.main:main",
            "pvacbind = pvactools.tools.pvacbind.main:main",
            "pvacfuse = pvactools.tools.pvacfuse.main:main",
            "pvacvector = pvactools.tools.pvacvector.main:main",
            "pvacview = pvactools.tools.pvacview.main:main",
            "pvacsplice = pvactools.tools.pvacsplice.main:main",
        ]
    },
    install_requires=[
        'vcfpy==0.13.8',
        'requests',
        'PyYAML>=5.1',
        'biopython==1.77',
        'networkx',
        'simanneal',
        'numpy==1.26.4',
        'pandas<2.1.0',
        'wget',
        'pysam',
        'Pillow',
        'pymp-pypi',
        'mock',
        'vaxrank>=1.1.0',
        'varcode>=1.1.0',
        'mhcnuggets==2.4.1',
        'mhcflurry==2.0.6',
        'testfixtures',
        'gtfparse==2.0.1',
        'pyfaidx>=0.7.1',
        'fsspec<=2025.3.0',
        'packaging',
        'pyarrow',
        'polars==0.16.18',
        'XlsxWriter',
        'openpyxl',
        'deepdiff',
        'scikit-learn>=1.6.0,<1.8.0',
        'imblearn',
    ],
    packages=find_namespace_packages(),
    include_package_data=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    author = "Susanna Kiwala, Jasreet Hundal, Joshua McMichael, Christopher A Miller, Luke Hendrickson, Alexander T Wollam, Huiming Xia, Miller Richters, Connor J Liu, Evelyn Schmidt, Jennie Yao, My Hoang, Sidi Zhao, Yang-Yang Feng, Aaron P Graubert, Jayden Zebrowski, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith",
    author_email = "help@pvactools.org",
    description = "A cancer immunotherapy tools suite",
    long_description = readme(),
    long_description_content_type="text/markdown",
    license = "BSD-3-Clause-Clear",
    keywords = "immunogenomics, immunotherapy, antigens, neoantigens, cancer, sequencing, variant, variants, fusion, fusions",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/pVACtools",
)
