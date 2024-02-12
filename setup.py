from setuptools import setup
import os

import sys
if (sys.version_info.major, sys.version_info.minor) < (3,7):
    print("This python version is not supported:")
    print(sys.version)
    print("pVACtools requires python 3.7 or greater")
    sys.exit(1)

try:
    import pypandoc
    pypandoc.download_pandoc()
    long_description = pypandoc.convert_file('README.md', 'rst')
except:
    long_description = ""

pvacseq_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacseq/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
pvacfuse_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacfuse/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacfuse_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
pvacvector_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacvector/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacvector_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
pvacbind_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacbind/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacbind_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
pvacview_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacview"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacview_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacseq/VEP_plugins"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacseq/iedb_alleles"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))

setup(
    name="pvactools",
    version="4.0.7",
    packages=[
        "pvactools.tools",
        "pvactools.tools.pvacbind",
        "pvactools.tools.pvacfuse",
        "pvactools.tools.pvacvector",
        "pvactools.tools.pvacseq",
        "pvactools.tools.pvacview",
        "pvactools.lib",
    ],
    entry_points={
        "console_scripts":[
            "pvactools = pvactools.tools.main:main",
            "pvacseq = pvactools.tools.pvacseq.main:main",
            "pvacbind = pvactools.tools.pvacbind.main:main",
            "pvacfuse = pvactools.tools.pvacfuse.main:main",
            "pvacvector = pvactools.tools.pvacvector.main:main",
            "pvacview = pvactools.tools.pvacview.main:main",
        ]
    },
    install_requires=[
        'vcfpy',
        'requests',
        'PyYAML>=5.1',
        'biopython==1.77',
        'networkx',
        'simanneal',
        'pandas',
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
        'polars==0.16.18',
        'bigmhc @ git+https://github.com/griffithlab/bigmhc.git#egg=bigmhc',
        'deepimmuno @ git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno',
    ],
    package_data={
        'pvactools.tools.pvacseq': pvacseq_data_files,
        'pvactools.tools.pvacfuse': pvacfuse_data_files,
        'pvactools.tools.pvacvector': pvacvector_data_files,
        'pvactools.tools.pvacbind': pvacbind_data_files,
        'pvactools.tools.pvacview': pvacview_data_files,
    },
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],

    author = "Jasreet Hundal, Susanna Kiwala, Joshua McMichael, Christopher A Miller, Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng, Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker, William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith",
    author_email = "help@pvactools.org",
    description = "A cancer immunotherapy tools suite",
    long_description = long_description,
    license = "BSD-3-Clause-Clear",
    keywords = "antigens neoantigens cancer sequencing variant variants fusion fusions",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/pVACtools",
)
