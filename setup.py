from setuptools import setup
import os

import sys
if sys.version_info < (3,6):
    print("This python version is not supported:")
    print(sys.version)
    print("pVACtools requires python 3.6 or greater")
    sys.exit(1)

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
    version="3.0.0rc",
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
        'PyVCF',
        'requests',
        'PyYAML>=5.1',
        'biopython==1.76',
        'networkx',
        'simanneal',
        'pandas',
        'wget',
        'pysam',
        'Pillow',
        'pymp-pypi',
        'mock',
        'vaxrank>=1.1.0',
        'keras==2.4.3',
        'tensorflow==2.2.2',
        'mhcnuggets==2.3.3',
        'mhcflurry==2.0.1',
        'testfixtures'
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

        "Programming Language :: Python :: 3.6"
    ],

    author = "Jasreet Hundal, Susanna Kiwala, Joshua McMichael, Yang-Yang Feng, Christopher A. Miller, Aaron Graubert, Amber Wollam, Connor Liu, Jonas Neichin, Megan Neveau, Jason Walker, Elaine R. Mardis, Obi L. Griffith, Malachi Griffith",
    author_email = "help@pvactools.org",
    description = "A cancer immunotherapy tools suite",
    license = "BSD-3-Clause-Clear",
    keywords = "antigens neoantigens cancer sequencing variant variants fusion fusions",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/pVACtools",
)
