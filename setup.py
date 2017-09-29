from setuptools import setup
import os

import sys
if sys.version_info < (3,5):
    print("This python version is not supported:")
    print(sys.version)
    print("pVAC-Seq requires python 3.5 or greater")
    sys.exit(1)

data_files = []
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            data_files.append(os.path.join('..', dirpath, filename))
for dirpath, dirnames, filenames in os.walk("tools/pvacfuse/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            data_files.append(os.path.join('..', dirpath, filename))
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/VEP_plugins"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            data_files.append(os.path.join('..', dirpath, filename))
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/iedb_alleles"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            data_files.append(os.path.join('..', dirpath, filename))
server_data = []
for dirpath, dirnames, filenames in os.walk("utils/pvacapi"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            server_data.append(os.path.join(
                os.path.relpath(
                    dirpath,
                    'utils/pvacapi'
                ),
                filename
            ))

setup(
    name="pvacseq",
    version="4.1.0b4",
    packages=["tools.pvacseq", "lib", "utils.pvacapi", "utils.pvacapi.controllers"],
    entry_points={
        "console_scripts":[
            "pvactools = tools.main:main",
            "pvacseq = tools.pvacseq.main:main",
            "pvacfuse = tools.pvacfuse.main:main",
            "pvacseq-api = utils.pvacapi.app:main [API]"
        ]
    },
    install_requires=[
        'PyVCF',
        'requests',
        'PyYAML',
        'connexion',
        'biopython',
        'networkx',
        'simanneal',
        'pandas'
    ],
    package_data={
        'tools.pvacseq' : data_files,
        'utils.pvacapi' : server_data,
    },
    extras_require={
        'API':[
            'connexion',
            'py-postgresql',
            'watchdog',
            'flask-cors',
            'bokeh',
            'pvacseq-client'
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5"
    ],

    author = "Jasreet Hundal, Susanna Kiwala, Aaron Graubert, Jason Walker, Chris Miller, Malachi Griffith and Elaine Mardis",
    author_email = "pvacseq-support@genome.wustl.edu",
    description = "Personalized Variant Antigens by Cancer Sequencing (pVAC-Seq)",
    long_description = "A cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVAC-Seq)",
    license = "NPOSL-3.0",
    keywords = "antigens neoantigens cancer sequencing variant variants",
    url = "https://github.com/griffithlab/pVAC-Seq",
)
