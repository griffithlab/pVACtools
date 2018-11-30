from setuptools import setup
import os

import sys
if sys.version_info < (3,5):
    print("This python version is not supported:")
    print(sys.version)
    print("pVACtools requires python 3.5 or greater")
    sys.exit(1)

pvacseq_data_files = []
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
pvacfuse_data_files = []
for dirpath, dirnames, filenames in os.walk("tools/pvacfuse/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacfuse_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
pvacvector_data_files = []
for dirpath, dirnames, filenames in os.walk("tools/pvacvector/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacvector_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/VEP_plugins"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("tools/pvacseq/iedb_alleles"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
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
client_data = []
for dirpath, dirnames, filenames in os.walk("utils/pvacviz/client"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            client_data.append(os.path.join(
                os.path.relpath(
                    dirpath,
                    'utils/pvacapi/client'
                ),
                filename
            ))

setup(
    name="pvactools",
    version="1.1.4",
    packages=[
        "tools",
        "tools.pvacfuse",
        "tools.pvacvector",
        "tools.pvacseq",
        "lib",
        "utils.pvacapi",
        "utils.pvacapi.controllers",
        "utils.pvacviz"
    ],
    entry_points={
        "console_scripts":[
            "pvactools = tools.main:main",
            "pvacseq = tools.pvacseq.main:main",
            "pvacfuse = tools.pvacfuse.main:main",
            "pvacvector = tools.pvacvector.main:main",
            "pvacapi = utils.pvacapi.app:main [API]",
            "pvacviz = utils.pvacviz.app:main [API]"
        ]
    },
    install_requires=[
        'PyVCF',
        'requests',
        'PyYAML',
        'connexion==1.4.2',
        'biopython',
        'networkx',
        'simanneal',
        'pandas',
        'wget',
        'mhcflurry',
        'mhcnuggets',
        'pysam',
        'tensorflow==1.8.0'
    ],
    package_data={
        'tools.pvacseq': pvacseq_data_files,
        'tools.pvacfuse': pvacfuse_data_files,
        'tools.pvacvector': pvacvector_data_files,
        'utils.pvacapi': server_data,
        'utils.pvacviz': client_data,
    },
    extras_require={
        'API':[
            'connexion==1.4.2',
            'py-postgresql',
            'watchdog',
            'flask-cors',
            'bokeh==0.13.0',
            'tornado==5.0.2',
            'swagger-spec-validator==2.1.0',
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5"
    ],

    author = "Jasreet Hundal, Susanna Kiwala, Joshua McMichael, Yang-Yang Feng, Christopher A. Miller, Aaron Graubert, Amber Wollam, Connor Liu, Jonas Neichin, Megan Neveau, Jason Walker, Elaine R. Mardis, Obi L. Griffith, Malachi Griffith",
    author_email = "help@pvactools.org",
    description = "A cancer immunotherapy tools suite",
    license = "NPOSL-3.0",
    keywords = "antigens neoantigens cancer sequencing variant variants fusion fusions",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/pVACtools",
)
