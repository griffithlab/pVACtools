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
pvacbind_data_files = []
for dirpath, dirnames, filenames in os.walk("tools/pvacbind/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacbind_data_files.append(os.path.join(os.path.relpath(dirpath, 'tools/pvacseq'), filename))
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
                    'utils/pvacapi'
                ),
                filename
            ))

setup(
    name="pvactools",
    version="2.0.1",
    packages=[
        "tools",
        "tools.pvacbind",
        "tools.pvacfuse",
        "tools.pvacvector",
        "tools.pvacseq",
        "tools.pvacview",
        "lib",
        "utils.pvacapi",
        "utils.pvacapi.controllers",
        "utils.pvacviz"
    ],
    entry_points={
        "console_scripts":[
            "pvactools = tools.main:main",
            "pvacseq = tools.pvacseq.main:main",
            "pvacbind = tools.pvacbind.main:main",
            "pvacfuse = tools.pvacfuse.main:main",
            "pvacvector = tools.pvacvector.main:main",
            "pvacview = tools.pvacview.main:main",
            "pvacapi = utils.pvacapi.main:main",
            "pvacviz = utils.pvacviz.app:main"
        ]
    },
    install_requires=[
        'PyVCF',
        'requests',
        'PyYAML>=5.1',
        'connexion==1.4.2',
        'biopython==1.76',
        'networkx',
        'simanneal',
        'pandas',
        'wget',
        'mhcflurry',
        'mhcnuggets',
        'pysam',
        'Pillow',
        'pymp-pypi',
        'connexion==1.4.2',
        'py-postgresql',
        'watchdog',
        'flask-cors',
        'bokeh==0.13.0',
        'tornado==5.0.2',
        'swagger-spec-validator==2.1.0',
        'jsonschema==2.6.0',
        'mock',
        'vaxrank>=1.1.0',
    ],
    package_data={
        'tools.pvacseq': pvacseq_data_files,
        'tools.pvacfuse': pvacfuse_data_files,
        'tools.pvacvector': pvacvector_data_files,
        'tools.pvacbind': pvacbind_data_files,
        'utils.pvacapi': server_data,
        'utils.pvacviz': client_data,
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
    license = "BSD-3-Clause-Clear",
    keywords = "antigens neoantigens cancer sequencing variant variants fusion fusions",
    #This needs to be the url where the code is being hosted
    url = "https://github.com/griffithlab/pVACtools",
)
