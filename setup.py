from setuptools import setup
import os

import sys
if (sys.version_info.major, sys.version_info.minor) < (3, 9) or (sys.version_info.major, sys.version_info.minor) > (3, 11):
    print("This Python version is not supported:")
    print(sys.version)
    print("pVACtools supports Python versions 3.9, 3.10, and 3.11")
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
            pvacfuse_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacfuse'), filename))
pvacvector_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacvector/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacvector_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacvector'), filename))
pvacbind_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacbind/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacbind_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacbind'), filename))
pvacsplice_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacsplice/example_data"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacbind_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacsplice'), filename))
pvacview_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacview"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacview_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacview'), filename))
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacseq/VEP_plugins"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvacseq/iedb_alleles"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
for dirpath, dirnames, filenames in os.walk("pvactools/supporting_files"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvacseq_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvacseq'), filename))
pvaccompare_data_files = []
for dirpath, dirnames, filenames in os.walk("pvactools/tools/pvaccompare/html_report"):
    for filename in filenames:
        if not (filename.endswith(".py") or filename.endswith(".pyc")):
            pvaccompare_data_files.append(os.path.join(os.path.relpath(dirpath, 'pvactools/tools/pvaccompare'), filename))

setup(
    name="pvactools",
    version="6.0.2",
    packages=[
        "pvactools.tools",
        "pvactools.tools.pvacbind",
        "pvactools.tools.pvacfuse",
        "pvactools.tools.pvacvector",
        "pvactools.tools.pvacseq",
        "pvactools.tools.pvacview",
        "pvactools.tools.pvacsplice",
        "pvactools.tools.pvaccompare",
        "pvactools.tools.pvaccompare.compare_tools",
        "pvactools.tools.pvaccompare.comparisons",
        "pvactools.tools.pvaccompare.runners",
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
    ],
    package_data={
        'pvactools.tools.pvacseq': pvacseq_data_files,
        'pvactools.tools.pvacfuse': pvacfuse_data_files,
        'pvactools.tools.pvacvector': pvacvector_data_files,
        'pvactools.tools.pvacbind': pvacbind_data_files,
        'pvactools.tools.pvacview': pvacview_data_files,
        'pvactools.tools.pvacsplice': pvacsplice_data_files,
        'pvactools.tools.pvaccompare': pvaccompare_data_files,
    },
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

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
