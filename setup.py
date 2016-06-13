from setuptools import setup

import sys
if sys.version_info < (3,5):
    print("This python version is not supported:")
    print(sys.version)
    print("pVAC-Seq requires python 3.5 or greater")
    sys.exit(1)

setup(
    name="pvacseq",
    version="1.0.2",
    packages=["pvacseq", "pvacseq.lib"],
    entry_points={
        "console_scripts":[
            "pvacseq = pvacseq.pvacseq:main"
        ]
    },
    install_requires=[
        'PyVCF',
    ],
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
    url = "https://github.com/griffithlab/pVAC-Seq",   #
)
