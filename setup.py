from setuptools import setup

setup(
    name="pvacseq",
    version="1.0.1",
    packages=["pvacseq", "pvacseq.lib"],
    entry_points={
        "console_scripts":[
            "pvacseq = pvacseq.pvacseq:main"
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5"
    ],

    author = "Griffith Lab",
    author_email = "pvacseq-support@genome.wustl.edu",
    description = "Personalized Variant Antigens by Cancer Sequencing (pVAC-Seq)",
    long_description = "A cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVAC-Seq)",
    license = "NPOSL-3.0",
    keywords = "antigens neoantigens cancer sequencing variant variants",
    url = "https://github.com/griffithlab/pVAC-Seq",   #
)
