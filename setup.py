from setuptools import setup

setup(
    name="pVAC-Seq",
    version="1.5",
    packages=["pvac_seq", "pvac_seq.lib"],
    entry_points={
        "console_scripts":[
            "pVAC-Seq = pvac_seq.pVAC_Seq:main"
        ]
    },

    author = "Griffith Lab",
    author_email = "abc@email.com",
    description = "A cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVAC-Seq)",
    license = "NPOSL-3.0",
    keywords = "antigens neoantigens cancer sequencing variant variants",
    url = "https://github.com/griffithlab/pVAC-Seq",   #
)
