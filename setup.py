from setuptools import setup

setup(
    name="pVAC-Seq",
    version="1.5",
    packages=["src", "src.pvac_seq"],
    entry_points={
        "console_scripts":[
            "pVAC-Seq = src.pVAC_Seq:main"
        ]
    }
)
