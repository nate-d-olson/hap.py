# Modern setup for Cython extension modules
import os

from setuptools import setup

ext_modules = []

setup(
    name="happy",
    version="0.4.0",
    description="Haplotype Comparison Tools",
    packages=["happy", "happy.Haplo", "happy.Tools"],
    package_dir={"happy": "src/python"},
    ext_modules=ext_modules,
    install_requires=[
        "numpy>=1.15.0",
        "pysam>=0.15.0",
        "scipy>=1.0.0",
        "pandas>=0.23.0",
        "cython>=0.29.0",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "hap.py = happy.hap:main",
            "qfy = happy.qfy:main",
            "pre = happy.pre:main",
            "ftx = happy.ftx:main",
            "cnx = happy.cnx:main",
            "ovc = happy.ovc:main",
        ],
    },
)
