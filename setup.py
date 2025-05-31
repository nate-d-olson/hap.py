# Modern setup for Cython extension modules
from setuptools import find_packages, setup

setup(
    name="happy",
    version="0.4.0",
    description="Haplotype Comparison Tools",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
    install_requires=[
        "numpy>=1.15.0",
        "pysam>=0.15.0",
        "scipy>=1.0.0",
        "pandas>=0.23.0",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "hap.py = happy.hap:main",
            "qfy = happy.qfy:main",
            "pre = happy.pre:main",
        ],
    },
)
