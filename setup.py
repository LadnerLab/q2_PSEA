#! /usr/bin/env python
import versioneer

from setuptools import setup, find_packages


setup(
    name="q2-PSEA",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_date={},
    author="Asa Henry",
    author_email="ajh728@nau.edu",
    url="https://github.com/LadnerLab/q2-PSEA",
    python_requires=">=3.8",
    install_requirements=[
        "numpy>=1.24.4",
        "pandas>=1.5.3",
        "rpy2=3.5.15",
        "r-base>=4.2.3",
        "r-essentials",
        "r-ggraph",
        "r-ggforce",
        "r-scatterpie",
        "r-igraph",
        "r-biocmanager",
        "bioconductor-clusterprofiler",
        "bioconductor-enrichplot",
        "bioconductor-ggtree",
        "pillow=9.4",
        "r-ggplot2=3.3"
    ],
    entry_points={
        "qiime2.plugins": ["q2-PSEA=q2_PSEA.plugin_setup:plugin"]
    },
    zip_safe=False,
    download_url="https://github.com/LadnerLab/q2-PSEA.git"
)
