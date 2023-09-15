#! /usr/bin/env python
import versioneer

from setuptools import setup, find_packages


setup(
    name = "q2-long-gsea",
    version = versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    packages = find_packages(),
    package_date = {},
    author = "Asa Henry",
    author_email = "ajh728@nau.edu",
    description = "",
    license = "",
    url = "https://github.com/LadnerLab/q2-long-gsea",
    entry_points = {
        "qiime2.plugins": ["q2-long-gsea=q2_long_gsea.plugin_setup:plugin"]
    },
    zip_safe = False
)