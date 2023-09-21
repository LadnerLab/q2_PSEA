#! /usr/bin/env python
import versioneer

from setuptools import setup, find_packages


setup(
    name = "q2-PSEA",
    version = versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    packages = find_packages(),
    package_date = {},
    author = "Asa Henry",
    author_email = "ajh728@nau.edu",
    description = "",
    license = "",
    url = "https://github.com/LadnerLab/q2-PSEA",
    entry_points = {
        "qiime2.plugins": ["q2-PSEA=q2_PSEA.plugin_setup:plugin"]
    },
    zip_safe = False
)