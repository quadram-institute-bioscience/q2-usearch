# ----------------------------------------------------------------------------
# Copyright (c) 2024, Thanh Le Viet.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = ("A template QIIME 2 plugin.")

setup(
    name="q2-usearch",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Thanh Le Viet",
    author_email="thanh@cloudbinfies.com",
    description=description,
    url="https://github.com/quadram-institute-bioscience/q2-usearch",
    entry_points={
        "qiime2.plugins": [
            "q2_usearch="
            "q2_usearch"
            ".plugin_setup:plugin"]
    },
    package_data={
        "q2_usearch": ["citations.bib"],
        "q2_usearch.tests": ["data/*"],
    },
    zip_safe=False,
)
