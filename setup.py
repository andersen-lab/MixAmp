# ----------------------------------------------------------------------------
# Copyright (c) 2021-, Andersen Lab development team.
#
# Distributed under the terms of the BSD 2-Clause License
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()


description = ("Amplicon read simualtor")


setup(
    name="amp-seq-sim",
    version="2024.08",
    packages=find_packages(include=['pywgsim*']),
    author="Maryam Ahmadi Jeshvaghane",
    license='BSD 2-Clause',
    author_email="mahmadi@scripps.edu",
    url="https://github.com/mariaelf97/Amp-seq-sim",
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'amp-seq-sim=amp_seq_sim._cli:cli',  # Use underscores here
        ],
    },
    package_data={
        'amp-seq-sim': ['data/*', ],
    },
    install_requires=[
        "click", "pandas", "biopython",
        "regex", "numpy"
        ]
)