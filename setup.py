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
    name="mixamp",
    version="2024.10",
    packages=find_packages(include=['mixamp']),
    author="Maryam Ahmadi Jeshvaghane",
    license='BSD 2-Clause',
    author_email="mahmadi@scripps.edu",
    url="https://github.com/andersen-lab/MixAmp",
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points='''
        [console_scripts]
        mixamp=mixamp._cli:cli
        ''',
    package_data={
        'mixamp': ['data/*', ],
    },
    install_requires=[
        "click", "pandas", "biopython",
        "regex", "numpy"
        ]
)