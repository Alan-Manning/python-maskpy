#!/usr/bin/env python
import re
from pathlib import Path

from setuptools import find_packages, setup


def read(*names, **kwargs):
    with Path(__file__).parent.joinpath(*names).open(encoding=kwargs.get("encoding", "utf8")) as fh:
        return fh.read()


setup(
    name="maskpy",
    version="1.0.0",
    license="LGPL-3.0-or-later",
    description="package for creating gds maks in python",
    long_description="{}\n{}".format(
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub("", read("README.rst")),
        re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst")),
    ),
    author="Alan Manning",
    author_email="Alan_Manning@Live.co.uk",
    url="https://github.com/Alan-Manning/python-maskpy",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[path.stem for path in Path("src").glob("*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)" "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Utilities",
    ],
    project_urls={
        "Documentation": "https://python-maskpy.readthedocs.io/",
        "Changelog": "https://python-maskpy.readthedocs.io/en/latest/changelog.html",
        "Issue Tracker": "https://github.com/Alan-Manning/python-maskpy/issues",
    },
    keywords=[
        "maskpy",
        "gds",
        "gdsii",
        ".gds",
    ],
    python_requires=">=3.10",
    install_requires=[
        "click",
        "gdspy",
        "openpyxl",
        "pandas",
        "phidl",
        "pyyaml",
        "shapely",
        "tqdm",
    ],
    extras_require={},
    entry_points={
        "console_scripts": [
            "maskpy = maskpy.cli:main",
        ]
    },
)
