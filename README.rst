========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |github-actions|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/python-maskpy/badge/?style=flat
    :target: https://python-maskpy.readthedocs.io/
    :alt: Documentation Status

.. |github-actions| image:: https://github.com/Alan-Manning/python-maskpy/actions/workflows/github-actions.yml/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/Alan-Manning/python-maskpy/actions

.. |codecov| image:: https://codecov.io/gh/Alan-Manning/python-maskpy/branch/main/graphs/badge.svg?branch=main
    :alt: Coverage Status
    :target: https://app.codecov.io/github/Alan-Manning/python-maskpy

.. |version| image:: https://img.shields.io/pypi/v/maskpy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/maskpy

.. |wheel| image:: https://img.shields.io/pypi/wheel/maskpy.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/maskpy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/maskpy.svg
    :alt: Supported versions
    :target: https://pypi.org/project/maskpy

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/maskpy.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/maskpy

.. |commits-since| image:: https://img.shields.io/github/commits-since/Alan-Manning/python-maskpy/v1.0.1.svg
    :alt: Commits since latest release
    :target: https://github.com/Alan-Manning/python-maskpy/compare/v1.0.1...main



.. end-badges

package for creating gds maks in python

* Free software: GNU Lesser General Public License v3 or later (LGPLv3+)

Installation
============

::

    pip install maskpy

You can also install the in-development version with::

    pip install https://github.com/Alan-Manning/python-maskpy/archive/main.zip


Documentation
=============


https://python-maskpy.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
