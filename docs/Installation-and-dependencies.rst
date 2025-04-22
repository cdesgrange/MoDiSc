.. _installation:

Installation and dependencies
-----------------------------

First create a new environment dedicated to run MoDiSc.

.. code-block:: bash

  $ conda create -n modisc_env python=python3.9

Activate this environment

.. code-block:: bash
  
  $ conda activate modisc_env

Install MoDiSc via pip:

.. code-block:: bash

  $ pip install MoDiSc

The packages required by MoDiSc are indicated in MoDiSc/packages.py.
Most of the packages can be installed via 

.. code-block:: bash

  $ pip install packagename

**Exceptions:**

- opencv must be installed via 

.. code-block:: bash 

  pip install opencv-python

- emcee: the appropriate version of the emcee package to be used by MoDiSc is located in MoDiSc/emcee/


