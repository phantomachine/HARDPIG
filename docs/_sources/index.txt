.. ted documentation master file, created by
   sphinx-quickstart on Tue Sep 24 13:10:27 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
===============================

| **H** eterogeneous 
|   **A** gent 
| **R** ecursive
|   **D** ynamic 
| **P** ublic
|   **I** nsurance
| **G** ames.
|
|

This web documentation details the methods proposed in the paper by Kam and
Stauber, "Computing Dynamic Public Insurance Games with Endogenous Agent
Distributions".

:download:`This <../build/latex/doc.pdf>` is a downloadable PDF version of the documentation. 



.. toctree::
    :maxdepth: 5
    :numbered:

    statespace
    actions_transitions
    intersections
    sampling
    compute_statespace
    payoff_approximation
    payoff_concept
    payoff_bilinear
    payoff_bilinear_punish
    payoff_compute


Software
==============================

HARDPIG relies on the following software:
    
* `MATLAB <http://www.mathworks.com>`_ platform

* `YALMIP <http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.BMIBNB>`_ BMIBNB global optimization solver
    
    * Bilinear Matrix Inequality Branch-and-Bound solver

* `GNU GLPK <http://www.gnu.org/software/glpk/>`_ linear programming
  solver (in ANSI C)

* `SNOPT <http://www.sbsi-sol-optimize.com/asp/sol_products_snopt_desc.htm>`_ general-purpose local optimization solver (Fortran)

    * `MATLAB Executable binaries files <http://www.scicomp.ucsd.edu/~peg/Software.html>`_ for up to 300 variables and 300 constraints available freely from Phillip E. Gill.

Source codes (in MATLAB) and executables (in C/Fortran) are available from the authors via email.


Indices and tables
==============================

* :ref:`genindex`
* :ref:`search`

