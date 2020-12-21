# mextractmodel

--------------
Micro Extraction Diffusion Model
--------------

A robust and general numerical model for gas-diffusion micro extraction mechamisms based on one dimensional diffusion differential equations.

Main developer:
Pedro S. Peixoto - University of São Paulo (ppeixoto@usp.br)

Development is mainly based on the following article:
Zhang, Z. and Pawliszyn, J., 1993. Headspace solid-phase microextraction. Analytical chemistry, 65(14), pp.1843-1852.

--------------
 Installation
--------------
0) You can use the code without installing it via binder:

A generic example:
https://mybinder.org/v2/gh/pedrospeixoto/mextractmodel/master?filepath=model.ipynb

Especific examples in the examples folder from here:
https://mybinder.org/v2/gh/pedrospeixoto/mextractmodel/master

or you can run locally on your computer,

1) Install python-3.8.x or higher

2) Install required python packages (pip install pkgname):
- numpy
- scipy
- matplotlib
- pandas

--------------
Execution
--------------

- With Jupyter-notebook:
1) Copy one of the notebook examples (e.g. model.ipynb -> mymodel.ipynb)
2) Open the notebook with jupyter-notebook 
3) Edit the parameters in the notebook and run the cells.
4) Outputs will be placed in an output/ folder

- With a parameter file:
1) Edit parameters file: mex_params.py
2) Run execution file: python model.py
3) Outputs will be placed in an output/ folder



