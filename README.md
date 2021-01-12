# SimExtract


Simulation of extraction processes
------------------

Model for the extraction of volatile and semi-volatile compounds based on numerical simulations of diffusion processes. It shows some of the equilibria established along time. This model can be on-line simulated with the input of variables like diffusion coefficients, Henry constant, or the thickness of each phase, among others.

A robust and general numerical model to simulate closed system extraction processes for analytical purposes based on one dimensional diffusion differential equations.

Main developers:
- Pedro S. Peixoto - University of São Paulo (ppeixoto@usp.br)
- Luís Moreira Gonçalves - University of São Paulo (lmgoncalves@iq.usp.br)

The software is freely provided without restriction for general purpose usage. Please cite the software as:
- (TODO)

Suggestions for improvements can be mailed to developers and are more than welcome! 

--------------
 Installation
--------------

**Use the code without installing it, via binder:** 

- A generic example: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pedrospeixoto/SimExtract/HEAD?filepath=simex_example.ipynb)

- Especific examples in the examples folder from here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pedrospeixoto/SimExtract/HEAD)

**Run locally on your computer:**

1) Install python-3.8.x or higher

2) Install required python packages (pip install pkgname):
- numpy
- scipy
- matplotlib
- pandas

3) git clone the repository (or download the repository files).

--------------
Execution
--------------

**With Jupyter-notebook:**
1) Copy one of the notebook examples (e.g. simex_example.ipynb -> mysimex_example.ipynb)
2) Open the notebook with jupyter-notebook 
3) Edit the parameters in the notebook and run the cells.
4) Outputs will be placed in an output/ folder

**With a parameter file:**
1) Edit parameters file: simex_params.py
2) Run execution file: python simex_run.py
3) Outputs will be placed in an output/ folder



