# Binding Energy and Saturation Density in the Walecka Model

This project contains a numerical solution to the problem of computing the
binding energy per nucleon in the Walecka model at zero temperature for
different values of the baryon density. The notebook reproduces results
from specific references and validatesthe use of given coupling constants.

## Problem Description

The task involves solving Eq. (3.23) from a specific textbook at zero temperature
numerically forvarious values of baryon density. The solution is then used to
compute the binding energy per nucleon, checking against provided values for
coupling constants. Additionally, the notebook aims to reproduce results from
Fig. 3.2 of the reference material.

## Repository Contents
- `Walecka_model.ipynb`: The main Jupyter notebook containing the code for
  solving the equations, plotting results, and performing validations.
- `functions/`: A directory containing Python scripts with helper functions used in the notebook.
  
    - `constants.py`: Defines constants used in the calculations.
    - `arrays.py`: Contains functions for creating and managing arrays of data.
    - `constant_couplings.py`: Determines the model coupling constants.
    - `solved_equations.py`: Includes functions for solving the equations relevant to the model.

- `plots/`: Directory where generated plots are saved.

## Reference

SCHMITT, Andreas. **Dense matter in compact stars: A pedagogical introduction.** Springer, 2010.

