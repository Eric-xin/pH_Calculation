# 🧪 pH Calculation 🤖
### Comprehensive Analysis of Acid-Base Equilibria Using MBE, CBE, and PBE

[![Python](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![Jupyter](https://img.shields.io/badge/jupyter-1.0.0-blue.svg)](https://jupyter.org/)
[![Scipy](https://img.shields.io/badge/scipy-1.13.1-blue.svg)](https://www.scipy.org/)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

This repository contains the code and data used in the paper "Comprehensive Analysis of Acid-Base Equilibria Using MBE, CBE, and PBE". The code is written in Python and can be run in any Python environment with the required libraries installed. Jupyter notebooks are also provided for interactive exploration of the calculations and results.

The paper is available in the `documents` folder, and the compiled pKa data used in the calculations is available in the `pKa_data` folder.

Want to read the paper? [Click here](documents/doc.pdf).

Also there is an online version of the paper available on [Click here](https://ph_calc.ericxin.eu).

## Introduction

This project explores the calculation of acid-base equilibrium calculations, focusing on the Material Balance Equation (MBE), Charge Balance Equation (CBE), and Proton Balance Equation (PBE). The code provided in this repository allows for the calculation of pH values for various acid-base systems using the CBE and PBE methods. 

## Theory

### Acid-Base Reactions and Equilibria in Solution

Acid-base reactions in solution are governed by the activity of ions, which refers to their effective concentration in chemical reactions. Due to electrostatic interactions between ions, their effective concentration (activity) can differ from their actual concentration. The relationship between the activity $a_i$ and the concentration $c_i$ of an ion $i$ is given by:

$$
a_i = \gamma_i \, c_i 
$$

where $\gamma_i$ is the activity coefficient. For dilute solutions, the Debye-Hückel equation can be used to approximate the activity coefficient:

$$
- \log \gamma_i = \frac{0.51 z_i^2 \sqrt{I}}{1 + B a_i \sqrt{I}}
$$

where $z_i$ is the charge of ion $i$, $B$ is a constant, $a_i$ is the ion volume parameter, and $I$ is the ionic strength of the solution.

### Equilibrium Constants

The equilibrium constant $K$ for a reaction can be expressed in terms of activities (thermodynamic constant $K^\circ$) or concentrations (concentration constant $K$):

$$
K^\circ = \frac{a_{\text{A}^-} a_{\text{HB}^+}}{a_{\text{B}} a_{\text{HA}}}
$$
$$
K = \frac{[\text{A}^-][\text{HB}^+]}{[\text{B}][\text{HA}]}
$$

The relationship between $K^\circ$ and $K$ involves the activity coefficients:

$$
K^\circ = K \cdot \frac{\gamma_{\text{A}^-} \gamma_{\text{HB}^+}}{\gamma_{\text{B}} \gamma_{\text{HA}}}
$$

### Material Balance Equation (MBE)

The Material Balance Equation (MBE) states that the total concentration of a substance in a chemical equilibrium system is equal to the sum of the equilibrium concentrations of all relevant species. For example, for a solution with concentration $c$ of $\text{H}_3\text{PO}_4$:

$$
[\text{H}_3\text{PO}_4] + [\text{H}_2\text{PO}_4^-] + [\text{HPO}_4^{2-}] + [\text{PO}_4^{3-}] = c
$$

### Charge Balance Equation (CBE)

The Charge Balance Equation (CBE) ensures that the total positive charge in a solution equals the total negative charge. (As the solution is always neutral in charge) For example, in a solution with concentration $c$ of $\text{NaCN}$:

$$
[\text{H}^+] + [\text{Na}^+] = [\text{OH}^-] + [\text{CN}^-]
$$

Below is a diagram of the CBE curve for a solution of $\text{(NH4)}_3\text{(PO4)}$:

![CBE Curve](assets/img/CBE.png)

### Proton Balance Equation (PBE)

The Proton Balance Equation (PBE) states that the amount of protons gained by substances (bases) equals the amount of protons lost by substances (acids). For example, in an aqueous solution of $(\text{NH}_4)\text{HPO}_4$:

$$
[\text{H}^+] + [\text{H}_2\text{PO}_4^-] + 2[\text{H}_3\text{PO}_4] = [\text{OH}^-] + [\text{PO}_4^{3-}] + [\text{NH}_3]
$$

These equations are essential tools for understanding and predicting the behavior of acid-base systems in solution.

Here is an example of the PBE curve for a solution of $\text{(NH4)}_3\text{(PO4)}$:

![PBE Curve](assets/img/PBE.png)

> For more details, please refer to the paper at [documents/doc.pdf](documents/doc.pdf).

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/Eric-xin/pH_Calculation.git
    cd pH_Calculation
    ```

2. Create a virtual environment and activate it:
    ```sh
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. Install the required libraries:
    ```sh
    pip install -r requirements.txt
    ```

## Running the Code

To run the code, you can use the provided Jupyter notebooks or run the Python scripts directly. For example, to run the `CBE.ipynb` notebook:

1. Start Jupyter Notebook:
    ```sh
    jupyter notebook
    ```

2. Open `CBE.ipynb` in the Jupyter interface and run the cells.

## Code Structure

The code is organized as follows:

```
.
├── CBE.ipynb # Jupyter notebook for Charge Balance Equation (CBE) based calculations
├── LICENSE # MIT License
├── PBE_investigation.ipynb # Jupyter notebook for Proton Balance Equation (PBE) based calculations (investigation)
├── assets
│   └── img # Image files used in the readme
├── calculation # Code for MBE, CBE, and PBE calculations
│   ├── Basics.py # Basic classes and functions for acid-base calculations, based on MBE
│   ├── CBE_calc.py
│   ├── PBE_calc.py
├── documents
│   ├── doc.pdf # Paper "Comprehensive Analysis of Acid-Base Equilibria Using MBE, CBE, and PBE"
├── pKa_data # Compiled pKa data for various acids, used in the calculations
│   ├── pka-compilation-williams.pdf
│   └── readme.md
└── readme.md # This file
```

## Examples

**Note: currently only CBE package is available for general use. PBE package is still under development. Code in PBE_investigation.ipynb is for demonstration on Ammonium Phosphate system only.**

### Example 1: CBE calculation of 0.01 M HCl

```python
from calculation import *

HCl = Inert(charge=-1, conc=0.01)
pH = CBE_calc(HCl)
pH.pH_calc()

print(pH.pH) # the result should be 1.9999977111816385
```

### Example 2: CBE calculation of 0.01 M (NH4)2(HPO4)

```python
NH4 = Acid(charge=1, conc=0.01*3, pKa=9.25)
pKa = [1.97, 6.82, 12.5]
P = Acid(charge=0, conc=0.01, pKa=pKa)

pH = CBE_calc(NH4, P)
pH.pH_calc()

print(pH.pH) # the result should be 8.952952575683597
```

For more examples, see the Jupyter notebooks.

## Contact

If you have any questions or suggestions, please feel free to contact me at [me@ericxin.eu](mailto:me@ericxin.eu).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.