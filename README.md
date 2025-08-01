# Non-Linear Finite Element Analysis of Euler-Bernoulli Beams

This is a MATLAB full implementation of an 1D analysis of Beams using the Euler-Bernoulli model and with the Finite Element Method solver for ODEs.

The purpose of this code is merely academic and its use isn't recommended or allowed with commerce purposes.

## How to use

The code was implemented with the aim that a carefull reader shoulb be able to understand and modified for your personal goals. Here, I'll provide some key points for how the analysed structure and algorithim's info should be provided in the code.

The main way to enter all the needed info is by an Excel or Google Sheet file with the required extension `.xlsx`. The name of the file should be placed at the variable `nome` (line 17). An example file is provided in this repo under the name **_Vigas.xlsx_**. This spreadsheet has 6 differents sheets that should provide the following info:

1. `Nós`: The beam's geometric info, as nodes, coordinates, bound constraints and prescribed displacements;

2. `Elementos (Vigas)`: The beam's material and section description (without values);

3. `Forças`: The external load configuration. It's only available the analysis with point load;

4. `Materiais`: The material's properties listed in the second spreadsheet should be provided here: Young Modulus (E), Poisson Coefficient $\nu$ and Mass Density;

5. `Seções`: The beam's section listed in the second spreadsheet should be informed here: base $b$ and height $h$, area of the section, moment of inertia and the neutral axis $y$ coordinate;

6. `Parâmetros Newton`: Adjustment of the non-linear system algorithm. Here we use the Newthon method for solving the Non-Linear system of equations.
