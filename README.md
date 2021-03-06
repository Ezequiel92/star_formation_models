<div align="center">
    <h1>✨ Star formation models</h1>
</div>

<p align="center">
    <a href="https://github1s.com/ezequiel92/star_formation_models"><img src="http://forthebadge.com/images/badges/built-with-science.svg"></a>
    <a href="https://julialang.org"><img src="https://forthebadge.com/images/badges/made-with-julia.svg"></a>
</p>

<p align="center">
    <a href="https://github.com/ezequiel92/star_formation_models/blob/main/LICENSE"><img src="https://img.shields.io/github/license/ezequiel92/star_formation_models?style=flat&logo=GNU&labelColor=2B2D2F"></a>
</p>

Explicit implementation of several star formation models.

The goal of the models is to improve the realism of the star formation rate (SFR) in simulations of galaxy formation and evolution.

- Models in `archived/models` are old models which served as experiments.

  - Each model has its folder `archived/models/model_XXX`, where the file `models/model_XXX/notebook_XXX.jl` is a [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook showcasing the model.
  - The files inside the folders `archived/models/model_XXX/alternative_implementations` are scripts implementing the same model in different programming languages, namely [Julia](https://julialang.org), [Python](https://www.python.org/), [C](http://www.open-std.org/jtc1/sc22/wg14/www/standards), [C++](https://isocpp.org/), and [Wolfram Language](https://www.wolfram.com/language/).
  - Each one of those scripts outputs a `.dat` file with the evolution of the variables after integrating the system for 1 Gyr. These files are saved in `archived/models/model_XXX/alternative_implementations/plots`, where the Mathematica notebook `plots.nb` will plot all the data together to compare the implementations. The goal is to test that all the implementations give the same results.
  - The notebooks and scripts in `archived/analysis` were use to study and understand the old models.

- The main model is described in `model/model_019/model_019.jl`, which can be view online with [![Binder](https://mybinder.org/badge_logo.svg)](https://binder.plutojl.org/v0.19.9/open?url=https%253A%252F%252Fgithub.com%252Fezequiel92%252Fstar_formation_models%252Fblob%252Fmain%252Fmodel_019%252Fmodel_019.jl%253Fraw%253Dtrue)
- The dependencies are given by the `Project.toml`, and `requirements.txt` files.

## 🧮 Solvers

The ODE solvers used in each implementation are:

- **Python**: [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
- **C**: [GSL](https://www.gnu.org/software/gsl/)
- **C++**: [Boost odeint](https://headmyshoulder.github.io/odeint-v2/index.html)
- **Mathematica**: [NDSolve](https://reference.wolfram.com/language/ref/NDSolve.html)
- **Julia**: [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/)

## 📘 Documentation

Although every implementation fully describes each system of ODEs, the main documentation for each model is the corresponding Pluto notebook. They show the system of equations, the values of the constants and plot the evolution of the main variables for a mock set of the initial conditions.

The Wolfram Language notebooks inside each `archived/models/model_XXX/alternative_implementations` plots the evolution of the variables too. They also generate the Jacobian functions needed for the C and C++ implementation, which are saved as `.../alternative_implementations/jacobian_gsl.txt` and `.../alternative_implementations/jacobian_boost.txt` respectively.

## 👨‍💻 Run Locally

- **Julia and Python**: The Julia and Python scripts in `alternative_implementations` can be run directly with the corresponding interpreters. To install the dependencies:

  - Clone the project

    ```bash
    git clone https://github.com/ezequiel92/star_formation_models.git
    ```

  - Go to the project directory

    ```bash
    cd path/to/star_formation_models
    ```

  - Install the dependencies specified by the `Project.toml` file (for Julia within the Julia REPL)

    ```julia
    (@v1.7) pkg> activate .

    (star_formation_models) pkg> instantiate
    ```

  - Or the ones specified by the `requirements.txt` file (for Python)

    ```python
    pip install -r requirements.txt
    ```
  
  - Run (for Julia within the Julia REPL)

    ```julia
    julia> include("path/to/script.jl")
    ```

  - Or run (for Python)

    ```python
    python path/to/script.py
    ```

- **C**: To run the C scripts you have to link at compile time the [GSL](https://www.gnu.org/software/gsl/) library.
- **C++**: To run the C++ scripts you have to include at compile time the [Boost](https://www.boost.org/) library.

## ⚠️ Warning

These scripts are written for my personal use and may break at any moment. So, use them at your own risk.
