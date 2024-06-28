# CollocInfer Method Implementation for Lotka-Volterra System

## Overview

This repository contains an implementation of the **CollocInfer** method for parameter estimation in ordinary differential equations (ODEs), specifically applied to the **Lotka-Volterra predator-prey model**. The CollocInfer method combines data fitting with differential equation constraints using basis function approximations, facilitating the estimation of both system parameters and initial conditions from noisy observational data.

## Mathematical Formulation

### Lotka-Volterra Model

The Lotka-Volterra equations describe the dynamics of a predator-prey system:

\[
\begin{aligned}
\frac{dC}{dt} &= \alpha C - \beta_C C B, \\
\frac{dB}{dt} &= \beta_B C B - \delta B.
\end{aligned}
\]

Where:
- \( C(t) \): Prey population (e.g., Chorella algae)
- \( B(t) \): Predator population (e.g., Brachionis rotifers)
- \( \alpha \): Birth rate of prey
- \( \beta_C \): Predation rate coefficient affecting prey mortality
- \( \beta_B \): Growth rate coefficient for predators due to predation
- \( \delta \): Death rate of predators

### CollocInfer Method

The CollocInfer method aims to estimate the parameters \( \theta = [\alpha, \beta_C, \beta_B, \delta] \) and initial conditions \( C(0) \) and \( B(0) \) from noisy observations of \( C(t) \) and \( B(t) \).

#### Basis Function Approximation

Each state variable is approximated using a basis function expansion:

\[
x(t, \mathbf{c}(\theta, \lambda)) = \sum_{k} c_k(\theta, \lambda) \phi_k(t)
\]

Where:
- \( \phi_k(t) \): B-spline basis functions
- \( c_k(\theta, \lambda) \): Coefficients dependent on parameters \( \theta \) and smoothing parameter \( \lambda \)

#### Differential Equations

The Lotka-Volterra system is expressed in terms of the logarithms of the populations:

\[
\begin{aligned}
\frac{d \ln C}{dt} &= \alpha - \beta_C B, \\
\frac{d \ln B}{dt} &= \beta_B C - \delta.
\end{aligned}
\]

#### Inner Criterion

The inner criterion \( J(\mathbf{c}, \theta, \lambda) \) combines data fitting and equation fitting:

\[
\begin{aligned}
J(\mathbf{c}, \theta, \lambda) &= \sum_{j} \left[ (\ln C_j - \ln x_C(t_j, \mathbf{c}_C))^2 + (\ln B_j - \ln x_B(t_j, \mathbf{c}_B))^2 \right] \\
&\quad + \lambda \int \left[ \left( \frac{d}{dt} \ln x_C(t, \mathbf{c}_C) - (\alpha - \beta_C x_B(t, \mathbf{c}_B)) \right)^2 \right. \\
&\quad \left. + \left( \frac{d}{dt} \ln x_B(t, \mathbf{c}_B) - (\beta_B x_C(t, \mathbf{c}_C) - \delta) \right)^2 \right] dt
\end{aligned}
\]

*Explanation:* The first term measures the discrepancy between the observed data and the model predictions. The second term enforces the differential equation constraints, weighted by the smoothing parameter \( \lambda \).

#### Outer Criterion

The outer criterion \( H(\theta, \lambda) \) is defined as:

\[
H(\theta, \lambda) = J(\mathbf{c}^*(\theta, \lambda), \theta, \lambda)
\]

Where \( \mathbf{c}^*(\theta, \lambda) = \arg \min_{\mathbf{c}} J(\mathbf{c}, \theta, \lambda) \) represents the optimal coefficients for given \( \theta \) and \( \lambda \).

#### Optimization Process

1. **Inner Optimization:** For fixed \( \theta \) and \( \lambda \), find \( \mathbf{c}^*(\theta, \lambda) \) by minimizing \( J(\mathbf{c}, \theta, \lambda) \) with respect to \( \mathbf{c} \).

2. **Outer Optimization:** Update \( \theta \) by minimizing \( H(\theta, \lambda) \).

3. **Lambda Progression:**
   - Start with a small \( \lambda_0 \).
   - Progressively increase \( \lambda \) using a logarithmic sequence.
   - For each \( \lambda_i \):
     - Use \( \theta^* \) from \( \lambda_{i-1} \) as the initial guess.
     - Perform the optimization steps until convergence.

*Explanation:* This approach balances the trade-off between fitting the data and satisfying the differential equation constraints.

### Lotka-Volterra Specifics

For the Lotka-Volterra model, we have:

\[
\begin{aligned}
x_C(t, \mathbf{c}_C(\theta, \lambda)) &= \exp\left( \sum_k c_{C,k}(\theta, \lambda) \phi_k(t) \right), \\
x_B(t, \mathbf{c}_B(\theta, \lambda)) &= \exp\left( \sum_k c_{B,k}(\theta, \lambda) \phi_k(t) \right), \\
f_C(x_C, x_B, \theta) &= \alpha - \beta_C x_B, \\
f_B(x_C, x_B, \theta) &= \beta_B x_C - \delta.
\end{aligned}
\]

**Inner Criterion for Lotka-Volterra:**

\[
\begin{aligned}
J(\mathbf{c}_C, \mathbf{c}_B, \theta, \lambda) &= \sum_j \left[ (\ln C_j - \ln x_C(t_j, \mathbf{c}_C))^2 + (\ln B_j - \ln x_B(t_j, \mathbf{c}_B))^2 \right] \\
&\quad + \lambda \int \left[ \left( \frac{d}{dt} \ln x_C(t, \mathbf{c}_C) - (\alpha - \beta_C x_B(t, \mathbf{c}_B)) \right)^2 \right. \\
&\quad \left. + \left( \frac{d}{dt} \ln x_B(t, \mathbf{c}_B) - (\beta_B x_C(t, \mathbf{c}_C) - \delta) \right)^2 \right] dt
\end{aligned}
\]

## Implementation Details

### Basis Functions

- **Type:** B-splines of order 5 (degree 4 polynomials)
- **Number of Basis Functions:** 34 per variable
- **Knots:** Equally spaced between times 0 and 15

### Quadrature

- **Method:** Simpson's Rule
- **Quadrature Points:** 151 equally spaced points over the interval [0, 15]

### Optimization

- **Inner Optimization:** Minimization of \( J(\mathbf{c}_C, \mathbf{c}_B, \theta, \lambda) \) with respect to \( \mathbf{c}_C \) and \( \mathbf{c}_B \) using the `L-BFGS-B` method.
- **Outer Optimization:** Minimization of \( H(\theta, \lambda) \) with respect to \( \theta \) using the `L-BFGS-B` method.
- **Gradient Computation:** Symbolic differentiation using SymPy to compute necessary partial derivatives.

### Handling Numerical Issues

- **Hessian Matrix:** To address potential negative eigenvalues in the Hessian matrix during optimization, a scoring-type Hessian (dropping second partial derivatives) is used.
- **Parameter Constraints:** Parameters are estimated in their logarithmic form to ensure positivity.

### Convergence Criteria

- **Inner and Outer Loops:** Optimization loops continue until convergence based on predefined criteria (e.g., tolerance levels on parameter updates).

### Visualization

- **Plots:** At each iteration, plots of the current model fit against the data are generated for both prey and predator populations.
- **Parameter Trends:** The evolution of parameter estimates across different \( \lambda \) values is visualized to assess stability.

## Usage

### Prerequisites

Ensure you have the following Python libraries installed:

- NumPy
- Matplotlib
- SciPy
- SymPy

You can install them using pip:

```bash
pip install numpy matplotlib scipy sympy
```

### Running the Code

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/...
   ```

2. **Execute the Notebook:**

   Open the Jupyter Notebook file in Jupyter and run all cells to perform the analysis and visualize the results.


### Outputs

- **Synthetic Data Plots:** Visualizations of the true trajectories and noisy observations for both prey and predator populations.
- **Basis Function Plots:** Display of selected B-spline basis functions used in the approximation.
- **Model Fit Plots:** Iterative plots showing the estimated trajectories against the noisy data for each \( \lambda \) value.
- **Parameter Estimates Plot:** Evolution of parameter estimates across different \( \lambda \) values.

## Reference

Ramsay, James, and Giles Hooker. *Dynamic Data Analysis: Modeling Data with Differential Equations*. 1st ed., Springer Publishing Company, Incorporated, 2017. ISBN: 1493971883.

*Abstract:* This text focuses on the use of smoothing methods for developing and estimating differential equations following recent developments in functional data analysis and building on techniques described in Ramsay and Silverman (2005) Functional Data Analysis. The central concept of a dynamical system as a buffer that translates sudden changes in input into smooth controlled output responses has led to applications of previously analyzed data, opening up entirely new opportunities for dynamical systems. The technical level has been kept low so that those with little or no exposure to differential equations as modeling objects can be brought into this data analysis landscape. There are already many texts on the mathematical properties of ordinary differential equations, or dynamic models, and there is a large literature distributed over many fields on models for real world processes consisting of differential equations. However, a researcher interested in fitting such a model to data, or a statistician interested in the properties of differential equations estimated from data will find rather less to work with. This book fills that gap.

## Acknowledgements

This implementation is inspired by the methodologies described in *Dynamic Data Analysis: Modeling Data with Differential Equations* by Ramsay and Hooker (2017).

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
```