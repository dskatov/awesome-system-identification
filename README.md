# README

## Efficient Sampling of the Pareto Frontier for ODE Parameter Fitting Using Sensitivity Variables and Gradient Matching

---

### **Overview**

This repository aim to present a novel approach to efficiently sample the Pareto frontier for multi-objective optimization in the context of parameter fitting for multivariable ODE systems. The method is specifically designed for large systems where computational efficiency is critical. By leveraging explicit sensitivity analysis via sensitivity variables method, linear approximations and modeling trajectories and derivatives data with splines, we avoid the computationally expensive task of continuously solving the ODE system inside an optimization loop, providing a practical solution for exploring parameter space and trade-offs between fitting different variables.

---

### **Motivation**

Fitting complex ODE systems to experimental data can be highly computationally demanding, especially when the system has a large number of parameters and state variables. Traditional methods that treat the ODE system as a black box require a large number of numerical integrations, making them infeasible for high-dimensional problems. This method addresses that challenge by:

- Utilizing explicit (symbolic) sensitivity analysis to linearize the relationship between parameters and state variables.
- Precomputing matrices that allow efficient exploration of parameter adjustments.
- Sampling the Pareto frontier using weighted sums, exploring the trade-offs between fitting different variables.

This approach is particularly well-suited for applications in physical and biochemical systems, where fitting accuracy varies across different state variables.

---

### **Methodology**

#### **1. Initial Estimation of Parameters Using Gradient Matching**

The process begins with an initial estimate of the parameters, using a **gradient matching** technique:

- **Spline Approximation:** Experimental data is interpolated using monotonic splines (such as PCHIP or Akima).
- **Derivative Estimation:** Derivatives are computed from the splines to estimate the time derivatives of the state variables.
- **Nonlinear System:** A nonlinear system is constructed based on the original ODE system, substituting spline values for state variables and their derivatives.
- **Solution:** This system is solved efficiently using nonlinear least squares, providing a good initial estimate of parameters vector.

#### **2. Formulating the Multi-Objective Optimization Problem**

Once the initial estimate of parameters vector is obtained, we formulate the multi-objective optimization problem. The objective is to minimize the discrepancies between the model and data for each state variable:

- For each variable, we compute the discrepancy between the experimental data and the model output.
- Using the sensitivity matrices, we approximate how parameter adjustments impact these discrepancies.
- The optimization problem is multi-objective, aiming to minimize the discrepancy for each variable independently, leading to trade-offs between variables.

#### **3. Precomputing Matrices for Efficient Computation**

To avoid repeatedly integrating the ODE system, we precompute key matrices for each state variable:

- **Sensitivity Matrix:** This captures how each state variable depends on the parameters.
- **Discrepancy Vector:** This measures the difference between the data and the model output for each variable.
- **Precomputed Matrices:** These matrices are used to set up a system of linear equations, which can be solved efficiently.

#### **4. Sampling the Pareto Frontier Using Weighted Sums**

The Pareto frontier is sampled using a weighted sum method:

- **Weighted Sum Objective:** We create a weighted sum of the objective functions, where each objective represents the discrepancy for a particular state variable.
- **Varying Weights:** By varying the weights applied to each objective, we generate different solutions, each corresponding to a different trade-off between variables.

#### **5. Efficient Computation of Solutions for Trade-Offs**

For each set of weights:

- We compute the weighted sensitivity and discrepancy matrices.
- A system of linear equations is solved to find the optimal parameter adjustments for the given weight configuration.
- This avoids the need for repeated numerical integration of the ODE system, significantly reducing computational cost.

#### **6. Analyzing and Visualizing the Pareto Frontier**

After generating a set of solutions, we analyze and visualize the Pareto frontier:

- **Objective Vectors:** The solutions are evaluated in terms of the discrepancies for each variable.
- **Dimensionality Reduction:** Techniques like UMAP are used to project high-dimensional objective vectors into 2D for visualization.
- **Clustering:** Clustering algorithms are applied to identify groups of solutions representing similar trade-offs.
- **Interpretation:** For each cluster, we analyze the parameter adjustments and identify patterns in the trade-offs between variables.

### **Key Features**

- **Computational Efficiency:** The method leverages precomputed sensitivity matrices and avoids repeated ODE integrations, making it computationally efficient even for large systems.
- **Trade-Off Exploration:** By systematically varying the weights in the objective function, we can explore different trade-offs between fitting various state variables, offering a comprehensive view of parameter space.
- **Scalability:** This approach is well-suited for high-dimensional systems, handling large numbers of state variables and parameters with ease.

---

### **Extensions and Enhancements**

- **Regularization:** To address ill-conditioning, regularization terms can be added to the system of equations, ensuring stable solutions.
- **Iterative Refinement:** For highly nonlinear systems, the method can be applied iteratively, updating the parameter estimates and recomputing the sensitivity matrices at each step.
- **Constraints:** The method can be extended to incorporate parameter constraints, such as bounds on parameter values or additional constraints based on prior knowledge.

### **Applications**

This method is particularly useful for applications in:

- **Biochemical Systems:** Where different components may fit data with varying degrees of accuracy due to complex interactions and measurement errors.
- **Physical Models:** Where certain variables may be more sensitive to parameter changes than others, leading to different trade-offs in model fitting.