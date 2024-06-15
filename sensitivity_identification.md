# Efficient Sampling of the Pareto Frontier for ODE Parameter Fitting Using Sensitivity Analysis

---

## **Introduction**

In the modeling of complex systems using ordinary differential equations (ODEs), parameter estimation is a critical task. Often, not all state variables can be fitted equally well to experimental data due to unmodeled phenomena, measurement errors, or system complexity. This leads to trade-offs between fitting different variables. To address this, we seek an efficient method to explore these trade-offs by sampling the Pareto frontier of a multi-objective optimization problem, where the objectives are the discrepancies between the model and data for each variable.

To make this approach feasible, we leverage sensitivity analysis and linear approximations, avoiding the repeated numerical integration of the ODE system for every parameter adjustment. This method offers a computationally efficient way to explore parameter space and sample the Pareto frontier, ultimately allowing the user to examine different trade-offs in fitting various state variables to data.

The proposed method avoids the need for continuous numerical integration of the ODE system within optimization loops, which is computationally expensive. Instead, an initial estimate of the parameter vector \( A_0 \) is obtained by approximating the desired trajectories and their derivatives from data using monotonic splines (such as PCHIP, Akima, or similar). This estimate is found by solving a nonlinear system derived from the original ODEs with respect to the parameters \( A \), using the spline values and their derivatives in place of the state variables and their time derivatives. In physical and biochemical applications, the resulting system is typically linear or polynomial in the parameters, allowing it to be efficiently solved using nonlinear least squares techniques.

---

## **Motivation**

Given the complexity of fitting high-dimensional ODE systems to data, traditional methods that treat the ODE system as a black box require a large number of numerical integrations, making them computationally expensive. To avoid this, we:

1. **Leverage Sensitivity Analysis:** Precompute the sensitivity matrices for each variable with respect to the parameters, allowing us to formulate the fitting process as a linear system.
   
2. **Efficiently Sample the Pareto Frontier:** By using weighted sums of the fitting objectives, we can efficiently explore trade-offs between variables and sample the Pareto frontier of parameter sets that provide different trade-offs.

---

## **Objective**

Develop a method that efficiently samples the Pareto frontier in the context of ODE parameter fitting using sensitivity analysis, linear approximations, and a precomputed system. This method will avoid the need for continuous ODE integrations, providing a practical solution for large-scale systems.

---

## **Overview of the Method**

1. **Initial Estimation of Parameters \( A_0 \) Using Gradient Matching**
2. **Formulate the Multi-Objective Optimization Problem Using Sensitivities**
3. **Precompute Matrices for Efficient Computation**
4. **Sample the Pareto Frontier Using the Weighted Sum Method**
5. **Efficiently Compute Solutions for Different Trade-offs**
6. **Analyze and Visualize the Pareto Frontier**

---

## **Step 1: Initial Estimation of Parameters \( A_0 \) Using Gradient Matching**

Before solving the multi-objective problem, we first compute an initial guess \( A_0 \) for the parameter set using a **gradient matching** approach.

### **Process:**

1. **Spline Approximation of Data:**
   - Use **PCHIP splines** to interpolate the noisy data points for each variable. This gives smooth approximations of the data.
   
2. **Estimate Derivatives from Splines:**
   - Compute the derivatives of the splines, providing estimates for the time derivatives of the variables \( \frac{d\text{data}_i}{dt} \).
   
3. **Minimize the Discrepancy Between the Model and Data Derivatives:**
   - Solve a least-squares optimization problem, where the objective is to match the time derivatives computed from the splines to the model's ODE right-hand sides.
   - The optimization problem is:
   
     \[
     \min_{A} \sum_{i=1}^N \| \dot{\text{data}}_i(t_k) - f_i(X, A, t_k) \|^2
     \]
   
   - The resulting parameter set \( A_0 \) is used as a starting point for further optimization.

### **Key Advantage:**
- This provides a **reasonable initial guess** for the parameters, ensuring that the subsequent exploration of the parameter space starts from a point close to a good fit.

---

## **Step 2: Formulate the Multi-Objective Optimization Problem Using Sensitivities**

Once \( A_0 \) is obtained, we formulate the multi-objective optimization problem using sensitivity analysis to efficiently explore the trade-offs between variables.

### **Given:**

- **ODE System:**

  \[
  \frac{dX}{dt} = f(X, A, t)
  \]

  where \( X(t) = [x_1(t), x_2(t), \dots, x_N(t)]^\top \) are the state variables, and \( A = [a_1, a_2, \dots, a_m]^\top \) is the parameter vector.

- **Data:**

  Experimental data for each variable \( \text{data}_i(t_k) \) at time points \( t_k \), for \( k = 1, \dots, \tau \).

- **Sensitivities:**

  \[
  S_{ij}(t) = \frac{\partial x_i(t)}{\partial a_j}
  \]

  These measure how changes in parameter \( a_j \) affect variable \( x_i \) at time \( t \).

### **Discrepancy for Variable \( i \):**

We compute the discrepancy between the data and the model for variable \( i \) as:

\[
D_i = [D_i(t_1), D_i(t_2), \dots, D_i(t_\tau)]^\top = \text{data}_i(t_k) - x_i(t_k; A_0)
\]

### **Linear Approximation Using Sensitivities:**

Using the sensitivities \( S_{ij}(t_k) \), we approximate the change in the discrepancy when adjusting the parameters \( A \) from \( A_0 \):

\[
D_i \approx S_i \Delta A
\]

where \( S_i \) is the sensitivity matrix for variable \( i \), defined as:

\[
S_i = \begin{bmatrix}
S_{i1}(t_1) & S_{i2}(t_1) & \dots & S_{im}(t_1) \\
S_{i1}(t_2) & S_{i2}(t_2) & \dots & S_{im}(t_2) \\
\vdots & \vdots & \ddots & \vdots \\
S_{i1}(t_\tau) & S_{i2}(t_\tau) & \dots & S_{im}(t_\tau)
\end{bmatrix}
\]

### **Objective for Variable \( i \):**

The objective for each variable is to minimize the squared discrepancy:

\[
f_i(\Delta A) = \| D_i - S_i \Delta A \|^2
\]

### **Multi-Objective Optimization Problem:**

The overall optimization problem is to minimize the discrepancies for all variables simultaneously:

\[
\min_{\Delta A} \quad [f_1(\Delta A), f_2(\Delta A), \dots, f_N(\Delta A)]
\]

Our goal is to explore trade-offs between different objectives (i.e., fitting different variables) by sampling the Pareto frontier.

---

## **Step 3: Precompute Matrices for Efficient Computation**

We precompute matrices that will allow us to efficiently solve the multi-objective optimization problem without having to re-integrate the ODE system.

For each variable \( i \):

1. **Sensitivity Matrix \( S_i \):**
   This is already defined as:

   \[
   S_i = \begin{bmatrix}
   S_{i1}(t_1) & S_{i2}(t_1) & \dots & S_{im}(t_1) \\
   \dots & \dots & \dots & \dots \\
   S_{i1}(t_\tau) & S_{i2}(t_\tau) & \dots & S_{im}(t_\tau)
   \end{bmatrix}
   \]

2. **Discrepancy Vector \( D_i \):**

   \[
   D_i = \text{data}_i(t_k) - x_i(t_k; A_0)
   \]

3. **Precompute the Following Matrices:**

   - **Matrix \( M_i \):**

     \[
     M_i = S_i^\top S_i \quad \text{(size } m \times m)
     \]

   - **Vector \( b_i \):**

     \[
     b_i = S_i^\top D_i \quad \text{(size } m \times 1)
     \]

By precomputing these matrices, we significantly reduce the computational cost of subsequent steps. We only need to solve linear systems involving these precomputed matrices, avoiding the need for continuous numerical integration of the ODE system.

---

## **Step 4: Sample the Pareto Frontier Using the Weighted Sum Method**

We use a **weighted sum scalarization** approach to sample the Pareto frontier. This method allows us to explore different trade-offs between fitting the various state variables by adjusting the weights in the scalarization.

### **Weighted Sum Objective:**

\[
\Phi(\Delta A; w) = \sum_{i=1}^N w_i f_i(\Delta A)
\]

where \( w = [w_1, w_2, \dots, w_N]^\top \) are non-negative weights (\( w_i \geq 0 \)).

By varying the weights \( w_i \), we can prioritize the fitting of certain variables over others and generate different solutions along the Pareto frontier.

Use the Dirichlet distribution for uniform sampling.

---

## **Step 5: Efficiently Compute Solutions for Different Trade-offs**

For each set of weights \( w \), we can efficiently compute the corresponding parameter adjustment \( \Delta A \) by solving a linear system.

### **1. Compute Weighted Matrices:**

- **Weighted \( M \) Matrix:**

 

 \[
  M(w) = \sum_{i=1}^N w_i M_i = \sum_{i=1}^N w_i S_i^\top S_i
  \]

- **Weighted \( b \) Vector:**

  \[
  b(w) = \sum_{i=1}^N w_i b_i = \sum_{i=1}^N w_i S_i^\top D_i
  \]

### **2. Solve the Linear System:**

Solve the following system to find \( \Delta A \):

\[
M(w) \Delta A = b(w)
\]

This system is of size \( m \times m \) and can be solved efficiently using standard linear algebra techniques.

### **3. Compute Residuals and Objectives:**

For each variable \( i \):

- **Residuals:**

  \[
  r_i = D_i - S_i \Delta A
  \]

- **Objective Value:**

  \[
  f_i(\Delta A) = \| r_i \|^2
  \]

### **4. Collect Solutions:**

For each set of weights \( w \), store the parameter adjustment \( \Delta A \) and the corresponding objective values \( [f_1(\Delta A), f_2(\Delta A), \dots, f_N(\Delta A)] \).

---

## **Step 6: Analyze and Visualize the Pareto Frontier**

Once a set of solutions has been computed, we analyze and visualize the Pareto frontier.

### **1. Collect Objective Vectors:**

For each solution \( \Delta A \), obtain the corresponding objective vector:

\[
F(\Delta A) = [f_1(\Delta A), f_2(\Delta A), \dots, f_N(\Delta A)]
\]

### **2. Visualization:**

- **Pareto Front Plotting:** In 2D or 3D (for 2 or 3 objectives), plot \( f_i(\Delta A) \) directly to visualize trade-offs.

- **Dimensionality Reduction:** Use techniques like **UMAP** (Uniform Manifold Approximation and Projection) to project higher-dimensional objective vectors into 2D for visualization.

### **3. Clustering:**

Apply clustering algorithms (e.g., **K-Means**, **DBSCAN**) to the objective vectors or reduced-dimensional projections to identify clusters of solutions representing similar trade-offs.

### **4. Interpretation:**

For each cluster:

- Analyze the typical parameter adjustments \( \Delta A \).
- Identify which variables are fitted well or poorly in each cluster.
- Present typical or average trajectories for each cluster.

### **5. Present Results:**

Provide representative parameter sets and trajectories for each cluster, offering the user a clear view of the different trade-offs present in the system.

---

## **Mathematical Justification**

### **Weighted Least Squares Solution**

Given the weighted sum objective:

\[
\Phi(\Delta A; w) = \sum_{i=1}^N w_i \| D_i - S_i \Delta A \|^2
\]

The normal equations for minimizing this objective are:

\[
\left( \sum_{i=1}^N w_i S_i^\top S_i \right) \Delta A = \sum_{i=1}^N w_i S_i^\top D_i
\]

This simplifies to:

\[
M(w) \Delta A = b(w)
\]

Thus, for any set of weights \( w \), we can efficiently compute the parameter adjustment \( \Delta A \) by solving this linear system.

---

## **Efficiency Considerations**

### **Matrix Computations:**

- Matrix additions and scalar multiplications are computationally inexpensive.
- Solving linear systems is efficient, especially with modern linear algebra libraries.

### **Avoiding ODE Integration:**

- The sensitivity matrices and discrepancy vectors are computed once using the initial parameter estimate \( A_0 \).
- This method relies on linear approximations, allowing us to avoid repeated numerical integrations of the ODE system.

---

## **Extensions and Enhancements**

### **Regularization:**

To handle ill-conditioned systems, add a regularization term:

\[
M_{\text{reg}}(w) = M(w) + \lambda I
\]

where \( \lambda \) is a small positive constant and \( I \) is the identity matrix.

### **Constraints Handling:**

If parameter bounds or constraints are present, use constrained optimization methods, such as quadratic programming or constrained least squares solvers.

### **Iterative Refinement:**

For highly nonlinear systems, iterative refinement can be employed:

1. **Update Parameters:**

   \[
   A_{\text{new}} = A_0 + \Delta A
   \]

2. **Recompute Sensitivities:** 

   Recompute the sensitivity matrices \( S_i \) at the new parameter set \( A_{\text{new}} \).

3. **Repeat the Process:** 

   Iterate until convergence criteria are met.

---

## **Conclusion**

By leveraging sensitivity analysis and linear approximations, we have developed an efficient method for sampling the Pareto frontier in ODE parameter fitting. This approach avoids the computational expense of continuously solving the ODE system and provides an efficient way to explore trade-offs between fitting different state variables.

---

## **Key Advantages**

- **Computational Efficiency:** 

  The method avoids repeated numerical integration, focusing on solving linear systems instead, which is computationally efficient.

- **Exploration of Trade-offs:** 

  By systematically varying the weights in the objective function, we can explore how improving the fit for one variable impacts the fit of others.

- **Scalability:** 

  This method is suitable for high-dimensional systems with many variables and parameters.

---

## **Next Steps**

1. **Implementation:**
   - Develop the full implementation of the method in a Python environment.
   - Ensure that the method is numerically stable and handles large-scale systems efficiently.

2. **Validation:**
   - Validate the method with synthetic and real data.
   - Compare the results to traditional methods of parameter fitting.

3. **Extension:**
   - Incorporate additional constraints and regularization techniques as needed.
   - Explore iterative refinement for nonlinear systems.
   - Extend the method to include uncertainty quantification.

---

## **References**

1. Miettinen, K. (1999). *Nonlinear Multiobjective Optimization*. Springer.
2. Saltelli, A., et al. (2008). *Global Sensitivity Analysis: The Primer*. Wiley.
3. Das, I., & Dennis, J. E. (1998). Normal-Boundary Intersection: A New Method for Generating the Pareto Surface in Nonlinear Multicriteria Optimization Problems. *SIAM Journal on Optimization*, 8(3), 631-657.

---

**Note:** This mathematical formulation provides the foundation for an efficient method of sampling the Pareto frontier in ODE parameter fitting using sensitivity analysis. While this method is designed to minimize computational effort, care must be taken to ensure the linear approximation remains valid.

Copyright (c) 2024: Dan Skatov.