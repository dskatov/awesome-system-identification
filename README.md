# Awesome System Identification: Parameters Fitting, System Identification and State Estimation in Python 

This repository is purposed to studying and developing both state-of-the-art and previously unseen methods of working with dynamical systems defined as systems of non-linear ordinary differential equations. The current state of affairs here is centered around the estimation of ODE systems' parameters, including large-scale systems of ODEs met in modern biochemistry and systems biology.

We've started with the studies in a set of Jupyter Notebooks, and are planning to expand this into a `pip` package. If you have an interest in the topic, and a word or two to tell about the subject, especially if you're coming from the angles of non-parametric methods, Bayessian processes, Gaussian processes, kernel methods, non-linear filtration, functional analysis etc, we'd like to invite you for a contribution.

## Considered Methods

In the task of estimating parameters of an ODE system `d(x_i)/dt = f_i(X|θ)`, `i=1..n`, a "common wisdom" approach is to formulate an objective in terms of ODE solution trajectories' `x_i(t)` proximity to a given known trajectory `x_i_0(t)`, usually given at some finite set of values of `t`. This assumes using a non-linear optimizer (L-BFGS-B, Nelder-Mead, genetic algorithm and the likes) with such objective function of proximity, which could surface as weighted RMSE, MLE, etc. This entails solving the ODE at each objective's call, and the objective's landscape could be highly non-regular and full of local minima. This sounds like a highly laborious path, and it in fact is. In literature, this approach is known as "trajectory matching" [1, 2].

However, it's desirable to be able to do better than that, as the right-hand sides of the equations are known in advance and explicitly given as formulas. As an immediate idea, one can formulate adjoint equations for so-called sensitivity variables `s_i_j = ∂(x_i)/∂(θ_j)`, `i=1..n`, `j=1..m` and solve an expended ODE with `m*n` additional equations of the form `d(s_i_j)/dt = ∂(f_i)/∂(x_j) * s_i_j + ∂(f_i)/∂(θ_j)`. As a result, additional high-precision and explicit gradient information becomes available to the chosen non-linear optimizer. A diagonal of the Hessian matrix for sensitivity variables may be computed as well for the second-order optimization methods. However, this doesn't improve the complexity of the objective's landscape and still requires solving of the much bigger-now system of ODEs at each call to the objective function. The difficulty of the objective's landscape requires it to run an optimizer for so-many starting points.

Another idea of a white-box use of the ODE system definition is to (1) build  splines through the points of known trajectories `x_i_0(t)`, `i=1..n` (preferably a monotone spline preserving the direction of the known trajectory points in the continuous spline between the points, i.e. preserving signs of derivatives, which is an idea applicable to both smoothing splines and interpolation splines), and (2) estimate the left-hand side derivatives' values over the span of `t=0..T_max` (say, at `T` distinct time points) in the defined ODE system. This yields a system of non-ODE equations in the amount of `n*T`, which solutions (or best-residual solutions) w.r.t. `θ_j`, `j=1..m`, aim to give a vector `θ` delivering a reasonable discrepancy between the r.h.s. and the l.h.s. of the original ODE system. This direction of research is called in literature as "gradient matching" [1], and it paves a path of estimating parameters of an ODE system without actually solving this ODE system. What is remarkable about this approach is that the resulting equational (non-ODE) system, in many of the biochemical and biological applications, is, at worst, polynomial (multinomial) of low degree w.r.t. to `θ_j`, thus the functional landscapes arising in the problems being solved with gradient matching are expected to be "not as bad" as the ones in trajectory matching approach.

### 1. Combined trajectory and gradient matching ("Collocation inference")

Development of the gradient matching approach leads us to combining trajectory matching and derivative matching in a trade-off style into a single objective, then exploring the (Pareto-frontier of) this trade-off. The use of splines, in terms of which the sought trajectories would be re-formulated, as a carrier of estimates to the values of `θ_j`, leads us to an optimization problem which is in some way recursive and, seemingly, can't be directly solved. However, the split of this problem into two and solving them iteratively one by another in a sequence (somewhat similar to the Expectation-Maximization algorithm or to the Forward-Backward algorithm of Baum and Welch), gives us a workable way towards the desired estimation of `θ_j`. We discovered that such idea is mostly developed in [1] (Ramsay, 2017) and may be further improved in our work. The presented approach exhibits fast and high-quality convergence to plausible matches of the fits `θ_j`, while being 2-3 orders of magnitude faster than the "common wisdom" raw state space search approach of trajectory matching to fits of somewhat comparable quality (this becomes worse with the growth of dimentionality). The trade-off between the matching of trajectory and the matching of the ODE equations is explored in the fullest, thus allowing to explore the entire space of fits in a continuous fashion, but not just by a few "lucky" fits. Another feature is that the method doesn't appear to be highly sensitive to the starting point of `θ` values, thus conceptually requiring less multi-starts in the process of sampling best parameter fits.

The aspect of convergence guarantees in the proposed approach is an active area of our research.

This method is also known as CollocInfer in the R language ecosystem [2].

As for the practical expansions of the method, a few are considered:

* Require the higher-order derivatives to match the equation, start with second-order gradient matching, so that both first and second derivatives of the state variables are matched, thus possibly improving the convergence speed and quality  
* Understand the influence of using the Hessian matrix and the selection of optimization method in the adjoint optimization problem
* Understand the influeuce of spline type selection on the method efficiency and convergence
* Expand the method to estimating the `θ` as a function of time `θ(t)`, as in [1]
* The method where the continuous ODE is solved via representing the solution as a finite sum over the basis functions is known in literature as orthogonal collocation. The aim could be to reconcile the method proposed here, which has some traits of the collocation method, with more advanced methods of orthogonal collocation (as seen below)

Below is a picture of the ODE trajectories for the set of parameters closely matching the given (known) trajectory. This fit was attained using the combined matching in 1 minute and 1 CPU core on i7-12800H for a random initial point. For the reference, a comparable match using only trajectory matching (and a special treatment of the estimated objective surface, which we don't discuss further here and which is aimed at improving the speed to convergence) was obtained in 10 minutes of run time on 40 CPU cores in parallel using compiled version of the ODE system and compiled version of the ODE solver (via Julia bridge to Python using `diffeqpy`):

![](CollocInfer.png)

Note how both the final ODE solution trajectory and spline shapes are plotted. It is typical to get the splines fitted very fast but that is not a guarantee of the ODE trajectories fitted all that fast too. Therefore, this is actually one point where the ODE system with the currently estimated parameters needs to be solved, and that is in between each of the two steps, merely for the validation of resulting trajectories. We would like to investigate further how selection of the spline basis could help get considerably closer to faster matching of the actual ODE trajectories to the tightness of the spline fitting, maybe without even needing to solve ODE for validation.

### 2. Pareto-sampling of the ODEs system's trajectories

Unlike demo setups with "known" reference trajectories pre-generated through the ODE system solution itself, it is hard to expect near-equal fit of all variables of the system to real-world known trajectories of the mechanistic ODE models of natural phenomena. Often it happens that, while some of the variables are describing the real process well, other variables are not, and this is because their corresponding right-hand sides in the ODEs system are not describing the real process well enough. In this case, it is informative to explore the Pareto-space of vectors of partial per-variable fits: i.e., we know that, while we can't attain tight fits on all variables, we at least can estimate how well we can fit by some sub-groups of these variables while sacrificing the others.

One approach is to fit not against the entire set of `n` variables, but against all or some of the `2^n` sub-groups of variables, while leaving the other variables away from the fitness objective. For systems with 10 or so variables, combined with a high-performance method like in p.1 above, it doesn't sound like a bad idea.

Another approach is to weight the terms of the objective function with a one-normed vector of weights sampled from a Dirichlet distribution. With doing so, we effectively sample around the Pareto-frontier instead of selecting it from the solution trajectories obtained by randomly sampled initial points.

Lastly, we could try incorporating the weighted Dirichlet-vector right into the combined matching method from p.1 above ("Collocation inference"), incorporating the weights right into the equations rather than just/only in the objective functions of the iterative method. This, in theory, could allow us estimate the whole shape of Pareto-frontier of fitness from the resulting equational system. This demands future research.

### 3. Parameter selection via nearest neighbors of trajectories

The idea is to generate many trajectories for different starting points (of initial conditions and parameters) and then match the desired trajectory against the collected obtained trajectories via nearest neighbors search in a reasonable similarity metric for time series. The starting points are selected from a pre-defined hypercube using Sobol sequences to sample the space in a more uniform fashion, not relying on random generator prior run-to-run. The method would rely on the ability to massively compute ODE solutions trajectories on GPUs.

### 4. Surrogate ML models for the estimation of parameters

In this approach (e.g. as in [3]), we train an ML model, such as a decision forest, on the resulting trajectories for some representative set of sampled parameters (and initial conditions). Then we utilize this ML model to give us an idea about a shape of the parameter space where the fits are best, and possibly sample a Pareto-frontier of the fits.

### 5. Simultaneous solving of the ODE and estimating of parameters

GEKKO [6] uses a simultaneous approach for solving dynamic optimization problems, including parameter estimation for ODE systems. This method discretizes the differential equations and solves the resulting nonlinear programming (NLP) problem. The key aspects of this approach are: (1) Discretization: The continuous-time ODE system is converted into a set of algebraic equations using methods like orthogonal collocation. (2) Nonlinear Programming: The discretized system, along with the parameter estimation objective, is formulated as a large-scale NLP problem. (3) Solver Integration: GEKKO leverages built-in solvers, particularly IPOPT (Interior Point Optimizer) and APOPT (Advanced Process Optimizer), to solve the resulting NLP problem.

While GEKKO itself doesn't introduce a novel parameter estimation method, it builds upon well-established techniques in the field of dynamic optimization and parameter estimation. Some key concepts and methods that form the theoretical foundation of GEKKO's approach include: (1) Simultaneous Discretization: This method, also known as direct transcription, is widely used in optimal control and parameter estimation problems. (2) Interior Point Methods: The IPOPT solver used by GEKKO implements interior point optimization algorithms. (3) Sequential Quadratic Programming: The APOPT solver in GEKKO uses SQP methods for solving NLP problems.

We've intensively tried the GEKKO setting for parameter estimation, which worked fairly well for simple, vanilla ODE systems like Lotka-Volterra examples, but struggled with real-world kinetic models of systems biology and biochemistry, whatever we tried. The main problem was that IPOPT solver (which was the only one that worked reliably and was capable of parameter estimation), which was simply get stuck in the solving without any way for getting it unstuck. We therefore proceeded to developing our own methods, which was a success.

The CollocInfer method for parameter estimation in ODEs uses a nested optimization approach, separating the problem into inner and outer optimizations. It approximates state variables with basis functions and incorporates ODEs as soft constraints in the objective function. The method balances data fit and equation fit using a smoothing parameter λ, which is progressively increased. In contrast, the simultaneous discretization approach with IPOPT formulates the entire problem as a single large-scale NLP, discretizing the time horizon and treating ODEs as hard constraints. While IPOPT can handle larger systems and is more general-purpose, it may struggle with parameter estimation in large ODE systems due to the increased problem size and complexity resulting from discretization. CollocInfer's success in these cases likely stems from its continuous representation of state variables, which reduces the problem dimensionality, and its nested optimization structure, which allows for more focused parameter updates. Additionally, the progressive increase of λ in CollocInfer may provide a more gradual path to convergence, avoiding the potential pitfalls of getting stuck in local optima that IPOPT might encounter when dealing with highly nonlinear, large-scale discretized systems.

### 6. Incorporating Methods from Dynamical Systems Theory and Exploring the Phase Space

While the above methods were focusing on exploitation of the state space, we should question what's available in the phase space of the ODE system. We could explore advanced analytical techniques to deepen the understanding of the ODE system's qualitative behavior. Bifurcation Analysis is a method to identify critical parameter values where the system's stability and dynamics change, involving steps like finding equilibrium points, linearizing the system, and analyzing eigenvalues using SymPy for symbolic computation. Poincaré Maps are introduced as a tool to reduce continuous dynamics to discrete maps, facilitating the study of periodic orbits and chaos by defining a Poincaré section, simulating trajectories, and recording their intersections. Additionally, the computation of Lyapunov Exponents may be used as a means to quantify the divergence of nearby trajectories, with practical implementation using the nolds library to detect chaotic behavior. The section also emphasizes leveraging SymPy for symbolic analysis, such as deriving the Jacobian matrix and performing parameter sweeps to find stability criteria.

### 7. Incremental parameter estimation of kinetic metabolic network models

We aim to adopt the approach outlined in [7] (2012). A robust and efficient technique for parameter estimation is crucial when constructing biological models based on ordinary differential equations (ODE). Conventional estimation methods typically involve simultaneously locating the global minimum of data fitting residuals across the entire parameter space. However, this approach often leads to prohibitively high computational demands due to the vast number of parameters and the challenge of complete parameter identifiability (i.e., the inability to uniquely identify all parameters).

The study [7] introduced an incremental strategy for estimating ODE model parameters from concentration time profiles. This method was specifically designed to address a common scenario in metabolic network modeling, where the number of metabolic fluxes (reaction rates) surpasses that of metabolites (chemical species). In this approach, the minimization of model residuals is conducted over a subset of the parameter space, focusing on the degrees of freedom in dynamic flux estimation derived from concentration time-slopes.

The effectiveness of this technique was showcased using two generalized mass action (GMA) models, demonstrating significant improvements over traditional single-step estimation methods. Furthermore, the authors presented an extension of their estimation approach to accommodate missing data.

This incremental estimation method effectively addresses the issue of incomplete parameter identifiability and substantially reduces the computational burden associated with model parameter estimation. These advancements are expected to facilitate the kinetic modeling of genome-scale cellular metabolism in future research endeavors.

### 8. Parameter estimation using Differential Variational Inequalities

The paper [5] (2006) presents a framework for parameter estimation in metabolic flux balance models of batch fermentation. The key idea is to formulate the fermentation dynamics as a system of differential variational inequalities (DVIs) that capture the interaction between cell metabolism and the changing environment. The DVIs are discretized and reformulated as a mathematical program with complementarity constraints (MPCC), which is then solved using an interior point algorithm. The authors apply this approach to estimate biomass composition parameters for a yeast fermentation model using both simulated and experimental data. The results show the method can accurately estimate parameters and fit fermentation profiles when using simulated data. For experimental data, the model provides reasonable fits to glucose and biomass measurements, though with some limitations in capturing biomass yields. Overall, the paper demonstrates the potential of this DVI-based framework for parameter estimation in complex fermentation models.

## Literature

1. Ramsay, James, Giles Hooker. *Dynamic Data Analysis: Modeling Data with Differential Equations*. 1st ed., Springer Publishing Company, Incorporated, 2017. ISBN: 1493971883. ([Link](https://www.amazon.com/Dynamic-Data-Analysis-Differential-Statistics/dp/1493971883))

2. Hooker, G., Ramsay, J. O., & Xiao, L. (2016). CollocInfer: Collocation Inference in Differential Equation Models. Journal of Statistical Software, 75(2), 1–52. https://doi.org/10.18637/jss.v075.i02 ([Link - Paper](https://www.jstatsoft.org/article/view/v075i02), [Link - R](https://cran.r-project.org/web/packages/CollocInfer/index.html))

3. David de Oliveira, Rafael & Procópio, Dielle & Basso, Thiago & Carrillo Le Roux, Galo. (2022). Parameter estimation in dynamic metabolic models applying a surrogate approximation. 10.1016/B978-0-323-95879-0.50036-9. ([Link](https://www.sciencedirect.com/science/article/abs/pii/B9780323958790500369))

4. Leppävuori, J., Domach, M. M., & Biegler, L. T. (2011). Parameter estimation in batch bioreactor simulation using metabolic models: Sequential solution with direct sensitivities. Industrial & Engineering Chemistry Research, 50(21), 12080-12091. https://doi.org/10.1021/ie201020g ([Link](https://pubs.acs.org/doi/10.1021/ie201020g))

5. Raghunathan, A.U., PÉRez-Correa, J.R., Agosin, E. et al. Parameter estimation in metabolic flux balance models for batch fermentation—Formulation & Solution using Differential Variational Inequalities (DVIs). Ann Oper Res 148, 251–270 (2006). https://doi.org/10.1007/s10479-006-0086-8 ([Link](https://link.springer.com/article/10.1007/s10479-006-0086-8))

6. Beal, L.D.R., Hill, D., Martin, R.A., and Hedengren, J. D., GEKKO Optimization Suite, Processes, Volume 6, Number 8, 2018, doi: 10.3390/pr6080106. Article - BibTeX - RIS

7. Jia, G., Stephanopoulos, G. & Gunawan, R. Incremental parameter estimation of kinetic metabolic network models. BMC Syst Biol 6, 142 (2012). https://doi.org/10.1186/1752-0509-6-142 ([Link](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-142))

8. https://en.wikipedia.org/wiki/Trajectory_optimization


## License

The MIT License. Copyright (c) 2024: Dan Skatov.
