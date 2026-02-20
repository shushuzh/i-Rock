# i-Rock: Expected Shortfall Regression via Optimization

This repository contains the R and C++ implementation of the **i-Rock** methodology, a framework for robust statistical modeling of large-scale datasets via Expected Shortfall (ES) estimation.

## Repository Structure

* **`src/` (Core Methodology):**
    * `mRock_KNN_disjoint.cpp`: C++ source for computationally efficient kNN-based non-parametric ES estimation.
    * `INT_ES_disjoint.R`: Primary R implementation of the i-Rock formulations. The function `mRock_KNN_q0_both` provides the specific kNN-based ES estimation utilized in the manuscript.
    * `2step.R`: Implementation of the weighted two-step linear ES regression approach (Barendse, 2020) used for comparison.
* **`simulations/` (Replication):**
    * Contains three examples demonstrating the framework’s performance across various data-generating processes. These scripts also serve as templates for using the `mRock_KNN_q0_both` function.
    * **Replication:** To replicate the manuscript results, use parallel computing to call the functions across the sample sizes $n$ and tail levels $\tau$ specified in the text, using 500 replications (seeds 1–500).
* **`data_applications/` (Empirical Analysis):**
    * `data_application_processing.R`: Parses raw NCHS Natality data, cleans/recodes variables (e.g., birth weight, WIC, smoking status), and outputs the design matrix.
    * `data_application_final.R`: Executes the i-Rock approach on the processed dataset to obtain final ES regression results.
    * `data_application_bootstrap.R`: Use bootstrap to get standard error of the method. 

## Workflow for Reproducibility

*   **Environment Setup:** Ensure R and a C++ compiler are installed. Install dependencies: `Rcpp`, `RcppArmadillo`, `conquer`, `quantreg`, and `FNN`.
*   **Simulation:** 
    * Run the three R scripts in `/simulations` to generate the performance metrics for the three synthetic examples.
    * Use `Process_data_final.ipynb` to generate figures and tables in the main manuscript and supplementary materials. 
*   **Data Application:**
    * Run `data_application_processing.R` to prepare the data.
    * Run `data_application_final.R` to produce the regression coefficients.
    * Run `data_application_bootstrap.R` in parallel to estimate the standard error.
    * Use `Process_data_application.ipynb` to generate figures and tables in the main manuscript and supplementary materials. 

**References**

Barendse, S. (2020). Efficiently weighted estimation of tail and interquantile expectations. *Tinbergen Institute Discussion Paper 2017-034/III*.
