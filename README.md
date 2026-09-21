<!-- PROJECT LOGO -->
<br />
<h3 align="center">Quantifying Long-Range Dependence in Object-Valued Time Series</h3>

</div>

<!-- ABOUT THE PROJECT -->
## Abstract
Object-valued time series, including distributions, covariance matrices, networks, and compositions, lack the linear structure required by conventional autocovariance-based definitions of long memory. We develop an intrinsic framework for defining and estimating long memory in metric spaces of negative type. An isometric embedding yields a centered Hilbert-valued process, and the trace of its lag-covariance operator provides a signed, additive measure of temporal dependence. Crucially, this trace equals the difference between the marginal mean pairwise distance and expected lagged distance, so the framework and its estimators use only distances between the original objects. We define the memory parameter through the nonsummable decay of this trace and show that it coincides with the usual parameter in Hilbert-valued settings. Estimation uses Bartlett aggregates of distance-based lag measures. Estimating the common marginal distance from the same dependent sample induces a common-centering bias in these aggregates. We derive its finite-sample form and propose iterated block-difference corrections, log-ratio and multi-bandwidth log-slope estimators, and a localized self-consistency refinement. We establish consistency, document substantial bias reduction in simulations, and find long-memory evidence in foreign-exchange return distributions and U.S. electricity-generation compositions.

## Setup Instructions for Reproducibility

To ensure full reproducibility of the results presented in the paper, please configure your environment as follows:
### Software Requirements
The analysis was performed using **R version 4.5.1** and **Python 3.13.13**. The following R and Python packages and their specific versions are required to replicate the results:

* **R Packages:** ggplot2 (version 3.5.2), plot3D (version 1.4.1), transport (version 0.15-4), tidyr (version 1.3.1), parallel (version 4.5.1), foreach (version 1.5.2), doSNOW (version 1.0.20), readr (version 2.2.0), dplyr (version 1.1.4), tidyverse (version 2.0.0), patchwork (version 1.3.2).
* **Python Packages:** Standard library modules for Python 3.13.13 --- os, glob. Third-party packages --- numpy (version 2.4.4), pandas (version 3.0.2), scipy (version 1.17.1), pyfinancialdata (available via running the command 'pip3 install https://github.com/FutureSharks/financial-data/archive/master.zip
').



### Code
1. main_function.R --- basic R functions for numerical studies. 
2. simulation.R --- R functions to generate required simulation results.
3. data_download.ipynb --- Python functions to download and prepare the log return distribution data, with raw data obtained following <https://www.kaggle.com/datasets/arashnic/stock-data-intraday-minute-bar/data>.
4. real_data.R --- R functions to generate required real data analysis results.

### Data
1. US_energy.csv --- US energy generation composition dataset downloaded from <https://www.eia.gov/electricity/data/browser/>.
2. energy_final_residual.RData --- Preprocessed US energy generation composition data to remove trend and seasonality. For details, see Section 5.2 of the paper.

### Procedure to reproduce simulation results
Step 1. Run simulation.R to reproduce simulation results.

### Procedure to reproduce real data results
Step 1. Run data_download.ipynb to generate the required dataset for log return distributions.

Step 2. Run real_data.R to reproduce real data analysis results.




