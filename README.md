# BP Regression Residuals

This repository contains the R codes used to reproduce the simulation study and real data applications presented in the manuscript

> Justino, M. E. C., Pereira, T. L., & Souza, T. C. *Residual-Based Diagnostics for Assessing the Adequacy of Beta Prime Regression Models*. (Under review).

## Authors

- Maria Eduarda da Cruz Justino
- Tarciana Liberal Pereira
- Tatiene Correia de Souza

## Programming language

- R

## Repository structure

### `scripts/`

Simulation studies presented in the manuscript.

| File | Description |
|------|-------------|
| `residuals simulation - correct specification.R` | Simulation study under the correctly specified beta prime regression model. |
| `fixed precision error.R` | Simulation study under misspecification of the beta prime regression model due to fixed precision. |
| `error in the covariate of the mean submodel.R` | Simulation study under misspecification of the mean submodel due to omission of a relevant covariate. |

### `source codes/`

Supporting functions used by the simulation scripts and applications.

| File | Description |
|------|-------------|
| `Residual_H_Log_Like_BP.R` | Functions for residuals, hat matrix, log-likelihood, diagnostic measures, and diagnostic plots. |
| `gamlss_BP.R` | Beta prime family for the `gamlss` framework. |

### `applications/`

R scripts reproducing the real data applications presented in the manuscript.

| File | Description |
|------|-------------|
| `landrent application.R` | Application to the land rent data from Weisberg (2014). |
| `cherry application.R` | Application to the black cherry trees data from Ryan et al. (1976). |
| `cherry.dat` | Black cherry trees dataset used in `cherry application.R`. |

### Root directory

| File | Description |
|------|-------------|
| `README.md` | Repository description and usage information. |
| `reproducibility_info.txt` | R session information, package versions, and RNG settings used to ensure reproducibility. |
