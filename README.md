<img src="Figures/Figure2.png" align="left" width="150" height="230" alt="A21CSCASIS"/>

# Overview 
This repository contains the `R` source code and data to reproduce analyses in the paper:

"A Twenty-First Century Structural Change in Antarctica’s Sea ice System"

by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock. It will appear in *Nature-Communications Earth & Environment* in 2025.

This repository is live and will be maintained to ensure the code runs correctly. A static version of this code that coincides with the publication of the paper is available at [a zenodo repository](https://doi.org/10.5281/zenodo.14741173).

# System Requirements
The code is all in the `R` language for statistical computation. It will run on any system that supports `R` (e.g., Windows, Macintosh, Linux).
Version 4.0 or higher of `R` is required to be installed on your
local machine, as well as some standard `R` packages from `CRAN`. For an introduction to `R`, see [here](https://www.r-project.org/). 

# Installation
To explore the code, the easy way to clone the repository on your local machine and run the `R` code locally on your machine.
So first, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).
The install time should be a few seconds. 

# What is in the repository?
The repository is comprised of four components: At the top level are  separate files of `R` code.

* **Data** This folder contains the data files:
   * `nsidcV4.RData`: A version of the recorded sea ice extents from the National Snow and Ice Data Center (NSIDC).
   * `reconstructions`: An ensemble of reconstructions of the monthly sea ice extent for each of the five sectors. There are 2500 in number. The reconstruction file is large (110Mb) so is read in from the cloud using a URL in the `R` files.
   * `reconstructions_total`: An ensemble of reconstructions of the monthly sea ice extent for all of Antarctica. The reconstruction file is read in from the cloud using a URL in the `R` files.
   * `persistence.RDS`: Summary persistence measures for each member of the ensemble of reconstructions of the monthly sea ice extent for the total Antarctica. These are computed by `Nature-COMMSENV-23-1836A-dynamic-ARIMA.R` and are provided here to speed the production of the plots as are separately of interest.
   * `persistence-[X].RDS`: Summary persistence measures for each member of the ensemble of reconstructions of the monthly sea ice extent for the Xth sector of Antarctica. Here X is one of 1, 2, 3, 4, 5, corresponding to "King Haakon VII", "Ross-Admundsen", "East Antarctica", "Weddell", and "Bellingshausen Amundsen", respectively. These are computed by `Nature-COMMSENV-23-1836A-dynamic-ARIMA.R` and are provided here to speed the production of the plots as are separately of interest.
    * `pos_pred_fits`: Posterior predictive distributions for the satellite observed period. These are provided here to speed the production of the plots as are separately of interest. The file is large so is read in from the cloud using a URL in the `R` files.
    * `pos_pred_fits_brief.RData`: Summary statistics of the posterior predictive distributions for the satellite observed period. These are provided here to speed the production of the plots as are separately of interest.
* **Figures**: When the code is run the figures are saved in this folder as PDF and PNG files.
* **ExpectedOutput**: This folder contains example outputs resulting from code runs. They are here to validate your runs.
* **R**:  At the top level are separate files of `R` code.
   * `Nature-COMMSENV-23-1836A-Plots.R`: Running this `R` code will reconstruct the plots in the paper and save them to the folder `Figures`. This takes about two minutes to run.
   * `Nature-COMMSENV-23-1836A-Events.R`: Running this `R` code will recompute the probabilities of the events analyzed in Section 3 of the paper and also create Supplementary Figure 3. This takes about 250 seconds to run.
   * `Nature-COMMSENV-23-1836A-dynamic-ARIMA.R`: Running this `R` code will fit the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper. It will also produce various Markov Chain Monte Carlo (MCMC) diagnostics, save the cores results and save plots of them to the folder `Figures`. The example takes about 250 seconds to run, but the one run in the paper uses all 2500 ensemble members and takes a day to run. There is a setting in the file that takes about a day to run and produces the results in the paper. The amount of computing is excessive, but the overall computational cost is minor.
   * `Nature-COMMSENV-23-1836A-dynamic-ARIMA-analysis.R`: Running this `R` code will summarize the results of the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper (and the file `Nature-COMMSENV-23-1836A-dynamic-ARIMA-analysis.R`). It will also produce various diagnostics and save plots of them to the folder `Figures`.

See the following papers for more information and examples:

#### Methodology

* Maierhofer, Thomas J, Raphael, Marilyn N., Fogt, Ryan L., and Handcock, Mark S. (2024). [A Bayesian model for 20th century Antarctic sea ice extent reconstruction](https://doi.org/10.1029/2024EA003577). Earth and Space Science, 11, e2024EA003577. https://doi.org/10.1029/2024EA003577.
* Maierhofer, Thomas J (2023) [Statistical Reconstruction of 20th Century Antarctic Sea Ice](https://escholarship.org/uc/item/33m3c3mn)
  Dissertation, University of California at Los Angeles. ProQuest ID: Maierhofer_ucla_0031D_21798. Merritt ID: ark:/13030/m5z68gkk.
* Raphael, Marilyn N., Maierhofer, Thomas J, Fogt, Ryan L., Hobbs, William R. and Handcock, Mark S. (2024). [A Twenty-First Century Structural Change in Antarctica's Sea Ice System: Data and Code Repository](https://doi.org/10.5281/zenodo.14741173). Zenodo, version 1.0 [https://doi.org/10.1029/2024EA003577](https://doi.org/10.5281/zenodo.14741173).
* R Core Team (2021). [R: A language and environment for statistical computing](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
* R Core Team (2021). [The Comprehensive R Archive Network (CRAN)](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
