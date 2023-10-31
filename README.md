<img src="Figures/Fig2color-rel-recon.png" align="left" width="250" height="450" alt="A21CSCASIS"/>

#Overview 
This repository contains the `R` code and data to reproduce analyses in the paper:

"A Twenty-First Century Structural Change in Antarctica’s Sea ice System"

by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock

# Requirements
The code is all in the `R` language for statistical computation. To run it you will need `R` installed on your
local machine as well as some standard `R` packages from `CRAN`. For an introduction to `R`, see [here](https://www.r-project.org/). 

# Installation
To explore the code, the easy way to clone the repository on your local machine and run the `R` code locally on your machine.
So first, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

# What is in the repository?
The repository is comprised of four components: At the top level are  separate files of `R` code.

* **Data** This folder contains the two data files.
   * `nsidcV4.RData`: A version of the recorded sea ice extents from the National Snow and Ice Daat Center (NSIDC).
   * `reconstructions_regional.RDS`: An ensemble of reconstructions of the monthly sea ice extent. There are 2500 in number. The reconstruction file is large (110Mb).
* **Figures**: When the code is run the figures are saved in this folder as PDF and PNG files.
* **ExampleOutput**: This folder contains example outputs resulting from code runs. They are here to validate your runs.
* **R**:  At the top level are separate files of `R` code.
   * `Nature-2023-07-12845A-Plots.R`: Running this `R` code will reconstruct the plots in the paper and save them to the folder `Figures`. This
   * takes about 13 seconds to run.
   * `Nature-2023-07-12845A-Events.R`: Running this `R` code will recompute the probabilities of the events analyzed in Section 3 of the paper and also create Supplementary Figure 3. This takes about 250 seconds to run.
   * `Nature-2023-07-12845A-dynamic-ARIMA.R`: Running this `R` code will fit the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper. It will also produce various Markov Chain Monte Carlo (MCMC) diagnostics and save plots of them to the folder `Figures`. This takes about 64 seconds to run. There is a setting in the file that takes about 2 hours to run and produces the results in the paper. The amount of computing was excessive, but the computational cost is minor.

See the following papers for more information and examples:

#### Methodology

* Maierhofer, Thomas J (2023) [Statistical Reconstruction of 20th Century Antarctic Sea Ice](https://escholarship.org/uc/item/33m3c3mn)
  Dissertation, University of California at Los Angeles. ProQuest ID: Maierhofer_ucla_0031D_21798. Merritt ID: ark:/13030/m5z68gkk.
