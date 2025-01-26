<img src="../Figures/Fig2color-abs-with-legend_R2.png" align="left" width="150" height="230" alt="A21CSCASIS"/>

# Overview 
This repository contains the `R` source code and data to reproduce analyses in the paper:

"A Twenty-First Century Structural Change in Antarctica’s Sea ice System"

by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock. It will appear in *Nature-Communications Earth & Environment* in 2025.


The directory contains separate files of `R` code.
   * `Nature-COMMSENV-23-1836A-Plots.R`: Running this `R` code will reconstruct the plots in the paper and save them to the folder `Figures`. This takes about two minutes to run.
   * `Nature-COMMSENV-23-1836A-Events.R`: Running this `R` code will recompute the probabilities of the events analyzed in Section 3 of the paper and also create Supplementary Figure 3. This takes about 250 seconds to run.
   * `Nature-COMMSENV-23-1836A-dynamic-ARIMA.R`: Running this `R` code will fit the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper. It will also produce various Markov Chain Monte Carlo (MCMC) diagnostics, save the cores results and save plots of them to the folder `Figures`. The example takes about 250 seconds to run, but the one run in the paper uses all 2500 ensemble members and takes a day to run. There is a setting in the file that takes about a day to run and produces the results in the paper. The amount of computing is excessive, but the overall computational cost is minor.
   * `Nature-COMMSENV-23-1836A-dynamic-ARIMA-analysis.R`: Running this `R` code will summarize the results of the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper (and the file `Nature-COMMSENV-23-1836A-dynamic-ARIMA-analysis.R`). It will also produce various diagnostics and save plots of them to the folder `Figures`.

See the following papers for more information and examples:

#### Methodology

* Maierhofer, Thomas J, Raphael, Marilyn N., Fogt, Ryan L., and Handcock, Mark S. (2024). [A Bayesian model for 20th century Antarctic sea ice extent reconstruction](https://doi.org/10.1029/2024EA003577). Earth and Space Science, 11, e2024EA003577. https://doi.org/10.1029/2024EA003577.
* Maierhofer, Thomas J (2023) [Statistical Reconstruction of 20th Century Antarctic Sea Ice](https://escholarship.org/uc/item/33m3c3mn)
  Dissertation, University of California at Los Angeles. ProQuest ID: Maierhofer_ucla_0031D_21798. Merritt ID: ark:/13030/m5z68gkk.
* R Core Team (2021). [R: A language and environment for statistical computing](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
* R Core Team (2021). [The Comprehensive R Archive Network (CRAN)](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
