<img src="Figures/Figure2.png" align="left" width="150" height="230" alt="A21CSCASIS"/>

# Overview 
This repository contains the `R` source code to reproduce analyses in the paper:

"[A Twenty-First Century Structural Change in Antarctica’s Sea ice
System](https://doi.org/10.1038/s43247-025-02107-5)"

by Marilyn N. Raphael, Thomas J. Maierhofer, Ryan L. Fogt, William R. Hobbs, and Mark S. Handcock. It appears in *Nature-Communications Earth & Environment*, 2025.

This repository is live and will be maintained to ensure the code runs correctly. A static version of this code that coincides with the publication of the paper is available at [a zenodo repository](https://doi.org/10.5281/zenodo.14741173).

* At the top level are separate files of `R` code.
   * `Nature-s43247-025-02107-5-Plots.R`: Running this `R` code will reconstruct the plots in the paper and save them to the folder `Figures`. This takes about two minutes to run.
   * `Nature-s43247-025-02107-5-Events.R`: Running this `R` code will recompute the probabilities of the events analyzed in Section 3 of the paper and also create Supplementary Figure 3. This takes about 250 seconds to run.
   * `Nature-s43247-025-02107-5-dynamic-ARIMA.R`: Running this `R` code will fit the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper. It will also produce various Markov Chain Monte Carlo (MCMC) diagnostics, save the cores results and save plots of them to the folder `Figures`. The example takes about 250 seconds to run, but the one run in the paper uses all 2500 ensemble members and takes a day to run. There is a setting in the file that takes about a day to run and produces the results in the paper. The amount of computing is excessive, but the overall computational cost is minor.
   * `Nature-s43247-025-02107-5-dynamic-ARIMA-analysis.R`: Running this `R` code will summarize the results of the dynamic auto-regressive integrated moving average model (ARIMA) analyzed in Section 4 of the paper (and the file `Nature-s43247-025-02107-5-dynamic-ARIMA-analysis.R`). It will also produce various diagnostics and save plots of them to the folder `Figures`.

See the following papers for more information and examples:

#### Methodology

* Raphael, Marilyn N., Maierhofer, Thomas J, Fogt, Ryan L., Hobbs, William R. and Handcock, Mark S. (2025). [A Twenty-First Century Structural Change in Antarctica's Sea Ice System](https://doi.org/10.1038/s43247-025-02107-5), Nature-Communications Earth & Environment, DOI:[10.1038/s43247-025-02107-5](https://doi.org/10.1038/s43247-025-02107-5).
* Maierhofer, Thomas J, Raphael, Marilyn N., Fogt, Ryan L., and Handcock, Mark S. (2024). [A Bayesian model for 20th century Antarctic sea ice extent reconstruction](https://doi.org/10.1029/2024EA003577). Earth and Space Science, 11, e2024EA003577. https://doi.org/10.1029/2024EA003577.
* Maierhofer, Thomas J (2023) [Statistical Reconstruction of 20th Century Antarctic Sea Ice](https://escholarship.org/uc/item/33m3c3mn)
  Dissertation, University of California at Los Angeles. ProQuest ID: Maierhofer_ucla_0031D_21798. Merritt ID: ark:/13030/m5z68gkk.
* Raphael, Marilyn N., Maierhofer, Thomas J, Fogt, Ryan L., Hobbs, William R. and Handcock, Mark S. (2024). [A Twenty-First Century Structural Change in Antarctica's Sea Ice System: Data and Code Repository](https://doi.org/10.5281/zenodo.14741173). Zenodo, version 1.0 [https://doi.org/10.1029/2024EA003577](https://doi.org/10.5281/zenodo.14741173).
* R Core Team (2021). [R: A language and environment for statistical computing](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
* R Core Team (2021). [The Comprehensive R Archive Network (CRAN)](https://www.R-project.org/). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
