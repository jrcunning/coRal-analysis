This repository contains code to accompany the manuscript titled

### A dynamic bioenergetic model for coral-*Symbiodinium* symbioses and coral bleaching as an alternate stable state

by **Ross Cunning, Erik B. Muller, Ruth D. Gates, and Roger M. Nisbet**

which is available as a [pre-print at bioRxiv](https://doi.org/10.1101/120733), and has been submitted for peer-review at the *Journal of Theoretical Biology*.

In this manuscript, we describe and evaluate a model of coral growth and symbiosis dynamics across gradients of light, nutrients, and prey availability, with a module of photo-oxidative stress that allows for simulation of coral bleaching and recovery. The code for this model has been developed into an R package called **coRal** ([github.com/jrcunning/coRal](http://github.com/jrcunning/coRal)) to facilitate its further development and application by the scientific community. This repository contains the code used to describe the behavior of the model that is presented in the above manuscript. To use this code, you will need to install the coRal package, which can be accomplished by ```devtools::install_github("jrcunning/coRal")```.

#### Repository contents:

* **R/:** Contains R scripts to produce each figure in the manuscript (Fig*.R), and accessory scripts for computation.
* **data/:** Contains data used for model parameter estimation (described in the Supplementary Information).
* **img/:** Contains the figures presented in the manuscript.
* **ms/:** Contains the manuscript in R Markdown (ms.Rmd), and the rendered .tex and .pdf output.
* **supp/:** Contains the manuscript's Supplementary Information in R Markdown (supp.Rmd), and the rendered HTML output.
* **Makefile:** A script to reproduce the analysis and all components of this project. Executing make -B in this repository will re-run all simulations, generate figures, and knit the manuscript pdf and Supplementary Info HTML documents. Warning: many of these simulations are time consuming. The code has been optimized for use with multiple processors, and the entire build completes using 24 CPUs in a few hours.