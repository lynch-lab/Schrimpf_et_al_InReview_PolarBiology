# Schrimpf_etal_InReview_PolarBiology
Source code for the multi-state occupancy models of Antarctic Peninsula breeding birds described in the paper

This folder contains code for the manuscript titled "Regional breeding bird assessment of the Antarctic Peninsula" which has been submitted to Polar Biology.

The raw data for the projected are in the file Raw_ASI_Data.csv.

The R script Running_model.R contains the code for running the occupancy model for a single species on JAGS, and the script Compiling_Results.R allows the user to re-package the results from the whole suite of species.

Basic output from our model runs can be found in the Model_Output folder, and in several of the .Rdata files in the root folder.

A copy of the electronic supplemental file ESM_5.csv contains both metadata about each of the 196 sites, as well as the basic breeding probability results from the model.

The remaining .R files allow the user to make maps of the results (Species_Maps.R), create the graphs from the posterior predictive check (Posterior_Pred_Check.R), and the species accumulation curves (Accumulation_Curves.R).

Any questions should be directed to Michael Schrimpf (michael.schrimpf@stonybrook.edu).
