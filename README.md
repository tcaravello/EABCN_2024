# EABCN_2024
Slides and codes for practice sessions. Contents:

_data: contains all data used by the codes for Sessions 1 and 2 (Session 3 does not require data).

1_estim_nonlin: the /estim subfolder contains code to run VAR, LPs (standard and bias corrected), Bayesian VAR with prior selection as in Giannone, Lenza & Primiceri (2015), and Smoothed Local Projections (Barnichon & Brownlees, 2019)  
 the /sign_size subfolder contains replication codes to test for sign and size non-linearities, compute impulse responses and implicit weight schemes as in Caravello & Martinez-Bruera (2024)

2_mp_modelcnfctls: contains all necessary code to run the procedure in Caravello, McKay and Wolf (2024). Specifically:
 - /var_inputs: estimates Wold representation, forecasts and the IRF for the Aruoba-Drechsel monetary policy shock using VARs.
 - /irf_matching: obtains draws of parameters from the posterior, as well as the implied monetary shock causal effect matrices. Also computes posterior model probabilities.
 - /applications: runs the three applications in Section 5 of the paper.
 - /suff_stats: stores the draws for the monetary causal effect matrices. Due to space, we share a small number of draws for each case (50, whereas we use 1000 in the paper).
- /auxiliary_functions: contains auxiliary files.

3_micro_macro/olg : Computes the C_y mapping of an OLG model for different parameter values. Also computes different transformations corresponding to different behavioral models following Auclert, Rognilie & Straub (2020). 
