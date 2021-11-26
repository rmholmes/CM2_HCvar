# Plotting routines and processed data for ACCESS-CM2 adibatic/diabatic ocean temperature variability analysis

Contains data and code to plot the figures in the following manuscript:

Holmes, R.M., Sohail, T. and Zika, J.D. (2021): Adiabatic and diabatic signatures of ocean temperature variability, Journal of Climate, **XXX XXX**

Process_ACCESS_CM2.m contains the code to process raw ACCESS-CM2 output on NCI into .mat files in temperature, depth and latitude coordinates

Collate_ACCESS_CM2.m collates the output from Process_ACCESS_CM2 into a single .mat file.

PostProcess_ACCESS_CM2.m processes the data produced by Collate_ACCESS-CM2.m including percentile binning, de-drifting and removing the climatology.

ACCESS_SpecificHeat_PIcontrol_SWP.mat contains processed data from the ACCESS-CM2 PI control run. 

Plot_ACCESS_CM2.m contains the plotting code itself.

plot_budget_ACCESS_CM2.m contains the plotting code for the mean state heat budgets.

The raw ACCESS-CM2 data is available on NCI. Processed data is available on Zenodo at [https://dx.doi.org/10.5281/zenodo.5728574](https://dx.doi.org/10.5281/zenodo.5728574)

If you have any questions, please contact me at r.holmes@sydney.edu.au
