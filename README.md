# **Data and code for: Model projections of the epidemiological benefit of pyrethroid-pyrrole insecticide treated nets against malaria **
## **Authors: Thomas S. Churcher*, Isaac J. Stopard*, Arran Hamlet1, Dominic P. Dee, Antoine Sanou, Mark Rowland, Moussa W. Guelbeogo, Basiliana Emidi, Jacklin F. Mosha, Joseph D. Challenger, Adrian Denz, Andrew Glover,  Giovanni D. Charles, Emma L. Russell, Rich Fitzjohn, Pete Winskill, Christen Fornadel, Tom Mclean, Peder Digre, Joe Wagman, Frank Mosha, Jackie Cook, Martin C Akogbéto, Luc S. Djogbenou, Hilary Ranson, Philip McCall, Alphaxard Manjurano, Sagnon N’Falé, Natacha Protopopoff, Manfred Accrombessi, Corine Ngufor, Geraldine Foster, Ellie Sherrard-Smith **

This GitHub repository provides all code necessary to run the model validation for this paper.

:one: read_data.R - reads in the trial data and values needed to parameterise the model. 
:two: run_model.R - runs and saves the simulations.
:three: run_model_functions.R helper functions to run the model simulations.
:four: calc_net_top_up.R - estimates the net coverages for all the trial arms.
:five: retention_fit_top_up.R generic Stan model and code to estimate the coverages for a single trial arm.
:six: functions.R - helper functions to estimate the net coverages.
:seven: plot_model.R - generates the plots and extracts the required values from all simulations.
