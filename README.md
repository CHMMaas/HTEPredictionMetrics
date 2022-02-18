# HTEPredictionMetrics package
This package includes various metrics that can be used to evaluate heterogeneous treatment effect predictions.

# Install 

library(remotes)

remotes::install_github("CHMMaas/HTEPredictionMetrics")

# E.for.Benefit
This function calculates the calibration metrics, i.e. the E-for-benefit metrics: Eavg-for-benefit, E50-for-benefit, and E90-for-benefit, as proposed by C.C.H.M. Maas et al. (2022) by matching patients based on patient characterstics or individualized treatment effect predictions. 
Calibration can be assessed by a smoothed calibration curve obtained by a local regression, with default values for span and degree of polynomials, of the observed treatment effect of matched patient pairs on predicted treatment effect of matched patient pairs.
The Eavg-for-benefit, E50-for-benefit, and E90-for-benefit were defined as the weighted average, median, and 90th percentile of the absolute differences between the smoothed calibration curve and the diagonal line of perfect calibration.

For an example of how to use the "E.for.Benefit" function, have look at "Examples" under "?E.for.Benefit".

# C.for.Benefit
This function calculates the C-for-benefit as proposed by D. van Klaveren et al. (2018) by matching patients based on patient characterstics or individualized treatment effect predictions.

The C-for-benefit measures the discriminative ability of a model that estimates heterogeneous treatment effect. 
The C-for-benefit is the probability that from two randomly chosen matched patient pairs with unequal observed treatment effect, the pair with greater observed treatment effect also has a higher predicted treatment effect. 
Observed treatment effect was defined as the difference between outcomes in pairs of patients matched on patient characteristics or on individualized treatment effect predictions. 
Predicted treatment effect of a matched pair was defined as the difference between the predicted outcome probability of the untreated patient minus the predicted outcome probability of the treated patient.
The C-for-benefit was calculated by the number of concordant pairs divided by the number of concordant and discordant pairs. 
Two patient pairs are concordant if the pair with the greater observed benefit also has higher predicted benefit. 
Two patient pairs are discordant if the pair with greater observed benefit has lower predicted benefit. Two patient pairs are uninformative if the pairs have the same observed benefit.

For an example of how to use the "C.for.Benefit" function, have look at "Examples" under "?C.for.Benefit".

# OP.for.Benefit
This function calculates the overall performance metrics, i.e. the logistic-loss-for-benefit and Brier-for-benefit, as proposed by C.C.H.M. Maas et al. (2022) by matching patients based on patient characterstics or individualized treatment effect predictions.
The log-loss-for-benefit was defined as the average logarithmic distance between predicted and observed treatment effect of matched patient pairs. 
The Brier-for-benefit was defined as the average squared distance between predicted and observed treatment effect of matched patient pairs.

For an example of how to use the "OP.for.Benefit" function, have look at "Examples" under "?OP.for.Benefit".

# Other functions

# val.surv.mi
This function can calculate the intercept, slope, and C-index for survival outcomes for complete case analysis and predictions from multiple imputed data sets.

# rubin.combine
This function combines multiple estimates with multiples standard errors using the Rubin's rule.

# References

van Klaveren D, Steyerberg EW, Serruys PW, Kent DM. The proposed 'concordance-statistic for benefit' provided a useful metric when modeling heterogeneous treatment effects. J Clin Epidemiol. 2018 Feb;94:59-68. doi: 10.1016/j.jclinepi.2017.10.021. Epub 2017 Nov 11. PMID: 29132832; PMCID: PMC7448760.

Maas C.C.H.M., van Klaveren D. Performance metrics for models predicting individualized treatment effect of patients. 2022. Not published yet.
