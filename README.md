# HTEPredictionMetrics package
This package includes various metrics that can be used to evaluate heterogeneous treatment effect predictions.

# Install 

library(remotes)

remotes::install_github("CHMMaas/HTEPredictionMetrics")

Note, if you want to update the package you can also use the above call. If this doesn't work, use remove the package using remove.packages() and reopen R and run the above call.

If you receive the error: 'installation of package ‘~/HTEPredictionMetrics_1.0.0.tar.gz’ had non-zero exit status', please first install and update the necessary packages that HTEPredictionMetrics depends on.

# Version

Most up-to-date version: 1.0.0

# Load

library(HTEPredictionMetrics)

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
Some other useful functions included in this package.

# plot.calibration
This function enables to plot the calibration plot of treatment effect by plotting the observed versus the predicted treatment effect of matched pairs.

# Rubin.combine
This function combines multiple estimates with multiples standard errors using the Rubin's rule.

# References

van Klaveren D, Steyerberg EW, Serruys PW, Kent DM. The proposed 'concordance-statistic for benefit' provided a useful metric when modeling heterogeneous treatment effects. J Clin Epidemiol. 2018 Feb;94:59-68. doi: 10.1016/j.jclinepi.2017.10.021. Epub 2017 Nov 11. PMID: 29132832; PMCID: PMC7448760.

Maas C.C.H.M., Kent D.M., Hughes M.C., Dekker R., Lingsma H.F., van Klaveren D. Performance metrics for models designed to predict treatment effect. 2022. <a href="https://www.medrxiv.org/content/10.1101/2022.06.14.22276387v1">Preprint</a>.
