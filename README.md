# C.for.Benefit from the HTEPredictionMetrics package
This function calculates the C-for-benefit as proposed by D. van Klaveren et al. (2018) by matching patients based on patient characterstics or individualized treatment effect predictions.

The C-for-benefit measures the discriminative ability of a model that estimates heterogeneous treatment effect. 
The C-for-benefit is the probability that from two randomly chosen matched patient pairs with unequal observed treatment effect, the pair with greater observed treatment effect also has a higher predicted treatment effect. 
Observed treatment effect was defined as the difference between outcomes in pairs of patients matched on patient characteristics or on individualized treatment effect predictions. 
Predicted treatment effect of a matched pair was defined as the difference between the predicted outcome probability of the untreated patient minus the predicted outcome probability of the treated patient.
The C-for-benefit was calculated by the number of concordant pairs divided by the number of concordant and discordant pairs. 
Two patient pairs are concordant if the pair with the greater observed benefit also has higher predicted benefit. 
Two patient pairs are discordant if the pair with greater observed benefit has lower predicted benefit. Two patient pairs are uninformative if the pairs have the same observed benefit.

<b>Install using</b> 

library(remotes)

remotes::install_github("CHMMaas/HTEPredictionMetrics")

For an example of how to use the "C.for.Benefit" function, have look at "Examples" under "?C.for.Benefit".

<b>Reference</b>

van Klaveren D, Steyerberg EW, Serruys PW, Kent DM. The proposed 'concordance-statistic for benefit' provided a useful metric when modeling heterogeneous treatment effects. J Clin Epidemiol. 2018 Feb;94:59-68. doi: 10.1016/j.jclinepi.2017.10.021. Epub 2017 Nov 11. PMID: 29132832; PMCID: PMC7448760.
