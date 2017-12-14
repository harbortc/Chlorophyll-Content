# Chlorophyll-Content
This script calculates Chlorophyll content (Total, A, B, and ratio) in mg/gFW from absorbance values of a DMSO extraction, plots the results in box plots by Genotype and performs ANOVA with Tukey's post-hoc test.

Input file is a tab-delim with . as decimal separator
Columns: "Genotype" "FW" "A652" "A663" "A645"
No dashes allowed in Genotype, though line included to change Col-0 to Col0 in case

Makes plots of input FW (used for chlorophyll extraction), Total, A, B, and ratio A/B with stats.
FW is not Shoot Fresh Weight as growth measurement (Different script)
