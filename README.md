# Chlorophyll-Content
This script calculates Chlorophyll content (Total, A, B, and ratio) in mg/gFW from absorbance values of a DMSO extraction, plots the results in box plots by Genotype and performs ANOVA with Tukey's post-hoc test.

Input file is a tab-delim with . as decimal separator
Columns: "Genotype" "FW" "A652" "A663" "A645"
No dashes allowed in Genotype, though line included to change Col-0 to Col0 in case

Makes plots of input FW (used for chlorophyll extraction), Total, A, B, and ratio A/B with stats.
FW is not Shoot Fresh Weight as growth measurement (Different script)
Performs anova and labels groups based on 95% confidence, then saves plots in the working directory

If extraction was done in a volume of DMSO besides 1ml, can adjust by defining volume in the functions, if left blank default it 1ml

Used internet example for stats on my data without understanding it
https://stackoverflow.com/questions/18771516/is-there-a-function-to-add-aov-post-hoc-testing-results-to-ggplot2-boxplot
This script gives me what I want, though I don't quite understand how yet. Works for now...

Doesn't yet take different media types to make a grouped boxplot
