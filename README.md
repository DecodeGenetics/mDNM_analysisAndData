# mDNM_analysisAndData
Scripts for mDNM analysis and the data the analysis was run on.

attributeRegression.R
Poisson regression testing effect of reference length, purity, GC motif content and motif length on mDNM rate.
Also splits up by motif length and tests effect of remaining attributes on mDNM rate.
Last, splits by motif length and reference length to test effect of purity
Input file: mutRateDataAll

mutRateCIs.R
Get CI for mDNM rate and CIs for motif length specific mDNM rates.
Input file: mutRateDataAll

em.R
Test effect of maternal and paternal age on number of transmitted mDNMs and interpolate to genome wide age effect.
Input file: totalPerPn
