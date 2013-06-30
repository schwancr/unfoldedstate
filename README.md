SAXS
====

Goals / Questions to attempt to answer:
---------------------------------------

1. Does MD unfolded state fit the experimental saxs data?
2. Can we re-weight an MSM to fit the experimental data?
  a. If so, how does the distribution of structures change?
  b. What do the re-weighted structures look like? Are they extended? Or compact?
  c. Are there many distributions that can fit the experiment? If so how do these differ?
  d. If a compact distribution and an extended distribution can both fit the experiment well, this is interesting.
  e. But it is also interesting if only one fits.
3. Can we re-weight a kinetic model to fit the experiment?
4. How do explicit vs. implicit predicted saxs differ (for NuG2)? 
  a. Problem here is we really only have native state experimental data, unless we want denatured experiments
  b. We could instead do protein L or ubiquitin since Diwakar has explicit simulations 
  c. Or we should just use DESRES ubiquitin



Things we need to bug James about:
----------------------------------
1. What is the range of q? Can we have more data?
2. We need all of the data that made the error bars (Because we need to calculate covariance matrix for the likelihood)
3. Ubiquitin / Protein L data as well (since we have explicit and implicit)
