# CBCTrendAnalysis
This is the official repository for the standard analysis of Audubon Christmas Bird Count data. The goal of the analysis is to produce annual abundance indices and temporal trends in those indices for more than 500 bird species that winter in Canada and the continental USA for use in studies of avian conservation and ecology.

Abundance indices and trends are estimated for spatial units (analytical strata) that are intersections of provinces, states, and Bird Conservation Regions (BCRs) using the spatial first-difference model described in Smith et al. (2024). Abundance indices are then scaled to larger areas and trends are estimated per province, state, BCR, country, and for the entire survey area. 

Note that abundance estimates and trends produced using this code will be different from those produced using previous methods described by Soykan et al. (2016). This will be especially true at the edges of time windows and geographic ranges, where we expect the current methods to produce more accurate and precise results (Smith et al. 2024).

This code base evolved from ideas and code developed and shared by Bill Link, John Sauer, and Dan Niven (Link et al. 2006); Candan Soykan (Soykan et al. 2016); Adam Smith and Brandon Edwards (Smith and Edwards 2021); and others. Special thanks to Adam Smith for contributing much of the Stan modeling and post processing code.

## Citations

Link, W. A., Sauer, J. R. and Niven, D. K. 2006. A hierarchical model for regional analysis of population change using Christmas Bird Count data, with application to the American Black Duck. The Condor, 108(1), pp.13-24.

Smith, A. C. and Edwards, B. P. 2021. North American Breeding Bird Survey status and trend estimates to inform a wide range of conservation needs, using a flexible Bayesian hierarchical generalized additive model. The Condor, 123(1), p.duaa065.

Smith, A. C., D. Binley, A., Daly, L., Edwards, B. P., Ethier, D., Frei, B., Iles, D., Meehan, T. D., Michel, N. L. and Smith, P. A. 2024. Spatially explicit Bayesian hierarchical models improve estimates of avian population status and trends. Ornithological Applications, 126(1), p.duad056.

Soykan, C. U., Sauer, J., Schuetz, J. G., LeBaron, G. S., Dale, K. and Langham, G. M. 2016. Population trends for North American winter birds based on hierarchical models. Ecosphere, 7(5), p.e01351.
