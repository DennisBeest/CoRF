# CoRF
This repository contains code for fitting a Co-data moderated randomForest (CoRF). 
Link to paper:
https://link.springer.com/article/10.1186/s12859-017-1993-1

For more information on what is "co-data" see:
https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6732


Install with:

library(devtools)
install_github("DennisBeest/CoRF")

See help files for intructions and examples (?CoRF)

See the help file on CoRF for the example of the paper. There are instructions on how to use the package for fitting, but there are also instructions on how to directly use underlying procedures (scam and rfsrc). The latter option gives more freedom in what scam model you want to fit. The CorRF is a wrapper function.

