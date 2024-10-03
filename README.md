# **Unsupervised linear discrimination using skewness**

This project contains code to reproduce the results of the paper with the same name by Radojicic, U., Nordhausen, K. and Virta, J.

The goal of the project was to investigate properties of skewness based estimators of the linear discriminant in an unsupervised framework in a setting for two component Gaussian location mixture.

**Main functions**

All code here is written in R and requires the packages MASS, ICtest, parallel and mvtnorm.

The main functions are:

-   ThetaM: the moment based estimator.

-   ThetaR: the affine equivariant moment estimator (Loperfido 2015).

-   ThetaL: third moment based estimator introduced by Loperfido (2013).

-   ThetaJ: the 3-JADE estimator.

-   ThetaP: the estimator based on projection pursuit using skewness (Virta et al. 2015). Wrapper around the function NGPP from the ICtest package.

-   ThetaLDA: the supervised benchmark. Wrapper around the lda function from the MASS package for linear discriminant analysis.

In addition the project contains several helper functions for data generation, computing the constants C for the affine equivariant estimators and for the simulations performed in the paper.

**Authors**

Radojicic, U., Nordhausen, K. and Virta, J.

**License**

GNU GPLv3

**References**

Loperfido, N. (2013): Skewness and the linear discriminant function. Statistics & Probability Letters, 83, 93-99.\
Loperfido, N. (2015): Vector-valued skewness for model-based clustering. Statistics & Probability Letters, 99, 230-237.\
Virta, J., Nordhausen, K. and Oja, H. (2015): Joint use of third and fourth cumulants in independent component analysis, arXiv:1505.02613.
