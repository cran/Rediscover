# version 0.3.1:

remove dependency of 'TCGAbiolinks'.

# version 0.3.0:  

getMutex and getMutexAB computes the exact p-value using the ShiftConvolvePoibin R package (Peres, N., Lee, A., and Keich, U.,2020). Mixed and exact have been enhanced in terms of speed and memory usage.

# version 0.2.0:

getMutex and getMutexAB have three different approximations of the Poison-Binomial distribution function: a "RefinedNormal"" Approximation (Volkova, 1996), "Binomial"" with two parameters (Cam, 1960), and a "Shifted Binomial" with three parameters (Peköz, Shwartz, Christiansen, & Berlowitz, 2010). This last is the approximation by default for both functions. Moreover, some bugs have been fixed.
