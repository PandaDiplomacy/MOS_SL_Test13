MOS_SL_Test13
=============

Python 2.7 code by Ben Richardson (2013) as part of a Part III Astrophysics project. 

Predicts the radial velocity accuracy of a Sandage-Loeb campaign using emission line galaxies.

Uses the numpy, scipy and cosmolopy (http://roban.github.io/CosmoloPy/) packages.

Primary Modules:

SL_Test.py: The main program. Options to predit RV errors, megagalaxy errors and Monte Carlo from a pickled
z,L,F list. 

  Sky_Background.py: Fitted sky background data from the E-ELT reference mission

LZ_Dist.py: Returns an array N(z,L) giving the number of galaxies in a series of z,L bins.

  LF_Data.py: Fitted Schechter profiles for H-alpha and [OII] emission line galaxies

Secondary Modules:

Useful_Functions.py: A series of small, utility functions. Details within.
Constants.py: Scientific constants (including cosmological density parameters)


