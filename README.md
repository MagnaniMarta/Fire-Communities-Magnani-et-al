# Fire-Communities (Magnani et al. 2023)
This repository provides the code to model plant-fire relationship in the Mediterranean, tropical and boreal communities, as in Magnani et al. 

# Summary

This program models the dynamics of communities in fire prone environments. 
It includes 3 vegetation types that have different flammability (L), fire response (R),
growth rate (c) and mortality rate (m). The model represents the time series of fractional
vegetation cover (b) of the 3 plant types. The deterministic succession of plant 
(Tilman - Ecology, 1994) is perturbed by stochastic fires. Fire is represented as 
a non-stationary Poisson process, with average fire return time, 'tf'. 
The average fire return time depends on plant cover and flammability, leading to
a fire-vegetation feedback.

# Content

2 folders:

TimeSeries/

Codes in this folder produce only one time series for set parameter values (e.g. Fig. 2 in Magnani et al., 2023). Contains:

  mainfire.f90 : master code simulating the time series of vegetation cover and fires.
               Subroutines reported at the end of the file.

  parafire.f90 : module containing parameter values.
               Change the parameter values to represent different biomes.

  Three examples are enclosed:
  - Mediterranean basin 'parafireMed.f90'
  - tropical forest and savannas 'parafireTrop.f90'
  - boreal forest 'parafireBor.f90'

SensitivityAnalysis/

Codes in this folder explore a section of the parameter space identified by 2 parameters (e.g. Fig. 4 in Magnani et al., 2023).
Change 'idummy' and 'MSEED' in run3() several times to obtain the parameter plane section of the paper. Contains:

  mainfire.f90 : master code iterating over 2 parameters and simulating the time series of vegetation cover and fires.
               Subroutines reported at the end of the file.

  parafire.f90 : module containing parameter values. Iteration parameters must not be fixed.
               Change the parameter values to represent different biomes as described above.


See the paper for plant type and parameter estimation.

# How to run

The code is provided in Fortran language.

Example of code compilation:
gfortran parafireMed.f90 mainfire.f90 -o runFire


# References and Contacts

This code has been developed within the manuscript:
Magnani M., Diaz-Sierra R., Sweneey L., Provenzale A., Baudena M.(2023)'Fire responses shape plant communities in a minimal model for fire ecosystems across the world'. The American Naturalist.

Please refer to the paper for the underlying assumptions and equations.

Authors: Marta Magnani and Mara Baudena

For more information please contact
Dr. Marta Magnani ( marta.magnani@igg.cnr.it )
