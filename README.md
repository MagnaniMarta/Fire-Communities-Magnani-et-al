# Fire-Communities-Magnani-et-al
This repository provides the code to model plant-fire relationship in the Mediterranean, tropical and boreal communities, as in Magnani et al. 

# Summay

This program models the dynamics of communities in fire prone environments. 
It includes 3 vegetation types that have differet flammability (L), fire response (R),
growth rate (c) and mortality rate (m). The model represents the time series of fractional
vegetation cover (b) of the 3 plant types. The deterministic succession of plant 
(Tilman,ecology 1994) is perturbed by stochastic fires. Fire is represented as 
a non-stationary Poisson process, with average fire return time, 'tf'. 
The average fire return time depends on plant cover and flammability, leading to
a fire-vegetation feedback.

# Content

2 folders:

TimeSeries/

Codes in this foder produce only one time series for set parameter values. Contains:

  mainfire.f90 : master code simulating the time series of vegetation cover and fires.
               Subrutines reported at the end of the file.

  parafire.f90 : module containing parameter values.
               Change the parameter values to represent different biomes.

  Three explaes are enclosed:
  - Mediterranean bassin 'parafireMed.f90'
  - tropical forest and savannas 'parafireTrop.f90'
  - boreal forest 'parafireBor.f90'

SensitivityAnalisis/

Codes in this folder explore a section of the parameter space identified by 2 parameters.
Change 'idummy' and 'MSEED' in run3() several times to obtain the parameter plane section of the paper. Contains:

  mainfire.f90 : master code iterating over 2 parameters and simulating the time series of vegetation cover and fires.
               Subrutines reported at the end of the file.

  parafire.f90 : module containing parameter values. Itaration parameters must not be fixed.
               Change the parameter values to represent different biomes as described above.


See the paper for plant type and parameter estimation.

# How to run

The code is provided in Fortran language.

Example of code compilation:
gfortran parafireMed.f90 mainfire.f90 -o runFire


# References and Contacts

This code has been developed within the manuscript:
'Fire responses shape plant communities in a minimal model for fire ecosystems worldwide'
(Magnani M.,Diaz-Sierra R., Sweneey L., Provenzale A., Baudena M.). TEMPORARY REFERENCE.

Please refer to the paper for the underlying assumptions and equations.

Authors: Marta Magnani and Mara Baudena

For more information please contact
Dr. Marta Magnani ( marta.magnani@igg.cnr.it )
