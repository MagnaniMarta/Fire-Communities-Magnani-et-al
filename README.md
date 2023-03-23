# Fire-Communities (Magnani et al. 2023)
This repository provides the code to model plant-fire relationship in the Mediterranean, tropical and boreal communities, as in Magnani et al. 2023.

DOI: 
10.5281/zenodo.7763275

Cite as: 
MagnaniMarta. (2023). MagnaniMarta/Fire-Communities-Magnani-et-al: Fire-Communities-Magnani-et-al-2023 (v1.0). Zenodo. 
https://doi.org/10.5281/zenodo.7763275

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

3 folders:

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
               
 CalibrationBoreal/
 
 Codes in this folder can be used for model calibration.
 
   main_nofire.f90 : master code simulating the time series of vegetation cover of one species without fires. Subroutines at the end of the file.

   parafireBor.f90 : module containing parameter values. Change 'c1' to find the value for which the cover attains its asymptotic values
                     within a time comparable to the value reported in 'boreal-species.csv' in 'time of max cover after burning'.    

   boreal-species.csv: values used for calibration in Magnani et al., 2023 (see the Online Supplement B).


See the paper for plant type and parameter estimation.

# How to run

The code is provided in Fortran language. Fortran compilers can for instance be downloaded from the webpage: 
https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#hpc-kit
For an introductory guide to fortran language see Perrin, C. L. (1997). Numerical Recipes in Fortran 90: The Art of Scientific Computing, Volume 2. 
By William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. Cambridge University Press: New York, 1996.  

Example of code compilation:
ifort parafireMed.f90 mainfire.f90 -o runFire


# References and Contacts

This code has been developed within the manuscript:
Magnani M., Diaz-Sierra R., Sweneey L., Provenzale A., Baudena M.(2023)'Fire responses shape plant communities in a minimal model for fire ecosystems across the world'. The American Naturalist.

Please refer to the paper for the underlying assumptions and equations.

Authors: Marta Magnani and Mara Baudena

For more information please contact
Dr. Marta Magnani ( marta.magnani@igg.cnr.it )
