## Endogenous_Health_Groups
In this repository, I explain the codes for replicating Endogenous Health Groups and Heterogeneous Dynamics of the Elderly, by Amengual, Bueren, and Crego.

The replication material is divided in three blocks:
1. **Data Preparation**: from the HRS to the files to be read by the estimation program.
1. **Estimation Program**: reads the data from the previous step and estimates the econometric model parameters.
1. **Health Classification**: from the econometric model parameter, it classifies individuals into health groups.

### Data Preparation

This folder contains a STATA do-file. This do file reads the data from the HRS rand contributed files, cleans the data, and produces a series of csv files:
1. data_all.csv has individual level information on I-ADLs for all the interviews.
1. ages_all.csv has individual level information on age at first and last interview.
2. gender_all.csv has individual level information on gender.
3. educ_all.csv has individual level information on education.

Before you run it, you need to change the directory in line 2 of the do file

### Estimation Program

This folder contains a set Fortran 90 files.

global_var.f90: includes all the set of global variables. In this file, you need to change the location and length of the string of the location of the different files:

The concatenation of the strings path and path_s_ini is where the set of initial conditions will be stored.

The concatenation of the strings path and path_s_fin is where the set of the posterior distribution of the estimated parameters will be stored.

In you main path you need to create a folder named "Data" where you include the csv files from Data Preparation.

main.f90 is the main script of the code.

It first calls charge_data.f90 which loads the csv files from Data Preparation. Then it call for the set of initial conditions using initial_conditions.f90 and finally runs the main estimation exercise full_posterior.f90

