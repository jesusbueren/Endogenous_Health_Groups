## Endogenous_Health_Groups
In this repository, I explain the codes for replicating Endogenous Health Groups and Heterogeneous Dynamics of the Elderly, by Amengual, Bueren, and Crego.

The replication material is divided in three blocks:
1. **Data Preparation**: from the HRS to the files to be read by the estimation program.
1. **Estimation Program**: reads the data from the previous step and estimates the econometric model parameters.
1. **Health Classification**: from the econometric model parameter, it classifies individuals into health groups.

### Data Preparation

This folder contains STATA do-file. This do file reads the data from the HRS rand contributed files, cleans the data, and produces a series of csv files:
1. data_all.csv has individual level information on I-ADLs for all the interviews.
1. ages_all.csv has individual level information on age at first and last interview.
2. gender_all.csv has individual level information on gender.
3. educ_all.csv has individual level information on education.
