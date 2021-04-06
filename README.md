## Endogenous_Health_Groups
In this repository, I explain the codes for replicating Endogenous Health Groups and Heterogeneous Dynamics of the Elderly, by Amengual, Bueren, and Crego.

The replication material is divided in three blocks:
1. **Data Preparation**: from the HRS to the files to be read by the estimation program.
1. **Estimation Program**: reads the data from the previous step and estimates the econometric model parameters.
1. **Health Classification**: from the econometric model parameter, it classifies individuals into health groups.

### 1. Data Preparation

This folder contains a STATA do-file. This do file reads the data from the HRS rand contributed files, cleans the data, and produces a series of csv files:
1. data_all.csv has individual level information on I-ADLs for all the interviews.
1. ages_all.csv has individual level information on age at first and last interview.
2. gender_all.csv has individual level information on gender.
3. educ_all.csv has individual level information on education.

Before you run it, you need to change the directory in line 2 of the do file

### 2. Estimation Program

This folder contains a set Fortran 90 files.

global_var.f90: includes all the set of global variables. In this file, you need to set the number of clusters (clusters=x) and change the location and length of the string of the location of the different files:

The concatenation of the strings path and path_s_ini is where the set of initial conditions will be stored.

The concatenation of the strings path and path_s_fin is where the set of the posterior distribution of the estimated parameters will be stored.

In you main path you need to create a folder named "Data" where you include the csv files from Data Preparation.

main.f90 is the main script of the code. It first calls charge_data.f90 which loads the csv files from Data Preparation. Then it call for the set of initial conditions using initial_conditions.f90 and finally runs the main estimation exercise full_posterior.f90

#### 2.1 Initial Conditions

Initial conditions are obtained in two blocks:
1. The first block estimates the initial conditions for the probability of I-ADLs in each group.
2. The second block estimatates the initial conditions for the parameters drinving the transition probabilities.

##### Initial Conditions for the Probability of I-ADLs in Each Group

This first step is done by estimating a mixture model by pooling all individuals with available information on I-ADLs. This initial model is estimated using an EM algorithm ignoring all the time series information from transitions. This step produces the probability that each interviewed individual belongs to each health group. Given these probabilities, we sample an individual health group at random.  A detailed exposition on how to estimate this class of models can be found in section 9.3.3 of "Pattern Recognition and Machine Learning" by Christopher M. Bishop.  

##### Initial Conditions for Transition Probabilities

The second step of the initial conditions is estimated taking as observed the previously assigned health groups. We thus perform a Bayesian estimation of a multinomial logit model using a Metropolis algorithm. In order to speed up the mixing in the proposal we make use of the adaptive metropolis algorithm proposed by Haario et al. (2001). We save the mean of the posterior distribution and the variance covariance matrix of the proposal. 

