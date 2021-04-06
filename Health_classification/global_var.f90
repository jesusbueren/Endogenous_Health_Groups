module global_var
    implicit none
    !Define number of clusters (same number in both variables)
    integer,parameter::clusters=4
    character(LEN=1)::s_c
    integer,parameter::variables=12, indv=26486, generations=20, initial_age=60,L_e=2,L_gender=2,covariates=6 !(indv=normal sample:26486, forecasting: 25063)
    integer::G=4000
    integer,parameter::K=6+(clusters-1)*2
    integer::sample_size
    integer,dimension(indv,variables,generations)::data
    integer,dimension(indv)::first_age,last_age,gender,educ,first_age2,last_age2
    character(LEN=34)::path="C:\Users\jbueren\Google Drive\ABC\"  
    character(LEN=71)::path_s_fin="C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\Results_ABC\"
    character(LEN=27)::path_s2="implied_cluster_pr\Results\"
    integer::sims
    double precision, parameter :: pi=3.14159265358979323846264338327950288419716939937510

end module