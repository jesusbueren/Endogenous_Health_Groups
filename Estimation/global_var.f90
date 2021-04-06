module global_var
    implicit none
    !Define number of clusters
    integer,parameter::clusters=4
    character(LEN=1)::s_c
    integer,parameter::variables=12, indv=26486, generations=20, initial_age=60,L_e=2,L_gender=2,covariates=6 !(indv=normal sample:26486, forecasting: 25063)
    integer,dimension(indv,variables,generations)::data
    integer,dimension(indv)::first_age,last_age,gender,educ
    character(LEN=32)::path="C:\Users\jesus\Google Drive\ABC\"  
    character(LEN=39)::path_s_ini="estimation\Results\initital_conditions\"
    !character(LEN=33)::path_s_fin="estimation\Results\final_results\"
    character(LEN=69)::path_s_fin="C:\Users\jesus\OneDrive - Istituto Universitario Europeo\Results_ABC\"
end module