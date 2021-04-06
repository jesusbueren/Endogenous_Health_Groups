program main
    use global_var
    implicit none
    integer,dimension(1)::seed=254  
    double precision,dimension(variables,clusters)::p
    double precision,dimension(covariates*clusters**2,1)::c_tr
    double precision,dimension(covariates*clusters*clusters,covariates*clusters*clusters)::cov_tr
    character::end_of_program
    
    call random_seed(PUT=seed)
    
    Write( s_c, '(I1)' )  clusters 
    
    !Load original data
    call charge_data()
    
    !Initial conditions
    !call initial_conditions(p,c_tr,cov_tr)
    
    call load_initial_conditions(p,c_tr,cov_tr)
    
    !Compute full posterior
    call full_posterior(p,c_tr,cov_tr)    
    
    PRINT*, 'End of program, press any key to continue'                                                             
    READ*, end_of_program 
    
end program