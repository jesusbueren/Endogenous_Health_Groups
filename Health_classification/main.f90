program main
    use global_var
    implicit none
    integer,dimension(1)::seed=789 
    double precision,dimension(variables,clusters)::p_star
    double precision,dimension(clusters**2*covariates,1)::c_tr_star
    double precision,dimension(covariates*clusters*clusters,covariates*clusters*clusters)::cov_tr
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender)::H
    double precision,dimension(clusters,generations,L_e,L_gender)::dist_init
    double precision,dimension(indv,clusters,generations)::likelihood
    double precision,dimension(indv,clusters,generations)::filtered_states
    double precision,dimension(indv,clusters+1,generations)::smoothed_states
    integer,dimension(indv,generations)::sample_k
    character::end_of_program
    integer::ns,burn,g_l,i_l,it
    double precision::counter,chg
    
    
    call random_seed(PUT=seed) 
    Write( s_c, '(I1)' )  clusters
    call charge_data()
    
    !Load high density point from estimated posterior distribution 
    ! and take it as true value
    call load_high_density(p_star,c_tr_star)
    
    call transitions(c_tr_star,H)
    !Initial guess for the smoothed states
    smoothed_states=1.0d0/dble(clusters)
    
      
    print*,'Compute filtered probability for the agent: assume that before seeing his health, he has the av health of his education and sex'
    
    call load_init_cond(H,smoothed_states)
    call likelihood_all(p_star,likelihood)
    call filtration(H,p_star,smoothed_states,likelihood,filtered_states) 
    
    open(unit=9,file=path//path_s2//'naive_pr_'//s_c//'.txt',status='replace')
    do i_l=1,indv; do g_l=1,generations;
        write(9,'(I5,I3,<clusters>F6.3,I3)'), i_l,g_l,filtered_states(i_l,:,g_l),maxloc(filtered_states(i_l,:,g_l))
    end do;end do
    close(9)

    read*,end_of_program
end
