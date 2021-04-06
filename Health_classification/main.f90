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
    
    !print*,'Compute smoothed probability/full info for econometritian'
    !call likelihood_all(p_star,likelihood)
    !do it=1,100 !I run it 100 so that initial conditions for filter don't matter
    !    print*,it
    !    call filtration(H,p_star,smoothed_states,likelihood,filtered_states)
    !    call smoothing(filtered_states,H,p_star,sample_k,smoothed_states)  
    !end do
    !
    !open(unit=9,file=path//path_s2//'smoothed_pr_'//s_c//'.txt',status='replace')
    !do i_l=1,indv; do g_l=1,generations;
    !    if (clusters==2) then
    !       ! write(9,'(I5,I3,F6.3,F6.3,I3)'), i_l,g_l,smoothed_states(i_l,1,g_l),smoothed_states(i_l,2,g_l),maxloc(smoothed_states(i_l,:,g_l))
    !    elseif (clusters==3) then
    !        !write(9,'(I5,I3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,smoothed_states(i_l,1,g_l),smoothed_states(i_l,2,g_l),smoothed_states(i_l,3,g_l),maxloc(smoothed_states(i_l,:,g_l))
    !    elseif (clusters==4) then
    !        write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,smoothed_states(i_l,1,g_l),smoothed_states(i_l,2,g_l),smoothed_states(i_l,3,g_l),smoothed_states(i_l,4,g_l),maxloc(smoothed_states(i_l,:,g_l))
    !    elseif (clusters==5) then
    !        !write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,smoothed_states(i_l,1,g_l),smoothed_states(i_l,2,g_l),smoothed_states(i_l,3,g_l),smoothed_states(i_l,4,g_l),smoothed_states(i_l,5,g_l),maxloc(smoothed_states(i_l,:,g_l))
    !    end if
    !end do;end do
    !close(9)
    !
    !!Compute average absolute deviation between filter and smoother
    !counter=0.0d0
    !chg=0.0d0
    !do i_l=1,indv; do g_l=1,generations;
    !    if (sample_k(i_l,g_l)/=-1 .and. sample_k(i_l,g_l)/=clusters+1) then
    !        if (maxloc(smoothed_states(i_l,:,g_l),1)/=maxloc(filtered_states(i_l,:,g_l),1)) then
    !            chg=chg+1
    !        end if
    !        counter=counter+1.0d0
    !    end if
    !end do;end do
    !
    !print*,'Percent of people that change status from filter to smoother'
    !print*,chg/counter*100
    
    print*,'Compute filtered probability for the agent: assume that before seeing his health, he has the av health of his education and sex'
    
    call load_init_cond(H,smoothed_states)
    call likelihood_all(p_star,likelihood)
    call filtration(H,p_star,smoothed_states,likelihood,filtered_states) 
    
    open(unit=9,file=path//path_s2//'naive_pr_'//s_c//'.txt',status='replace')
    do i_l=1,indv; do g_l=1,generations;
        if (clusters==2) then
            !write(9,'(I5,I3,F6.3,F6.3,I3)'), i_l,g_l,filtered_states(i_l,1,g_l),filtered_states(i_l,2,g_l),maxloc(filtered_states(i_l,:,g_l))
        elseif (clusters==3) then
            !write(9,'(I5,I3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,filtered_states(i_l,1,g_l),filtered_states(i_l,2,g_l),filtered_states(i_l,3,g_l),maxloc(filtered_states(i_l,:,g_l))
        elseif (clusters==4) then
            write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,filtered_states(i_l,1,g_l),filtered_states(i_l,2,g_l),filtered_states(i_l,3,g_l),filtered_states(i_l,4,g_l),maxloc(filtered_states(i_l,:,g_l))
        elseif (clusters==5) then
            !write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,filtered_states(i_l,1,g_l),filtered_states(i_l,2,g_l),filtered_states(i_l,3,g_l),filtered_states(i_l,4,g_l),filtered_states(i_l,5,g_l),maxloc(filtered_states(i_l,:,g_l))
        end if
    end do;end do
    close(9)

    read*,end_of_program
end