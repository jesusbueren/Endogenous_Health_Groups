subroutine posterior_fct(H,p,sample_k,posterior_mu,posterior_beta,prior)
    use global_var
    implicit none
    double precision,intent(out)::posterior_mu,posterior_beta
    integer,dimension(indv,generations),intent(in)::sample_k
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender),intent(in)::H
    double precision,dimension(variables,clusters),intent(in)::p
    double precision,dimension(indv)::posterior_mu_i,posterior_beta_i
    integer::i_l,g_l
    double precision::likelihood,prior
    character::end_of_program
    
    posterior_mu_i=0.0d0
    posterior_beta_i=0.0d0
    !$OMP PARALLEL num_threads(3) DEFAULT(SHARED) PRIVATE(i_l,g_l,likelihood)
    !$OMP DO SCHEDULE(GUIDED)
    do i_l=1,indv; do g_l=first_age(i_l),last_age(i_l)
        if (data(i_l,1,g_l)/=-1) then
            call likelihood_i(i_l,sample_k(i_l,g_l),g_l,p,likelihood)
            posterior_mu_i(i_l)=posterior_mu_i(i_l)+log(likelihood)
        end if
        if (sample_k(i_l,g_l)>-1 .and. g_l>first_age(i_l)) then
            posterior_beta_i(i_l)=posterior_beta_i(i_l)+log(H(sample_k(i_l,g_l-1),sample_k(i_l,g_l),g_l-1,educ(i_l),gender(i_l)))
        end if
        if (isnan(posterior_beta_i(i_l)))then
            print*,'Error in posterior'
            print*,'middle_age',H(sample_k(i_l,g_l-1),sample_k(i_l,g_l),g_l-1,educ(i_l),gender(i_l))          
            READ*, end_of_program 
        end if
    end do; end do
    !$OMP END DO
    !$OMP END PARALLEL
    posterior_mu=sum(posterior_mu_i)
    posterior_beta=sum(posterior_beta_i)+prior
end subroutine