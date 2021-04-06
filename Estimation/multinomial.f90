subroutine multinomial(c_tr,posterior_beta)
    use global_ini;use global_var
    implicit none
    double precision,dimension(covariates*clusters**2,1),intent(in)::c_tr
    integer::c_l,c_l2,g_l,e_d,ge_l,age,ge_d,it,i_l,e_l
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender)::H
    double precision,dimension(clusters,clusters+1)::betas,alphas,c_e,c_gender,c_axg,c_axe
    double precision,dimension(clusters,generations,L_e,L_gender)::dist_init
    double precision,dimension(clusters,clusters,generations,L_e,L_gender)::H_new
    double precision,dimension(indv)::posterior_beta_i
    double precision,intent(out)::posterior_beta
    
    call c2beta(c_tr,betas,alphas,c_e,c_gender,c_axg,c_axe)
    
    do e_l=1,L_e
        e_d=e_l-1
        do c_l=1,clusters ; do c_l2=1,clusters+1; do g_l=1,generations; do ge_l=1,2
            age=initial_age+(g_l-1)*2
            ge_d=ge_l-1
            H(c_l,c_l2,g_l,e_l,ge_l)=1 &
                                /sum(exp(-betas(c_l,:)-alphas(c_l,:)*age-c_e(c_l,:)*e_d-c_gender(c_l,:)*ge_d &
                                         -c_axg(c_l,:)*age*ge_d-c_axe(c_l,:)*age*e_d- &
                                    (-betas(c_l,c_l2)-alphas(c_l,c_l2)*age-c_e(c_l,c_l2)*e_d-c_gender(c_l,c_l2)*ge_d &
                                         -c_axg(c_l,c_l2)*age*ge_d-c_axe(c_l,c_l2)*age*e_d)))
        end do;end do; end do; end do
    end do
    
    !No resurection
    H(clusters+1,1:clusters,:,:,:)=0 
    H(clusters+1,clusters+1,:,:,:)=1 
    
    do c_l=1,clusters ; do c_l2=1,clusters; do g_l=1,generations;  do e_l=1,L_e; do ge_l=1,2
            H_new(c_l,c_l2,g_l,e_l,ge_l)=H(c_l,c_l2,g_l,e_l,ge_l)/sum(H(c_l,1:clusters,g_l,e_l,ge_l))
    end do; end do; end do; end do; end do  
    
    dist_init(:,:,:,:)=1/dble(clusters)
    do g_l=1,generations; do e_l=1,L_e; do ge_l=1,2
        if (g_l==1) then
            do it=1,40
                dist_init(:,g_l,e_l,ge_l)=matmul(dist_init(:,g_l,e_l,ge_l),H_new(:,:,g_l,e_l,ge_l))     
            end do
        else
            dist_init(:,g_l,e_l,ge_l)=matmul(dist_init(:,g_l-1,e_l,ge_l),H_new(:,:,g_l,e_l,ge_l))
        end if 
    end do; end do; end do
    
    !$OMP PARALLEL num_threads(3) DEFAULT(SHARED) PRIVATE(i_l,g_l)
    !$OMP DO 
    do i_l=1,indv; do g_l=first_age(i_l),last_age(i_l)
        if (g_l==first_age(i_l)) then
          posterior_beta_i(i_l)=log(dist_init(sample_k_ini(i_l,g_l),g_l,educ(i_l),gender(i_l)))
        elseif (sample_k_ini(i_l,g_l)>-1 .and. g_l>first_age(i_l)) then
            posterior_beta_i(i_l)=posterior_beta_i(i_l)+log(H(sample_k_ini(i_l,g_l-1),sample_k_ini(i_l,g_l),g_l-1,educ(i_l),gender(i_l)))
        end if
    end do; end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    posterior_beta=sum(posterior_beta_i)
    
    end subroutine
    