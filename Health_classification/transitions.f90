   subroutine transitions(c_tr,H)
    use global_var
    implicit none
    double precision,dimension(covariates*clusters**2,1),intent(in)::c_tr
    integer::c_l,c_l2,g_l,ge_l,it,e_l
    double precision,dimension(clusters,clusters+1)::betas,alphas,c_e,c_gender,c_axg,c_axe
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender),intent(out)::H
    double precision::age,e_d,ge_d
    double precision,dimension(clusters,clusters,generations,L_e,L_gender)::H_new
    character::continue_program
    
    call c2beta(c_tr,betas,alphas,c_e,c_gender,c_axg,c_axe)
    do e_l=1,L_e
        e_d=dble(e_l-1)
        do c_l=1,clusters ; do c_l2=1,clusters+1; do g_l=1,generations; do ge_l=1,2
            age=dble(initial_age+(g_l-1)*2)
            ge_d=dble(ge_l-1)
            H(c_l,c_l2,g_l,e_l,ge_l)=1 &
                                /sum(exp(-betas(c_l,:)-alphas(c_l,:)*age-c_e(c_l,:)*e_d-c_gender(c_l,:)*ge_d &
                                         -c_axg(c_l,:)*age*ge_d-c_axe(c_l,:)*age*e_d- &
                                    (-betas(c_l,c_l2)-alphas(c_l,c_l2)*age-c_e(c_l,c_l2)*e_d-c_gender(c_l,c_l2)*ge_d &
                                         -c_axg(c_l,c_l2)*age*ge_d-c_axe(c_l,c_l2)*age*e_d)))
        end do; end do; end do; end do
    end do
    H(clusters+1,1:clusters,:,:,:)=0.0d0
    H(clusters+1,clusters+1,:,:,:)=1.0d0
        
    do c_l=1,clusters ; do c_l2=1,clusters; do g_l=1,generations;  do e_l=1,L_e; do ge_l=1,2
            H_new(c_l,c_l2,g_l,e_l,ge_l)=H(c_l,c_l2,g_l,e_l,ge_l)/sum(H(c_l,1:clusters,g_l,e_l,ge_l))
    end do; end do; end do; end do; end do

end subroutine
