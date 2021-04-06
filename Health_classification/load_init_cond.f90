subroutine load_init_cond(H,smoothed_states)
    use global_var
    implicit none
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender),intent(in)::H
    double precision,dimension(indv,clusters+1,generations),intent(out)::smoothed_states    
    integer,parameter,dimension(5)::iterations=(/0,30268,119492,82748,225796/)
    double precision,dimension(clusters*L_e*L_gender,iterations(clusters))::init_cond_v
    double precision,dimension(clusters,L_e,L_gender,iterations(clusters))::init_cond
    double precision,dimension(clusters,L_e,L_gender)::init_cond_star
    double precision,dimension(clusters,L_e,L_gender,generations)::init_cond_t
    integer::e_l,ge_l,t_l,i_l
    
    open(unit=9,file=path_s_fin//'init_cond'//s_c//'.txt')
            read(9,*) init_cond_v
    close(9)
    
    init_cond=reshape(init_cond_v,(/clusters,L_e,L_gender,iterations(clusters)/))
    init_cond_star=sum(init_cond(:,:,:,iterations(clusters)-10000:iterations(clusters)),4)/(10001)
    
    !Average health across time
    do e_l=1,2;do ge_l=1,L_gender; do t_l=1,generations
        if (t_l==1)then
            init_cond_t(:,e_l,ge_l,t_l)=init_cond_star(:,e_l,ge_l)
        else
            init_cond_t(:,e_l,ge_l,t_l)=matmul(transpose(H(1:clusters,1:clusters,t_l,e_l,ge_l)),init_cond_t(:,e_l,ge_l,t_l-1))
            init_cond_t(:,e_l,ge_l,t_l)=init_cond_t(:,e_l,ge_l,t_l)/sum(init_cond_t(:,e_l,ge_l,t_l))
        end if
    end do; end do;end do
    
    smoothed_states=0.0d0
    do i_l=1,indv;do t_l=1,generations
        smoothed_states(i_l,1:clusters,t_l)=init_cond_t(:,educ(i_l),gender(i_l),t_l)
    end do;end do
end subroutine