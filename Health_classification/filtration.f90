subroutine filtration(H,p,smoothed_states,L_i,filtered_states)
    use global_var
    implicit none
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender),intent(in)::H
    double precision,dimension(variables,clusters),intent(in)::p
    double precision,dimension(indv,clusters+1,generations),intent(in)::smoothed_states
    double precision,dimension(indv,clusters,generations),intent(in)::L_i
    double precision,dimension(indv,clusters,generations),intent(out)::filtered_states
    double precision,dimension(clusters,clusters)::joint_states
    double precision::A
    integer::g_l,i_l,c_l1,c_l2
    
    filtered_states=0.0d0
    !$OMP PARALLEL num_threads(4)  DEFAULT(SHARED) PRIVATE(g_l,i_l,c_l1,c_l2,A,joint_states)
    !$OMP DO 
    do i_l=1,indv
        do g_l=first_age(i_l),last_age(i_l)
            if (data(i_l,1,g_l)/=-1) then
                do c_l1=1,clusters; do c_l2=1,clusters
                    if (g_l==first_age(i_l)) then
                        joint_states(c_l2,c_l1)=smoothed_states(i_l,c_l2,g_l)*L_i(i_l,c_l2,g_l)
                    elseif (g_l>first_age(i_l) .and. g_l<=last_age(i_l)) then
                        joint_states(c_l2,c_l1)=filtered_states(i_l,c_l1,g_l-1)*H(c_l1,c_l2,g_l-1,educ(i_l),gender(i_l))*L_i(i_l,c_l2,g_l)
                    end if
                end do;end do
                A=sum(joint_states)
                do c_l2=1,clusters
                    filtered_states(i_l,c_l2,g_l)=sum(joint_states(c_l2,:))/A
                end do
            end if
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

end subroutine