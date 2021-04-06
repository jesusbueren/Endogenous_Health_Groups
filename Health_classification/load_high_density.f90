subroutine load_high_density(p,c_mean)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(out)::p
    integer,parameter,dimension(5)::burn=(/0,1000,60000,60000,130000/)
    integer,parameter,dimension(5)::iterations=(/0,30268,119492,82748,225796/)
    double precision,dimension(variables,clusters,iterations(clusters))::p_all
    double precision,dimension(clusters**2*covariates,iterations(clusters))::c_vec
    double precision,dimension(clusters**2*covariates,1),intent(out)::c_mean
    integer::c_l,c_l2,v_l,i,j
    
    open(unit=9,file=path_s_fin//'p'//s_c//'.txt')
        read(9,'(F20.5)') p_all
    close(9)
    open(unit=9,file=path_s_fin//'c_tr'//s_c//'.txt')
        read(9,'(F20.5)') c_vec
    close(9)
    
    !Compute high density point as mean of posterior density
    do v_l=1,variables; do c_l=1,clusters
        p(v_l,c_l)=sum(p_all(v_l,c_l,burn(clusters):iterations(clusters)))/dble(iterations(clusters)-burn(clusters)+1)
    end do; end do
    
    do i=1,covariates*clusters**2
        c_mean(i,1)=sum(c_vec(i,burn(clusters):iterations(clusters)))/dble(iterations(clusters)-burn(clusters)+1)
    end do
        
end subroutine