subroutine load_initial_conditions(p,c_tr,cov_tr)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(out)::p
    double precision,dimension(covariates*clusters**2,1)::c_tr
    double precision,dimension(covariates*clusters**2,covariates*clusters**2),intent(out)::cov_tr
    
    open(unit=9,file=path//path_s_ini//'p_ini'//s_c//'.txt')
        read(9,'(F14.10)') p
    close(9)
    open(unit=9,file=path//path_s_ini//'c_tr_ini'//s_c//'.txt')
        read(9,'(F14.10)') c_tr
    close(9)
    open(unit=9,file=path//path_s_ini//'cov_tr_ini'//s_c//'.txt')
        read(9,'(F14.10)') cov_tr
    close(9)
    
end subroutine