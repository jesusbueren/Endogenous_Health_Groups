subroutine save_results(c_tr,p,posterior_beta,posterior_mu,init_cond,it)
    use global_var
    implicit none
    double precision,dimension(clusters**2*covariates,1),intent(in)::c_tr
    double precision,dimension(variables,clusters),intent(in)::p
    double precision,intent(in)::posterior_beta,posterior_mu
    double precision,dimension(clusters,L_e,L_gender),intent(in)::init_cond
    integer,intent(in)::it
    
    if (it==0) then
        open(unit=10,file=path_s_fin//'c_tr'//s_c//'_new.txt',status='replace')!path//
            write(10,'(F20.5)') c_tr
        close(10)
        open(unit=10,file=path_s_fin//'p'//s_c//'_new.txt',status='replace')!path//
            write(10,'(F20.5)') p
        close(10)
        open(unit=10,file=path_s_fin//'init_cond'//s_c//'_new.txt',status='replace')!path//
            write(10,'(F20.5)') init_cond
        close(10)
        open(unit=10,file=path_s_fin//'posterior'//s_c//'_new.txt',status='replace')!path//
            write(10,'(F14.2,F14.2)') posterior_beta,posterior_mu
        close(10)
    else
        open(unit=10,file=path_s_fin//'c_tr'//s_c//'_new.txt',access='append')!path//
            write(10,'(F20.5)') c_tr
        close(10)
        open(unit=10,file=path_s_fin//'p'//s_c//'_new.txt',access='append')!path//
            write(10,'(F20.5)') p
        close(10)
        open(unit=10,file=path_s_fin//'init_cond'//s_c//'_new.txt',access='append')!path//
            write(10,'(F20.5)') init_cond
        close(10)
        open(unit=10,file=path_s_fin//'posterior'//s_c//'_new.txt',access='append')!path//
            write(10,'(F14.2,F14.2)') posterior_beta,posterior_mu
        close(10)
    end if
        
end subroutine