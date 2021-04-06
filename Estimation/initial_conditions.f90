module global_ini
    use global_var
    implicit none
    integer,dimension(indv,generations)::sample_k_ini
end module
    
subroutine initial_conditions(p,c_tr,cov_tr)
    use global_var;use global_ini
    implicit none
    double precision,dimension(variables,clusters),intent(out)::p
    double precision,dimension(covariates*clusters**2,1),intent(out)::c_tr
    double precision,dimension(covariates*clusters*clusters,covariates*clusters*clusters),intent(out)::cov_tr
    double precision,dimension(indv,generations,clusters)::gamma_ik
    integer::i_l,g_l
    character::pause_k
    !Initial conditions for p (MLE)
    call initial_conditions_p(p,gamma_ik)
    
    print*,p
    !Save results
    open(unit=9,file=path//path_s_ini//'p_ini'//s_c//'.txt')
        write(9,'(F14.10)') p
    close(9)
    
    open(unit=9,file=path//path_s_ini//'no_dyn_pr_'//s_c//'.txt',status='replace')
    do i_l=1,indv; do g_l=1,generations;
        if (clusters==2) then
           ! write(9,'(I5,I3,F6.3,F6.3,I3)'), i_l,g_l,gamma_ik(i_l,g_l,1),gamma_ik(i_l,g_l,2),maxloc(gamma_ik(i_l,g_l,:))
        elseif (clusters==3) then
            !write(9,'(I5,I3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,gamma_ik(i_l,g_l,1),gamma_ik(i_l,g_l,2),gamma_ik(i_l,g_l,3),maxloc(gamma_ik(i_l,g_l,:))
        elseif (clusters==4) then
            !write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,I3)'), i_l,g_l,gamma_ik(i_l,g_l,1),gamma_ik(i_l,g_l,2),gamma_ik(i_l,g_l,3),gamma_ik(i_l,g_l,4),maxloc(gamma_ik(i_l,g_l,:))
        elseif (clusters==5) then
            write(9,'(I5,I3,F6.3,F6.3,F6.3,F6.3,F6.3,I3)'),i_l,g_l,gamma_ik(i_l,g_l,1),gamma_ik(i_l,g_l,2),gamma_ik(i_l,g_l,3),gamma_ik(i_l,g_l,4),gamma_ik(i_l,g_l,5),maxloc(gamma_ik(i_l,g_l,:))
        end if
    end do;end do
    close(9)   
    
    !Sample states based on the MLE of p
    call cluster_rand_assign(gamma_ik,sample_k_ini)
    
    !Posterior distribution of transition parameters given the sampled states
    call initial_conditions_tr(c_tr,cov_tr)
    
    open(unit=9,file=path//path_s_ini//'c_tr_ini'//s_c//'.txt',status='replace')
        write(9,'(F14.10)') c_tr
    close(9)
    open(unit=9,file=path//path_s_ini//'cov_tr_ini'//s_c//'.txt',status='replace')
        write(9,'(F14.10)') cov_tr
    close(9)  
end subroutine
    
subroutine initial_conditions_p(p,gamma_ik)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(out)::p
    double precision,dimension(indv,generations,clusters),intent(out)::gamma_ik
    double precision::indv_all,E_L_new,E_L_old,crit,stop_rule=1.0e-12
    integer::i_l,g_l,c_l,it
    double precision,dimension(clusters)::pi_k
    
    !EM algortihm for finding the MLE of p
    indv_all=0.0d0
    do i_l=1,indv; do g_l=first_age(i_l),last_age(i_l)
        if (data(i_l,1,g_l)/=-1) then
            indv_all=indv_all+1.0d0
        end if
    end do; end do
    
    do c_l=1,clusters
        p(:,c_l)=1/(dble(clusters)+1)*c_l
    end do
        
    pi_k=1/(dble(clusters))
    
    E_L_new=-999999999
    crit=1.0
    it=0
    do while (crit>stop_rule)
        it=it+1
        E_L_old=E_L_new
        call E_part(p,pi_k,gamma_ik,E_L_new)  
        call M_part(indv_all,gamma_ik,p,pi_k)
        crit=abs(E_L_old-E_L_new)/abs(E_L_new)
        print*,'crit',crit
        print*,'E_L_new',E_L_new
        print*,''
        if (E_L_old>E_L_new) then
            print*,'error'
        end if
    end do
    
end subroutine
    
subroutine cluster_rand_assign(gamma_ik,sample_k)
    use global_var
    implicit none
    double precision,dimension(indv,generations,clusters),intent(in)::gamma_ik
    integer,dimension(indv,generations),intent(out)::sample_k
    double precision::u
    integer::i_l,g_l,ind
    
    sample_k=0
    do i_l=1,indv; do g_l=generations,1,-1
        if (data(i_l,1,g_l)==-1 .and. g_l>=first_age(i_l) .and. g_l<=last_age(i_l)) then  !Dead
            sample_k(i_l,g_l)=clusters+1
        elseif (g_l<first_age(i_l) .or. g_l>last_age(i_l)) then !Not interviewed
            sample_k(i_l,g_l)=-1
        else
            call RANDOM_NUMBER(u)
            ind=1
            do while (sample_k(i_l,g_l)==0)
            if (u<sum(gamma_ik(i_l,g_l,1:ind))) then
                sample_k(i_l,g_l)=ind
            end if
            ind=ind+1
            end do
        end if
    end do; end do
    
end subroutine

subroutine initial_conditions_tr(c_tr,cov_tr)
use global_var
implicit none
double precision,dimension(covariates*clusters*clusters,covariates*clusters*clusters),intent(out)::cov_tr
double precision,dimension(covariates*clusters**2,1),intent(out)::c_tr
double precision,dimension(covariates*clusters**2,1)::mean_tr,new_mean_tr,factor
double precision::posterior_beta,posterior_beta_g,u,s_d=2.4**2/(covariates*clusters**2),prior
integer::chg_beta,ind,it,burn=200000
double precision,dimension(covariates*clusters*clusters,covariates*clusters*clusters)::C_rnd_W
double precision,dimension(covariates*clusters*clusters,1)::c_tr_g
double precision::factor_A=1

!Initial guess
c_tr(:,1)=0.0d0

mean_tr=c_tr
cov_tr=0.0d0

call multinomial(c_tr,posterior_beta)
call prior_tr(c_tr,prior)
posterior_beta=posterior_beta+prior

chg_beta=0
ind=-1

factor(1:clusters*clusters,1)=1/200.0d0
factor(clusters*clusters+1:2*clusters*clusters,1)=1/10000.0d0
factor(2*clusters*clusters+1:3*clusters*clusters,1)=1/600.0d0
factor(3*clusters*clusters+1:4*clusters*clusters,1)=1/600.0d0
factor(4*clusters*clusters+1:5*clusters*clusters,1)=1/10000.0d0
factor(5*clusters*clusters+1:6*clusters*clusters,1)=1/10000.0d0

call C_random_walk(C_rnd_W,covariates*clusters*clusters,factor)

do it=0,1000000
    ind=ind+1
    
    if (it>=burn) then
        if (it==burn) then
            mean_tr=c_tr
        end if
        call compute_mean(mean_tr,c_tr,it-burn,clusters**2*covariates,new_mean_tr)
        call compute_cov(cov_tr,it-burn,s_d,mean_tr,new_mean_tr,c_tr,clusters**2*covariates)
        mean_tr=new_mean_tr
    end if
    
    !Adaptive Metropolis-Hastings
    if (it>burn+1000) then
        C_rnd_W=cov_tr*factor_A
        call choldc(C_rnd_W,clusters**2*covariates)
    end if
        
    call proposal(c_tr,c_tr_g,C_rnd_W,covariates*clusters**2)

    call multinomial(c_tr_g,posterior_beta_g)
    call prior_tr(c_tr_g,prior)
    posterior_beta_g=posterior_beta_g+prior
    
    call RANDOM_NUMBER(u)  
    if (u<min(exp(posterior_beta_g-posterior_beta),1.0d0)) then
        c_tr=c_tr_g
        posterior_beta=posterior_beta_g
        chg_beta=chg_beta+1
    end if
    if (ind==1000 .or. it==0) then
        print*,'iteration',it
        print*,'acceptance rate:',real(chg_beta)/1000*100,'%'
        print*,'posterior_beta:',posterior_beta
        print*,' '
        
        !if (real(chg_beta)/1000*100<15) then
        !    print*,'got here'
        !    factor_A=factor_A/1.5
        !end if
        
        chg_beta=0
        ind=0
        
        !call save_results(c_tr,0.0d0,posterior_beta,0.0d0,it)
    end if
    
end do


end subroutine
    
SUBROUTINE choldc(a,n)
    IMPLICIT NONE
    integer,intent(in)::n
    DOUBLE PRECISION, DIMENSION(n,n), INTENT(INOUT) :: a
    DOUBLE PRECISION, DIMENSION(n) :: p
    INTEGER :: i,j
    DOUBLE PRECISION :: summ
    do i=1,n
	    summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
	    p(i)=sqrt(summ)
	    a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
    end do
    
    do i=1,n
        do j=1,n
            if (i==j) then
                a(i,i)=p(i)
            elseif (i<j) then
                a(i,j)=0
            end if
        end do
    end do
    
END SUBROUTINE choldc