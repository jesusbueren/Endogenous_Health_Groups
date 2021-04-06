subroutine full_posterior(p,c_tr,C_rnd_W_tr_or)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(inout)::p
    double precision,dimension(variables,clusters)::p_g
    double precision,dimension(covariates*clusters**2,covariates*clusters**2)::C_rnd_W_tr_or,C_rnd_W_tr,C_A_tr
    double precision,dimension(covariates*clusters**2,1)::c_tr,c_tr_g,mean_tr,new_mean_tr
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender)::H,H_g
    double precision,dimension(indv,clusters,generations)::likelihood,filtered_states
    double precision,dimension(indv,clusters+1,generations)::smoothed_states
    integer,dimension(indv,generations)::sample_k
    double precision::posterior_mu,posterior_beta,posterior_mu_g,posterior_beta_g,factor_p=1.0D-4,u,s_d=2.4**2/(covariates*clusters**2),s_d2=2.4**2/(variables*clusters)
    double precision,dimension(variables*clusters,1)::c_p,c_p_g,mean_p,new_mean_p
    double precision,dimension(variables*clusters,variables*clusters)::C_rnd_W_p,C_A_p
    integer::c_l,c_l2,it,ind,v_l,chg_mu,chg_beta,it2,burn=200000,burn2=200000,it3
    double precision::factor_beta,factor_mu,prior,prior_g
    double precision,dimension(clusters,L_e,L_gender)::init_cond
    character::continue_program
    !Timer variables
    real:: Times1,Times2
    INTEGER:: I,J, iTimes1,iTimes2, rate
    
1   C_rnd_W_tr=C_rnd_W_tr_or*2.5d0
    C_rnd_W_p=0.0d0
    do c_l=1,variables*clusters; do c_l2=1,variables*clusters
        if (c_l==c_l2) then
            C_rnd_W_p(c_l,c_l2)=factor_p
        end if
    end do;end do
     
    call transitions(c_tr,H)
    !Initial guess for the smoothed states
    smoothed_states=0.0d0
    smoothed_states(:,1:clusters,:)=1.0d0/dble(clusters)
    
    call likelihood_all(p,likelihood)
    call filtration(H,p,smoothed_states,likelihood,filtered_states)
    
    it2=-1
    it3=-1
    chg_mu=0
    chg_beta=0
    factor_beta=1.0d0
    factor_mu=1.0d0
    
    !Timer
    !CALL system_clock(count_rate=rate)
    !call CPU_TIME(Times1)
    !call SYSTEM_CLOCK(iTimes1)
    
    do it=0,100000000
        it2=it2+1  
        it3=it3+1
        call smoothing(filtered_states,H,p,sample_k,smoothed_states,init_cond)
        call prior_tr(c_tr,prior)
        call posterior_fct(H,p,sample_k,posterior_mu,posterior_beta,prior)
        !Proposal transitions
            !Adaptive
            if (it>=burn) then
                if (it==burn) then
                    mean_tr=c_tr
                end if
                call compute_mean(mean_tr,c_tr,it-burn,clusters**2*covariates,new_mean_tr)
                call compute_cov(C_A_tr,it-burn,s_d,mean_tr,new_mean_tr,c_tr,clusters**2*covariates)
                mean_tr=new_mean_tr
            end if
            if (it>burn+1000) then
                C_rnd_W_tr=C_A_tr*factor_beta
                call choldc(C_rnd_W_tr,clusters**2*covariates)
                do c_l=1,clusters**2*covariates; do c_l2=1,clusters**2*covariates
                    if (isnan(C_rnd_W_tr(c_l,c_l2)))then
                        print*,'c_l',c_l,'c_l2',c_l2
                        print*,C_A_tr
                        READ*, continue_program
                    end if            
                end do; end do
            end if
            if (it==burn+1000) then
                factor_beta=1.0d0
            end if
        call proposal(c_tr,c_tr_g,C_rnd_W_tr*factor_beta,covariates*clusters**2)
        
        !Proposal mixture of Bernouilli parameter
        c_p(:,:)=reshape(p(:,:),(/variables*clusters,1/))
            !Adaptive
            if (it>=burn2) then
                if (it==burn2) then
                    mean_p=c_p
                end if
                call compute_mean(mean_p,c_p,it-burn2,clusters*variables,new_mean_p)
                call compute_cov(C_A_p,it-burn2,s_d2,mean_p,new_mean_p,c_p,clusters*variables)
                mean_p=new_mean_p
            end if
            if (it>burn2+1000) then
                C_rnd_W_p=C_A_p*factor_mu
                call choldc(C_rnd_W_p,clusters*variables)
            end if
        if (it==burn2+1000) then
            factor_mu=1.0d0
        end if  
        call proposal(c_p,c_p_g,C_rnd_W_p*factor_mu,variables*clusters)
        do c_l=1,clusters;do v_l=1,variables
            ind=v_l+(c_l-1)*variables
            p_g(v_l,c_l)=min(max(c_p_g(ind,1),1.0d-6),0.9999d0)
        end do;end do
        
        !Evaluate posterior at new guess
        call transitions(c_tr_g,H_g)
        call prior_tr(c_tr_g,prior_g)
        call posterior_fct(H_g,p_g,sample_k,posterior_mu_g,posterior_beta_g,prior_g)
        if (isnan(posterior_beta_g)) then
            print*,'isnan posterior_beta_g'
        end if
        
        !Acceptance or rejection of new guess
        call RANDOM_NUMBER(u)
        ind=0
        if (u<min(exp(posterior_mu_g-posterior_mu),1.0d0)) then
            p=p_g
            call likelihood_all(p,likelihood) 
            ind=1
            chg_mu=chg_mu+1
        end if     
        call RANDOM_NUMBER(u)
        if (u<min(exp(posterior_beta_g-posterior_beta),1.0d0)) then
            c_tr=c_tr_g
            H=H_g
            ind=1
            chg_beta=chg_beta+1
        end if
        if (ind==1)then
            call filtration(H,p,smoothed_states,likelihood,filtered_states)
        end if
        !Print acceptance rate
        if (it2==50 .or. it==0) then
            it2=0
            call save_results(c_tr,p,posterior_beta,posterior_mu,init_cond,it)
        end if
        !Adapt the variance of the proposal
        if ((it3==50 .and. it<100000) .or. (it3==500 .and. it>=100000) ) then
            print*,'completed ',it,' iterations'
            if (it>0) then
                print*,'acceptance rate:'
                print*,'mu:',real(chg_mu)/real(it3)*100,'%'
                print*,'beta:',real(chg_beta)/real(it3)*100,'%' 
            end if
            print*,'posterior:',real(posterior_beta)+real(posterior_mu)
            print*,' '
            
            if (real(chg_mu)/real(it3)*100<=15) then
                factor_mu=factor_mu/2.0d0
                print*,'got here A','factor: ',factor_mu
            end if
            if (real(chg_mu)/real(it3)*100>=50) then
                factor_mu=factor_mu*1.5
                print*,'got here B','factor: ',factor_mu
            end if
            if (real(chg_beta)/real(it3)*100<=15) then
                factor_beta=factor_beta/2.0d0
                print*,'got here C','factor: ',factor_beta
            end if
            if (real(chg_beta)/real(it3)*100>=50) then
                factor_beta=factor_beta*1.5
                print*,'got here D','factor: ',factor_beta
            end if
            chg_mu=0
            chg_beta=0
            it3=0
        end if
    
    end do
           
end     
