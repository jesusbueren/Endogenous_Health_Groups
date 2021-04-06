subroutine c2beta(c_tr,betas,alphas,c_e,c_gender,c_axg,c_axe)
    use global_var
    implicit none
    double precision,dimension(clusters**2*covariates,1),intent(in)::c_tr
    double precision,dimension(clusters,clusters+1),intent(out)::betas,alphas,c_e,c_gender,c_axg,c_axe
    integer::c_l,c_l2,ind
    
    do c_l2=1,clusters;do c_l=1,clusters; 
        ind=c_l+(c_l2-1)*clusters
        betas(c_l,c_l2)=c_tr(ind,1)
        alphas(c_l,c_l2)=c_tr(ind+clusters**2,1)
        c_e(c_l,c_l2)=c_tr(ind+clusters**2*2,1)
        c_gender(c_l,c_l2)=c_tr(ind+clusters**2*3,1)
        c_axg(c_l,c_l2)=c_tr(ind+clusters**2*4,1)
        c_axe(c_l,c_l2)=c_tr(ind+clusters**2*5,1)
    end do;end do
    
    betas(:,clusters+1)=0.0d0
    alphas(:,clusters+1)=0.0d0
    c_e(:,clusters+1)=0.0d0
    c_gender(:,clusters+1)=0.0d0    
    c_axg(:,clusters+1)=0.0d0
    c_axe(:,clusters+1)=0.0d0
    
end subroutine
    
subroutine beta2c(betas,alphas,c_e,c_gender,c_axg,c_axe,c_tr)
    use global_var
    implicit none
    double precision,dimension(clusters**2*covariates,1),intent(out)::c_tr
    double precision,dimension(clusters,clusters+1),intent(in)::betas,alphas,c_e,c_gender,c_axg,c_axe
    
    c_tr(:,1)=[ reshape(betas(:,1:clusters),(/clusters**2,1/)), reshape(alphas(:,1:clusters),(/clusters**2,1/)) &
              , reshape(c_e(:,1:clusters),(/clusters**2,1/)), reshape(c_gender(:,1:clusters),(/clusters**2,1/)) &
              , reshape(c_axg(:,1:clusters),(/clusters**2,1/)), reshape(c_axe(:,1:clusters),(/clusters**2,1/)) ]    
    
end subroutine