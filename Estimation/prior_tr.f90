subroutine prior_tr(c_tr,prior)
    use global_var
    double precision,dimension(covariates*clusters**2,1),intent(in)::c_tr
    double precision,intent(out)::prior
    integer::i
    double precision::sigma2,pi
    sigma2=100.0d0
    pi=4.0d0*atan(1.0d0)
    
    prior=0.0d0
    do i=1,covariates*clusters**2
        prior=prior-1.0d0/2.0d0*log(2.0d0*pi*sigma2)-1.0d0/2.0d0*(c_tr(i,1)**2.0d0/sigma2)
    end do
    
    
end 