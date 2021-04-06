subroutine likelihood_all(p,likelihood)
    use global_var
    implicit none
    double precision,dimension(indv,clusters,generations),intent(out)::likelihood
    integer::g_l,i_l,v_l,c_l
    double precision,dimension(variables,clusters),intent(in)::p
    
    likelihood(:,:,:)=1.0d0
    !$OMP PARALLEL num_threads(3) DEFAULT(SHARED) PRIVATE(g_l,i_l,v_l,c_l)
    !$OMP DO SCHEDULE(GUIDED)
    do i_l=1,indv;do g_l=first_age(i_l),last_age(i_l)    
        if (data(i_l,1,g_l)/=-1) then
          do v_l=1,variables; do c_l=1,clusters
            if (data(i_l,v_l,g_l)==1) then
                likelihood(i_l,c_l,g_l)=likelihood(i_l,c_l,g_l)*p(v_l,c_l)
            elseif (data(i_l,v_l,g_l)==0) then
                likelihood(i_l,c_l,g_l)=likelihood(i_l,c_l,g_l)*(1-p(v_l,c_l))
            end if 
          end do; end do      
        end if
    end do; end do
    !$OMP END DO
    !$OMP END PARALLEL
    end subroutine
    
subroutine likelihood_i(i_l,c_l,g_l,p,likelihood)
    use global_var
    implicit none
    double precision,intent(out)::likelihood
    integer,intent(in)::i_l,c_l,g_l
    double precision,dimension(variables,clusters),intent(in)::p
    integer::v_l
    
    likelihood=1.0d0
    do v_l=1,variables
        if (data(i_l,v_l,g_l)==1) then
            likelihood=likelihood*p(v_l,c_l)
        elseif (data(i_l,v_l,g_l)==0) then
            likelihood=likelihood*(1.0d0-p(v_l,c_l))
        end if
    end do
end subroutine
    
    