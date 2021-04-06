subroutine E_part(p,pi_k,gamma_ik,E_L_new)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(in)::p
    double precision,dimension(clusters),intent(in)::pi_k
    double precision,dimension(indv,generations,clusters),intent(out)::gamma_ik
    double precision,intent(out)::E_L_new
    double precision,dimension(indv,clusters,generations)::L_i
    double precision,dimension(indv,generations,clusters)::E_L_i
    integer::g_l,i_l,c_l
    double precision::A
    
    call likelihood_all(p,L_i)
    E_L_i=0.0d0
    gamma_ik=0.0d0
    !$OMP PARALLEL num_threads(3) DEFAULT(SHARED) PRIVATE(g_l,i_l,A,c_l)
    !$OMP DO
    do i_l=1,indv; do g_l=1,generations
      if (data(i_l,1,g_l) /= -1) then
        A=sum(pi_k*L_i(i_l,:,g_l))
        do c_l=1,clusters
            gamma_ik(i_l,g_l,c_l)=(pi_k(c_l)*L_i(i_l,c_l,g_l))/A
            E_L_i(i_l,g_l,c_l)=gamma_ik(i_l,g_l,c_l)*(log(pi_k(c_l))+log(L_i(i_l,c_l,g_l)))
            if (isnan(E_L_i(i_l,g_l,c_l))) then
                print*,'error'
            end if
        end do
      end if
    end do; end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    E_L_new=sum(E_L_i)
    
    end subroutine
    
    subroutine M_part(indv_all,gamma_ik,p,pi_k)
    use global_var
    implicit none
    double precision,dimension(variables,clusters),intent(out)::p
    double precision,dimension(clusters),intent(out)::pi_k
    double precision,dimension(indv,generations,clusters),intent(in)::gamma_ik
    double precision,intent(in)::indv_all
    double precision,dimension(clusters)::N_k
    double precision,dimension(indv,generations)::p_i
    integer::c_l,v_l,g_l,i_l
    
    N_k=0
    do c_l=1,clusters
      N_k(c_l)=sum(gamma_ik(:,:,c_l))
      pi_k(c_l)=N_k(c_l)/indv_all
    end do

    do c_l=1,clusters; do v_l=1,variables
      p_i=0
      do i_l=1,indv;do g_l=1,generations
        if (data(i_l,1,g_l)/=-1 .and. data(i_l,v_l,g_l)/=-9) then
          p_i(i_l,g_l)=gamma_ik(i_l,g_l,c_l)*data(i_l,v_l,g_l)
        elseif (data(i_l,v_l,g_l)==-9) then
          p_i(i_l,g_l)=gamma_ik(i_l,g_l,c_l)*1/dble(clusters)
        end if
      end do; end do
      p(v_l,c_l)=1/N_k(c_l)*sum(p_i(:,:))
    end do; end do
    
end subroutine
    