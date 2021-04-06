subroutine smoothing(filtered_states,H,p,sample_k,smoothed_states,init_cond)
    use global_var
    implicit none
    integer,dimension(indv,generations),intent(out)::sample_k
    double precision,dimension(indv,clusters,generations),intent(in)::filtered_states
    double precision,dimension(indv,clusters+1,generations),intent(out)::smoothed_states
    double precision,dimension(clusters+1,clusters+1,generations,L_e,L_gender),intent(in)::H
    double precision,dimension(variables,clusters),intent(in)::p
    double precision,dimension(clusters,L_e,L_gender),intent(out)::init_cond
    double precision,dimension(clusters)::smoothed_states_cond
    double precision,dimension(clusters+1,clusters)::smoothed_states_uncond
    integer::g_l,i_l,ind,c_l2,e_l
    double precision::u
    integer,dimension(L_e,L_gender)::counter
    character::pause_key
    
    !Initialize smoother
    smoothed_states=0.0d0
    sample_k(:,:)=0
    counter=0
    init_cond=0.0d0
    
    !$OMP PARALLEL num_threads(3) DEFAULT(SHARED) PRIVATE(g_l,i_l,u,ind,smoothed_states_cond,smoothed_states_uncond)
    !$OMP DO SCHEDULE(GUIDED)
    do i_l=1,indv;do g_l=last_age(i_l),1,-1;
        if (data(i_l,1,g_l)==-1 .and. g_l>=first_age(i_l) .and. g_l<=last_age(i_l)) then  !Dead
            sample_k(i_l,g_l)=clusters+1
            smoothed_states(i_l,:,g_l)=0.0d0
            smoothed_states(i_l,clusters+1,g_l)=1.0d0
        elseif (g_l>last_age(i_l)) then !Not interviewed
            sample_k(i_l,g_l)=-1
        elseif (g_l==last_age(i_l))  then !Last interview when alive if always alive 
            smoothed_states(i_l,:,g_l)=0.0d0
            smoothed_states(i_l,1:clusters,g_l)=filtered_states(i_l,1:clusters,g_l)
            call RANDOM_NUMBER(u)
            ind=1
            do while (sample_k(i_l,g_l)==0)
                if (ind<clusters) then
                    if (u<sum(filtered_states(i_l,1:ind,g_l))) then 
                        sample_k(i_l,g_l)=ind
                    end if
                else
                    sample_k(i_l,g_l)=ind
                end if
                ind=ind+1
            end do 
        elseif (g_l<last_age(i_l) .or. (data(i_l,1,g_l)/=-1 .and. sample_k(i_l,g_l+1)==clusters+1)) then !not last interview and alive or last alive interview
            if (sample_k(i_l,g_l+1)==0) then
                print*,'Error'
            end if
            !Kim smoother
            smoothed_states_cond(:)=(H(1:clusters,sample_k(i_l,g_l+1),g_l,educ(i_l),gender(i_l))*filtered_states(i_l,1:clusters,g_l))/sum((H(1:clusters,sample_k(i_l,g_l+1),g_l,educ(i_l),gender(i_l))*filtered_states(i_l,1:clusters,g_l))) !revisa
            call RANDOM_NUMBER(u)
            ind=1 
            do while (sample_k(i_l,g_l)==0)
                if (ind<clusters) then
                    if (u<sum(smoothed_states_cond(1:ind))) then
                        sample_k(i_l,g_l)=ind
                    end if
                else
                    sample_k(i_l,g_l)=ind
                end if
                ind=ind+1
            end do
            !Hamilton smoother
            do c_l2=1,clusters+1
                if (smoothed_states(i_l,c_l2,g_l+1)==0.0d0)then
                    smoothed_states_uncond(c_l2,:)=0.0d0
                else
                    smoothed_states_uncond(c_l2,:)=smoothed_states(i_l,c_l2,g_l+1)* &
                                                (H(1:clusters,c_l2,g_l,educ(i_l),gender(i_l))*filtered_states(i_l,1:clusters,g_l))/&
                    sum((H(1:clusters,c_l2,g_l,educ(i_l),gender(i_l))*filtered_states(i_l,1:clusters,g_l)))
                end if
            end do
            smoothed_states(i_l,1:clusters,g_l)=sum(smoothed_states_uncond,1)
        end if
    end do;
    counter(educ(i_l),gender(i_l))=counter(educ(i_l),gender(i_l))+1
    init_cond(:,educ(i_l),gender(i_l))=init_cond(:,educ(i_l),gender(i_l))+smoothed_states(i_l,1:clusters,1)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    do e_l=1,L_e; do g_l=1,L_gender
        init_cond(:,e_l,g_l)=init_cond(:,e_l,g_l)/dble(counter(e_l,g_l))
    end do;end do
    
    
end subroutine
    
