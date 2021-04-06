subroutine charge_data()
    use global_var
    implicit none 
    integer,dimension(indv*variables*generations)::data_vec
    integer,dimension(indv*2)::ages
    integer::i_l,g_l,index
    
    open(unit=10,file=path//"Data\ages_all.csv")
        read(10,*) ages
    close(10)
    open(unit=10,file=path//"Data\gender_all.csv")
        read(10,*) gender
    close(10)
    open(unit=10,file=path//"Data\educ_all.csv")
        read(10,*) educ
    close(10) 
    open(unit=10,file=path//"Data\data_all.csv")
        read(10,*) data_vec
    close(10)
    
    do i_l=1,indv
        do g_l=1,generations 
            index=(i_l-1)*(generations*variables)+(g_l-1)*variables+1
            data(i_l,:,g_l)=data_vec(index:index+variables-1)
        end do
    end do
    
    do i_l=1,indv
        index=(i_l-1)*2+1
        first_age(i_l)=ages(index)
        last_age(i_l)=ages(index+1)
        if (educ(i_l)>2) then
            educ(i_l)=2
        end if
    end do
end subroutine
    