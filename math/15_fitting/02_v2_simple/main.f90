program main
    use module_common
    implicit none
    !real(8)          :: test(2*m,2*m)
    real(8)          :: test(2*m)
    
    !*****************  get data
    filename = "data.txt"
    open(10,file = filename)
    do i = 1,m
        a_para(i) = dble(i)
        b_para(i) = dble(m + 1 - i)
    end do
    do  i = 1,n
        x(i) = dble(i)*dt
    end do
    y = funx(x,a_para,b_para)
    do i = 1,n
        write(10,*)x(i),y(i)
    end do
    
    close(10)
    !  get parameter
    !call getpara(anew,bnew,x,y)
    !print "('a = ',<m>f18.5)",anew(:)
    !print "('b = ',<m>f18.5)",bnew(:)


    !  test 
    !a_para = 3.0
    test = equations(a_para,b_para,x,y)
    print *,test
    !test = getpartial(a_para,b_para,x,y)
    


    

end program main
