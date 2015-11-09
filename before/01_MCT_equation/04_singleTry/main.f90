program main
    !use module_fst 
    use module_mct
    integer               :: itmp
    integer               :: jtmp

    ! get dk
    do itmp = 1,n
        dk(itmp) = dble(itmp)*h
        !print *,dk(itmp)
        !pause
    end do
    ! get sk
    filename = "sk.txt"
    open(10,file = filename,status = "old",iostat = ierror)

    do itmp  = 1,ncut
        read(10,*,iostat = ierror)xtmp, sk(itmp)
        sk(itmp) = sk(itmp)/smul
        !print *,sk(itmp)
        !pause
        if(ierror /= 0)exit
    end do
    close(10)
    write(filename,"('fq',i2.2,'.txt')")int(smul)
    !filename = "data.txt"
    open(10,file = filename)
    
    ! get final f
    finalf = getfinal()
    lambda = judge(finalf,sk(1:ncut),ncut)
    print *,"lambda = ",lambda

    do itmp = 1,ncut
        write(10,"(4f18.6)")dk(itmp),finalf(itmp)/sk(itmp)
    end do
    close(10)

end program main
!code before{{{
    !do q = 1,300
    !        
    !    wri
    !    open(10,file = filename)
    !    do i = 1,n
    !        write(10,"(4f12.5)")(i - 1)*h,f(1,1,i,3),f(1,2,i,3),f(2,2,i,3)
    !    end do
    !    close(10)
    !end do

    !filename = "data.dat"
    !open(10,file = filename)
    !do i = 1,ncut
    !    tmp = (i*deltar)**2.0
    !    write(10,*)i*deltar,tmp
    !end do
    !close(10)



!}}}
