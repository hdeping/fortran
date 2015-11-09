program main
    !use mkl95_lapack
    use lapack95
    
    real*8        ::    M(3,3),Z(3,3),A(3,3),B(3,3)
    
    M(1,:)=(/2.,-1.,0./)
    M(2,:)=(/-1.,2.,-1./)
    M(3,:)=(/0.,-1.,3./)
    
    A=M
    call potrf(A,uplo='L')
    
    do i=1,2
        A(i,i+1:)=0.
    enddo
    Z=A
    call potri(Z,uplo='L')
    do i=1,2
        Z(i,i+1:)=Z(i+1:,i)
    enddo
    print*, 'M:'
    do i=1,3
        print*,M(i,:)
    enddo
    print*, 'inv(M):'
    do i=1,3
        print*,Z(i,:)
    enddo
    print*, 'inv(M)*M:'
    do i=1,3
        print*,sum(M(i,:)*Z(:,1)),sum(M(i,:)*Z(:,2)),sum(M(i,:)*Z(:,3))
    enddo
    print*, 'A:'
    do i=1,3
        print*,A(i,:)
    enddo
    print*, 'A*AT:'
    do i=1,3
        print*,sum(A(i,:)*A(1,:)),sum(A(i,:)*A(2,:)),sum(A(i,:)*A(3,:))
    enddo
    
    B=Z
    call potrf(B,uplo='L')
    do i=1,2
        B(i,i+1:)=0.
    enddo
    print*, 'B:'
    do i=1,3
        print*,B(i,:)
    enddo
    print*, 'B*BT:'
    do i=1,3
        print*,sum(B(i,:)*B(1,:)),sum(B(i,:)*B(2,:)),sum(B(i,:)*B(3,:))
    enddo
    
    
end program
