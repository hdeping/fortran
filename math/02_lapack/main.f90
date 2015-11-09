program main
    use mkl95_lapack

  implicit none
  integer,parameter     :: n=3
  integer,parameter     :: lda=n
  integer,parameter     :: lwork=2*n-1
  character             :: jobz
  character             :: uplo
  integer               :: info=0
  real(8)               :: rwork(3*n-2)
  real(8)               :: w(n)
  complex(16)           :: a(lda,n)
  complex(16)           :: work(lwork)
  complex(16)           :: b(5,3)
  a(1,:) = (/1,2,3/)
  a(1,:) = (/2,4,5/)
  a(1,:) = (/3,5,6/)
  b(2:4,:) = a

 !输出中a的每一列是一个本征矢。
 !call zheev('v', 'u', n, a, lda, w, work, lwork, rwork, info )

  call heev(a, w, 'V','U',info )
  print*, w
  print*, '--------------------   --------------------    ------------------'
  print*, real(a(1,:))
  print*, real(a(2,:))
  print*, real(a(3,:))
  print*
  !call zheev('v', 'u', n, b(2:5,:), 4, w, work, lwork, rwork, info )
  print*, w
  print*, '--------------------   --------------------    ------------------'
  print '(3f25.15)', real(b(1,:))
  print '(3f25.15)', real(b(2,:))
  print '(3f25.15)', real(b(3,:))
  print '(3f25.15)', real(b(4,:))
  print '(3f25.15)', real(b(5,:))

end program main
