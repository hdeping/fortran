module parameter_module
include 'mpif.h'
!include 'omp_lib.h'
integer, parameter                    ::                n=128, sn=400, ini=40, stn=2, wstep=1000
integer                                ::                istep
real*8, parameter                    ::                pi=3.1415926
real*8                                ::                dk, dr, h, dt, dro
real*8, dimension(0:n-1)            ::                sk, ck, gammak, ak
real*8                                ::                df, di

!--**********************并行参数*************************--!

integer                                ::                rank, ierror
integer, parameter                    ::                node=2
integer, parameter                    ::                nnode=n/node
integer                                            ::                thread=1
integer                                            ::                sprocs, senode
integer, dimension(mpi_status_size)                ::                status
integer, dimension(0:node-1, 2)                    ::                inode

real*8, dimension(:, :), allocatable            ::            dphi, ke_me
real*8, dimension(0:n-1)                        ::            phi

real*8                                            ::            rho
real*8                                            ::            det_time, beta
real*8                                            ::            v0



end module parameter_module
