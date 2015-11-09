program main
use common_module
use Go_Structure
use GO_force
use dynamic_module

real*8    :: Q
real*8    :: x_tmp

CALL RANDOM_SEED()

do jF=1,1
!    F=dble(jF-20)*0.2
!    temp=dble(jF*50)+1000.
    call initial(jF)
    write(fName,"('ts_',I2.2,'.dat')") jF
    OPEN(3, FILE=fName)
    write(fName,"('cfg_',I2.2,'.xyz')") jF
    open(1,file=fName)
    call Save2xyz( 1 )
    it=0;cnt_it=0
    call getNonNeighbor
    Q=sum((R_nonNei-R0_nonNei)**2)
    Q=sqrt(2.*Q/dble(num**2-3*num+2))
    print*,temp,F,Q
    print*,matrix(1,4),matrix(4,28)
    print*,EPS_bond,EPS_theta,EPS_DiTheta(1,4),EPS_nonNei(1,4),EPS_nonNei(4,28)
!    pause
    DO it=1,nt
!        if(it==20000) F=0.
        call RunDynamic()
!        f_printV=0
        call getOP(op)
        if(mod(it,20)==0) then
!            f_printV=1
            call Save2xyz( 1 )
            call getNonNeighbor
            Q=sum((R_nonNei-R0_nonNei)**2)
            Q=sqrt(2.*Q/dble(num**2-3*num+2))
            write(3,"(10G18.9)") dble(it*nT_jump)*dt,Q,op,R_bond(1),R0_bond(1),Theta(1),Theta0(1),DiTheta(1),DiTheta01(1)
            PRINT*, dble(it*nT_jump)*dt,Q,op,F
        endif
    enddo 
    pause
    write(fName,"('initial/initialX.dat')")
    OPEN(10, FILE=fName)
    write(fName,"('initial/initialV.dat')")
    OPEN(11, FILE=fName)
    write(fName,"('initial/initialF.dat')")
    OPEN(12, FILE=fName)
    do it=1,Num
        write(10,'(3G18.9)') x(it,:)    
        write(11,'(3G18.9)') V(it,:)
        write(12,'(3G18.9)') F_tot(it,:)
    enddo !it    
    close(10)
    close(11)
    close(12)
    write(fName,"('cfg_',I2.2,'.dat')") jf
    OPEN(4, FILE=fName)
    DO I=1, NUM
        WRITE(4, '(A4, I7, A4, A5, I6, F12.3, 2F8.3, A22 )') 'ATOM', I, 'CA', 'ALA', I, X(I,:), AF
    ENDDO
    close(4)
    close(1)
    close(3)
enddo !irun
end program main
