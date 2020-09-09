!     algorithm 419 collected algorithms from acm.
!     algorithm appeared in comm. acm, vol. 15, no. 02,
!     p. 097.
!
! added changes from remark on algorithm 419 by david h. withers
! cacm (march 1974) vol 17 no 3 p. 157 
program cpolydr
!
!driver to test cpoly
!
!
    dimension mult(50)
    logical flag(50)
    real*8 fzr(50),fzi(50),bnd(50)
    logical fail
    double precision  p(50),pi(50),zr(50),zi(50)
    write(6,100)
  100 format('1example 1.  polynomial with zeros 1,2,...,10.')
    p(1) = 1
    p(2) = -55
    p(3) = 1320
    p(4) = -18150
    p(5) = 157773
    p(6) = -902055
    p(7)  =  3416930
    p(8) = -8409500
    p(9) = 12753576
    p(10) = -10628640
    p(11) = 3628800
    do 10 i = 1,11
   10 pi(i) = 0
    call prtc(11,p,pi)
    call cpoly(p,pi,10,zr,zi,fail)
    if(fail) go to 95
    call prtz (10,zr,zi)
    2 write(6,101)
  101 format('1example 2. zeros on imaginary axis degree 3.')
    p(1) = 1
    p(2) = 0
    p(3) = -10001.0001d0
    p(4) = 0
    pi(1) = 0
    pi(2) = -10001.0001d0
    pi(3) = 0
    pi(4) = 1
    call prtc(4,p,pi)
    call cpoly(p,pi,3,zr,zi,fail)
    if (fail) go to 96
    call prtz (3,zr,zi)
    3 write(6,102)
  102 format('1example 3. zeros at 1+i,1/2*(1+i)....1/(2**-9)*(1+i)')
    p(1) = 1.0
    p(2) = -1.998046875
    p(3) = 0.0
    p(4) = .7567065954208374d0
    p(5) = -.2002119533717632d0
    p(6) = 1.271507365163416d-2
    p(7) = 0
    p(8) = -1.154642632172909d-5
    p(9) = 1.584803612786345d-7
    p(10) = -4.652065399568528d-10
    p(11) = 0
    pi(1) = 0
    pi(2) = p(2)
    pi(3) = 2.658859252929688d0
    pi(4) = -7.567065954208374d-1
    pi(5) = 0
    pi(6) = p(6)
    pi(7) = -7.820779428584501d-4
    pi(8) = -p(8)
    pi(9) = 0
    pi(10) = p(10)
    pi(11) = 9.094947017729282d-13
    call prtc(11,p,pi)
    call cpoly(p,pi,10,zr,zi,fail)
    if (fail) go to 97
    call prtz(10,zr,zi)
    4 write(6,103)
  103 format('1example 4. multiple zeros')
    p(1) = 1
    p(2) = -10
    p(3) = 3
    p(4) = 284
    p(5) = -1293
    p(6) = 2374
    p(7) = -1587
    p(8) = -920
    p(9) = 2204
    p(10) = -1344
    p(11) = 288
    pi(1) = 0
    pi(2) = -10
    pi(3) = 100
    pi(4) = -334
    pi(5) = 200
    pi(6) = 1394
    pi(7)  = -3836
    pi(8) = 4334
    pi(9) = -2352
    pi(10) = 504
    pi(11) = 0
    call prtc(11,p,pi)
    call cpoly(p,pi,10,zr,zi,fail)
    if (fail) go to 98
    call prtz(10,zr,zi)
    5 write(6,104)
  104 format('1example 5. 12 zeros evenly distribute on a circle of radi
     -us 1. centered at 0+2i.')
    p(1) = 1
    p(2) = 0
    p(3) = -264
    p(4) = 0
    p(5) = 7920
    p(6) = 0
    p(7) = -59136
    p(8) = 0
    p(9) = 126720
    p(10) = 0
    p(11) = -67584
    p(12) = 0
    p(13) = 4095
    pi(1) = 0
    pi(2) = -24
    pi(3) = 0
    pi(4) = 1760
    pi(5) = 0
    pi(6) = -25344
    pi(7) = 0
    pi(8) = 101376
    pi(9) = 0
    pi(10) = -112640
    pi(11) = 0
    pi(12) = 24576
    pi(13) = 0
    call prtc(13,p,pi)
    call cpoly(p,pi,12,zr,zi,fail)
    if(fail) go to 99
    call prtz(12,zr,zi)
    return
   95 write(6,105)
    go to 2
   96 write(6,105)
    go to 3
   97 write(6,105)
    go to 4
   98 write(6,105)
    go to 5
   99 write(6,105)
    return
  105 format(//' cpoly has failed on this example')
    end program cpolydr
