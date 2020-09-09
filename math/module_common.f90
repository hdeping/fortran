module module_common
    !implicit none
    !integer,parameter            :: n = 100
    !integer                      :: i
    !integer                      :: j
    !integer                      :: k
    !character(10)                :: filename 
    contains
!subroutine prtc{{{
subroutine prtc(n,p,q)
        
    double precision p(50),q(50)
    write(6,10) (p(i),q(i) ,i = 1,n)
   10 format(//' coefficients' /50(2d26.16/))
    return
    end subroutine prtc
!}}}
!subroutine prtz{{{
subroutine prtz(n,zr,zi)
    double precision zr(50),zi(50)
    write(6,10) (zr(i),zi(i) ,i = 1,n)
   10 format(//' zeros'/ 50(2d26.16/))
    return
end subroutine prtz
!}}}
!subroutine cpoly{{{
subroutine cpoly(opr,opi,degree,zeror,zeroi,fail)                 cpol  10
!
! added changes from remark on algorithm 419 by david h. withers
! cacm (march 1974) vol 17 no 3 p. 157 
!
! finds the zeros of a complex polynomial.
! opr, opi  -  double precision vectors of real and
! imaginary parts of the coefficients in
! order of decreasing powers.
! degree    -  integer degree of polynomial.
! zeror, zeroi  -  output double precision vectors of
! real and imaginary parts of the zeros.
! fail    -  output logical parameter,  true  only if
! leading coefficient is zero or if cpoly
! has found fewer than degree zeros.
! the program has been written to reduce the chance of overflow
! occurring. if it does occur, there is still a possibility that
! the zerofinder will work provided the overflowed quantity is
! replaced by a large number.
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
! to change the size of polynomials which can be solved, replace
! the dimension of the arrays in the common area.
    double precision xx,yy,cosr,sinr,smalno,base,xxx,zr,zi,bnd,
     *    opr(1),opi(1),zeror(1),zeroi(1),
     *    cmod,scale,cauchy,dsqrt
    logical fail,conv
    integer degree,cnt1,cnt2
! initialization of constants
    call mcon(eta,infin,smalno,base)
    are  =  eta
    mre  =  2.0d0*dsqrt(2.0d0)*eta
    xx  =  .70710678
    yy  =  -xx
    cosr  =  -.069756474
    sinr  =  .99756405
    fail  =  .false.
    nn  =  degree+1
! algorithm fails if the leading coefficient is zero.
    if (opr(1) .ne. 0.0d0 .or. opi(1) .ne. 0.0d0) go to 10
        fail  =  .true.
        return
! remove the zeros at the origin if any.
   10 if (opr(nn) .ne. 0.0d0 .or. opi(nn) .ne. 0.0d0) go to 20
        idnn2  =  degree-nn+2
        zeror(idnn2)  =  0.0d0
        zeroi(idnn2)  =  0.0d0
        nn  =  nn-1
        go to 10
! make a copy of the coefficients.
   20 do 30 i  =  1,nn
        pr(i)  =  opr(i)
        pi(i)  =  opi(i)
        shr(i)  =  cmod(pr(i),pi(i))
   30 continue
! scale the polynomial.
    bnd  =  scale (nn,shr,eta,infin,smalno,base)
    if (bnd .eq. 1.0d0) go to 40
    do 35 i  =  1,nn
        pr(i)  =  bnd*pr(i)
        pi(i)  =  bnd*pi(i)
   35 continue
! start the algorithm for one zero .
   40 if (nn.gt. 2) go to 50
! calculate the final zero and return.
        call cdivid(-pr(2),-pi(2),pr(1),pi(1),zeror(degree),
     *    zeroi(degree))
        return
! calculate bnd, a lower bound on the modulus of the zeros.
   50 do 60 i  =  1,nn
        shr(i)  =  cmod(pr(i),pi(i))
   60 continue
    bnd  =  cauchy(nn,shr,shi)
! outer loop to control 2 major passes with different sequences
! of shifts.
    do 100 cnt1  =  1,2
! first stage calculation, no shift.
        call noshft(5)
! inner loop to select a shift.
        do 90 cnt2  =  1,9
! shift is chosen with modulus bnd and amplitude rotated by
! 94 degrees from the previous shift
             xxx  =  cosr*xx-sinr*yy
             yy  =  sinr*xx+cosr*yy
             xx  =  xxx
             sr  =  bnd*xx
             si  =  bnd*yy
! second stage calculation, fixed shift.
             call fxshft(10*cnt2,zr,zi,conv)
             if (.not. conv) go to 80
! the second stage jumps directly to the third stage iteration.
! if successful the zero is stored and the polynomial deflated.
                  idnn2  =  degree-nn+2
                  zeror(idnn2)  =  zr
                  zeroi(idnn2)  =  zi
                  nn  =  nn-1
                  do 70 i  =  1,nn
                       pr(i)  =  qpr(i)
                       pi(i)  =  qpi(i)
   70             continue
                  go to 40
   80        continue
! if the iteration is unsuccessful another shift is chosen.
   90     continue
! if 9 shifts fail, the outer loop is repeated with another
! sequence of shifts.
  100 continue
! the zerofinder has failed on two major passes.
! return empty handed.
    fail  =  .true.
    return
end subroutine cpoly
!}}}
!subroutine  noshft{{{
subroutine  noshft(l1)                                            nosh1130
! computes  the derivative  polynomial as the initial h
! polynomial and computes l1 no-shift h polynomials.
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
    double precision xni,t1,t2,cmod
    n  =  nn-1
    nm1  =  n-1
    do 10 i  =  1,n
        xni  =  nn-i
        hr(i)  =  xni*pr(i)/float(n)
        hi(i)  =  xni*pi(i)/float(n)
   10 continue
    do 50 jj  =  1,l1
        if (cmod(hr(n),hi(n)) .le. eta*10.0d0*cmod(pr(n),pi(n)))
     *    go to 30
        call cdivid(-pr(nn),-pi(nn),hr(n),hi(n),tr,ti)
        do 20 i  =  1,nm1
             j  =  nn-i
             t1  =  hr(j-1)
             t2  =  hi(j-1)
             hr(j)  =  tr*t1-ti*t2+pr(j)
             hi(j)  =  tr*t2+ti*t1+pi(j)
   20     continue
        hr(1)  =  pr(1)
        hi(1)  =  pi(1)
        go to 50
! if the constant term is essentially zero, shift h coefficients.
   30     do 40 i  =  1,nm1
             j  =  nn-i
             hr(j)  =  hr(j-1)
             hi(j)  =  hi(j-1)
   40     continue
        hr(1)  =  0.0d0
        hi(1)  =  0.0d0
   50 continue
    return
end subroutine noshft
!}}}
!subroutine fxshft{{{
subroutine fxshft(l2,zr,zi,conv)                                  fxsh1550
! computes l2 fixed-shift h polynomials and tests for
! convergence.
! initiates a variable-shift iteration and returns with the
! approximate zero if successful.
! l2 - limit of fixed shift steps
! zr,zi - approximate zero if conv is .true.
! conv  - logical indicating convergence of stage 3 iteration
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
    double precision zr,zi,otr,oti,svsr,svsi,cmod
        logical conv,test,pasd,bool
    n  =  nn-1
! evaluate p at s.
    call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
    test  =  .true.
    pasd  =  .false.
! calculate first t  =  -p(s)/h(s).
    call calct(bool)
! main loop for one second stage step.
    do 50 j  =  1,l2
        otr  =  tr
        oti  =  ti
! compute next h polynomial and new t.
        call nexth(bool)
        call calct(bool)
        zr  =  sr+tr
        zi  =  si+ti
! test for convergence unless stage 3 has failed once or this
! is the last h polynomial .
        if ( bool .or. .not. test .or. j .eq. l2) go to 50
        if (cmod(tr-otr,ti-oti) .ge. .5d0*cmod(zr,zi)) go to 40
             if (.not. pasd) go to 30
! the weak convergence test has been passed twice, start the
! third stage iteration, after saving the current h polynomial
! and shift.
                  do 10 i  =  1,n
                       shr(i)  =  hr(i)
                       shi(i)  =  hi(i)
   10             continue
                  svsr  =  sr
                  svsi  =  si
                  call vrshft(10,zr,zi,conv)
                  if (conv) return
! the iteration failed to converge. turn off testing and restore
! h,s,pv and t.
                  test  =  .false.
                  do 20 i  =  1,n
                       hr(i)  =  shr(i)
                       hi(i)  =  shi(i)
   20             continue
                  sr  =  svsr
                  si  =  svsi
                  call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
                  call calct(bool)
                  go to 50
   30        pasd  =  .true.
             go to 50
   40     pasd  =  .false.
   50 continue
! attempt an iteration with final h polynomial from second stage.
    call vrshft(10,zr,zi,conv)
    return
end subroutine fxshft
!}}}
!subroutine vrshft{{{
subroutine vrshft(l3,zr,zi,conv)                                  vrsh2230
! carries out the third stage iteration.
! l3 - limit of steps in stage 3.
! zr,zi   - on entry contains the initial iterate, if the
! iteration converges it contains the final iterate
! on exit.
! conv    -  .true. if iteration converges
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
    double precision zr,zi,mp,ms,omp,relstp,r1,r2,cmod,dsqrt,errev,tp
    logical conv,b,bool
    conv  =  .false.
    b  =  .false.
    sr  =  zr
    si  =  zi
! main loop for stage three
    do 60 i  =  1,l3
! evaluate p at s and test for convergence.
        call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
        mp  =  cmod(pvr,pvi)
        ms  =  cmod(sr,si)
        if (mp .gt. 20.0d0*errev(nn,qpr,qpi,ms,mp,are,mre))
     *     go to 10
! polynomial value is smaller in value than a bound on the error
! in evaluating p, terminate the iteration.
             conv  =  .true.
             zr  =  sr
             zi  =  si
             return
   10     if (i .eq. 1) go to 40
             if (b .or. mp .lt.omp .or. relstp .ge. .05d0)
     *          go to 30
! iteration has stalled. probably a cluster of zeros. do 5 fixed
! shift steps into the cluster to force one zero to dominate.
                  tp  =  relstp
                  b  =  .true.
                  if (relstp .lt. eta) tp  =  eta
                  r1  =  dsqrt(tp)
                  r2  =  sr*(1.0d0+r1)-si*r1
                  si  =  sr*r1+si*(1.0d0+r1)
                  sr  =  r2
                  call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
                  do 20 j  =  1,5
                       call calct(bool)
                       call nexth(bool)
   20             continue
    omp  =  infin
                  go to 50
! exit if polynomial value increases significantly.
   30        if (mp*.1d0 .gt. omp) return
   40     omp  =  mp
! calculate next iterate.
   50     call calct(bool)
        call nexth(bool)
        call calct(bool)
        if (bool) go to 60
        relstp  =  cmod(tr,ti)/cmod(sr,si)
        sr  =  sr+tr
        si  =  si+ti
   60 continue
    return
end subroutine vrshft
!}}}
!subroutine calct{{{
subroutine calct(bool)                                            calc2890
! computes  t  =  -p(s)/h(s).
! bool   - logical, set true if h(s) is essentially zero.
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
    double precision hvr,hvi,cmod
    logical bool
    n  =  nn-1
! evaluate h(s).
    call polyev(n,sr,si,hr,hi,qhr,qhi,hvr,hvi)
    bool  =  cmod(hvr,hvi) .le. are*10.0d0*cmod(hr(n),hi(n))
    if (bool) go to 10
        call cdivid(-pvr,-pvi,hvr,hvi,tr,ti)
        return
   10 tr  =  0.0d0
    ti  =  0.0d0
    return
end subroutine calct
!}}}
!subroutine nexth{{{
subroutine nexth(bool)                                            next3110
! calculates the next shifted h polynomial.
! bool   -  logical, if .true. h(s) is essentially zero
! common area
    common/global/pr,pi,hr,hi,qpr,qpi,qhr,qhi,shr,shi,
     *    sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,nn
    double precision sr,si,tr,ti,pvr,pvi,are,mre,eta,infin,
     *    pr(50),pi(50),hr(50),hi(50),qpr(50),qpi(50),qhr(50),
     *    qhi(50),shr(50),shi(50)
    double precision t1,t2
    logical bool
    n  =  nn-1
    nm1  =  n-1
    if (bool) go to 20
        do 10 j  =  2,n
             t1  =  qhr(j-1)
             t2  =  qhi(j-1)
             hr(j)  =  tr*t1-ti*t2+qpr(j)
             hi(j)  =  tr*t2+ti*t1+qpi(j)
   10     continue
        hr(1)  =  qpr(1)
        hi(1)  =  qpi(1)
        return
! if h(s) is zero replace h with qh.
   20 do 30 j  =  2,n
        hr(j)  =  qhr(j-1)
        hi(j)  =  qhi(j-1)
   30 continue
    hr(1)  =  0.0d0
    hi(1)  =  0.0d0
    return
end subroutine nexth
!}}}
!subroutine polyev{{{
subroutine polyev(nn,sr,si,pr,pi,qr,qi,pvr,pvi)                   poly3430
! evaluates a polynomial  p  at  s  by the horner recurrence
! placing the partial sums in q and the computed value in pv.
    double precision pr(nn),pi(nn),qr(nn),qi(nn),
     *    sr,si,pvr,pvi,t
    qr(1)  =  pr(1)
    qi(1)  =  pi(1)
    pvr  =  qr(1)
    pvi  =  qi(1)
    do 10 i  =  2,nn
        t  =  pvr*sr-pvi*si+pr(i)
        pvi  =  pvr*si+pvi*sr+pi(i)
        pvr  =  t
        qr(i)  =  pvr
        qi(i)  =  pvi
   10 continue
    return
end subroutine polyev
!}}}
!function errev{{{
double precision function errev(nn,qr,qi,ms,mp,are,mre)           erre3610
! bounds the error in evaluating the polynomial by the horner
! recurrence.
! qr,qi - the partial sums
! ms    -modulus of the point
! mp    -modulus of polynomial value
! are, mre -error bounds on complex addition and multiplication
    double precision qr(nn),qi(nn),ms,mp,are,mre,e,cmod
    e  =  cmod(qr(1),qi(1))*mre/(are+mre)
    do 10 i  =  1,nn
        e  =  e*ms+cmod(qr(i),qi(i))
   10 continue
    errev  =  e*(are+mre)-mp*mre
    return
    end function error    
!}}}
!function cauchy{{{
    double precision function cauchy(nn,pt,q)                         cauc3760
! cauchy computes a lower bound on the moduli of the zeros of a
! polynomial - pt is the modulus of the coefficients.
    double precision q(nn),pt(nn),x,xm,f,dx,df,
     *   dabs,dexp,dlog
    pt(nn)  =  -pt(nn)
! compute upper estimate of bound.
    n  =  nn-1
    x  =  dexp( (dlog(-pt(nn)) - dlog(pt(1)))/float(n) )
    if (pt(n).eq.0.0d0) go to 20
! if newton step at the origin is better, use it.
        xm  =  -pt(nn)/pt(n)
        if (xm.lt.x) x = xm
! chop the interval (0,x) unitl f le 0.
   20 xm  =  x*.1d0
    f  =  pt(1)
    do 30 i  =  2,nn
        f  =  f*xm+pt(i)
   30 continue
    if (f.le. 0.0d0) go to 40
        x  =  xm
        go to 20
   40 dx  =  x
! do newton iteration until x converges to two decimal places.
   50 if (dabs(dx/x) .le. .005d0) go to 70
        q(1)  =  pt(1)
        do 60 i  =  2,nn
             q(i)  =  q(i-1)*x+pt(i)
   60     continue
        f  =  q(nn)
        df  =  q(1)
        do 65 i  =  2,n
             df  =  df*x+q(i)
   65     continue
        dx  =  f/df
        x  =  x-dx
        go to 50
   70 cauchy  =  x
    return
        end function cauchy
!}}}
!function scale{{{
    double precision function scale(nn,pt,eta,infin,smalno,base)      scal4160
! returns a scale factor to multiply the coefficients of the
! polynomial. the scaling is done to avoid overflow and to avoid
! undetected underflow interfering with the convergence
! criterion.  the factor is a power of the base.
! pt - modulus of coefficients of p
! eta,infin,smalno,base - constants describing the
! floating point arithmetic.
    double precision pt(nn),eta,infin,smalno,base,hi,lo,
     *    max,min,x,sc,dsqrt,dlog
! find largest and smallest moduli of coefficients.
    hi  =  dsqrt(infin)
    lo  =  smalno/eta
    max  =  0.0d0
    min  =  infin
    do 10 i  =  1,nn
        x  =  pt(i)
        if (x .gt. max) max  =  x
        if (x .ne. 0.0d0 .and. x.lt.min) min  =  x
   10 continue
! scale only if there are very large or very small components.
    scale  =  1.0d0
    if (min .ge. lo .and. max .le. hi) return
    x  =  lo/min
    if (x .gt. 1.0d0) go to 20
        sc  =  1.0d0/(dsqrt(max)*dsqrt(min))
        go to 30
   20 sc  =  x
    if (infin/sc .gt. max) sc  =  1.0d0
   30 l  =  dlog(sc)/dlog(base) + .500
    scale  =  base**l
    return
    end
end function scale
!}}}
!subroutine cdivid{{{
subroutine cdivid(ar,ai,br,bi,cr,ci)                              cdiv4490
! complex division c  =  a/b, avoiding overflow.
    double precision ar,ai,br,bi,cr,ci,r,d,t,infin,dabs
    if (br .ne. 0.0d0  .or. bi .ne. 0.0d0) go to 10
! division by zero, c  =  infinity.
        call mcon (t,infin,t,t)
        cr  =  infin
        ci  =  infin
        return
   10 if (dabs(br) .ge. dabs(bi)) go to 20
        r  =  br/bi
        d  =  bi+r*br
        cr  =  (ar*r+ai)/d
        ci  =  (ai*r-ar)/d
        return
   20 r  =  bi/br
    d  =  br+r*bi
    cr  =  (ar+ai*r)/d
    ci  =  (ai-ar*r)/d
    return
end subroutine cdivid
!}}}
!function cmod{{{
double precision function cmod(r,i)                               cmod4700
! modulus of a complex number avoiding overflow.
    double precision r,i,ar,ai,dabs,dsqrt
    ar  =  dabs(r)
    ai  =  dabs(i)
    if (ar .ge. ai) go to 10
        cmod  =  ai*dsqrt(1.0d0+(ar/ai)**2)
        return
   10 if (ar .le. ai) go to 20
        cmod  =  ar*dsqrt(1.0d0+(ai/ar)**2)
        return
   20 cmod  =  ar*dsqrt(2.0d0)
    return
end function cmod
!}}}
!subroutine mcon{{{
subroutine mcon(eta,infiny,smalno,base)                           mcon4840
! mcon provides machine constants used in various parts of the
! program. the user may either set them directly or use the
! statements below to compute them. the meaning of the four
! constants are -
! eta     the maximum relative representation error
! which can be described as the smallest positive
! floating-point number such that 1.0d0 + eta is
! greater than 1.0d0.
! infiny    the largest floating-point number
! smalno    the smallest positive floating-point number
! base    the base of the floating-point number system used
! let t be the number of base-digits in each floating-point
! number(double precision). then eta is either .5*b**(1-t)
! or b**(1-t) depending on whether rounding or truncation
! is used.
! let m be the largest exponent and n the smallest exponent
! in the number system. then infiny is (1-base**(-t))*base**m
! and smalno is base**n.
! the values for base,t,m,n below correspond to the ibm/360.
    double precision eta,infiny,smalno,base
    integer m,n,t
    base    =  16.0d0
    t       =  14
    m       =  63
    n       =  -65
    eta     =  base**(1-t)
    infiny  =  base*(1.0d0-base**(-t))*base**(m-1)
    smalno  =  (base**(n+3))/base**3
    return
    end subroutine mcon
!}}}
end module module_common
