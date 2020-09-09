module module_common
    implicit none
    integer,parameter            :: n = 100
    integer                      :: i
    integer                      :: j
    integer                      :: k
    real(8)                      :: x1
    character(10)                :: filename 
    !  common variables for laplace
    real(8)                      :: ta
    real(8)                      :: tb
    real(8)                      :: t0
    real(8)                      :: conopt
    real(8)                      :: absf
    real(8)                      :: lval
    integer                      :: hmono
    
    contains
!subroutine lapin{{{
subroutine lapin(ti,tn,n,iman,ilapin,ikonv,&
                ns1,ns2,icon,ikor,con,h,e,ier,iout)
    real(8)             ::    ti
    real(8)             ::    tn
    integer             ::    n
    integer             ::    iman
    integer             ::    ilapin
    integer             ::    ikonv
    integer             ::    ns1
    integer             ::    ns2
    integer             ::    icon
    integer             ::    ikor
    integer             ::    con
    integer             ::    ier
    integer             ::    iout
    real(8)             ::    h(6,n)
    real(8)             ::    e(3,ns1)
    real(8)             ::    abrn
    real(8)             ::    absf
    real(8)             ::    conopt
    real(8)             ::    con1
    real(8)             ::    con2
    real(8)             ::    del
    real(8)             ::    fn
    real(8)             ::    fns1
    real(8)             ::    freal
    real(8)             ::    fimagjh
    real(8)             ::    pi
    real(8)             ::    racc
    real(8)             ::    rnsum
    real(8)             ::    rnsumk
    real(8)             ::    t
    real(8)             ::    ta
    real(8)             ::    tb
    real(8)             ::    tk
    real(8)             ::    t0
    real(8)             ::    v1
    real(8)             ::    v2
    real(8)             ::    w
    integer             ::    j3
    integer             ::    l3
    integer             ::    kor1
    integer             ::    jump
    !initialize array h
    do i = 1,n
        do j = 1,6
            h(j,i)=0.D0
        enddo
    enddo
    
    ier=0
    if(tn<ti)ier=1
    if(n<1)ier=ier+10
    if(iman==0)goto 100
    if(iman==i)goto 110
    ier=ier+100
    goto 120
    
    !parametersforiman=0
    100 if(ns1/=60)ier=ier+100000
    if(ier/=0)goto 120
    ilapin = 1
    ikonv  = 2
    icon   = 1
    ikor   = 0
    ns2    = 0
    goto 200
    iman=1
    110 if(ilapin<1.or.ilapin>2)ier=ier+1000
    if(ikonv<1.or.ikonv>2)ier=ier+10000
    if(ns1<i.or.(ns2<1.and.ikor==1))ier=ier+100000
    if(icon<0.or.icon>2.or.(icon==2.and.ilapin==2))then
    ier=ier+1000000
    endif
    if(ikor<0.or.ikor>1)ier=ier+10000000
    if(icon==0.and.con<=0.D0)ier=ier+100000000
    if(ier==0)goto 200
    !120 call  error(ier,ti,tn,n,iman,ilapin,ikonv,ns1,&
    !ns2,icon,ikor,con,iout)
    return
    200 pi = 4.D0*datan(1.D0)
    con1   = 20.D0
    con2   = con1 - 2.D0
    absf   = 0.D0
    j3     = (3-icon)/2+n*icon/2
    ta     = ti
    tb     = tn
    
    do 830 l3 = 1,j3
    lval   = l3
    hmono  = 0
    kor1   = ikor
    !computation of the optimal parameters
    210 if(icon-l)215,230,220
    215 jump=0
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,ikor,con,h,e,jump)
    goto 830
    220 ta=ti+float(l3)*(tn-ti)/float(n+1)
    tb=ta
    230 nh=n/2
    t=ta+(tb-ta)*float(nh)/float(n+1)
    tk=float(2-ilapin)*t+float(ilapin-1)*tb
    !computation of the truncation error(rnsum)
    to   = t
    con  = con1
    jump = 1
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,&
                icon,ikor,con,h,e,jump)
    240 fn = h(1,l3)
    fns1   = e(1,ns1)
    con    = con2
    jump   = 2
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,&
                icon,ikor,con,h,e,jump)
    250 if(fn/=h(1,l3).and.fns1/=e(1,ns1))goto 255
    conopt = con
    absf   = 0.D0
    goto 320
    255 rnsum=tk*(fn-h(1,l3))/(dexp(con1) - dexp(con2))
    if(ilapin.eq.2)goto 260
    !computation of the acceleration factor(del)
    racc = t*(fns1-e(i,ns1))/(dexp(con1) - dexp(con2))
    del  = rnsum/racc
    260 if(ikor == 1)goto 280
    !optimal parameters(method A)
    to   = 2.D0*tk+t
    con  = con1/4.D0
    jump = 3
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,ikor,con,h,e,jump)
    270 fn = h(i,l3)
    conopt = -tk/(2.D0*tk+t)*dlog(dabs(rnsum/(tk*fn)))
    goto 310
    !optimal parameters for the Korrektur method(method A)
    280 to = 4.D0*tk+t
    con    = con1/4.D0
    jump   = 4
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,ikor,con,h,e,jump)
    290 fn = h(1,l3)
    to     = 8.D0*tk+t
    jump   = 5
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,ikor,con,h,e,jump)
    300 fn = fn-h(i,l3)
    conopt = -tk/(4.D0*tk+t)*dlog(dabs(rnsum/(tk*fn)))
    !optimal parameters(method B)
    310 if(ilapin.eq.1) goto 315
    absf=dabs(dexp(conopt)*rnsum*2.D0/tk)
    goto 320
    315 v1 = con1/t
    v2     = conopt/t
    w      = pi*float(ns1)/t
    !call fft(v2,w,freal,fimag)
    rnsumk =rnsum*freal
    !call fft(v1,w,freal,fimag)
    rnsumk = rnsumk/freal
    abrn   = (rnsumk-rnsum)/(v2-v1)
    conopt = - dlog(dabs((abrn/t+rnsumk)/(t*float(ikor*2+2)*fn)))/&
            float(3+2*ikor)
    vl     = v2
    v2     = conopt/t
    !call fft(v2,w,freal,fimag)
    rnsumk = rnsumk*freal
    !call fft(vi,w,freal,fimag)
    rnsumk = rnsumk/freal
    absf   = dexp(conopt)/t*dabs(rnsumk)+&
             dabs(dexp(-2.D0*conopt)*fn*float(ikor-i)&
             +dexp(-4.D0*conopt)*fn*float(ikor))
    320 if(conopt<=0.D0)conopt = 1.D0
    jump = 6
    call lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,ikor,con,h,e,jump)
    830 continue
    return
end subroutine lapin
!}}}
!subroutine lapin2{{{
subroutine lapin2(ti,tn,n,ilapin,ikonv,ns1,ns2,icon,&
           ikor,con,h,e,jump)
    real(8)         ::   ti
    real(8)         ::   tn
    integer         ::   n
    integer         ::   ilapin
    integer         ::   ikonv
    integer         ::   ns1
    integer         ::   ns2
    integer         ::   icon
    integer         ::   ikor
    integer         ::   con
    integer         ::   jump
    real(8)         ::   h(6,n)
    real(8)         ::   e(3,ns1)
    real(8)         ::   e3(3)
    real(8)         ::   a
    real(8)         ::   absf
    real(8)         ::   b
    real(8)         ::   conopt
    real(8)         ::   delta
    real(8)         ::   divi
    real(8)         ::   e1
    real(8)         ::   e2
    real(8)         ::   eins
    real(8)         ::   faktor
    real(8)         ::   fimag(n)
    real(8)         ::   freal(n)
    real(8)         ::   pi
    real(8)         ::   pit
    real(8)         ::   pite
    real(8)         ::   ral
    real(8)         ::   suim
    real(8)         ::   sure
    real(8)         ::   ta
    real(8)         ::   tb
    real(8)         ::   te
    real(8)         ::   tl
    real(8)         ::   tm
    real(8)         ::   tm1
    real(8)         ::   tt
    real(8)         ::   t0
    real(8)         ::   x1
    real(8)         ::   x2
    real(8)         ::   x3
    real(8)         ::   v
    real(8)         ::   w
    real(8)         ::   y1
    real(8)         ::   y2
    real(8)         ::   y3
    integer         ::   richt
    integer         ::   richta
    integer         ::   hmono
    integer         ::   nh
    common  /clapin/ ta,tb,t0,conopt,absf,l3,hmono

    pi  = 4.D0 *datan(1.D0)
    if(jump == 0)goto 360
    if(jump <  6)goto 370 
    to   = ta
    tt   = tb
    i1   = l3
    j1   = (2-icon)*n+(icon-i)*l3
    con  = conopt
    kor1 = ikor
    goto 380
    360 to = ti
    tt     = tn
    i1     = 1
    j1     = n
    kori   = ikor
    jump   = 6
    goto  380
    370 kor1 = 0
    tt       = to
    ii       = l3
    j1       = l3
    380 delta= (tt-t0)/float(n+1)
    nsum     = nsl
    !computation of the t-values from t0,tt
    do k2=1,2
        if(ilapin.eq.1) goto 420
        te = float(k2)*tt
        v  =con/te
        !call fft(v,0.D0,freal,fimag)
        ral = - 0.5D0*freal
        pite=pi/te
    405 do l=i,nsum
            w      = float(l-1)*pite
            !call fft(v,w,freal,fimag)
            e(2,l) = freal
            e(3,l) = fimag
        enddo  ! l
    420 do ki = i1,j1
            tl=t0+float(k1)*delta
            if(ilapin.eq.2)goto 440
            if(k2.eq.2)tl=3.D0*tl
            v      = con/tl
            faktor = dexp(v*tl)/tl
            !call fft(v,0.D0,freal,fimag)
            ral  = - 0.5D0*freal
            pit  = pi/tl
            eins = 1.D0
            sure = 0.D0
!method of Durbin
      425   do l=i,nsum
            w   =float(l-1)*pit
            !call fft(v,w,freal,fimag)
            sure   = sure+freal*eins
            eins   = -eins
            e(i,l) = faktor*(ral+sure)
            enddo
            goto 460
            
       440  if(k2.eq.2) tl=tl+te
            faktor =dexp(v*tl)/te
            sure=0.D0 
            suim=0.D0
            do l = i,nsum
                w=float(l-i)*pite
                sure   = sure+e(2,l)*dcos(w*tl)
                suim   = suim+e(3,l)*dsin(w*tl)
                e(i,l) = faktor*(ral+sure-suim)
            enddo
            !search for stationary values
       460  nmax   = nsum*2/3
            monoto = 0
            k      = 0
            richta = dsign(1.5D0,(e(i,nsum)-e(i,nsum-1)))
            do l=i,nmax
                j     = nsum-l
                richt = dsign(1.5D0,(e(i,j)-e(1,j-1)))
                if(richt.eq.richta)cycle
                k      = k+i
                e3(k)  = e(1,j)
                richta = richt
                if(k.eq.3) goto 510
            enddo
            if(k.eq.0)goto 700
            h(k2,k1) =  e(1,nsum)
            if(k2.eq.i) h(4,k1)=0
            goto 790
            510 ke = 2
            if((e3(ke)-e(1,j))*float(richta)>0.D0) goto 560
            jmin = nsum/3
            jmax = j-1
            do jj=jmin,jmax
                j=j-1
                richt=dsign(1.5D0,(e(i,j)-e(1,j-1)))
                if(richt==richta)cycle
                richta = richt
                ke     = 3-ke
                if((e3(ke)-e(i,j))*float(richta)>0.D0)goto 560
            enddo
      550   monoto=1
            if(ikonv==2)goto 630
            !minimum-max1mum method(minimax)
      560   h(k2,ki)=(e3(1)+e3(3))/4.D0+e3(2)/2.D0
            if(k2.eq.i) h(4,ki)=1
            goto 790
            !epsilon algorithm (epal)
      630   k = 0
            nsummi = nsum-1
            e2     = e(1,1)
            do l = 1,nsumm1
                e1  = e(1,1)
                tm  = 0.D0
                lp1 = l+1
                do m = 1,l
                    mm   = lp1-m
                    tmi  = e(i,mm)
                    divi = e(1,mm+1)-e(1,mm)
                    if(dabs(divi)>1.D-20) goto 640
                    k=l
                    goto 670
                640 e(i,mm)=tm+1.D0/divi
                    tm = tm1
                enddo
                e2=e1
            enddo
        670 if(dabs(e1)>dabs(e2))el=e2
            if(dabs(e(1,1))>dabs(e1)) e(1,1) = e1
            h(k2,k1) = e(1,1)
            if(k2==1) h(4,k1)=k+2
            goto 790
            !curve fitting(cfm)
       700  xl = float(nsum)- 2.D0 
            x2 = float(nsum)- 1.D0
            x3 = float(nsum)- 0.D0
            y1 = e(1,nsum-2)
            y2 = e(1,nsum-1)
            y3 = e(1,nsum)
            b  = ((y3-y1)*x3*x3*(xl+x2)/(xl-x3) &
                 -(y2_y1)*x2,x2x(x1+x3)/(x1_x2))/(x3_x2)
            a  =((y2-y1)-b*(x1-x2)/(xl*x2))*(x2*x2*x1*x1)/(x1*x1-x2*x2)
            h(k2,k1) = y1-(a/x1+b)/x1
            monoto   = 1
            if(k2.eq.1)h(4,k1)=3
        790 h(3,k1)= con
            h(5,k1)= absf
            hmono  = hmono+k2*monoto*10**(6-jump)
            if(jump<6)cycle
            h(6,k1)=float((2-k2)*hmono)+ float(k2-1)*&
                    (float(2*monoto)+h(6,k1))
            hmono = hmono/(10*10)
        enddo
        if(kor1.eq.0)return
        if(k2==2)goto 810
        nsum=ns2
    enddo    ! k2
    !korrekturmethod
810 faktor = - dexp(- 2.D0*con)
    do  k  = i1,j1
        h(1,k) = h(1,k)+faktor*h(2,k)
    end do
    return
end subroutine lapin2
!}}}
!  fft
!subroutine fft{{{
subroutine fft(x_re,x_im,y_re,y_im)
    real(8),intent(in)    :: x_re(n)
    real(8),intent(in)    :: x_im(n)
    real(8),intent(out)   :: y_re(n)
    real(8),intent(out)   :: y_im(n)
    integer               :: a(n)    ! original order
    integer               :: b(l)    ! 2 power 
    integer               :: c(l)    ! n/b(l) 
    integer               :: ia 
    integer               :: ib 
    integer               :: k 
    integer               :: num 
    integer               :: cy

    b(1) = 2
    do i = 1,l-1
       b(i+1) = b(i)*2 
    end do
    do i = 1,l-1
       c(i) = b(l - i)
    end do
    c(n) = 1
    ! calculate a(n)
    call inverse1(a,n,l)

    !  fft
    do i = 1,l
        num = b(i)
        do j = 1,c(i) 
            ia = (j-1)*num
            ib = (j-1)*num + num/2
            call butterfly(y_re(ia:ib),y_im(ia:ib),num)
        end do
    end do
end subroutine fft
!}}}
!  butterfly for butterfly algorithm
!subroutine butterfly{{{
subroutine butterfly(y_re,y_im,n)
    integer,intent(in)       :: n
    real(8),intent(inout)    :: y_re(n)
    real(8),intent(inout)    :: y_im(n)
    real(8)                  :: x_re(n)
    real(8)                  :: x_im(n)
    real(8)                  :: w(n)
    
    x_re = y_re
    x_im = y_im
    do i = 1,n/2
        j = i + n/2 
        w(i) = cos(2*pi*(i-1)/dble(n))
        w(j) = sin(2*pi*(i-1)/dble(n))
        y_re(i) = x_re(i) + w(i)*x_re(j) + w(j)*x_im(j)
        y_im(i) = x_im(i) + w(i)*x_im(j) - w(j)*x_re(j)
        y_re(j) = x_re(i) - w(i)*x_re(j) - w(j)*x_im(j)
        y_im(j) = x_im(i) - w(i)*x_im(j) + w(j)*x_re(j)
    end do
end subroutine butterfly
!}}}
!subroutine error{{{
subroutine error(ier,ti,tn,n,iman,ilapin,ikonvpns1,&
                     ns2,icon,ikor, con,iout)
    double precision   ti,tn,con
    write(iout,*)ier,ti,tn,n,iman,ilapin,ikonv,ns1,ns2,icon,ikor,con
end subroutine error
!}}}
end module module_common
