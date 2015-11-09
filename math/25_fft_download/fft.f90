!subroutine cfftb1{{{
subroutine cfftb1(n,c,ch,wa,ifac)
    dimension       ch(1)      ,c(1)       ,wa(1)      ,ifac(1)
    nf = ifac(2)
    na = 0
    l1 = 1
    iw = 1
    do 116 k1=1,nf
       ip = ifac(k1+2)
       l2 = ip*l1
       ido = n/l2
       idot = ido+ido
       idl1 = idot*l1
       if (ip .ne. 4) go to 103
       ix2 = iw+idot
       ix3 = ix2+idot
       if (na .ne. 0) go to 101
       call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
       go to 102
  101    call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
       go to 115
  103    if (ip .ne. 2) go to 106
       if (na .ne. 0) go to 104
       call passb2 (idot,l1,c,ch,wa(iw))
       go to 105
  104    call passb2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
       go to 115
  106    if (ip .ne. 3) go to 109
       ix2 = iw+idot
       if (na .ne. 0) go to 107
       call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
       go to 108
  107    call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
       go to 115
  109    if (ip .ne. 5) go to 112
       ix2 = iw+idot
       ix3 = ix2+idot
       ix4 = ix3+idot
       if (na .ne. 0) go to 110
       call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
       go to 111
  110    call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
       go to 115
  112    if (na .ne. 0) go to 113
       call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
       go to 114
  113    call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
       iw = iw+(ip-1)*idot
  116 continue
    if (na .eq. 0) return
    n2 = n+n
    do 117 i=1,n2
       c(i) = ch(i)
  117 continue
    return
    end
!}}}
!subroutine cfftb {{{
subroutine cfftb (n,c,wsave)
    dimension       c(1)       ,wsave(1)
    if (n .eq. 1) return
    iw1 = n+n+1
    iw2 = iw1+n+n
    call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
    return
    end
!}}}
!subroutine cfftf1 {{{
subroutine cfftf1 (n,c,ch,wa,ifac)
    dimension       ch(1)      ,c(1)       ,wa(1)      ,ifac(1)
    nf = ifac(2)
    na = 0
    l1 = 1
    iw = 1
    do 116 k1=1,nf
       ip = ifac(k1+2)
       l2 = ip*l1
       ido = n/l2
       idot = ido+ido
       idl1 = idot*l1
       if (ip .ne. 4) go to 103
       ix2 = iw+idot
       ix3 = ix2+idot
       if (na .ne. 0) go to 101
       call passf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
       go to 102
  101    call passf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
       go to 115
  103    if (ip .ne. 2) go to 106
       if (na .ne. 0) go to 104
       call passf2 (idot,l1,c,ch,wa(iw))
       go to 105
  104    call passf2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
       go to 115
  106    if (ip .ne. 3) go to 109
       ix2 = iw+idot
       if (na .ne. 0) go to 107
       call passf3 (idot,l1,c,ch,wa(iw),wa(ix2))
       go to 108
  107    call passf3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
       go to 115
  109    if (ip .ne. 5) go to 112
       ix2 = iw+idot
       ix3 = ix2+idot
       ix4 = ix3+idot
       if (na .ne. 0) go to 110
       call passf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
       go to 111
  110    call passf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
       go to 115
  112    if (na .ne. 0) go to 113
       call passf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
       go to 114
  113    call passf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
  115    l1 = l2
       iw = iw+(ip-1)*idot
  116 continue
    if (na .eq. 0) return
    n2 = n+n
    do 117 i=1,n2
       c(i) = ch(i)
  117 continue
    return
    end
!}}}
!subroutine cfftf {{{
subroutine cfftf (n,c,wsave)
    dimension       c(1)       ,wsave(1)
    if (n .eq. 1) return
    iw1 = n+n+1
    iw2 = iw1+n+n
    call cfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
    return
    end
!}}}
!subroutine cffti1 {{{
subroutine cffti1 (n,wa,ifac)
    dimension       wa(1)      ,ifac(1)    ,ntryh(4)
    data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
    nl = n
    nf = 0
    j = 0
  101 j = j+1
    if (j-4) 102,102,103
  102 ntry = ntryh(j)
    go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
    nr = nl-ntry*nq
    if (nr) 101,105,101
  105 nf = nf+1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry .ne. 2) go to 107
    if (nf .eq. 1) go to 107
    do 106 i=2,nf
       ib = nf-i+2
       ifac(ib+2) = ifac(ib+1)
  106 continue
    ifac(3) = 2
  107 if (nl .ne. 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    tpi = 6.28318530717959
    argh = tpi/float(n)
    i = 2
    l1 = 1
    do 110 k1=1,nf
       ip = ifac(k1+2)
       ld = 0
       l2 = l1*ip
       ido = n/l2
       idot = ido+ido+2
       ipm = ip-1
       do 109 j=1,ipm
          i1 = i
          wa(i-1) = 1.
          wa(i) = 0.
          ld = ld+l1
          fi = 0.
          argld = float(ld)*argh
          do 108 ii=4,idot,2
             i = i+2
             fi = fi+1.
             arg = fi*argld
             wa(i-1) = cos(arg)
             wa(i) = sin(arg)
  108     continue
          if (ip .le. 5) go to 109
          wa(i1-1) = wa(i-1)
          wa(i1) = wa(i)
  109    continue
       l1 = l2
  110 continue
    return
    end
!}}}
!subroutine cffti {{{
subroutine cffti (n,wsave)
    dimension       wsave(1)
    if (n .eq. 1) return
    iw1 = n+n+1
    iw2 = iw1+n+n
    call cffti1 (n,wsave(iw1),wsave(iw2))
    return
    end
!}}}
!subroutine cosqb1 {{{
subroutine cosqb1 (n,x,w,xh)
    dimension       x(1)       ,w(1)       ,xh(1)
    ns2 = (n+1)/2
    np2 = n+2
    do 101 i=3,n,2
       xim1 = x(i-1)+x(i)
       x(i) = x(i)-x(i-1)
       x(i-1) = xim1
  101 continue
    x(1) = x(1)+x(1)
    modn = mod(n,2)
    if (modn .eq. 0) x(n) = x(n)+x(n)
    call rfftb (n,x,xh)
    do 102 k=2,ns2
       kc = np2-k
       xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
       xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
    if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
    do 103 k=2,ns2
       kc = np2-k
       x(k) = xh(k)+xh(kc)
       x(kc) = xh(k)-xh(kc)
  103 continue
    x(1) = x(1)+x(1)
    return
    end
!}}}
!subroutine cosqb {{{
subroutine cosqb (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    data tsqrt2 /2.82842712474619/
    if (n-2) 101,102,103
  101 x(1) = 4.*x(1)
    return
  102 x1 = 4.*(x(1)+x(2))
    x(2) = tsqrt2*(x(1)-x(2))
    x(1) = x1
    return
  103 call cosqb1 (n,x,wsave,wsave(n+1))
    return
    end
!}}}
!subroutine cosqf1 {{{
subroutine cosqf1 (n,x,w,xh)
    dimension       x(1)       ,w(1)       ,xh(1)
    ns2 = (n+1)/2
    np2 = n+2
    do 101 k=2,ns2
       kc = np2-k
       xh(k) = x(k)+x(kc)
       xh(kc) = x(k)-x(kc)
  101 continue
    modn = mod(n,2)
    if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
    do 102 k=2,ns2
       kc = np2-k
       x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
       x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
    if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
    call rfftf (n,x,xh)
    do 103 i=3,n,2
       xim1 = x(i-1)-x(i)
       x(i) = x(i-1)+x(i)
       x(i-1) = xim1
  103 continue
    return
    end
!}}}
!subroutine cosqf {{{
subroutine cosqf (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    data sqrt2 /1.4142135623731/
    if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
    x(2) = x(1)-tsqx
    x(1) = x(1)+tsqx
  102 return
  103 call cosqf1 (n,x,wsave,wsave(n+1))
    return
    end
!}}}
!subroutine cosqi {{{
subroutine cosqi (n,wsave)
    dimension       wsave(1)
    data pih /1.57079632679491/
    dt = pih/float(n)
    fk = 0.
    do 101 k=1,n
       fk = fk+1.
       wsave(k) = cos(fk*dt)
  101 continue
    call rffti (n,wsave(n+1))
    return
    end
!}}}
!subroutine cost {{{
subroutine cost (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    nm1 = n-1
    np1 = n+1
    ns2 = n/2
    if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
    x(2) = x(1)-x(2)
    x(1) = x1h
    return
  102 if (n .gt. 3) go to 103
    x1p3 = x(1)+x(3)
    tx2 = x(2)+x(2)
    x(2) = x(1)-x(3)
    x(1) = x1p3+tx2
    x(3) = x1p3-tx2
    return
  103 c1 = x(1)-x(n)
    x(1) = x(1)+x(n)
    do 104 k=2,ns2
       kc = np1-k
       t1 = x(k)+x(kc)
       t2 = x(k)-x(kc)
       c1 = c1+wsave(kc)*t2
       t2 = wsave(k)*t2
       x(k) = t1-t2
       x(kc) = t1+t2
  104 continue
    modn = mod(n,2)
    if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
    call rfftf (nm1,x,wsave(n+1))
    xim2 = x(2)
    x(2) = c1
    do 105 i=4,n,2
       xi = x(i)
       x(i) = x(i-2)-x(i-1)
       x(i-1) = xim2
       xim2 = xi
  105 continue
    if (modn .ne. 0) x(n) = xim2
  106 return
    end
!}}}
!subroutine costi {{{
subroutine costi (n,wsave)
    dimension       wsave(1)
    data pi /3.14159265358979/
    if (n .le. 3) return
    nm1 = n-1
    np1 = n+1
    ns2 = n/2
    dt = pi/float(nm1)
    fk = 0.
    do 101 k=2,ns2
       kc = np1-k
       fk = fk+1.
       wsave(k) = 2.*sin(fk*dt)
       wsave(kc) = 2.*cos(fk*dt)
  101 continue
    call rffti (nm1,wsave(n+1))
    return
    end
!}}}
!subroutine ezfft1 {{{
subroutine ezfft1 (n,wa,ifac)
    dimension       wa(1)      ,ifac(1)    ,ntryh(4)
    data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/ ,tpi/6.28318530717959/
    nl = n
    nf = 0
    j = 0
  101 j = j+1
    if (j-4) 102,102,103
  102 ntry = ntryh(j)
    go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
    nr = nl-ntry*nq
    if (nr) 101,105,101
  105 nf = nf+1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry .ne. 2) go to 107
    if (nf .eq. 1) go to 107
    do 106 i=2,nf
       ib = nf-i+2
       ifac(ib+2) = ifac(ib+1)
  106 continue
    ifac(3) = 2
  107 if (nl .ne. 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    argh = tpi/float(n)
    is = 0
    nfm1 = nf-1
    l1 = 1
    if (nfm1 .eq. 0) return
    do 111 k1=1,nfm1
       ip = ifac(k1+2)
       l2 = l1*ip
       ido = n/l2
       ipm = ip-1
       arg1 = float(l1)*argh
       ch1 = 1.
       sh1 = 0.
       dch1 = cos(arg1)
       dsh1 = sin(arg1)
       do 110 j=1,ipm
          ch1h = dch1*ch1-dsh1*sh1
          sh1 = dch1*sh1+dsh1*ch1
          ch1 = ch1h
          i = is+2
          wa(i-1) = ch1
          wa(i) = sh1
          if (ido .lt. 5) go to 109
          do 108 ii=5,ido,2
             i = i+2
             wa(i-1) = ch1*wa(i-3)-sh1*wa(i-2)
             wa(i) = ch1*wa(i-2)+sh1*wa(i-3)
  108     continue
  109     is = is+ido
  110    continue
       l1 = l2
  111 continue
    return
    end
!}}}
!subroutine ezfftb {{{
subroutine ezfftb (n,r,azero,a,b,wsave)
    dimension       r(1)       ,a(1)       ,b(1)       ,wsave(1)
    if (n-2) 101,102,103
  101 r(1) = azero
    return
  102 r(1) = azero+a(1)
    r(2) = azero-a(1)
    return
  103 ns2 = (n-1)/2
    do 104 i=1,ns2
       r(2*i) = .5*a(i)
       r(2*i+1) = -.5*b(i)
  104 continue
    r(1) = azero
    if (mod(n,2) .eq. 0) r(n) = a(ns2+1)
    call rfftb (n,r,wsave(n+1))
    return
    end
!}}}
!subroutine ezfftf {{{
subroutine ezfftf (n,r,azero,a,b,wsave)
!
!                     version 3  june 1979
!
    dimension       r(1)       ,a(1)       ,b(1)       ,wsave(1)
    if (n-2) 101,102,103
  101 azero = r(1)
    return
  102 azero = .5*(r(1)+r(2))
    a(1) = .5*(r(1)-r(2))
    return
  103 do 104 i=1,n
       wsave(i) = r(i)
  104 continue
    call rfftf (n,wsave,wsave(n+1))
    cf = 2./float(n)
    cfm = -cf
    azero = .5*cf*wsave(1)
    ns2 = (n+1)/2
    ns2m = ns2-1
    do 105 i=1,ns2m
       a(i) = cf*wsave(2*i)
       b(i) = cfm*wsave(2*i+1)
  105 continue
    if (mod(n,2) .eq. 1) return
    a(ns2) = .5*cf*wsave(n)
    b(ns2) = 0.
    return
    end
!}}}
!subroutine ezffti {{{
subroutine ezffti (n,wsave)
    dimension       wsave(1)
    if (n .eq. 1) return
    call ezfft1 (n,wsave(2*n+1),wsave(3*n+1))
    return
    end
!}}}
!subroutine passb2 {{{
subroutine passb2 (ido,l1,cc,ch,wa1)
    dimension       cc(ido,2,l1),ch(ido,l1,2), wa1(1)
    if (ido .gt. 2) go to 102
    do 101 k=1,l1
       ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
       ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
       ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
       ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
          tr2 = cc(i-1,1,k)-cc(i-1,2,k)
          ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
          ti2 = cc(i,1,k)-cc(i,2,k)
          ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
          ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passb3 {{{
subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
    dimension cc(ido,3,l1),ch(ido,l1,3),wa1(1),wa2(1)
    data taur,taui/-.5,.866025403784439/
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       tr2 = cc(1,2,k)+cc(1,3,k)
       cr2 = cc(1,1,k)+taur*tr2
       ch(1,k,1) = cc(1,1,k)+tr2
       ti2 = cc(2,2,k)+cc(2,3,k)
       ci2 = cc(2,1,k)+taur*ti2
       ch(2,k,1) = cc(2,1,k)+ti2
       cr3 = taui*(cc(1,2,k)-cc(1,3,k))
       ci3 = taui*(cc(2,2,k)-cc(2,3,k))
       ch(1,k,2) = cr2-ci3
       ch(1,k,3) = cr2+ci3
       ch(2,k,2) = ci2+cr3
       ch(2,k,3) = ci2-cr3
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          tr2 = cc(i-1,2,k)+cc(i-1,3,k)
          cr2 = cc(i-1,1,k)+taur*tr2
          ch(i-1,k,1) = cc(i-1,1,k)+tr2
          ti2 = cc(i,2,k)+cc(i,3,k)
          ci2 = cc(i,1,k)+taur*ti2
          ch(i,k,1) = cc(i,1,k)+ti2
          cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
          ci3 = taui*(cc(i,2,k)-cc(i,3,k))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
          ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
          ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
          ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passb4 {{{
subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
    dimension cc(ido,4,l1),ch(ido,l1,4),wa1(1),wa2(1),wa3(1)
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       ti1 = cc(2,1,k)-cc(2,3,k)
       ti2 = cc(2,1,k)+cc(2,3,k)
       tr4 = cc(2,4,k)-cc(2,2,k)
       ti3 = cc(2,2,k)+cc(2,4,k)
       tr1 = cc(1,1,k)-cc(1,3,k)
       tr2 = cc(1,1,k)+cc(1,3,k)
       ti4 = cc(1,2,k)-cc(1,4,k)
       tr3 = cc(1,2,k)+cc(1,4,k)
       ch(1,k,1) = tr2+tr3
       ch(1,k,3) = tr2-tr3
       ch(2,k,1) = ti2+ti3
       ch(2,k,3) = ti2-ti3
       ch(1,k,2) = tr1+tr4
       ch(1,k,4) = tr1-tr4
       ch(2,k,2) = ti1+ti4
       ch(2,k,4) = ti1-ti4
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ti1 = cc(i,1,k)-cc(i,3,k)
          ti2 = cc(i,1,k)+cc(i,3,k)
          ti3 = cc(i,2,k)+cc(i,4,k)
          tr4 = cc(i,4,k)-cc(i,2,k)
          tr1 = cc(i-1,1,k)-cc(i-1,3,k)
          tr2 = cc(i-1,1,k)+cc(i-1,3,k)
          ti4 = cc(i-1,2,k)-cc(i-1,4,k)
          tr3 = cc(i-1,2,k)+cc(i-1,4,k)
          ch(i-1,k,1) = tr2+tr3
          cr3 = tr2-tr3
          ch(i,k,1) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
          ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
          ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
          ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
          ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
          ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passb5 {{{
subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
    dimension cc(ido,5,l1),ch(ido,l1,5),wa1(1),wa2(1),wa3(1),wa4(1)
    data tr11,ti11,tr12,ti12 /0.309016994374947,0.951056516295154,&
    0.190983005625053, 0.587785252292473/
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       ti5 = cc(2,2,k)-cc(2,5,k)
       ti2 = cc(2,2,k)+cc(2,5,k)
       ti4 = cc(2,3,k)-cc(2,4,k)
       ti3 = cc(2,3,k)+cc(2,4,k)
       tr5 = cc(1,2,k)-cc(1,5,k)
       tr2 = cc(1,2,k)+cc(1,5,k)
       tr4 = cc(1,3,k)-cc(1,4,k)
       tr3 = cc(1,3,k)+cc(1,4,k)
       ch(1,k,1) = cc(1,1,k)+tr2+tr3
       ch(2,k,1) = cc(2,1,k)+ti2+ti3
       cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
       ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
       cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
       ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
       cr5 = ti11*tr5+ti12*tr4
       ci5 = ti11*ti5+ti12*ti4
       cr4 = ti12*tr5-ti11*tr4
       ci4 = ti12*ti5-ti11*ti4
       ch(1,k,2) = cr2-ci5
       ch(1,k,5) = cr2+ci5
       ch(2,k,2) = ci2+cr5
       ch(2,k,3) = ci3+cr4
       ch(1,k,3) = cr3-ci4
       ch(1,k,4) = cr3+ci4
       ch(2,k,4) = ci3-cr4
       ch(2,k,5) = ci2-cr5
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ti5 = cc(i,2,k)-cc(i,5,k)
          ti2 = cc(i,2,k)+cc(i,5,k)
          ti4 = cc(i,3,k)-cc(i,4,k)
          ti3 = cc(i,3,k)+cc(i,4,k)
          tr5 = cc(i-1,2,k)-cc(i-1,5,k)
          tr2 = cc(i-1,2,k)+cc(i-1,5,k)
          tr4 = cc(i-1,3,k)-cc(i-1,4,k)
          tr3 = cc(i-1,3,k)+cc(i-1,4,k)
          ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
          ch(i,k,1) = cc(i,1,k)+ti2+ti3
          cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
          ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
          cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
          ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
          ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
          ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
          ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
          ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
          ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
          ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
          ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passb {{{
subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
    dimension ch(ido,l1,ip),cc(ido,ip,l1),c1(ido,l1,ip),wa(1),c2(idl1,ip),ch2(idl1,ip)
    idot = ido/2
    nt = ip*idl1
    ipp2 = ip+2
    ipph = (ip+1)/2
    idp = ip*ido
    if (ido .lt. l1) go to 106
    do 103 j=2,ipph
       jc = ipp2-j
       do 102 k=1,l1
          do 101 i=1,ido
             ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
             ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101     continue
  102    continue
  103 continue
    do 105 k=1,l1
       do 104 i=1,ido
          ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
    go to 112
  106 do 109 j=2,ipph
       jc = ipp2-j
       do 108 i=1,ido
          do 107 k=1,l1
             ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
             ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107     continue
  108    continue
  109 continue
    do 111 i=1,ido
       do 110 k=1,l1
          ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
    inc = 0
    do 116 l=2,ipph
       lc = ipp2-l
       idl = idl+ido
       do 113 ik=1,idl1
          c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
          c2(ik,lc) = wa(idl)*ch2(ik,ip)
  113    continue
       idlj = idl
       inc = inc+ido
       do 115 j=3,ipph
          jc = ipp2-j
          idlj = idlj+inc
          if (idlj .gt. idp) idlj = idlj-idp
          war = wa(idlj-1)
          wai = wa(idlj)
          do 114 ik=1,idl1
             c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
             c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
  114     continue
  115    continue
  116 continue
    do 118 j=2,ipph
       do 117 ik=1,idl1
          ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
    do 120 j=2,ipph
       jc = ipp2-j
       do 119 ik=2,idl1,2
          ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
          ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
          ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
          ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
    nac = 1
    if (ido .eq. 2) return
    nac = 0
    do 121 ik=1,idl1
       c2(ik,1) = ch2(ik,1)
  121 continue
    do 123 j=2,ip
       do 122 k=1,l1
          c1(1,k,j) = ch(1,k,j)
          c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
    if (idot .gt. l1) go to 127
    idij = 0
    do 126 j=2,ip
       idij = idij+2
       do 125 i=4,ido,2
          idij = idij+2
          do 124 k=1,l1
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  124     continue
  125    continue
  126 continue
    return
  127 idj = 2-ido
    do 130 j=2,ip
       idj = idj+ido
       do 129 k=1,l1
          idij = idj
          do 128 i=4,ido,2
             idij = idij+2
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  128     continue
  129    continue
  130 continue
    return
    end
!}}}
!subroutine passf2 {{{
subroutine passf2 (ido,l1,cc,ch,wa1)
    dimension cc(ido,2,l1),ch(ido,l1,2),wa1(1)
    if (ido .gt. 2) go to 102
    do 101 k=1,l1
       ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
       ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
       ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
       ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
          tr2 = cc(i-1,1,k)-cc(i-1,2,k)
          ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
          ti2 = cc(i,1,k)-cc(i,2,k)
          ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
          ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passf3 {{{
subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
    dimension cc(ido,3,l1),ch(ido,l1,3),wa1(1),wa2(1)
    data taur,taui /- 0.5,- 0.866025403784439/
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       tr2 = cc(1,2,k)+cc(1,3,k)
       cr2 = cc(1,1,k)+taur*tr2
       ch(1,k,1) = cc(1,1,k)+tr2
       ti2 = cc(2,2,k)+cc(2,3,k)
       ci2 = cc(2,1,k)+taur*ti2
       ch(2,k,1) = cc(2,1,k)+ti2
       cr3 = taui*(cc(1,2,k)-cc(1,3,k))
       ci3 = taui*(cc(2,2,k)-cc(2,3,k))
       ch(1,k,2) = cr2-ci3
       ch(1,k,3) = cr2+ci3
       ch(2,k,2) = ci2+cr3
       ch(2,k,3) = ci2-cr3
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          tr2 = cc(i-1,2,k)+cc(i-1,3,k)
          cr2 = cc(i-1,1,k)+taur*tr2
          ch(i-1,k,1) = cc(i-1,1,k)+tr2
          ti2 = cc(i,2,k)+cc(i,3,k)
          ci2 = cc(i,1,k)+taur*ti2
          ch(i,k,1) = cc(i,1,k)+ti2
          cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
          ci3 = taui*(cc(i,2,k)-cc(i,3,k))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
          ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
          ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
          ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passf4 {{{
subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
    dimension cc(ido,4,l1),ch(ido,l1,4),wa1(1),wa2(1),wa3(1)
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       ti1 = cc(2,1,k)-cc(2,3,k)
       ti2 = cc(2,1,k)+cc(2,3,k)
       tr4 = cc(2,2,k)-cc(2,4,k)
       ti3 = cc(2,2,k)+cc(2,4,k)
       tr1 = cc(1,1,k)-cc(1,3,k)
       tr2 = cc(1,1,k)+cc(1,3,k)
       ti4 = cc(1,4,k)-cc(1,2,k)
       tr3 = cc(1,2,k)+cc(1,4,k)
       ch(1,k,1) = tr2+tr3
       ch(1,k,3) = tr2-tr3
       ch(2,k,1) = ti2+ti3
       ch(2,k,3) = ti2-ti3
       ch(1,k,2) = tr1+tr4
       ch(1,k,4) = tr1-tr4
       ch(2,k,2) = ti1+ti4
       ch(2,k,4) = ti1-ti4
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ti1 = cc(i,1,k)-cc(i,3,k)
          ti2 = cc(i,1,k)+cc(i,3,k)
          ti3 = cc(i,2,k)+cc(i,4,k)
          tr4 = cc(i,2,k)-cc(i,4,k)
          tr1 = cc(i-1,1,k)-cc(i-1,3,k)
          tr2 = cc(i-1,1,k)+cc(i-1,3,k)
          ti4 = cc(i-1,4,k)-cc(i-1,2,k)
          tr3 = cc(i-1,2,k)+cc(i-1,4,k)
          ch(i-1,k,1) = tr2+tr3
          cr3 = tr2-tr3
          ch(i,k,1) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1+tr4
          cr4 = tr1-tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
          ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
          ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
          ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
          ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
          ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passf5 {{{
subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
    dimension cc(ido,5,l1),ch(ido,l1,5),wa1(1),wa2(1),wa3(1),wa4(1)
    data tr11,ti11,tr12,ti12 /0.309016994374947,- 0.951056516295154,0.190983005625053,- 0.587785252292473/
    if (ido .ne. 2) go to 102
    do 101 k=1,l1
       ti5 = cc(2,2,k)-cc(2,5,k)
       ti2 = cc(2,2,k)+cc(2,5,k)
       ti4 = cc(2,3,k)-cc(2,4,k)
       ti3 = cc(2,3,k)+cc(2,4,k)
       tr5 = cc(1,2,k)-cc(1,5,k)
       tr2 = cc(1,2,k)+cc(1,5,k)
       tr4 = cc(1,3,k)-cc(1,4,k)
       tr3 = cc(1,3,k)+cc(1,4,k)
       ch(1,k,1) = cc(1,1,k)+tr2+tr3
       ch(2,k,1) = cc(2,1,k)+ti2+ti3
       cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
       ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
       cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
       ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
       cr5 = ti11*tr5+ti12*tr4
       ci5 = ti11*ti5+ti12*ti4
       cr4 = ti12*tr5-ti11*tr4
       ci4 = ti12*ti5-ti11*ti4
       ch(1,k,2) = cr2-ci5
       ch(1,k,5) = cr2+ci5
       ch(2,k,2) = ci2+cr5
       ch(2,k,3) = ci3+cr4
       ch(1,k,3) = cr3-ci4
       ch(1,k,4) = cr3+ci4
       ch(2,k,4) = ci3-cr4
       ch(2,k,5) = ci2-cr5
  101 continue
    return
  102 do 104 k=1,l1
       do 103 i=2,ido,2
          ti5 = cc(i,2,k)-cc(i,5,k)
          ti2 = cc(i,2,k)+cc(i,5,k)
          ti4 = cc(i,3,k)-cc(i,4,k)
          ti3 = cc(i,3,k)+cc(i,4,k)
          tr5 = cc(i-1,2,k)-cc(i-1,5,k)
          tr2 = cc(i-1,2,k)+cc(i-1,5,k)
          tr4 = cc(i-1,3,k)-cc(i-1,4,k)
          tr3 = cc(i-1,3,k)+cc(i-1,4,k)
          ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
          ch(i,k,1) = cc(i,1,k)+ti2+ti3
          cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
          ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
          cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
          ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
          ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
          ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
          ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
          ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
          ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
          ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
          ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
  103    continue
  104 continue
    return
    end
!}}}
!subroutine passf {{{
subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
    dimension ch(ido,l1,ip),cc(ido,ip,l1),c1(ido,l1,ip),wa(1),c2(idl1,ip),ch2(idl1,ip)
    idot = ido/2
    nt = ip*idl1
    ipp2 = ip+2
    ipph = (ip+1)/2
    idp = ip*ido
    if (ido .lt. l1) go to 106
    do 103 j=2,ipph
       jc = ipp2-j
       do 102 k=1,l1
          do 101 i=1,ido
             ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
             ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101     continue
  102    continue
  103 continue
    do 105 k=1,l1
       do 104 i=1,ido
          ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
    go to 112
  106 do 109 j=2,ipph
       jc = ipp2-j
       do 108 i=1,ido
          do 107 k=1,l1
             ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
             ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107     continue
  108    continue
  109 continue
    do 111 i=1,ido
       do 110 k=1,l1
          ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
  112 idl = 2-ido
    inc = 0
    do 116 l=2,ipph
       lc = ipp2-l
       idl = idl+ido
       do 113 ik=1,idl1
          c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
          c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  113    continue
       idlj = idl
       inc = inc+ido
       do 115 j=3,ipph
          jc = ipp2-j
          idlj = idlj+inc
          if (idlj .gt. idp) idlj = idlj-idp
          war = wa(idlj-1)
          wai = wa(idlj)
          do 114 ik=1,idl1
             c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
             c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
  114     continue
  115    continue
  116 continue
    do 118 j=2,ipph
       do 117 ik=1,idl1
          ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
    do 120 j=2,ipph
       jc = ipp2-j
       do 119 ik=2,idl1,2
          ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
          ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
          ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
          ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
    nac = 1
    if (ido .eq. 2) return
    nac = 0
    do 121 ik=1,idl1
       c2(ik,1) = ch2(ik,1)
  121 continue
    do 123 j=2,ip
       do 122 k=1,l1
          c1(1,k,j) = ch(1,k,j)
          c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
    if (idot .gt. l1) go to 127
    idij = 0
    do 126 j=2,ip
       idij = idij+2
       do 125 i=4,ido,2
          idij = idij+2
          do 124 k=1,l1
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  124     continue
  125    continue
  126 continue
    return
  127 idj = 2-ido
    do 130 j=2,ip
       idj = idj+ido
       do 129 k=1,l1
          idij = idj
          do 128 i=4,ido,2
             idij = idij+2
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  128     continue
  129    continue
  130 continue
    return
    end
!}}}
!subroutine radb2 {{{
subroutine radb2 (ido,l1,cc,ch,wa1)
    dimension cc(ido,2,l1),ch(ido,l1,2),wa1(1)
    do 101 k=1,l1
       ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
       ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
    if (ido-2) 107,105,102
  102 idp2 = ido+2
    do 104 k=1,l1
       do 103 i=3,ido,2
          ic = idp2-i
          ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
          tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
          ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
          ti2 = cc(i,1,k)+cc(ic,2,k)
          ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
          ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
    if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
       ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
       ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
  107 return
    end
!}}}
!subroutine radb3 {{{
subroutine radb3 (ido,l1,cc,ch,wa1,wa2)
    dimension cc(ido,3,l1),ch(ido,l1,3),wa1(1),wa2(1)
    data taur,taui /-.5,.866025403784439/
    do 101 k=1,l1
       tr2 = cc(ido,2,k)+cc(ido,2,k)
       cr2 = cc(1,1,k)+taur*tr2
       ch(1,k,1) = cc(1,1,k)+tr2
       ci3 = taui*(cc(1,3,k)+cc(1,3,k))
       ch(1,k,2) = cr2-ci3
       ch(1,k,3) = cr2+ci3
  101 continue
    if (ido .eq. 1) return
    idp2 = ido+2
    do 103 k=1,l1
       do 102 i=3,ido,2
          ic = idp2-i
          tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
          cr2 = cc(i-1,1,k)+taur*tr2
          ch(i-1,k,1) = cc(i-1,1,k)+tr2
          ti2 = cc(i,3,k)-cc(ic,2,k)
          ci2 = cc(i,1,k)+taur*ti2
          ch(i,k,1) = cc(i,1,k)+ti2
          cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
          ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
          dr2 = cr2-ci3
          dr3 = cr2+ci3
          di2 = ci2+cr3
          di3 = ci2-cr3
          ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
          ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
          ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
          ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
    return
    end
!}}}
!subroutine radb4 {{{
subroutine radb4 (ido,l1,cc,ch,wa1,wa2,wa3)
    dimension cc(ido,4,l1),ch(ido,l1,4),wa1(1),wa2(1),wa3(1)
    data sqrt2 /1.414213562373095/
    do 101 k=1,l1
       tr1 = cc(1,1,k)-cc(ido,4,k)
       tr2 = cc(1,1,k)+cc(ido,4,k)
       tr3 = cc(ido,2,k)+cc(ido,2,k)
       tr4 = cc(1,3,k)+cc(1,3,k)
       ch(1,k,1) = tr2+tr3
       ch(1,k,2) = tr1-tr4
       ch(1,k,3) = tr2-tr3
       ch(1,k,4) = tr1+tr4
  101 continue
    if (ido-2) 107,105,102
  102 idp2 = ido+2
    do 104 k=1,l1
       do 103 i=3,ido,2
          ic = idp2-i
          ti1 = cc(i,1,k)+cc(ic,4,k)
          ti2 = cc(i,1,k)-cc(ic,4,k)
          ti3 = cc(i,3,k)-cc(ic,2,k)
          tr4 = cc(i,3,k)+cc(ic,2,k)
          tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
          tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
          ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
          tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
          ch(i-1,k,1) = tr2+tr3
          cr3 = tr2-tr3
          ch(i,k,1) = ti2+ti3
          ci3 = ti2-ti3
          cr2 = tr1-tr4
          cr4 = tr1+tr4
          ci2 = ti1+ti4
          ci4 = ti1-ti4
          ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
          ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
          ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
          ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
          ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
          ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
    if (mod(ido,2) .eq. 1) return
  105 continue
    do 106 k=1,l1
       ti1 = cc(1,2,k)+cc(1,4,k)
       ti2 = cc(1,4,k)-cc(1,2,k)
       tr1 = cc(ido,1,k)-cc(ido,3,k)
       tr2 = cc(ido,1,k)+cc(ido,3,k)
       ch(ido,k,1) = tr2+tr2
       ch(ido,k,2) = sqrt2*(tr1-ti1)
       ch(ido,k,3) = ti2+ti2
       ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
  107 return
    end
!}}}
!subroutine radb5 {{{
subroutine radb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
    dimension cc(ido,5,l1),ch(ido,l1,5),wa1(1),wa2(1),wa3(1),wa4(1)
    data tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,0.190983005625053,0.587785252292473/
    do 101 k=1,l1
       ti5 = cc(1,3,k)+cc(1,3,k)
       ti4 = cc(1,5,k)+cc(1,5,k)
       tr2 = cc(ido,2,k)+cc(ido,2,k)
       tr3 = cc(ido,4,k)+cc(ido,4,k)
       ch(1,k,1) = cc(1,1,k)+tr2+tr3
       cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
       cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
       ci5 = ti11*ti5+ti12*ti4
       ci4 = ti12*ti5-ti11*ti4
       ch(1,k,2) = cr2-ci5
       ch(1,k,3) = cr3-ci4
       ch(1,k,4) = cr3+ci4
       ch(1,k,5) = cr2+ci5
  101 continue
    if (ido .eq. 1) return
    idp2 = ido+2
    do 103 k=1,l1
       do 102 i=3,ido,2
          ic = idp2-i
          ti5 = cc(i,3,k)+cc(ic,2,k)
          ti2 = cc(i,3,k)-cc(ic,2,k)
          ti4 = cc(i,5,k)+cc(ic,4,k)
          ti3 = cc(i,5,k)-cc(ic,4,k)
          tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
          tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
          tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
          tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
          ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
          ch(i,k,1) = cc(i,1,k)+ti2+ti3
          cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
          ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
          cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
          ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
          cr5 = ti11*tr5+ti12*tr4
          ci5 = ti11*ti5+ti12*ti4
          cr4 = ti12*tr5-ti11*tr4
          ci4 = ti12*ti5-ti11*ti4
          dr3 = cr3-ci4
          dr4 = cr3+ci4
          di3 = ci3+cr4
          di4 = ci3-cr4
          dr5 = cr2+ci5
          dr2 = cr2-ci5
          di5 = ci2-cr5
          di2 = ci2+cr5
          ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
          ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
          ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
          ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
          ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
          ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
          ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
          ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
    return
    end
!}}}
!subroutine radbg {{{
subroutine radbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
    dimension ch(ido,l1,ip),cc(ido,ip,l1),c1(ido,l1,ip),c2(idl1,ip),ch2(idl1,ip),wa(1)
    data tpi/6.28318530717959/
    arg = tpi/float(ip)
    dcp = cos(arg)
    dsp = sin(arg)
    idp2 = ido+2
    nbd = (ido-1)/2
    ipp2 = ip+2
    ipph = (ip+1)/2
    if (ido .lt. l1) go to 103
    do 102 k=1,l1
       do 101 i=1,ido
          ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
    go to 106
  103 do 105 i=1,ido
       do 104 k=1,l1
          ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
  106 do 108 j=2,ipph
       jc = ipp2-j
       j2 = j+j
       do 107 k=1,l1
          ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
          ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
    if (ido .eq. 1) go to 116
    if (nbd .lt. l1) go to 112
    do 111 j=2,ipph
       jc = ipp2-j
       do 110 k=1,l1
          do 109 i=3,ido,2
             ic = idp2-i
             ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
             ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
             ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
             ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109     continue
  110    continue
  111 continue
    go to 116
  112 do 115 j=2,ipph
       jc = ipp2-j
       do 114 i=3,ido,2
          ic = idp2-i
          do 113 k=1,l1
             ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
             ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
             ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
             ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113     continue
  114    continue
  115 continue
  116 ar1 = 1.
    ai1 = 0.
    do 120 l=2,ipph
       lc = ipp2-l
       ar1h = dcp*ar1-dsp*ai1
       ai1 = dcp*ai1+dsp*ar1
       ar1 = ar1h
       do 117 ik=1,idl1
          c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
          c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
       dc2 = ar1
       ds2 = ai1
       ar2 = ar1
       ai2 = ai1
       do 119 j=3,ipph
          jc = ipp2-j
          ar2h = dc2*ar2-ds2*ai2
          ai2 = dc2*ai2+ds2*ar2
          ar2 = ar2h
          do 118 ik=1,idl1
             c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
             c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118     continue
  119    continue
  120 continue
    do 122 j=2,ipph
       do 121 ik=1,idl1
          ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
    do 124 j=2,ipph
       jc = ipp2-j
       do 123 k=1,l1
          ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
          ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
    if (ido .eq. 1) go to 132
    if (nbd .lt. l1) go to 128
    do 127 j=2,ipph
       jc = ipp2-j
       do 126 k=1,l1
          do 125 i=3,ido,2
             ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
             ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
             ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
             ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125     continue
  126    continue
  127 continue
    go to 132
  128 do 131 j=2,ipph
       jc = ipp2-j
       do 130 i=3,ido,2
          do 129 k=1,l1
             ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
             ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
             ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
             ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129     continue
  130    continue
  131 continue
  132 continue
    if (ido .eq. 1) return
    do 133 ik=1,idl1
       c2(ik,1) = ch2(ik,1)
  133 continue
    do 135 j=2,ip
       do 134 k=1,l1
          c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
    if (nbd .gt. l1) go to 139
    is = -ido
    do 138 j=2,ip
       is = is+ido
       idij = is
       do 137 i=3,ido,2
          idij = idij+2
          do 136 k=1,l1
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136     continue
  137    continue
  138 continue
    go to 143
  139 is = -ido
    do 142 j=2,ip
       is = is+ido
       do 141 k=1,l1
          idij = is
          do 140 i=3,ido,2
             idij = idij+2
             c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
             c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140     continue
  141    continue
  142 continue
  143 return
    end
!}}}
!subroutine radf2 {{{
subroutine radf2 (ido,l1,cc,ch,wa1)
    dimension ch(ido,2,l1),cc(ido,l1,2),wa1(1)
    do 101 k=1,l1
       ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
       ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
    if (ido-2) 107,105,102
  102 idp2 = ido+2
    do 104 k=1,l1
       do 103 i=3,ido,2
          ic = idp2-i
          tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
          ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
          ch(i,1,k) = cc(i,k,1)+ti2
          ch(ic,2,k) = ti2-cc(i,k,1)
          ch(i-1,1,k) = cc(i-1,k,1)+tr2
          ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
    if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
       ch(1,2,k) = -cc(ido,k,2)
       ch(ido,1,k) = cc(ido,k,1)
  106 continue
  107 return
    end
!}}}
!subroutine radf3 {{{
subroutine radf3 (ido,l1,cc,ch,wa1,wa2)
    dimension ch(ido,3,l1),cc(ido,l1,3),wa1(1),wa2(1)
    data taur,taui /-.5,.866025403784439/
    do 101 k=1,l1
       cr2 = cc(1,k,2)+cc(1,k,3)
       ch(1,1,k) = cc(1,k,1)+cr2
       ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
       ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
    if (ido .eq. 1) return
    idp2 = ido+2
    do 103 k=1,l1
       do 102 i=3,ido,2
          ic = idp2-i
          dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
          di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
          dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
          di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
          cr2 = dr2+dr3
          ci2 = di2+di3
          ch(i-1,1,k) = cc(i-1,k,1)+cr2
          ch(i,1,k) = cc(i,k,1)+ci2
          tr2 = cc(i-1,k,1)+taur*cr2
          ti2 = cc(i,k,1)+taur*ci2
          tr3 = taui*(di2-di3)
          ti3 = taui*(dr3-dr2)
          ch(i-1,3,k) = tr2+tr3
          ch(ic-1,2,k) = tr2-tr3
          ch(i,3,k) = ti2+ti3
          ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
    return
    end
!}}}
!subroutine radf4 {{{
subroutine radf4 (ido,l1,cc,ch,wa1,wa2,wa3)
    dimension cc(ido,l1,4),ch(ido,4,l1),wa1(1),wa2(1),wa3(1)
    data hsqt2 /.7071067811865475/
    do 101 k=1,l1
       tr1 = cc(1,k,2)+cc(1,k,4)
       tr2 = cc(1,k,1)+cc(1,k,3)
       ch(1,1,k) = tr1+tr2
       ch(ido,4,k) = tr2-tr1
       ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
       ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
    if (ido-2) 107,105,102
  102 idp2 = ido+2
    do 104 k=1,l1
       do 103 i=3,ido,2
          ic = idp2-i
          cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
          ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
          cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
          ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
          cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
          ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
          tr1 = cr2+cr4
          tr4 = cr4-cr2
          ti1 = ci2+ci4
          ti4 = ci2-ci4
          ti2 = cc(i,k,1)+ci3
          ti3 = cc(i,k,1)-ci3
          tr2 = cc(i-1,k,1)+cr3
          tr3 = cc(i-1,k,1)-cr3
          ch(i-1,1,k) = tr1+tr2
          ch(ic-1,4,k) = tr2-tr1
          ch(i,1,k) = ti1+ti2
          ch(ic,4,k) = ti1-ti2
          ch(i-1,3,k) = ti4+tr3
          ch(ic-1,2,k) = tr3-ti4
          ch(i,3,k) = tr4+ti3
          ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
    if (mod(ido,2) .eq. 1) return
  105 continue
    do 106 k=1,l1
       ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
       tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
       ch(ido,1,k) = tr1+cc(ido,k,1)
       ch(ido,3,k) = cc(ido,k,1)-tr1
       ch(1,2,k) = ti1-cc(ido,k,3)
       ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
  107 return
    end
!}}}
!subroutine radf5 {{{
subroutine radf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
    dimension cc(ido,l1,5),ch(ido,5,l1),wa1(1),wa2(1),wa3(1),wa4(1)
    data tr11,ti11,tr12,ti12 /.309016994374947,0.951056516295154,0.190983005625053,0.587785252292473/
    do 101 k=1,l1
       cr2 = cc(1,k,5)+cc(1,k,2)
       ci5 = cc(1,k,5)-cc(1,k,2)
       cr3 = cc(1,k,4)+cc(1,k,3)
       ci4 = cc(1,k,4)-cc(1,k,3)
       ch(1,1,k) = cc(1,k,1)+cr2+cr3
       ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
       ch(1,3,k) = ti11*ci5+ti12*ci4
       ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
       ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
    if (ido .eq. 1) return
    idp2 = ido+2
    do 103 k=1,l1
       do 102 i=3,ido,2
          ic = idp2-i
          dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
          di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
          dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
          di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
          dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
          di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
          dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
          di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
          cr2 = dr2+dr5
          ci5 = dr5-dr2
          cr5 = di2-di5
          ci2 = di2+di5
          cr3 = dr3+dr4
          ci4 = dr4-dr3
          cr4 = di3-di4
          ci3 = di3+di4
          ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
          ch(i,1,k) = cc(i,k,1)+ci2+ci3
          tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
          ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
          tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
          ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
          tr5 = ti11*cr5+ti12*cr4
          ti5 = ti11*ci5+ti12*ci4
          tr4 = ti12*cr5-ti11*cr4
          ti4 = ti12*ci5-ti11*ci4
          ch(i-1,3,k) = tr2+tr5
          ch(ic-1,2,k) = tr2-tr5
          ch(i,3,k) = ti2+ti5
          ch(ic,2,k) = ti5-ti2
          ch(i-1,5,k) = tr3+tr4
          ch(ic-1,4,k) = tr3-tr4
          ch(i,5,k) = ti3+ti4
          ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
    return
    end
!}}}
!subroutine radfg {{{
subroutine radfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
    dimension ch(ido,l1,ip),cc(ido,ip,l1),c1(ido,l1,ip),c2(idl1,ip),ch2(idl1,ip),wa(1)
    data tpi/6.28318530717959/
    arg = tpi/float(ip)
    dcp = cos(arg)
    dsp = sin(arg)
    ipph = (ip+1)/2
    ipp2 = ip+2
    idp2 = ido+2
    nbd = (ido-1)/2
    if (ido .eq. 1) go to 119
    do 101 ik=1,idl1
       ch2(ik,1) = c2(ik,1)
  101 continue
    do 103 j=2,ip
       do 102 k=1,l1
          ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
    if (nbd .gt. l1) go to 107
    is = -ido
    do 106 j=2,ip
       is = is+ido
       idij = is
       do 105 i=3,ido,2
          idij = idij+2
          do 104 k=1,l1
             ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
             ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104     continue
  105    continue
  106 continue
    go to 111
  107 is = -ido
    do 110 j=2,ip
       is = is+ido
       do 109 k=1,l1
          idij = is
          do 108 i=3,ido,2
             idij = idij+2
             ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
             ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108     continue
  109    continue
  110 continue
  111 if (nbd .lt. l1) go to 115
    do 114 j=2,ipph
       jc = ipp2-j
       do 113 k=1,l1
          do 112 i=3,ido,2
             c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
             c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
             c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
             c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112     continue
  113    continue
  114 continue
    go to 121
  115 do 118 j=2,ipph
       jc = ipp2-j
       do 117 i=3,ido,2
          do 116 k=1,l1
             c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
             c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
             c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
             c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116     continue
  117    continue
  118 continue
    go to 121
  119 do 120 ik=1,idl1
       c2(ik,1) = ch2(ik,1)
  120 continue
  121 do 123 j=2,ipph
       jc = ipp2-j
       do 122 k=1,l1
          c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
          c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
    ar1 = 1.
    ai1 = 0.
    do 127 l=2,ipph
       lc = ipp2-l
       ar1h = dcp*ar1-dsp*ai1
       ai1 = dcp*ai1+dsp*ar1
       ar1 = ar1h
       do 124 ik=1,idl1
          ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
          ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
       dc2 = ar1
       ds2 = ai1
       ar2 = ar1
       ai2 = ai1
       do 126 j=3,ipph
          jc = ipp2-j
          ar2h = dc2*ar2-ds2*ai2
          ai2 = dc2*ai2+ds2*ar2
          ar2 = ar2h
          do 125 ik=1,idl1
             ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
             ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125     continue
  126    continue
  127 continue
    do 129 j=2,ipph
       do 128 ik=1,idl1
          ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
    if (ido .lt. l1) go to 132
    do 131 k=1,l1
       do 130 i=1,ido
          cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
    go to 135
  132 do 134 i=1,ido
       do 133 k=1,l1
          cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
       jc = ipp2-j
       j2 = j+j
       do 136 k=1,l1
          cc(ido,j2-2,k) = ch(1,k,j)
          cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
    if (ido .eq. 1) return
    if (nbd .lt. l1) go to 141
    do 140 j=2,ipph
       jc = ipp2-j
       j2 = j+j
       do 139 k=1,l1
          do 138 i=3,ido,2
             ic = idp2-i
             cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
             cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
             cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
             cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138     continue
  139    continue
  140 continue
    return
  141 do 144 j=2,ipph
       jc = ipp2-j
       j2 = j+j
       do 143 i=3,ido,2
          ic = idp2-i
          do 142 k=1,l1
             cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
             cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
             cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
             cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142     continue
  143    continue
  144 continue
    return
    end
!}}}
!subroutine rfftb1 {{{
subroutine rfftb1 (n,c,ch,wa,ifac)
    dimension       ch(1)      ,c(1)       ,wa(1)      ,ifac(1)
    nf = ifac(2)
    na = 0
    l1 = 1
    iw = 1
    do 116 k1=1,nf
       ip = ifac(k1+2)
       l2 = ip*l1
       ido = n/l2
       idl1 = ido*l1
       if (ip .ne. 4) go to 103
       ix2 = iw+ido
       ix3 = ix2+ido
       if (na .ne. 0) go to 101
       call radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
       go to 102
  101    call radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
       go to 115
  103    if (ip .ne. 2) go to 106
       if (na .ne. 0) go to 104
       call radb2 (ido,l1,c,ch,wa(iw))
       go to 105
  104    call radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
       go to 115
  106    if (ip .ne. 3) go to 109
       ix2 = iw+ido
       if (na .ne. 0) go to 107
       call radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
       go to 108
  107    call radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
       go to 115
  109    if (ip .ne. 5) go to 112
       ix2 = iw+ido
       ix3 = ix2+ido
       ix4 = ix3+ido
       if (na .ne. 0) go to 110
       call radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
       go to 111
  110    call radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
       go to 115
  112    if (na .ne. 0) go to 113
       call radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
       go to 114
  113    call radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
       iw = iw+(ip-1)*ido
  116 continue
    if (na .eq. 0) return
    do 117 i=1,n
       c(i) = ch(i)
  117 continue
    return
    end
!}}}
!subroutine rfftb {{{
subroutine rfftb (n,r,wsave)
    dimension       r(1)       ,wsave(1)
    if (n .eq. 1) return
    call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
    return
    end
!}}}
!subroutine rfftf1 {{{
subroutine rfftf1 (n,c,ch,wa,ifac)
    dimension       ch(1)      ,c(1)       ,wa(1)      ,ifac(1)
    nf = ifac(2)
    na = 1
    l2 = n
    iw = n
    do 111 k1=1,nf
       kh = nf-k1
       ip = ifac(kh+3)
       l1 = l2/ip
       ido = n/l2
       idl1 = ido*l1
       iw = iw-(ip-1)*ido
       na = 1-na
       if (ip .ne. 4) go to 102
       ix2 = iw+ido
       ix3 = ix2+ido
       if (na .ne. 0) go to 101
       call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
       go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
       go to 110
  102    if (ip .ne. 2) go to 104
       if (na .ne. 0) go to 103
       call radf2 (ido,l1,c,ch,wa(iw))
       go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
       go to 110
  104    if (ip .ne. 3) go to 106
       ix2 = iw+ido
       if (na .ne. 0) go to 105
       call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
       go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
       go to 110
  106    if (ip .ne. 5) go to 108
       ix2 = iw+ido
       ix3 = ix2+ido
       ix4 = ix3+ido
       if (na .ne. 0) go to 107
       call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
       go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
       go to 110
  108    if (ido .eq. 1) na = 1-na
       if (na .ne. 0) go to 109
       call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
       na = 1
       go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
       na = 0
  110    l2 = l1
  111 continue
    if (na .eq. 1) return
    do 112 i=1,n
       c(i) = ch(i)
  112 continue
    return
    end
!}}}
!subroutine rfftf {{{
subroutine rfftf (n,r,wsave)
    dimension       r(1)       ,wsave(1)
    if (n .eq. 1) return
    call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
    return
    end
!}}}
!subroutine rffti1 {{{
subroutine rffti1 (n,wa,ifac)
    dimension       wa(1)      ,ifac(1)    ,ntryh(4)
    data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
    nl = n
    nf = 0
    j = 0
  101 j = j+1
    if (j-4) 102,102,103
  102 ntry = ntryh(j)
    go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
    nr = nl-ntry*nq
    if (nr) 101,105,101
  105 nf = nf+1
    ifac(nf+2) = ntry
    nl = nq
    if (ntry .ne. 2) go to 107
    if (nf .eq. 1) go to 107
    do 106 i=2,nf
       ib = nf-i+2
       ifac(ib+2) = ifac(ib+1)
  106 continue
    ifac(3) = 2
  107 if (nl .ne. 1) go to 104
    ifac(1) = n
    ifac(2) = nf
    tpi = 6.28318530717959
    argh = tpi/float(n)
    is = 0
    nfm1 = nf-1
    l1 = 1
    if (nfm1 .eq. 0) return
    do 110 k1=1,nfm1
       ip = ifac(k1+2)
       ld = 0
       l2 = l1*ip
       ido = n/l2
       ipm = ip-1
       do 109 j=1,ipm
          ld = ld+l1
          i = is
          argld = float(ld)*argh
          fi = 0.
          do 108 ii=3,ido,2
             i = i+2
             fi = fi+1.
             arg = fi*argld
             wa(i-1) = cos(arg)
             wa(i) = sin(arg)
  108     continue
          is = is+ido
  109    continue
       l1 = l2
  110 continue
    return
    end
!}}}
!subroutine rffti {{{
subroutine rffti (n,wsave)
    dimension       wsave(1)
    if (n .eq. 1) return
    call rffti1 (n,wsave(n+1),wsave(2*n+1))
    return
    end
!}}}
!subroutine sinqb {{{
subroutine sinqb (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    if (n .gt. 1) go to 101
    x(1) = 4.*x(1)
    return
  101 ns2 = n/2
    do 102 k=2,n,2
       x(k) = -x(k)
  102 continue
    call cosqb (n,x,wsave)
    do 103 k=1,ns2
       kc = n-k
       xhold = x(k)
       x(k) = x(kc+1)
       x(kc+1) = xhold
  103 continue
    return
    end
!}}}
!subroutine sinqf {{{
subroutine sinqf (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    if (n .eq. 1) return
    ns2 = n/2
    do 101 k=1,ns2
       kc = n-k
       xhold = x(k)
       x(k) = x(kc+1)
       x(kc+1) = xhold
  101 continue
    call cosqf (n,x,wsave)
    do 102 k=2,n,2
       x(k) = -x(k)
  102 continue
    return
    end
!}}}
!subroutine sinqi {{{
subroutine sinqi (n,wsave)
    dimension       wsave(1)
    call cosqi (n,wsave)
    return
    end
!}}}
!subroutine sint1{{{
subroutine sint1(n,war,was,xh,x,ifac)
    dimension war(1),was(1),x(1),xh(1),ifac(1)
    data sqrt3 /1.73205080756888/
    do 100 i=1,n
    xh(i) = war(i)
    war(i) = x(i)
  100 continue
    if (n-2) 101,102,103
  101 xh(1) = xh(1)+xh(1)
    go to 106
  102 xhold = sqrt3*(xh(1)+xh(2))
    xh(2) = sqrt3*(xh(1)-xh(2))
    xh(1) = xhold
    go to 106
  103 np1 = n+1
    ns2 = n/2
    x(1) = 0.
    do 104 k=1,ns2
       kc = np1-k
       t1 = xh(k)-xh(kc)
       t2 = was(k)*(xh(k)+xh(kc))
       x(k+1) = t1+t2
       x(kc+1) = t2-t1
  104 continue
    modn = mod(n,2)
    if (modn .ne. 0) x(ns2+2) = 4.*xh(ns2+1)
    call rfftf1 (np1,x,xh,war,ifac)
    xh(1) = .5*x(1)
    do 105 i=3,n,2
       xh(i-1) = -x(i)
       xh(i) = xh(i-2)+x(i-1)
  105 continue
    if (modn .ne. 0) go to 106
    xh(n) = -x(n+1)
  106 do 107 i=1,n
    x(i) = war(i)
    war(i) = xh(i)
  107 continue
    return
    end
!}}}
!subroutine sint {{{
subroutine sint (n,x,wsave)
    dimension       x(1)       ,wsave(1)
    np1 = n+1
    iw1 = n/2+1
    iw2 = iw1+np1
    iw3 = iw2+np1
    call sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
    return
    end
!}}}
!subroutine sinti {{{
subroutine sinti (n,wsave)
    dimension       wsave(1)
    data pi /3.14159265358979/
    if (n .le. 1) return
    ns2 = n/2
    np1 = n+1
    dt = pi/float(np1)
    do 101 k=1,ns2
       wsave(k) = 2.*sin(k*dt)
  101 continue
    call rffti (np1,wsave(ns2+1))
    return
    end
!}}}
!program tstfft{{{
program tstfft
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!                    version 4  april 1985
!
!                      a test driver for
!       a package of fortran subprograms for the fast fourier
!        transform of periodic and other symmetric sequences
!
!                           by
!
!                    paul n swarztrauber
!
!    national center for atmospheric research  boulder,colorado 80307
!
!     which is sponsored by the national science foundation
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!          this program tests the package of fast fourier
!    transforms for both complex and real periodic sequences and
!    certian other symmetric sequences that are listed below.
!
!    1.   rffti     initialize  rfftf and rfftb
!    2.   rfftf     forward transform of a real periodic sequence
!    3.   rfftb     backward transform of a real coefficient array
!
!    4.   ezffti    initialize ezfftf and ezfftb
!    5.   ezfftf    a simplified real periodic forward transform
!    6.   ezfftb    a simplified real periodic backward transform
!
!    7.   sinti     initialize sint
!    8.   sint    sine transform of a real odd sequence
!
!    9.   costi     initialize cost
!    10.  cost    cosine transform of a real even sequence
!
!    11.  sinqi     initialize sinqf and sinqb
!    12.  sinqf     forward sine transform with odd wave numbers
!    13.  sinqb     unnormalized inverse of sinqf
!
!    14.  cosqi     initialize cosqf and cosqb
!    15.  cosqf     forward cosine transform with odd wave numbers
!    16.  cosqb     unnormalized inverse of cosqf
!
!    17.  cffti     initialize cfftf and cfftb
!    18.  cfftf     forward transform of a complex periodic sequence
!    19.  cfftb     unnormalized inverse of cfftf
!
!
    dimension nd(10),x(200),y(200),w(2000),a(100),b(100),ah(100),bh(100),xh(200),cx(200),cy(200)
    complex         cx         ,cy
    data nd(1),nd(2),nd(3),nd(4),nd(5),nd(6),nd(7)/120,54,49,32,4,3,2/
    sqrt2 = sqrt(2.)
    nns = 7
    do 157 nz=1,nns
       n = nd(nz)
       modn = mod(n,2)
       fn = float(n)
       tfn = fn+fn
       np1 = n+1
       nm1 = n-1
       do 101 j=1,np1
          x(j) = sin(float(j)*sqrt2)
          y(j) = x(j)
          xh(j) = x(j)
  101    continue
!
!!    test subroutines rffti,rfftf and rfftb
!
       call rffti (n,w)
       pi = 3.14159265358979
       dt = (pi+pi)/fn
       ns2 = (n+1)/2
       if (ns2 .lt. 2) go to 104
       do 103 k=2,ns2
          sum1 = 0.
          sum2 = 0.
          arg = float(k-1)*dt
          do 102 i=1,n
             arg1 = float(i-1)*arg
             sum1 = sum1+x(i)*cos(arg1)
             sum2 = sum2+x(i)*sin(arg1)
  102     continue
          y(2*k-2) = sum1
          y(2*k-1) = -sum2
  103    continue
  104    sum1 = 0.
       sum2 = 0.
       do 105 i=1,nm1,2
          sum1 = sum1+x(i)
          sum2 = sum2+x(i+1)
  105    continue
       if (modn .eq. 1) sum1 = sum1+x(n)
       y(1) = sum1+sum2
       if (modn .eq. 0) y(n) = sum1-sum2
       call rfftf (n,x,w)
       rftf = 0.
       do 106 i=1,n
          rftf = amax1(rftf,abs(x(i)-y(i)))
          x(i) = xh(i)
  106    continue
       rftf = rftf/fn
       do 109 i=1,n
          sum = .5*x(1)
          arg = float(i-1)*dt
          if (ns2 .lt. 2) go to 108
          do 107 k=2,ns2
             arg1 = float(k-1)*arg
             sum = sum+x(2*k-2)*cos(arg1)-x(2*k-1)*sin(arg1)
  107     continue
  108     if (modn .eq. 0) sum = sum+.5*float((-1)**(i-1))*x(n)
          y(i) = sum+sum
  109    continue
       call rfftb (n,x,w)
       rftb = 0.
       do 110 i=1,n
          rftb = amax1(rftb,abs(x(i)-y(i)))
          x(i) = xh(i)
          y(i) = xh(i)
  110    continue
       call rfftb (n,y,w)
       call rfftf (n,y,w)
       cf = 1./fn
       rftfb = 0.
       do 111 i=1,n
          rftfb = amax1(rftfb,abs(cf*y(i)-x(i)))
  111    continue
!
!     test subroutines sinti and sint
!
       dt = pi/fn
       do 112 i=1,nm1
          x(i) = xh(i)
  112    continue
       do 114 i=1,nm1
          y(i) = 0.
          arg1 = float(i)*dt
          do 113 k=1,nm1
             y(i) = y(i)+x(k)*sin(float(k)*arg1)
  113     continue
          y(i) = y(i)+y(i)
  114    continue
       call sinti (nm1,w)
       call sint (nm1,x,w)
       cf = .5/fn
       sintt = 0.
       do 115 i=1,nm1
          sintt = amax1(sintt,abs(x(i)-y(i)))
          x(i) = xh(i)
          y(i) = x(i)
  115    continue
       sintt = cf*sintt
       call sint (nm1,x,w)
       call sint (nm1,x,w)
       sintfb = 0.
       do 116 i=1,nm1
          sintfb = amax1(sintfb,abs(cf*x(i)-y(i)))
  116    continue
!
!     test subroutines costi and cost
!
       do 117 i=1,np1
          x(i) = xh(i)
  117    continue
       do 119 i=1,np1
          y(i) = .5*(x(1)+float((-1)**(i+1))*x(n+1))
          arg = float(i-1)*dt
          do 118 k=2,n
             y(i) = y(i)+x(k)*cos(float(k-1)*arg)
  118     continue
          y(i) = y(i)+y(i)
  119    continue
       call costi (np1,w)
       call cost (np1,x,w)
       costt = 0.
       do 120 i=1,np1
          costt = amax1(costt,abs(x(i)-y(i)))
          x(i) = xh(i)
          y(i) = xh(i)
  120    continue
       costt = cf*costt
       call cost (np1,x,w)
       call cost (np1,x,w)
       costfb = 0.
       do 121 i=1,np1
          costfb = amax1(costfb,abs(cf*x(i)-y(i)))
  121    continue
!
!     test subroutines sinqi,sinqf and sinqb
!
       cf = .25/fn
       do 122 i=1,n
          y(i) = xh(i)
  122    continue
       dt = pi/(fn+fn)
       do 124 i=1,n
          x(i) = 0.
          arg = dt*float(i)
          do 123 k=1,n
             x(i) = x(i)+y(k)*sin(float(k+k-1)*arg)
  123     continue
          x(i) = 4.*x(i)
  124    continue
       call sinqi (n,w)
       call sinqb (n,y,w)
       sinqbt = 0.
       do 125 i=1,n
          sinqbt = amax1(sinqbt,abs(y(i)-x(i)))
          x(i) = xh(i)
  125    continue
       sinqbt = cf*sinqbt
       do 127 i=1,n
          arg = float(i+i-1)*dt
          y(i) = .5*float((-1)**(i+1))*x(n)
          do 126 k=1,nm1
             y(i) = y(i)+x(k)*sin(float(k)*arg)
  126     continue
          y(i) = y(i)+y(i)
  127    continue
       call sinqf (n,x,w)
       sinqft = 0.
       do 128 i=1,n
          sinqft = amax1(sinqft,abs(x(i)-y(i)))
          y(i) = xh(i)
          x(i) = xh(i)
  128    continue
       call sinqf (n,y,w)
       call sinqb (n,y,w)
       sinqfb = 0.
       do 129 i=1,n
          sinqfb = amax1(sinqfb,abs(cf*y(i)-x(i)))
  129    continue
!
!     test subroutines cosqi,cosqf and cosqb
!
       do 130 i=1,n
          y(i) = xh(i)
  130    continue
       do 132 i=1,n
          x(i) = 0.
          arg = float(i-1)*dt
          do 131 k=1,n
             x(i) = x(i)+y(k)*cos(float(k+k-1)*arg)
  131     continue
          x(i) = 4.*x(i)
  132    continue
       call cosqi (n,w)
       call cosqb (n,y,w)
       cosqbt = 0.
       do 133 i=1,n
          cosqbt = amax1(cosqbt,abs(x(i)-y(i)))
          x(i) = xh(i)
  133    continue
       cosqbt = cf*cosqbt
       do 135 i=1,n
          y(i) = .5*x(1)
          arg = float(i+i-1)*dt
          do 134 k=2,n
             y(i) = y(i)+x(k)*cos(float(k-1)*arg)
  134     continue
          y(i) = y(i)+y(i)
  135    continue
       call cosqf (n,x,w)
       cosqft = 0.
       do 136 i=1,n
          cosqft = amax1(cosqft,abs(y(i)-x(i)))
          x(i) = xh(i)
          y(i) = xh(i)
  136    continue
       cosqft = cf*cosqft
       call cosqb (n,x,w)
       call cosqf (n,x,w)
       cosqfb = 0.
       do 137 i=1,n
          cosqfb = amax1(cosqfb,abs(cf*x(i)-y(i)))
  137    continue
!
!     test programs ezffti,ezfftf,ezfftb
!
       call ezffti(n,w)
       do 138 i=1,n
          x(i) = xh(i)
  138    continue
       tpi = 8.*atan(1.)
       dt = tpi/float(n)
       ns2 = (n+1)/2
       cf = 2./float(n)
       ns2m = ns2-1
       if (ns2m .le. 0) go to 141
       do 140 k=1,ns2m
          sum1 = 0.
          sum2 = 0.
          arg = float(k)*dt
          do 139 i=1,n
             arg1 = float(i-1)*arg
             sum1 = sum1+x(i)*cos(arg1)
             sum2 = sum2+x(i)*sin(arg1)
  139     continue
          a(k) = cf*sum1
          b(k) = cf*sum2
  140    continue
  141    nm1 = n-1
       sum1 = 0.
       sum2 = 0.
       do 142 i=1,nm1,2
          sum1 = sum1+x(i)
          sum2 = sum2+x(i+1)
  142    continue
       if (modn .eq. 1) sum1 = sum1+x(n)
       azero = .5*cf*(sum1+sum2)
       if (modn .eq. 0) a(ns2) = .5*cf*(sum1-sum2)
       call ezfftf (n,x,azeroh,ah,bh,w)
       dezf1 = abs(azeroh-azero)
       if (modn .eq. 0) dezf1 = amax1(dezf1,abs(a(ns2)-ah(ns2)))
       if (ns2m .le. 0) go to 144
       do 143 i=1,ns2m
          dezf1 = amax1(dezf1,abs(ah(i)-a(i)),abs(bh(i)-b(i)))
  143    continue
  144    ns2 = n/2
       if (modn .eq. 0) b(ns2) = 0.
       do 146 i=1,n
          sum = azero
          arg1 = float(i-1)*dt
          do 145 k=1,ns2
             arg2 = float(k)*arg1
             sum = sum+a(k)*cos(arg2)+b(k)*sin(arg2)
  145     continue
          x(i) = sum
  146    continue
       call ezfftb (n,y,azero,a,b,w)
       dezb1 = 0.
       do 147 i=1,n
          dezb1 = amax1(dezb1,abs(x(i)-y(i)))
          x(i) = xh(i)
  147    continue
       call ezfftf (n,x,azero,a,b,w)
       call ezfftb (n,y,azero,a,b,w)
       dezfb = 0.
       do 148 i=1,n
          dezfb = amax1(dezfb,abs(x(i)-y(i)))
  148    continue
!
!     test  cffti,cfftf,cfftb
!
       do 149 i=1,n
          cx(i) = cmplx(cos(sqrt2*float(i)),sin(sqrt2*float(i*i)))
  149    continue
       dt = (pi+pi)/fn
       do 151 i=1,n
          arg1 = -float(i-1)*dt
          cy(i) = (0.,0.)
          do 150 k=1,n
             arg2 = float(k-1)*arg1
             cy(i) = cy(i)+cmplx(cos(arg2),sin(arg2))*cx(k)
  150     continue
  151    continue
       call cffti (n,w)
       call cfftf (n,cx,w)
       dcfftf = 0.
       do 152 i=1,n
          dcfftf = amax1(dcfftf,cabs(cx(i)-cy(i)))
          cx(i) = cx(i)/fn
  152    continue
       dcfftf = dcfftf/fn
       do 154 i=1,n
          arg1 = float(i-1)*dt
          cy(i) = (0.,0.)
          do 153 k=1,n
             arg2 = float(k-1)*arg1
             cy(i) = cy(i)+cmplx(cos(arg2),sin(arg2))*cx(k)
  153     continue
  154    continue
       call cfftb (n,cx,w)
       dcfftb = 0.
       do 155 i=1,n
          dcfftb = amax1(dcfftb,cabs(cx(i)-cy(i)))
          cx(i) = cy(i)
  155    continue
       cf = 1./fn
       call cfftf (n,cx,w)
       call cfftb (n,cx,w)
       dcfb = 0.
       do 156 i=1,n
          dcfb = amax1(dcfb,cabs(cf*cx(i)-cy(i)))
  156    continue
       write (6,1001) n,rftf,rftb,rftfb,sintt,sintfb,costt,costfb,&
                      sinqft,sinqbt,sinqfb,cosqft,cosqbt,cosqfb,dezf1, &
                      dezb1,dezfb,dcfftf,dcfftb,dcfb
  157 continue
!
!
!
 1001 format (2h0n,i5,8h rfftf  ,e10.3,8h rfftb  ,e10.3,8h rfftfb ,&
            e10.3,8h sint   ,e10.3,8h sintfb ,e10.3,8h cost   ,e10.3/&
            7x,8h costfb ,e10.3,8h sinqf  ,e10.3,8h sinqb  ,e10.3,&
            8h sinqfb ,e10.3,8h cosqf  ,e10.3,8h cosqb  ,e10.3/7x,&
            8h cosqfb ,e10.3,8h dezf   ,e10.3,8h dezb   ,e10.3,&
            8h dezfb  ,e10.3,8h cfftf  ,e10.3,8h cfftb  ,e10.3/&
            7x,8h cfftfb ,e10.3)
!
    end
!}}}
