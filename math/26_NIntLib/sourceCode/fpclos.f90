      subroutine fpclos(iopt,idim,m,u,mx,x,w,k,s,nest,tol,maxit,k1,k2,
     * n,t,nc,c,fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
!  ..
!  ..scalar arguments..
      real s,tol,fp
      integer iopt,idim,m,mx,k,nest,maxit,k1,k2,n,nc,ier
!  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),z(nc),a1(nest,k1),
     * a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),q(m,k1)
      integer nrdata(nest)
!  ..local scalars..
      real acc,cos,d1,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv,
     * p1,p2,p3,sin,store,term,ui,wi,rn,one,con1,con4,con9,half
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jj,jk,jper,j1,j2,kk,
     * kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1,
     * nrint,n10,n11,n7,n8
!  ..local arrays..
      real h(6),h1(7),h2(6),xi(10)
!  ..function references..
      real abs,fprati
      integer max0,min0
!  ..subroutine references..
!    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
!  ..
!  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 1: determination of the number of knots and their position     c
!  **************************************************************      c
!  given a set of knots we compute the least-squares closed curve      c
!  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c
!  if iopt=-1 sinf(u) is the requested curve                           c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!    if fp <=s we will continue with the current set of knots.         c
!    if fp > s we will increase the number of knots and compute the    c
!       corresponding least-squares curve until finally fp<=s.         c
!  the initial choice of knots depends on the value of s and iopt.     c
!    if s=0 we have spline interpolation; in that case the number of   c
!    knots equals nmax = m+2*k.                                        c
!    if s > 0 and                                                      c
!      iopt=0 we first compute the least-squares polynomial curve of   c
!      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c
!      find that s(u) reduces to a fixed point.                        c
!      iopt=1 we start with the set of knots found at the last         c
!      call of the routine, except for the case that s > fp0; then     c
!      we compute directly the least-squares polynomial curve.         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
!  determine the length of the period of the splines.
      per = u(m)-u(1)
      if(iopt.lt.0) go to 50
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
!  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s.gt.0. .or. nmax.eq.nmin) go to 30
!  if s=0, s(u) is an interpolating curve.
      n = nmax
!  test whether the required storage space exceeds the available one.
      if(n.gt.nest) go to 620
!  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 .eq.k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = u(i)
  10  continue
      if(s.gt.0.) go to 50
      kk = k-1
      kk1 = k
      if(kk.gt.0) go to 50
      t(1) = t(m)-per
      t(2) = u(1)
      t(m+1) = u(m)
      t(m+2) = t(3)+per
      jj = 0
      do 15 i=1,m1
        j = i
        do 12 j1=1,idim
          jj = jj+1
          c(j) = x(jj)
          j = j+n
  12    continue
  15  continue
      jj = 1
      j = m
      do 17 j1=1,idim
        c(j) = c(jj)
        j = j+n
        jj = jj+n
  17  continue
      fp = 0.
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
         j = i+k
         t(j) = (u(i)+u(i-1))*half
  25  continue
      go to 50
!  if s > 0 our initial choice depends on the value of iopt.
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial curve. (i.e. a constant point).
!  if iopt=1 and fp0>s we start computing the least-squares closed
!  curve according the set of knots found at the last call of the
!  routine.
  30  if(iopt.eq.0) go to 35
      if(n.eq.nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 50
!  the case that s(u) is a fixed point is treated separetely.
!  fp0 denotes the corresponding sum of squared residuals.
  35  fp0 = 0.
      d1 = 0.
      do 37 j=1,idim
        z(j) = 0.
  37  continue
      jj = 0
      do 45 it=1,m1
        wi = w(it)
        call fpgivs(wi,d1,cos,sin)
        do 40 j=1,idim
          jj = jj+1
          fac = wi*x(jj)
          call fprota(cos,sin,fac,z(j))
          fp0 = fp0+fac**2
  40    continue
  45  continue
      do 47 j=1,idim
        z(j) = z(j)/d1
  47  continue
!  test whether that fixed point is a solution of our problem.
      fpms = fp0-s
      if(fpms.lt.acc .or. nmax.eq.nmin) go to 640
      fpold = fp0
!  test whether the required storage space exceeds the available one.
      if(n.ge.nest) go to 620
!  start computing the least-squares closed curve with one
!  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = u(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
!  main loop for the different sets of knots. m is a save upper
!  bound for the number of trials.
  50  do 340 iter=1,m
!  find nrint, the number of knot intervals.
        nrint = n-nmin+1
!  find the position of the additional knots which are needed for
!  the b-spline representation of s(u). if we take
!      t(k+1) = u(1), t(n-k) = u(m)
!      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
!      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
!  then s(u) will be a smooth closed curve if the b-spline
!  coefficients satisfy the following conditions
!      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
!  with n7=n-2*k-1.
        t(k1) = u(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = u(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
!  compute the b-spline coefficients of the least-squares closed curve
!  sinf(u). the observation matrix a is built up row by row while
!  taking into account condition (**) and is reduced to triangular
!  form by givens transformations .
!  at the same time fp=f(p=inf) is computed.
!  the n7 x n7 triangularised upper matrix a has the form
!            ! a1 '    !
!        a = !    ' a2 !
!            ! 0  '    !
!  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
!  matrix of bandwith k+1 ( n10 = n7-k).
!  initialization.
        do 65 i=1,nc
          z(i) = 0.
  65    continue
        do 70 i=1,nk1
          do 70 j=1,kk1
            a1(i,j) = 0.
  70    continue
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = 0.
        l = k1
        jj = 0
        do 290 it=1,m1
!  fetch the current data point u(it),x(it)
          ui = u(it)
          wi = w(it)
          do 75 j=1,idim
            jj = jj+1
            xi(j) = x(jj)*wi
  75      continue
!  search for knot interval t(l) <= ui < t(l+1).
  80      if(ui.lt.t(l+1)) go to 85
          l = l+1
          go to 80
!  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  85      call fpbspl(t,n,k,ui,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
!  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
          if(l5.lt.n10) go to 285
          if(jper.ne.0) go to 160
!  initialize the matrix a2.
          do 95 i=1,n7
          do 95 j=1,kk
              a2(i,j) = 0.
  95      continue
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik.le.0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
!  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
!  we take account of condition (**) for setting up the new row
!  of the observation matrix a. this row is stored in the arrays h1
!  (the part with respect to a1) and h2 (the part with
!  respect to a2).
 160      do 170 i=1,kk
            h1(i) = 0.
            h2(i) = 0.
 170      continue
          h1(kk1) = 0.
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1.le.0) go to 200
            if(l1.le.n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
!  rotate the new row of the observation matrix into triangle
!  by givens transformations.
          if(n10.le.0) go to 250
!  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv.ne.0.) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.
            go to 240
!  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
!  transformation to the right hand side.
            j1 = j
            do 217 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 217        continue
!  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j.eq.n10) go to 250
            i2 = min0(n10-j,kk)
!  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = 0.
 240      continue
!  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij.le.0) go to 270
            piv = h2(j)
            if(piv.eq.0.) go to 270
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
!  transformations to right hand side.
            j1 = ij
            do 255 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 255        continue
            if(j.eq.kk) go to 280
            j1 = j+1
!  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
!  add contribution of this row to the sum of squares of residual
!  right hand sides.
 280      do 282 j2=1,idim
            fp = fp+xi(j2)**2
 282      continue
          go to 290
!  rotation of the new row of the observation matrix into
!  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
!  at ui.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 140
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
!  transformations to right hand side.
            j1 = j
            do 125 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 125        continue
            if(i.eq.kk1) go to 150
            i2 = 1
            i3 = i+1
!  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
!  add contribution of this row to the sum of squares of residual
!  right hand sides.
 150      do 155 j2=1,idim
            fp = fp+xi(j2)**2
 155      continue
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
!  backward substitution to obtain the b-spline coefficients .
        j1 = 1
        do 292 j2=1,idim
           call fpbacp(a1,a2,z(j1),n7,kk,c(j1),kk1,nest)
           j1 = j1+n
 292    continue
!  calculate from condition (**) the remaining coefficients.
        do 297 i=1,k
          j1 = i
          do 295 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 295      continue
 297    continue
        if(iopt.lt.0) go to 660
!  test whether the approximation sinf(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
!  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 350
!  if n=nmax, sinf(u) is an interpolating curve.
        if(n.eq.nmax) go to 630
!  increase the number of knots.
!  if n=nest we cannot increase the number of knots because of the
!  storage capacity limitation.
        if(n.eq.nest) go to 620
!  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
!  compute the sum of squared residuals for each knot interval
!  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k1
        jj = 0
        do 320 it=1,m1
          if(u(it).lt.t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.
          l0 = l-k2
          do 310 j2=1,idim
            fac = 0.
            j1 = l0
            do 305 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 305        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 310      continue
          fpart = fpart+term
          if(new.eq.0) go to 320
          if(l.gt.k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
!  add a new knot
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
!  if n=nmax we locate the knots as for interpolation
          if(n.eq.nmax) go to 5
!  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 340
 330    continue
!  restart the computations with the new set of knots.
 340  continue
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 2: determination of the smoothing closed curve sp(u).          c
!  **********************************************************          c
!  we have determined the number of knots and their position.          c
!  we now compute the b-spline coefficients of the smoothing curve     c
!  sp(u). the observation matrix a is extended by the rows of matrix   c
!  b expressing that the kth derivative discontinuities of sp(u) at    c
!  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
!  ponding weights of these additional rows are set to 1/p.            c
!  iteratively we then have to determine the value of p such that f(p),c
!  the sum of squared residuals be = s. we already know that the least-c
!  squares polynomial curve corresponds to p=0, and that the least-    c
!  squares periodic spline curve corresponds to p=infinity. the        c
!  iteration process which is proposed here, makes use of rational     c
!  interpolation. since f(p) is a convex and strictly decreasing       c
!  function of p, it can be approximated by a rational function        c
!  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
!  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
!  to calculate the new value of p such that r(p)=s. convergence is    c
!  guaranteed by taking f1>0 and f3<0.                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  evaluate the discontinuity jump of the kth derivative of the
!  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
!  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l.eq.0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
!  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
!  form the matrix g  as the matrix a extended by the rows of matrix b.
!  the rows of matrix b with weight 1/p are rotated into
!  the triangularised observation matrix a.
!  after triangularisation our n7 x n7 matrix g takes the form
!            ! g1 '    !
!        g = !    ' g2 !
!            ! 0  '    !
!  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
!  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
!  store matrix a into g
        do 358 i=1,nc
          c(i) = z(i)
 358    continue
        do 360 i=1,n7
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = 0.
          g2(i,1) = 0.
          do 360 j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
 360    continue
        l = n10
        do 370 j=1,k1
          if(l.le.0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
!  fetch a new row of matrix b and store it in the arrays h1 (the part
!  with respect to g1) and h2 (the part with respect to g2).
          do 380 j=1,idim
            xi(j) = 0.
 380      continue
          do 385 i=1,k1
            h1(i) = 0.
            h2(i) = 0.
 385      continue
          h1(k2) = 0.
          if(it.gt.n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0.eq.n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1.le.0) go to 450
            if(l1.le.n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11.le.0) go to 510
!  rotate this row into triangle by givens transformations
!  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
!  transformation to right hand side.
            j1 = j
            do 475 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 475        continue
!  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j.eq.n11) go to 510
            i2 = min0(n11-j,k1)
!  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.
 500      continue
!  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij.le.0) go to 530
            piv = h2(j)
!  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
!  transformation to the right hand side.
            j1 = ij
            do 515 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 515        continue
            if(j.eq.k1) go to 540
            j1 = j+1
!  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
!  backward substitution to obtain the b-spline coefficients
        j1 = 1
        do 542 j2=1,idim
          call fpbacp(g1,g2,c(j1),n7,k1,c(j1),k2,nest)
          j1 = j1+n
 542    continue
!  calculate from condition (**) the remaining b-spline coefficients.
        do 547 i=1,k
          j1 = i
          do 545 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 545      continue
 547    continue
!  computation of f(p).
        fp = 0.
        l = k1
        jj = 0
        do 570 it=1,m1
          if(u(it).lt.t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.
          do 565 j2=1,idim
            fac = 0.
            j1 = l0
            do 560 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 560        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 565      continue
          fp = fp+term*w(it)**2
 570    continue
!  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
!  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 600
!  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 580
        if((f2-f3) .gt. acc) go to 575
!  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2.lt.0.) ich3 = 1
 580    if(ich1.ne.0) go to 590
        if((f1-f2) .gt. acc) go to 585
!  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 595
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2.gt.0.) ich1 = 1
!  test whether the iteration process proceeds as theoretically
!  expected.
 590    if(f2.ge.f1 .or. f2.le.f3) go to 610
!  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 595  continue
!  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
!  the point (z(1),z(2),...,z(idim)) is a solution of our problem.
!  a constant function is a spline of degree k with all b-spline
!  coefficients equal to that constant.
      do 650 i=1,k1
        rn = k1-i
        t(i) = u(1)-rn*per
        j = i+k1
        rn = i-1
        t(j) = u(m)+rn*per
 650  continue
      n = nmin
      j1 = 0
      do 658 j=1,idim
        fac = z(j)
        j2 = j1
        do 654 i=1,k1
          j2 = j2+1
          c(j2) = fac
 654    continue
        j1 = j1+n
 658  continue
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
 660  return
      end
