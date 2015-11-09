      subroutine fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,
     * tu,nu,tv,nv,p,c,nc,fp,fpu,fpv,mm,mvnu,spu,spv,right,q,au,au1,
     * av,av1,bu,bv,nru,nrv)
!  ..
!  ..scalar arguments..
      real p,fp
      integer ifsu,ifsv,ifbu,ifbv,idim,mu,mv,mz,nu,nv,nc,mm,mvnu
!  ..array arguments..
      real u(mu),v(mv),z(mz*idim),tu(nu),tv(nv),c(nc*idim),fpu(nu),
     * fpv(nv),spu(mu,4),spv(mv,4),right(mm*idim),q(mvnu),au(nu,5),
     * au1(nu,4),av(nv,5),av1(nv,4),bu(nu,5),bv(nv,5)
      integer ipar(2),nru(mu),nrv(mv)
!  ..local scalars..
      real arg,fac,term,one,half,value
      integer i,id,ii,it,iz,i1,i2,j,jz,k,k1,k2,l,l1,l2,mvv,k0,muu,
     * ncof,nroldu,nroldv,number,nmd,numu,numu1,numv,numv1,nuu,nvv,
     * nu4,nu7,nu8,nv4,nv7,nv8
!  ..local arrays..
      real h(5)
!  ..subroutine references..
!    fpback,fpbspl,fpdisc,fpbacp,fptrnp,fptrpe
!  ..
!  let
!               |   (spu)    |            |   (spv)    |
!        (au) = | ---------- |     (av) = | ---------- |
!               | (1/p) (bu) |            | (1/p) (bv) |
!
!                                | z  ' 0 |
!                            q = | ------ |
!                                | 0  ' 0 |
!
!  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
!                coefficients.
!       z      : the mu x mv matrix which contains the function values.
!       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
!                according to the least-squares problems in the u-,resp.
!                v-direction.
!       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
!                containing the discontinuity jumps of the derivatives
!                of the b-splines in the u-,resp.v-variable at the knots
!  the b-spline coefficients of the smoothing spline are then calculated
!  as the least-squares solution of the following over-determined linear
!  system of equations
!
!    (1)  (av) c (au)' = q
!
!  subject to the constraints
!
!    (2)  c(nu-3+i,j) = c(i,j), i=1,2,3 ; j=1,2,...,nv-4
!            if(ipar(1).ne.0)
!
!    (3)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
!            if(ipar(2).ne.0)
!
!  set constants
      one = 1
      half = 0.5
!  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      muu = mu
      if(ipar(1).ne.0) muu = mu-1
      mvv = mv
      if(ipar(2).ne.0) mvv = mv-1
!  it depends on the value of the flags ifsu,ifsv,ifbu  and ibvand
!  on the value of p whether the matrices (spu), (spv), (bu) and (bv)
!  still must be determined.
      if(ifsu.ne.0) go to 50
!  calculate the non-zero elements of the matrix (spu) which is the ob-
!  servation matrix according to the least-squares spline approximation
!  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 40 it=1,muu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tu,nu,3,arg,l,h)
        do 30 i=1,4
          spu(it,i) = h(i)
  30    continue
        nru(it) = number
  40  continue
      ifsu = 1
!  calculate the non-zero elements of the matrix (spv) which is the ob-
!  servation matrix according to the least-squares spline approximation
!  problem in the v-direction.
  50  if(ifsv.ne.0) go to 100
      l = 4
      l1 = 5
      number = 0
      do 90 it=1,mvv
        arg = v(it)
  60    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(tv,nv,3,arg,l,h)
        do 80 i=1,4
          spv(it,i) = h(i)
  80    continue
        nrv(it) = number
  90  continue
      ifsv = 1
 100  if(p.le.0.) go to  150
!  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 110
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
!  calculate the non-zero elements of the matrix (bv).
 110  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
!  substituting (2)  and (3) into (1), we obtain the overdetermined
!  system
!         (4)  (avv) (cr) (auu)' = (qq)
!  from which the nuu*nvv remaining coefficients
!         c(i,j) , i=1,...,nu-4-3*ipar(1) ; j=1,...,nv-4-3*ipar(2) ,
!  the elements of (cr), are then determined in the least-squares sense.
!  we first determine the matrices (auu) and (qq). then we reduce the
!  matrix (auu) to upper triangular form (ru) using givens rotations.
!  we apply the same transformations to the rows of matrix qq to obtain
!  the (mv) x nuu matrix g.
!  we store matrix (ru) into au (and au1 if ipar(1)=1) and g into q.
 150  if(ipar(1).ne.0) go to 160
      nuu = nu4
      call fptrnp(mu,mv,idim,nu,nru,spu,p,bu,z,au,q,right)
      go to 180
 160  nuu = nu7
      call fptrpe(mu,mv,idim,nu,nru,spu,p,bu,z,au,au1,q,right)
!  we determine the matrix (avv) and then we reduce this matrix to
!  upper triangular form (rv) using givens rotations.
!  we apply the same transformations to the columns of matrix
!  g to obtain the (nvv) x (nuu) matrix h.
!  we store matrix (rv) into av (and av1 if ipar(2)=1) and h into c.
 180  if(ipar(2).ne.0) go to 190
      nvv = nv4
      call fptrnp(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,c,right)
      go to 200
 190  nvv = nv7
      call fptrpe(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,av1,c,right)
!  backward substitution to obtain the b-spline coefficients as the
!  solution of the linear system    (rv) (cr) (ru)' = h.
!  first step: solve the system  (rv) (c1) = h.
 200  ncof = nuu*nvv
      k = 1
      if(ipar(2).ne.0) go to 240
      do 220 ii=1,idim
      do 220 i=1,nuu
         call fpback(av,c(k),nvv,5,c(k),nv)
         k = k+nvv
 220  continue
      go to 300
 240  do 260 ii=1,idim
      do 260 i=1,nuu
         call fpbacp(av,av1,c(k),nvv,4,c(k),5,nv)
         k = k+nvv
 260  continue
!  second step: solve the system  (cr) (ru)' = (c1).
 300  if(ipar(1).ne.0) go to 400
      do 360 ii=1,idim
      k = (ii-1)*ncof
      do 360 j=1,nvv
        k = k+1
        l = k
        do 320 i=1,nuu
          right(i) = c(l)
          l = l+nvv
 320    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 340 i=1,nuu
           c(l) = right(i)
           l = l+nvv
 340    continue
 360  continue
      go to 500
 400  do 460 ii=1,idim
      k = (ii-1)*ncof
      do 460 j=1,nvv
        k = k+1
        l = k
        do 420 i=1,nuu
          right(i) = c(l)
          l = l+nvv
 420    continue
        call fpbacp(au,au1,right,nuu,4,right,5,nu)
        l = k
        do 440 i=1,nuu
           c(l) = right(i)
           l = l+nvv
 440    continue
 460  continue
!  calculate from the conditions (2)-(3), the remaining b-spline
!  coefficients.
 500  if(ipar(2).eq.0) go to 600
      i = 0
      j = 0
      do 560 id=1,idim
      do 560 l=1,nuu
         ii = i
         do 520 k=1,nvv
            i = i+1
            j = j+1
            q(i) = c(j)
 520     continue
         do 540 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 540     continue
 560  continue
      ncof = nv4*nuu
      nmd = ncof*idim
      do 580 i=1,nmd
         c(i) = q(i)
 580  continue
 600  if(ipar(1).eq.0) go to 700
      i = 0
      j = 0
      n33 = 3*nv4
      do 660 id=1,idim
         ii = i
         do 620 k=1,ncof
            i = i+1
            j = j+1
            q(i) = c(j)
 620     continue
         do 640 k=1,n33
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 640     continue
 660  continue
      ncof = nv4*nu4
      nmd = ncof*idim
      do 680 i=1,nmd
         c(i) = q(i)
 680  continue
!  calculate the quantities
!    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
!    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
!    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
!                  tu(r+3) <= u(i) <= tu(r+4)
!    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
!                  tv(r+3) <= v(j) <= tv(r+4)
 700  fp = 0.
      do 720 i=1,nu
        fpu(i) = 0.
 720  continue
      do 740 i=1,nv
        fpv(i) = 0.
 740  continue
      nroldu = 0
!  main loop for the different grid points.
      do 860 i1=1,muu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        iz = (i1-1)*mv
        do 840 i2=1,mvv
          numv = nrv(i2)
          numv1 = numv+1
          iz = iz+1
!  evaluate s(u,v) at the current grid point by making the sum of the
!  cross products of the non-zero b-splines at (u,v), multiplied with
!  the appropiate b-spline coefficients.
          term = 0.
          k0 = numu*nv4+numv
          jz = iz
          do 800 id=1,idim
          k1 = k0
          value = 0.
          do 780 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 760 l2=1,4
              k2 = k2+1
              value = value+fac*spv(i2,l2)*c(k2)
 760        continue
            k1 = k1+nv4
 780      continue
!  calculate the squared residual at the current grid point.
          term = term+(z(jz)-value)**2
          jz = jz+mz
          k0 = k0+ncof
 800      continue
!  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 820
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 820      nroldv = numv
          if(numu.eq.nroldu) go to 840
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 840    continue
        nroldu = numu
 860  continue
      return
      end
