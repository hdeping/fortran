program main
    use module_common
    implicit none
    real(8)        :: coltim(n)
    real(8)        :: rx(n)
    real(8)        :: ry(n)
    real(8)        :: timBig
    real(8)        :: rxij
    real(8)        :: ryij
    real(8)        :: diameter    

    call random_seed()
    
    do i = 1,n
        coltim(i) = timBig
    enddo !cycle ends

    siasq = a ** 2.0
    do i = 1,n - 1
        do j = i + 1,n
            rxij = rx(i) - ry(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)
            rxij = anInt(rxij)
            ryij = anInt(ryij)
            rzij = anInt(rzij)
            vxij = vx(i) - vy(j)
            vyij = vy(i) - vy(j)
            vzij = vz(i) - vz(j)
            bij  = rxij*vxij + ryij*vyij + rzij*vzij
            if ( bij < 0.0 )then
                rijsq = rxij ** 2.0 + ryij ** 2.0 +  rzij ** 2.0
                vijsq = vxij ** 2.0 + vyij ** 2.0 +  vzij ** 2.0
                discr = bij ** 2.0 - vijsq*(rijsq - sigsq)
                if ( discr > 0.0 )then
                    tij = (- bij - sqrt(discr))/vijsq
                    if ( tij < coltim(i) )then
                        coltim(i) = tij
                        partnr(i) = j
                    endif ! if ends
                endif ! if ends
            endif ! if ends
        enddo !cycle ends
    enddo !cycle ends

    tij = timbig
    do k = 1,n
        if ( coltim(k) < tij )then
            tij = coltim(k)
            i   = k
        endif ! if ends
        
    enddo !cycle ends
    do k = 1,n
        coltim(k)  = coltim(k) - tij
        rx(k)      = rx(k) + vk(k)*tij
        ry(k)      = ry(k) + vk(k)*tij
        rz(k)      = rz(k) + vk(k)*tij
        rx(k)      = rx(k) - anInt(rx(k)) 
        ry(k)      = ry(k) - anInt(ry(k)) 
        rz(k)      = rz(k) - anInt(rz(k)) 
    enddo !cycle ends
    
    ii = i
    jj = j

    do i = 1,n
        if ( i == ii .or. partnr(i) == ii .or. & 
            i == jj .or. partnr(i) == jj)then
            coltim(i) = timBig
            do j  = 1,n
                if ( j /= i )then
                    ! usual calculation for ij collision
                    if ( tij < coltim(i) )then
                        coltim(i) = tij
                        partnr(i) = j
                    endif ! if ends
                    if ( tij < coltim(j) )then
                        coltim(j) = tij
                        partnr(j) = i
                    endif ! if ends
                    
                endif ! if ends
                
            enddo !cycle ends
             
        endif ! if ends
    enddo !cycle ends
     

     
     

end program main
