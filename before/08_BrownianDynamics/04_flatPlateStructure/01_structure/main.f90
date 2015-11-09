program main
use module_common
    integer             :: oo
    integer             :: ikind
    real*8              :: packphi
    character(len=20)   :: filename

    call random_seed()
    
    call initial_parameter()
    
    
    do
        call evolution()
        if ( time > 5.0d0 )then
            exit
        endif ! if ends
    enddo

    filename = "init_55_18_18_3"
    open(10, file = filename)
    

    do ii = 1,pn
        write(10,"(i6,4f25.16)")ii,ra(ii),pp(ii,:)
    enddo !cycle ends
    open(10, file = filename)
