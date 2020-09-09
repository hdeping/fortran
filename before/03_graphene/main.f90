program main
    use module_graph
    use module_lattice
   
    implicit none
    integer   num(4),output(m)
    
    call RANDOM_SEED()
    call initgraphwindow(mycolor(17))
    call setviewport(ux,uy,dx,dy)
   state = 0
   call drawlines()
   call  latticeinterface()
    
    
    
    end program main