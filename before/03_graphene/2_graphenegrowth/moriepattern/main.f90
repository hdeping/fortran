program main
    use module_graph

   
    implicit none
   
    call initGraphWindow(mycolor(2))
   theta = 5
   
   !call setViewPort(50,50,800,800)
   
   do itmp = 1,900
    theta = itmp*0.1
     call lines_rot(150,theta)
    call lines(150)
    
    write(filename,'(i3.3,".bmp")')itmp
    jcolor = saveimage(filename,300,0,1100,800)
    call ClearScreen($Gviewport)
    end do 
   
  
    end program main