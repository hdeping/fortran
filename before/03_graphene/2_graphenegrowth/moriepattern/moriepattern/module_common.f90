 module module_common
        Use DFLIB
        use IFQWIN
        integer,parameter           ::  freq = int(1E5)
        integer*2,parameter       :: ux = 100,uy=50,dx=1400,dy=750
        real*8,parameter            :: x0 = 700,y0=400
         real*8,parameter           :: height = 4 , length=600, pi=3.1415927
        real*8   xcoor1,ycoor1,xcoor2,ycoor2
        real*8  x1,x2,x3,y1,y2,y3,theta
        character*20  filename
        integer  t1,t2,t3        
        integer  jcolor,ierror,i,j,itmp

        Type (WindowConfig) thisScreen
        Type (QWInfo)       theFrame, theChild 
        type(wxycoord) xy
        type(xycoord) intxy
    ! /* ------ Define my Own Color Palette for 16 Colors -------------- */ !
        Integer :: myColor(17)  =  (/#000000,   &        ! Black                /* --- 1
                                    #FFFFFF,   &        ! Bright White        /* --- 2
                                       #000080,   &        ! DULL Red            /* --- 3
                                    #0000FF,   &        ! Bright Red        /* --- 4
                                    #008000,   &        ! Dull Green        /* --- 5
                                    #00FF00,   &        ! Bright Green        /* --- 6
                                    #008080,   &        ! Dull Yellow        /* --- 7
                                    #00FFFF,   &        ! Bright Yellow        /* --- 8
                                    #800000,   &        ! Dull Blue            /* --- 9
                                    #FF0000,   &        ! Bright Blue        /* --- 10
                                    #800080,   &        ! Dull Magenta        /* --- 11
                                    #FF00FF,   &        ! Bright Magenta    /* --- 12
                                    #808000,   &        ! Dull Turquoise    /* --- 13
                                    #FFFF00,   &        ! Bright Turquoise    /* --- 14
                                    #808080,   &        ! Dark Gray            /* --- 15
                                    #C0C0C0,   &      ! Bright Gray        /* --- 16
                                    #C8E8C8 /)          ! Comfortable Green    /* ___17

!***************************************************************
    end module module_common

    