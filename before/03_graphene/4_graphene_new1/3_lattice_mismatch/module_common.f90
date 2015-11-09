 module module_common
        Use DFLIB
        use IFQWIN
        integer,parameter           :: ux = 200,uy=100,dx=1400,dy=800
        integer,parameter           :: n = 20,m=40
        integer,parameter           :: freq = int(1E5),k=5  
        real*8,parameter               :: r1 = 14,r2=5,pi=3.1415927
        real*8   hexa(n,m,2),tri(n-1,m-1,2),x0,y0
        integer  nei(2,2)
        integer  inew,jnew,snew,s
        integer  num(2),a(4),b(2,3),coor(2)
        integer  count,state(n,m),times
        integer  jcolor,ierror,i,j,icolor

        real*8  x1,x2,x3,x,y
        Type (WindowConfig) thisScreen
        Type (QWInfo)       theFrame, theChild 
        type(wxycoord) xy
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
                                    
    
    end module module_common
