 module module_common
        Use DFLIB
        use IFQWIN
        type node
            integer  i
            integer  j
            type(node),pointer   :: next
        end type node
        type(node),pointer    ::  head  ,p , last
        integer,parameter           :: n = 87,m=100,totaltimes=int(2E6)
        integer,parameter           ::  freq = int(1E3),k=5  
        integer*2,parameter       :: ux = 20,uy=20,dx=1400,dy=750
        real*8,parameter            :: r = 5,r0=12
        integer  state(n,m),nei1(4,2),nei2(6,2),num1,num2,i1,i2
        real*8   xcoor1,ycoor1,xcoor2,ycoor2
        real  x1,x2,x3,t1,t2
        character*20  filename
        integer  itmp,jtmp
        integer  junit
        integer  ::    border((m+n)*2,2)    ! for the border site
        real     desorb
        
        
    ! /* ------ Definitions of Globally Used Variables for Graphic Settings ------- */ !
        Integer g_xLeft, g_yLeft, g_xRight, g_yRight    ! g --> Global
        Integer g_jColor
        Integer g_xLeft1, g_yLeft1, g_xRight1, g_yRight1    ! g --> Global
        integer  jcolor,ierror,i,j

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
