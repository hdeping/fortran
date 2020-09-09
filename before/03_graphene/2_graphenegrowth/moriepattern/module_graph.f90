module graph_module
use msflib


contains

subroutine initial_window()
integer(4)            ::        results
integer            ::        result
type(qwinfo)    ::        winfo

winfo.type = qwin$max
results = setwsizeqq(qwin$framewindow, winfo)
results = setwsizeqq(0, winfo)

call setviewport(100, 100, 610, 610)
result = setwindow(.true., dble(-4.0), dble(-4.0), dble(4.0), dble(4.0))
result = setbkcolorrgb(#ffffff)
result = settextcolorrgb(#000000)
call clearscreen($gclearscreen)

end subroutine initial_window

!--*******************************--!

subroutine graph_boundary()
integer(2)            ::        status
type(wxycoord)        ::        wxy

status = setcolorrgb(#0000ff)


        call moveto_w(-4.0d0, -0.2d0, wxy)
        status = lineto_w(-3.8d0, -0.18d0)


end subroutine graph_boundary

!--*******************************--!


end module  graph_module
