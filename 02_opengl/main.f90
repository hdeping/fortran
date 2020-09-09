program main
    use module_common
    implicit none

    real   coor(3)



    call initializeAWindowPlease()
    call glclearcolor(0.0,0.0,0.0,0.0)
    call glclear(gl_color_buffer_bit)
    call glcolor3f(1.0,1.0,1.0)
    call glortho(0.0,1.0,0.0,1.0,- 1.0,1.0)
    call glbegin(gl_polygon)
        coor = (/0.25,0.75,0.00/)
        glvertex3f(coor(1),coor(2),coor(3))
        coor = (/0.75,0.25,0.00/)
        glvertex3f(coor(1),coor(2),coor(3))
        coor = (/0.75,0.75,0.00/)
        glvertex3f(coor(1),coor(2),coor(3))
        coor = (/0.25,0.75,0.00/)
        glvertex3f(coor(1),coor(2),coor(3))
    call glend()
    call glflush()
    call updateTheWindowAndCheckForEvents()
end program main
