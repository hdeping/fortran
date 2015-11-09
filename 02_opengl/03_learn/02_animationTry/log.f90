module_common.f90:24:42:

     call glrectf(- 25.0, - 25.0,25.0,25.0,)
                                          1
Error: Syntax error in argument list at (1)
module_common.f90:59:8:

     case glut_left_button
        1
Error: Syntax error in CASE specification at (1)
module_common.f90:60:37:

         if ( state == glut_down )then
                                     1
Error: Expected a CASE or END SELECT statement following SELECT CASE at (1)
module_common.f90:61:45:

             call glutidlefunc_gl(spindisplay)
                                             1
Error: Expected a CASE or END SELECT statement following SELECT CASE at (1)
module_common.f90:62:11:

         endif ! if ends
           1
Error: Expecting END SELECT statement at (1)
module_common.f90:63:8:

     case glut_left_button   
        1
Error: Syntax error in CASE specification at (1)
module_common.f90:64:37:

         if ( state == glut_down )then
                                     1
Error: Expected a CASE or END SELECT statement following SELECT CASE at (1)
module_common.f90:65:38:

             call glutidlefunc_gl(null)
                                      1
Error: Expected a CASE or END SELECT statement following SELECT CASE at (1)
module_common.f90:66:11:

         endif ! if ends
           1
Error: Expecting END SELECT statement at (1)
module_common.f90:31:8:

     spin = spin + 2.0
        1
Error: Symbol ‘spin’ at (1) has no IMPLICIT type
module_common.f90:22:23:

     call glrotatef(spin,0.0,0.0,1.0)
                       1
Error: Symbol ‘spin’ at (1) has no IMPLICIT type
module_common.f90:43:20:

     call glviewport(0.0,w,h)
                    1
Error: Type mismatch in argument ‘x’ at (1); passed REAL(4) to INTEGER(4)
module_common.f90:46:17:

     call glortho(- 50.0,50.0,- 50.0,50.0,- 1.0,1.0)
                 1
Error: Type mismatch in argument ‘left’ at (1); passed REAL(4) to REAL(8)
module_common.f90:34:38:

         call glutpostwindowredisplay()
                                      1
Error: Missing actual argument for argument ‘window’ at (1)
module_common.f90:20:18:

     call glclear()
                  1
Error: Missing actual argument for argument ‘mask’ at (1)
main.f90:2:8:

     use module_common
        1
Fatal Error: Can't open module file ‘module_common.mod’ for reading at (1): No such file or directory
compilation terminated.
make: *** [compile] Error 1
