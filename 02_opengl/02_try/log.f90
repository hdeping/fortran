module_graph.f90:35:28:

     call gllinestipple( 1 , #aaaa )            !...������
                            1
Error: Syntax error in argument list at (1)
module_graph.f90:145:29:

         drawpos(ip,:)=vk(ip)%p(:)
                             1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:148:31:

         call glvertex3f(vk(ip)%p(1),vk(ip)%p(2),0.0d0)
                               1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:152:34:

 end subroutine draw_vkticle_active
                                  1
Error: Expected label ‘draw_vkticle_simple’ for END SUBROUTINE statement at (1)
module_graph.f90:155:0:

 subroutine draw_vkticle_active( radius , nrd , lenv )
 1
Error: Unclassifiable statement at (1)
module_graph.f90:156:30:

     real               :: lenv
                              1
Error: Symbol ‘lenv’ at (1) already has basic type of REAL
module_graph.f90:157:32:

     real               :: radius
                                1
Error: Symbol ‘radius’ at (1) already has basic type of REAL
module_graph.f90:158:31:

     integer              :: nrd
                               1
Error: Symbol ‘nrd’ at (1) already has basic type of INTEGER
module_graph.f90:159:29:

     real               :: nsd
                             1
Error: Symbol ‘nsd’ at (1) already has basic type of REAL
module_graph.f90:160:32:

     real               :: thetav
                                1
Error: Symbol ‘thetav’ at (1) already has basic type of REAL
module_graph.f90:161:33:

     real               :: drawpos(nsum,2)
                                 1
Error: Symbol ‘drawpos’ at (1) already has basic type of REAL
module_graph.f90:162:35:

     real                 :: cnormal(3) , cnormal1(3)
                                   1
Error: Symbol ‘cnormal’ at (1) already has basic type of REAL
module_graph.f90:163:33:

     real                 :: ctest(3)  , ctest1(3)
                                 1
Error: Symbol ‘ctest’ at (1) already has basic type of REAL
module_graph.f90:164:30:

     real               :: vacc(3)
                              1
Error: Symbol ‘vacc’ at (1) already has basic type of REAL
module_graph.f90:178:29:

         drawpos(ip,:)=vk(ip)%p(:)
                             1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:186:56:

         call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                                                        1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:188:56:

         call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                                                        1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:190:67:

         call glvertex3f(drawpos(ip,1)+0.8*radius*cos(pi/15.+vk(ip)%ag),&
                                                                   1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:192:56:

         call glvertex3f(drawpos(ip,1)+radius*cos(vk(ip)%ag),&
                                                        1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:194:68:

         call glvertex3f(drawpos(ip,1)+0.8*radius*cos(-pi/15.+vk(ip)%ag),&
                                                                    1
Error: Symbol ‘vk’ at (1) has no IMPLICIT type
module_graph.f90:128:34:

     real               :: drawpos(nsum,2)
                                  1
Error: Variable ‘nsum’ cannot appear in the expression at (1)
module_graph.f90:209:24:

         call glvertex3f(lbox,0.0,0.0)
                        1
Error: Type mismatch in argument ‘x’ at (1); passed INTEGER(4) to REAL(4)
module_graph.f90:210:24:

         call glvertex3f(lbox,lbox,0.0)
                        1
Error: Type mismatch in argument ‘x’ at (1); passed INTEGER(4) to REAL(4)
module_graph.f90:211:28:

         call glvertex3f(0.0,lbox,0.0)
                            1
Error: Type mismatch in argument ‘y’ at (1); passed INTEGER(4) to REAL(4)
module_graph.f90:41:38:

         call glvertex3f(1.*ix*discell,lbox,0.0)
                                      1
Error: Type mismatch in argument ‘y’ at (1); passed INTEGER(4) to REAL(4)
module_graph.f90:48:24:

         call glvertex3f(lbox,1.*ix*discell,0.0)
                        1
Error: Type mismatch in argument ‘x’ at (1); passed INTEGER(4) to REAL(4)
module_graph.f90:83:30:

             call gltranslatef(dble(i*vstep), dble(j*vstep), 0.0)
                              1
Error: Type mismatch in argument ‘x’ at (1); passed REAL(8) to REAL(4)
module_graph.f90:223:20:

     call gluortho2d(-0.1, dble(lbox)+0.1, -0.1, dble(lbox)+0.1)
                    1
Error: Type mismatch in argument ‘left’ at (1); passed REAL(4) to REAL(8)
module_graph.f90:112:20:

     call glcolor3fv(loc(c1))
                    1
Error: Type mismatch in argument ‘v’ at (1); passed INTEGER(8) to REAL(4)
main.f90:2:8:

     use graph_module
        1
Fatal Error: Can't open module file ‘graph_module.mod’ for reading at (1): No such file or directory
compilation terminated.
make: *** [compile] Error 1
