PROGRAM randomspher_freeglut_main
   USE ISO_C_BINDING
   USE OpenGL_Example
   IMPLICIT NONE
   
   TYPE(SpinningSphere), TARGET :: sphere
   
   CALL TestGL(sphere)

END PROGRAM randomspher_freeglut_main
