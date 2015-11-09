Module Event_Module

Integer jfg_pause

Contains

  Subroutine PauseProgram( jUnit, jEvent, jKeyState, ix, iy )
  USE MSFLIB
  Integer    jUnit, jEvent, jKeyState, ix, iy

  j1  =  jUnit 
  j2  = jEvent 
  j3  = jKeyState 
  j4  = ix 
  j5  = iy
  
  jfg_pause  =  1
  
  End Subroutine PauseProgram

  Subroutine PauseProgram2( jUnit, jEvent, jKeyState, ix, iy )
  USE MSFLIB
  Integer    jUnit, jEvent, jKeyState, ix, iy

  j1  =  jUnit 
  j2  = jEvent 
  j3  = jKeyState 
  j4  = ix 
  j5  = iy
  
  jfg_pause  =  2
  
  End Subroutine PauseProgram2
    
    Subroutine SaveSpiral2File( NN, xyAry, nSpiral, jf_pause )
    Integer    NN, nSpiral, jf_Pause
    Real       xyAry( NN, NN )

    Character  fName * 20

    write( 6, * ) ' 0: Exit the program 1: Save the Spiral ? 2: Continue the Program'
    10  Read( 5, *, err = 10 ) jopt
    if( jopt < 0 .OR. jopt > 2 ) goto 10
    if( jopt  = = 0 ) STOP
    if( jopt  = = 2 ) goto 100
        
    nSpiral  =  nSpiral + 1
    write( fName, 1000 ) nSpiral 
    1000 Format( 'Spiral_', i2.2, '.dat' )

    open( 101, file = fname, status='unknown' )
    1100 Format( <NN>( f12.6 ) )
    Do ii  =  1, NN
       Write( 101, 1100 )  xyAry( ii,:) 
    Enddo
    close( 101 )
    
    100 jf_pause  =  0

    Return
    End Subroutine SaveSpiral2File
End Module Event_Module  
