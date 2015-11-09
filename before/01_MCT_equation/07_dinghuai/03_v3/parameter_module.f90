module parameter_module
    integer,parameter                     :: n   = 100
    integer,parameter                     :: sn  = 200
    integer,parameter                     :: ini = 50
    integer,parameter                     :: stn = 2    
    real(8),parameter                     :: pi  = 3.1415926
    real(8),dimension(0:n-1)              :: sk
    real(8),dimension(0:n-1)              :: ck
    real(8),dimension(0:n-1)              :: sk_neq            
    real(8),dimension(0:n-1)              :: au
    real(8),dimension(0:n-1)              :: bu
    real(8),dimension(0:n-1)              :: aus
    real(8),dimension(:,:),allocatable    :: phie
    real(8),dimension(:,:),allocatable    :: dphie
    real(8),dimension(:,:),allocatable    :: ke_mee
    real(8),dimension(0:n-1)              :: phia
    real(8)                               :: dk
    real(8)                               :: dr
    real(8)                               :: h
    real(8)                               :: v0
    real(8)                               :: di
    real(8)                               :: rho
    real(8)                               :: mf
    real(8)                               :: det_time
    real(8)                               :: df
    integer                               :: istep
    
end module parameter_module
