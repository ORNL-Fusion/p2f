module constants
    use read_namelist
    implicit none
    save
    
    real, parameter :: e_ = 1.602176e-19
    real :: mi
    real, parameter :: c = 3.0e8
    real :: q 
    real, parameter :: pi = 3.14159265
    integer, parameter :: dbl = selected_real_kind ( p=13, r = 200 ) 
   
contains
    subroutine set_constants 
        implicit none 
 
        mi = amu * 1.67262e-27
        q = AtomicZ * e_ 

    end subroutine set_constants

end module constants
