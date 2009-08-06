module constants
    use read_namelist
    implicit none
    save
    
    real, parameter :: e_ = 1.602176e-19
    real :: mi
    real, parameter :: c = 3.0e8
    real, parameter :: q = 1.0 * e_
    real, parameter :: pi = 3.14159265
   
contains
    subroutine set_constants 
        implicit none 
 
        mi = amu * 1.67262e-27

    end subroutine set_constants

end module constants
