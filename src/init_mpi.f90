module init_mpi
    use mpi 
    use read_namelist
    implicit none
    save
 
    integer :: mpi_iErr, mpi_nProcs, mpi_pId, mpi_nP
    integer :: mpi_start_, mpi_end_, mpi_fh
    integer :: mpi_status(MPI_STATUS_SIZE)
    !integer, parameter :: LONG = selected_int_kind(9)
    !integer(kind=LONG), parameter :: nP = 10000
    integer :: nP_wall, nP_bad, nP_off_vGrid, &
        nP_badWeight, nP_badEnergy
    
contains
    subroutine start_mpi ()
        implicit none
    
        call mpi_init ( mpi_iErr )
        call mpi_comm_size ( MPI_COMM_WORLD, mpi_nProcs, mpi_iErr )
        call mpi_comm_rank ( MPI_COMM_WORLD, mpi_pId, mpi_iErr )
    
        mpi_nP  = nP / mpi_nProcs
        mpi_start_   = mpi_nP * mpi_pId + 1
        mpi_end_ = mpi_start_ + mpi_nP - 1

        nP_wall = 0
        nP_bad  = 0
        nP_off_vGrid    = 0
        nP_badWeight    = 0
        nP_badEnergy    = 0

        write (*,*) 'Division: ', mpi_start_, mpi_end_
        
        !write(*,*) 'Proc id: ', mpi_pId, ' np range: ', mpi_start_, mpi_end_
     
    end subroutine start_mpi
    
    subroutine stop_mpi ()
        implicit none
        
        call mpi_finalize ( mpi_iErr )

    end subroutine stop_mpi 

end module init_mpi
