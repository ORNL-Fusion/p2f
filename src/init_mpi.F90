module init_mpi
#if PARALLEL==1
    use mpi 
#endif
    use read_namelist
    implicit none
    save
 
    integer :: mpi_iErr, mpi_nProcs, mpi_pId, mpi_nP
    integer :: mpi_start_, mpi_end_, mpi_fh
#if PARALLEL==1
    integer :: mpi_status_(MPI_STATUS_SIZE)
#endif
    integer :: nP_wall, nP_bad, nP_off_vGrid, &
        nP_badWeight, nP_badEnergy, nP_TookMaxStepsBeforeBounce
    
contains
    subroutine start_mpi ()
        implicit none
   
#if PARALLEL==1 
        call mpi_init ( mpi_iErr )
        call mpi_comm_size ( MPI_COMM_WORLD, mpi_nProcs, mpi_iErr )
        call mpi_comm_rank ( MPI_COMM_WORLD, mpi_pId, mpi_iErr )
    
        mpi_nP = nP / mpi_nProcs
        mpi_start_ = mpi_nP * mpi_pId + 1
        mpi_end_ = mpi_start_ + mpi_nP - 1
#else
        mpi_nP = nP
        mpi_start_ = 1
        mpi_end_ = nP
        mpi_pId = 1
#endif
        nP_wall = 0
        nP_bad  = 0
        nP_off_vGrid    = 0
        nP_badWeight    = 0
        nP_badEnergy    = 0
        nP_TookMaxStepsBeforeBounce = 0

        if(mpi_pId==1) write (*,*) 'Division: ', mpi_start_, mpi_end_
        
    end subroutine start_mpi
    
    subroutine stop_mpi ()
        implicit none
       
#if PARALLEL==1 
        call mpi_finalize ( mpi_iErr )
#endif

    end subroutine stop_mpi 

end module init_mpi