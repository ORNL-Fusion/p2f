!>  p2f is a program for converting particle lists from the ORBIT-rf code into
!!  continuous distribution functions for use in AORSA, amoung other things. 
!<

program p2f
    use eqdsk
    use gc_terms
    use interp
    use gc_integrate 
    use read_particle_list
    use rzvv_grid
    use init_mpi
    use constants
    !use plot
    use write_f_rzvv
    use read_namelist
 
    implicit none
    
    integer(kind=LONG) :: i,j,k
    real :: T1, T2
    integer :: mpi_count
    integer :: nP_wall_total, nP_bad_total, &
        nP_off_vGrid_total, nP_badWeight_total, nP_badEnergy_total
    real :: tmpDensity !< Test description

    call init_namelist () !< Read in namelist variables from p2f.nml
    call set_constants ()
    call read_geqdsk ( eqdsk_fileName, plot = .false. )
    call bCurvature ()
    call bGradient ()
    call init_interp ()
    call start_mpi ()
    call read_pl ()    
    call init_rzvv_grid ()

    call cpu_time (T1)

    do i=mpi_start_,mpi_end_
 
        !call gc_orbit ( 2.04280, -0.115082, 218696.93, -8985448.2 )

        if ( mod ( i, 100 ) == 0 .and. mpi_pId == 0 ) &
            write(*,*) nP, mpi_nP, i, p_R(i), p_z(i), p_vPerp(i), p_vPar(i)

        if ( p_weight(i) > 0 ) &
            call gc_orbit ( p_R(i), p_z(i), p_vPerp(i),&
                p_vPar(i), p_weight(i) , plot = plotOrbit )

    end do


    call mpi_barrier ( MPI_COMM_WORLD, mpi_iErr )

    mpi_count   = R_nBins * z_nBins * vPerp_nBins * vPar_nBins

    allocate ( f_rzvv_global ( R_nBins, z_nBins, vPerp_nBins, vPar_nBins ) )
    call mpi_reduce ( f_rzvv, f_rzvv_global, mpi_count, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )

    f_rzvv = f_rzvv_global
    deallocate(f_rzvv_global)

    call mpi_reduce ( nP_wall, nP_wall_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )
    call mpi_reduce ( nP_bad, nP_bad_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )
    call mpi_reduce ( nP_off_vGrid, nP_off_vGrid_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )
    call mpi_reduce ( nP_badWeight, nP_badWeight_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )
    call mpi_reduce ( nP_badEnergy, nP_badEnergy_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_iErr )

    ! Test the reduction operation
    write (*,*) 'Stopping MPI' 

    call stop_mpi () 

    write(*,*) 'Getting CPU time'   
    call cpu_time (T2)
 
    if ( mpi_pId == 0 ) then


        write (*,*) 'Time taken: ', T2-T1 
        write (*,'(a,f5.2,a)') 'Wall:      ', real ( nP_wall_total ) / real ( nP ) * 100.0, '%'
        write (*,'(a,f5.2,a)') 'Bad:       ', real ( nP_bad_total ) / real ( nP ) * 100.0, '%'
        write (*,'(a,f5.2,a,a)') 'off_vGrid: ', real ( nP_off_vGrid_total ) / real ( nP ) * 100.0, '%', &
            '   *** only applicable for gParticle = .false.'
        write (*,'(a,f5.2,a)') 'badWeight: ', real ( nP_badWeight_total ) / real ( nP ) * 100.0, '%'
        write (*,'(a,f5.2,a)') 'badEnergy: ', real ( nP_badEnergy_total ) / real ( nP ) * 100.0, '%'
        write (*,'(a,f5.1,a)') 'Suggested eNorm: ', &
            real ( maxVal ( vPerp_binCenters )**2 * 0.5 * mi / e_ / 1e3 ), ' keV'

        !   Integrate over velocity space to get a number for
        !   density for sanity check ;-)

        density = 0
        do i = 1, R_nBins
            do j = 1, z_nBins
    
                tmpDensity  = 0
                do k = 1, vPerp_nBins

                    tmpDensity = tmpDensity + &
                        sum ( f_rzvv(i,j,k,:) ) * vPerp_binCenters(k)

                end do
                
                density(i,j)    = tmpDensity * &
                    2.0 * pi * vPerp_binSize * vPar_binSize

            end do
        end do

        write(*,'(a,e8.2)') 'max ( density ): ', maxVal ( density )

        call write_f ()

    end if

end program p2f
