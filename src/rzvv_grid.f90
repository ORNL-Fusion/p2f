module rzvv_grid
    use constants

    implicit none

    save
    
    real :: R_range, z_Range, vPer_range, vPar_range
    real :: R_binSize, z_binSize, vPer_binSize, vPar_binSize
    real, allocatable :: R_binEdges (:), &
            z_binEdges (:), &
            vPer_binEdges (:), &
            vPar_binEdges (: )
    real, allocatable :: R_binCenters (:), &
            z_binCenters (:), &
            vPer_binCenters (:), &
            vPar_binCenters (: )
    real :: z_min, z_max, R_min, R_max

    real, allocatable :: f_rzvv (:,:,:,:), f_rzvv_global(:,:,:,:)
    real, allocatable :: density (:,:)
    real :: normFac_
    real(kind=dbl) :: v_sigma, x_sigma

contains
    subroutine init_rzvv_grid ()
        use read_particle_list
        use eqdsk
        use constants
        use read_namelist
        use init_mpi
        implicit none

        integer :: i
        i = 0

        allocate ( R_binEdges ( R_nBins + 1 ), &
                   z_binEdges ( z_nBins + 1 ), &
                   vPer_binEdges ( vPer_nBins + 1 ), &
                   vPar_binEdges ( vPar_nBins + 1 ), &
                   R_binCenters ( R_nBins ), &
                   z_binCenters ( z_nBins ), &
                   vPer_binCenters ( vPer_nBins ), &
                   vPar_binCenters ( vPar_nBins ) )
        allocate ( f_rzvv ( R_nBins, z_nBins, vPer_nBins, vPar_nBins ) )
        allocate ( density ( R_nBins, z_nBins ) )

        R_min   = minVal ( rbbbs )
        R_max   = maxVal ( rbbbs )    
        R_range = R_max - R_min
        R_binSize   = R_range / R_nBins
        R_binEdges  = (/ (i*R_binSize,i=0,R_nBins) /) + R_min
        R_binCenters    = R_binEdges(1:R_nBins) + R_binSize / 2.0

        z_min   = minVal ( zbbbs ) 
        z_max = maxVal ( zbbbs ) 
        z_range = ( z_max - z_min )
        z_binSize   = z_range / z_nBins
        z_binEdges  = (/ (i*z_binSize,i=0,z_nBins) /) + z_min
        z_binCenters    = z_binEdges(1:z_nBins) + z_binSize / 2.0

        if(mpi_pId==1) write(*,*) 'R_min: ', R_min
        if(mpi_pId==1) write(*,*) 'R_max: ', R_max
        if(mpi_pId==1) write(*,*) 'z_min: ', z_min
        if(mpi_pId==1) write(*,*) 'z_max: ', z_max

        vPer_range = vPerInRange / 100.0 * c 
        vPer_binSize   = vPer_range / vPer_nBins
        vPer_binEdges   = (/ (i*vPer_binSize,i=0,vPer_nBins) /)
        vPer_binCenters    = vPer_binEdges(1:vPer_nBins) &
            + vPer_binSize / 2.0

        vPar_range  = vParInRange / 100.0 * c
        vPar_binSize    = 2.0 * vPar_range / vPar_nBins
        vPar_binEdges   = (/ (i*vPar_binSize,i=0,vPar_nBins) /) - vPar_range
        vPar_binCenters = vPar_binEdges(1:vPar_nBins) + vPar_binSize / 2.0 

        f_rzvv  = 0.0
 
    end subroutine init_rzvv_grid

end module rzvv_grid
