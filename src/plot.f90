module plot

contains
    subroutine start_plot ()
        implicit none
    
        integer :: iStat, pgOpen
        iStat   = pgOpen ( 'p2f.ps/ps' )

    end subroutine start_plot 
    
    subroutine end_plot ()
        implicit none
    
        call pgClos () 

    end subroutine end_plot

    subroutine plot_f_vv ()
        use rzvv_grid 
        use constants
        implicit none
       
        integer :: nLevs, i=0
        real :: tr(6)
        real, allocatable :: levels(:)

        call pgEnv ( -1.0, 1.0, 0.0, 1.0, 1, 1 )
        call pgLab ( 'vPar [%c]', 'vPerp [%c]', 'f(vv)' )
        tr  = (/ ( minVal ( vPar_binCenters ) - vPar_binSize ) / c * 100.0, &
                vPar_binSize / c * 100.0, &
                0.0, &
                ( minVal ( vPerp_binCenters ) - vPerp_binSize ) / c * 100.0, &
                0.0, &
                vPerp_binSize / c * 100.0 /)

        nLevs   = 20
        allocate ( levels ( nLevs ) )
        levels  = (/ ((2.0**i)*1e-16,i=0,nLevs-1) /)

        call pgCont ( transpose ( f_vv ), vPar_nBins, vPerp_nBins, &
            1, vPar_nBins, 1, vPerp_nBins, &
            levels, nLevs, tr )

    end subroutine plot_f_vv

    subroutine plot_f_rzvv ( rPt, zPt )
        use rzvv_grid 
        use constants
        implicit none
       
        integer :: nLevs, i=0, rPt, zPt
        real :: tr(6)
        real, allocatable :: levels(:)

        call pgEnv ( -1.0, 1.0, 0.0, 1.0, 1, 1 )
        call pgLab ( 'vPar [%c]', 'vPerp [%c]', 'f(rzvv)' )
        tr  = (/ ( minVal ( vPar_binCenters ) - vPar_binSize ) / c * 100.0, &
                vPar_binSize / c * 100.0, &
                0.0, &
                ( minVal ( vPerp_binCenters ) - vPerp_binSize ) / c * 100.0, &
                0.0, &
                vPerp_binSize / c * 100.0 /)

        nLevs   = 20
        allocate ( levels ( nLevs ) )
        levels  = (/ ((2.0**i)*1e-15,i=0,nLevs-1) /)

        call pgCont ( transpose ( f_rzvv_(rPt,zPt,:,:) ), vPar_nBins, vPerp_nBins, &
            1, vPar_nBins, 1, vPerp_nBins, &
            levels, nLevs, tr )

    end subroutine plot_f_rzvv


end module plot
