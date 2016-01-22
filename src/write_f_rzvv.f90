module write_f_rzvv
    
contains
   
    subroutine write_f ()
        use rzvv_grid
        use netcdf
        use init_mpi
        use read_particle_list
        use read_namelist
        use dlg
        implicit none
    
        character ( len = 103 ) :: ncFileName
        integer :: nc_id, &
            R_nBins_id, &
            z_nBins_id, &
            vPer_nBins_id, &
            vPar_nBins_id, &
            f_rzvv_id, &
            R_binCenters_id, &
            z_binCenters_id, &
            vPer_binCenters_id, &
            vPar_binCenters_id, &
            R_nBins_p1_id, &
            z_nBins_p1_id, &
            vPer_nBins_p1_id, &
            vPar_nBins_p1_id, &
            R_binEdges_id, &
            z_binEdges_id, &
            vPer_binEdges_id, &
            vPar_binEdges_id, &
            scalar_id, &
            nP_id, &
            R_binSize_id, &
            z_binSize_id, &
            vPer_binSize_id, vPar_binSize_id, &
            vPer_range_id, vPar_range_id, &
            density_id

        integer :: plEnvLength, plEnvStatus 
        character(len=100) :: plEnvVar, DataDirChar
        logical :: DataDirExists

        plEnvVar = 'PL_FNAME'
        call get_environment_variable ( name = plEnvVar, &
            value = ncFileName, &
            length = plEnvLength, &
            status = plEnvStatus )

        inquire (file='data', direct=DataDirChar, exist=DataDirExists)
        if(.not.DataDirExists)then
            call system('mkdir data')
        endif
            
        ncFileName  = 'data/p2f_' // trim ( ncFileName(6:26) )

        if ( plEnvStatus .ne. 0 ) ncFileName = 'data/fdis.dav.nc' 

        !ncFileName = trim(pl_fileName) // '.nc'
        !ncFileName = 'data/fdis.dav.nc'

        call dlg_check ( nf90_create ( ncFileName, nf90_clobber, nc_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "R_nBins", R_nBins, R_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "z_nBins", z_nBins, z_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPer_nBins", vPer_nBins, vPer_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPar_nBins", vPar_nBins, vPar_nBins_id ) )

        call dlg_check ( nf90_def_dim ( nc_id, "R_nBins_p1", R_nBins+1, R_nBins_p1_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "z_nBins_p1", z_nBins+1, z_nBins_p1_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPer_nBins_p1", vPer_nBins+1, vPer_nBins_p1_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPar_nBins_p1", vPar_nBins+1, vPar_nBins_p1_id ) )

        call dlg_check ( nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )

        call dlg_check ( nf90_def_var ( nc_id, "f_rzvv", NF90_REAL, &
            (/ R_nBins_id, z_nBins_id, vPer_nBins_id, vPar_nBins_id /), &
            f_rzvv_id ) )
 
        call dlg_check ( nf90_def_var ( nc_id, "R_binCenters", NF90_REAL, &
            R_nBins_id, R_binCenters_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "z_binCenters", NF90_REAL, &
            z_nBins_id, z_binCenters_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPer_binCenters", NF90_REAL, &
            vPer_nBins_id, vPer_binCenters_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPar_binCenters", NF90_REAL, &
            vPar_nBins_id, vPar_binCenters_id ) ) 

        call dlg_check ( nf90_def_var ( nc_id, "R_binEdges", NF90_REAL, &
            R_nBins_p1_id, R_binEdges_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "z_binEdges", NF90_REAL, &
            z_nBins_p1_id, z_binEdges_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPer_binEdges", NF90_REAL, &
            vPer_nBins_p1_id, vPer_binEdges_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPar_binEdges", NF90_REAL, &
            vPar_nBins_p1_id, vPar_binEdges_id ) ) 

        call dlg_check ( nf90_def_var ( nc_id, "nP", NF90_INT, &
            scalar_id, nP_id ) )

        call dlg_check ( nf90_def_var ( nc_id, "R_binSize", NF90_REAL, &
            scalar_id, R_binSize_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "z_binSize", NF90_REAL, &
            scalar_id, z_binSize_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPer_binSize", NF90_REAL, &
            scalar_id, vPer_binSize_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPar_binSize", NF90_REAL, &
            scalar_id, vPar_binSize_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPar_range", NF90_REAL, &
            scalar_id, vPar_range_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "vPer_range", NF90_REAL, &
            scalar_id, vPer_range_id ) ) 
        call dlg_check ( nf90_def_var ( nc_id, "density", NF90_REAL, &
            (/ R_nBins_id, z_nBins_id /), density_id ) ) 


        call dlg_check ( nf90_enddef ( nc_id ) )

        call dlg_check ( nf90_put_var ( nc_id, nP_id, nP ) )
        call dlg_check ( nf90_put_var ( nc_id, f_rzvv_id, f_rzvv ) )
        call dlg_check ( nf90_put_var ( nc_id, R_binCenters_id, R_binCenters ) )
        call dlg_check ( nf90_put_var ( nc_id, z_binCenters_id, z_binCenters ) )
        call dlg_check ( nf90_put_var ( nc_id, vPer_binCenters_id, vPer_binCenters ) )
        call dlg_check ( nf90_put_var ( nc_id, vPar_binCenters_id, vPar_binCenters ) )
        call dlg_check ( nf90_put_var ( nc_id, R_binEdges_id, R_binEdges ) )
        call dlg_check ( nf90_put_var ( nc_id, z_binEdges_id, z_binEdges ) )
        call dlg_check ( nf90_put_var ( nc_id, vPer_binEdges_id, vPer_binEdges ) )
        call dlg_check ( nf90_put_var ( nc_id, vPar_binEdges_id, vPar_binEdges ) )
        call dlg_check ( nf90_put_var ( nc_id, R_binSize_id, R_binSize ) )
        call dlg_check ( nf90_put_var ( nc_id, z_binSize_id, z_binSize ) )
        call dlg_check ( nf90_put_var ( nc_id, vPer_binSize_id, vPer_binSize ) )
        call dlg_check ( nf90_put_var ( nc_id, vPar_binSize_id, vPar_binSize ) )
        call dlg_check ( nf90_put_var ( nc_id, vPar_range_id, vPar_range ) )
        call dlg_check ( nf90_put_var ( nc_id, vPer_range_id, vPer_range ) )
        call dlg_check ( nf90_put_var ( nc_id, density_id, density ) )

        call dlg_check ( nf90_close ( nc_id ) )

    end subroutine write_f

end module write_f_rzvv
