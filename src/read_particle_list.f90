module read_particle_list
    implicit none
    save

    real, allocatable :: p_R(:), p_z(:), &
        p_vPer(:), p_vPar(:), p_weight(:) 
   
contains    
    subroutine read_pl ()
        use constants
        use init_mpi
        !use mpi
        use read_namelist
        use netcdf
        use dlg

        integer(kind=LONG) :: i
        real, dimension(nP) :: p_v

        real, allocatable :: p_psi(:), p_theta(:), p_E(:), p_lambda(:)

        real, allocatable :: nc_R(:), nc_z(:), nc_vPer(:), nc_vPar(:), nc_weight(:) 
        integer :: nc_nP, nc_id, R_id, z_id, vPer_id, vPar_id, weight_id, R_dim_ids(1)
        integer :: plEnvLength, plEnvStatus 
        character(len=100) :: plEnvVar, plEnvFName

        allocate ( p_R(nP), p_z(nP), p_weight(nP), &
            p_vPer(nP), p_vPar(nP) )

            if ( pList_is_nCDF ) then 

            !   Check environment variable for scripting
            
            plEnvVar = 'PL_FNAME'
            call get_environment_variable ( name = plEnvVar, &
                value = plEnvFName, &
                length = plEnvLength, &
                status = plEnvStatus )

            if ( plEnvStatus .ne. 0 ) plEnvFName = pl_fileName 

            if ( mpi_pId .eq. 0 ) write (*,*) 'Using ', plEnvFName

            !   Read in variables from .nc particle list file

            call dlg_check ( nf90_open ( path = plEnvFName, mode = nf90_nowrite, ncid = nc_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'R', R_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'z', z_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'vPer', vPer_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'vPar', vPar_id ) )
            call dlg_check ( nf90_inq_varId ( nc_id, 'weight', weight_id ) )
            
            call dlg_check ( nf90_inquire_variable ( &
                nc_id, R_id, dimIds = R_dim_ids ) )
            
            call dlg_check ( nf90_inquire_dimension ( nc_id, R_dim_ids(1), &
                len = nc_nP ) ) 
            
            allocate ( nc_R ( nc_nP ), &
                nc_z ( nc_nP ), &
                nc_vPer ( nc_nP ), &
                nc_vPar ( nc_nP ), &
                nc_weight ( nc_nP ) )

            call dlg_check ( nf90_get_var ( nc_id, R_id, nc_R ) )
            call dlg_check ( nf90_get_var ( nc_id, z_id, nc_z ) )
            call dlg_check ( nf90_get_var ( nc_id, vPer_id, nc_vPer ) )
            call dlg_check ( nf90_get_var ( nc_id, vPar_id, nc_vPar ) )
            call dlg_check ( nf90_get_var ( nc_id, weight_id, nc_weight ) )
            
            call dlg_check ( nf90_close ( nc_id ) )

            p_R = nc_R(1:nP)
            p_z = nc_z(1:nP)
            p_vPer = nc_vPer(1:nP)
            p_vPar  = nc_vPar(1:nP)
            p_weight    = nc_weight(1:nP)

        else if ( pList_is_nubeam1 ) then
            
            !   read in variables from nubeam particle list

            allocate ( p_E(nP), p_lambda(nP) )
           
            open ( unit = 8, file = pl_fileName, status = 'OLD', action = 'READ' )

            do i=1,nP
               
                read (8,3000) p_R(i), p_z(i), p_E(i), p_lambda(i), p_weight(i)

            end do 

            !3000 format ( 7e14.6 )
            3000 format ( 7e16.6 )


            close ( unit = 8 )

            !   Convert R,z to [m] from [cm]

            p_R = p_R * 1e-2
            p_z = p_z * 1e-2

            !   Calculate particle vPer and vPar

            p_v = sqrt ( 2.0 * p_E * 1.0e3 * e_ / mi ) 
            p_vPar  = p_v * p_lambda
            p_vPer = sqrt ( p_v**2 - p_vPar**2 )

            deallocate ( p_E, p_lambda )

        else

            !   Read in variables from .dav particle list file

            allocate ( p_psi(nP), p_theta(nP), &
                p_E(nP), p_lambda(nP) )
           
            open ( unit = 8, file = pl_fileName, status = 'OLD', action = 'READ' )

            do i=1,nP
               
                read (8,2000) p_psi(i), p_theta(i), &
                    p_R(i), p_z(i), p_E(i), p_lambda(i), p_weight(i)

            end do 

            2000 format ( 7e16.6 )

            close ( unit = 8 )

            !   Convert R,z to [m] from [cm]

            p_R = p_R * 1e-2
            p_z = p_z * 1e-2

            !   Calculate particle vPer and vPar

            p_v = sqrt ( 2.0 * p_E * 1.0e3 * e_ / mi ) 
            p_vPar  = p_v * p_lambda
            p_vPer = sqrt ( p_v**2 - p_vPar**2 )

            deallocate ( p_psi, p_E, p_lambda, p_theta )

        end if

    end subroutine read_pl

end module read_particle_list
