module eqdsk
    implicit none
    save    
    character (len=10) :: case_ ( 6 )
    integer :: idum, nw, nh, nbbbs, limitr
    real :: rdim, zdim, rcentr, rleft, zmid, &
        rmaxis, zmaxis, simag, sibry, bcentr, & 
        current, xdum, rStep, zStep, fStep
    
    real, allocatable :: fpol (:), psizr (:,:), &
        pres (:), ffprim (:), pprime (:), qpsi (:), &
        rbbbs (:), zbbbs (:), rlim (:), zlim (:), &
        r (:), z (:), bR(:,:), bPhi(:,:), bz(:,:), &
        fluxGrid (:), fpolRZ(:,:), bMag(:,:)

contains
    subroutine read_geqdsk ( eqdsk_fileName, plot )
        !use dislin
        use dlg
 
        implicit none
        
        integer :: i, j, iErr
        character(len=*), intent(IN) :: eqdsk_fileName
        logical, intent(IN), optional :: plot
        real, allocatable :: yp(:), xp(:), temp(:), ss(:), yp_c(:)
        real :: spl1 = 0.0, spln = 0.0, t, curv2, sigma = 0.0

        !!   Plotting variables
    
        !integer :: nLevs
        !real :: levStep, lev

        !   Read in variables from geqdsk file

        open ( unit = 8, file = eqdsk_fileName, status = 'OLD' )

        read ( 8, 2000 ) ( case_ (i), i=1, 6 ), idum, nw, nh 
        read ( 8, 2020 ) rdim,zdim,rcentr,rleft,zmid 
        read ( 8, 2020 ) rmaxis,zmaxis,simag,sibry,bcentr 
        read ( 8, 2020 ) current,simag,xdum,rmaxis,xdum 
        read ( 8, 2020 ) zmaxis,xdum,sibry,xdum,xdum 
        
        allocate ( fpol ( nw ), pres ( nw ), ffprim ( nw ), &
            pprime ( nw ), psizr ( nw, nh ), qpsi ( nw ), &
            r ( nw ), z ( nh ), fluxGrid ( nw ) )
        
        read ( 8, 2020 ) ( fpol (i), i=1, nw ) 
        read ( 8, 2020 ) ( pres (i), i=1, nw ) 
        read ( 8, 2020 ) ( ffprim (i), i=1, nw ) 
        read ( 8, 2020 ) ( pprime (i), i=1, nw ) 
        read ( 8, 2020 ) ( ( psizr (i,j), i=1, nw ), j=1, nh ) 
        read ( 8, 2020 ) ( qpsi (i), i=1, nw ) 
        
        read ( 8, 2022 ) nbbbs,limitr 
        
        allocate ( rbbbs ( nbbbs ), zbbbs ( nbbbs ), &
            rlim ( limitr ), zlim ( limitr ) )
        
        read ( 8, 2020 ) ( rbbbs (i), zbbbs (i), i=1,nbbbs ) 
        read ( 8, 2020 ) ( rlim (i), zlim (i), i=1,limitr ) 
         
        2000 format (6a8,3i4) 
        2020 format (5e16.9)
        2022 format (2i5) 
        
        close ( unit = 8 )
        
        !   Calculate other required variables
       
        rStep   = rdim / ( nw - 1 )
        zStep   = zdim / ( nh - 1 )
        fStep   = ( sibry - simag ) / ( nw - 1 )

        r   = (/ (i,i=0,nw-1) /) * rStep + rleft
        z   = (/ (i,i=0,nh-1) /) * zStep + zmid - zdim / 2.0

        fluxGrid    = (/ (i,i=0,nW-1) /) * fStep + simag

        allocate ( bR ( nw, nh ), bz ( nw, nh ), &
            bPhi ( nw, nh ), bMag(nw,nh) )

        bR  = dlg_pDeriv ( psizr, 2, zStep ) * 1.0
        bz  = -1.0 * dlg_pDeriv ( psizr, 1, rStep ) * 1.0

        !   Remember psi = - R * A

        do i=1,nw
            do j=1,nh
                bR(i,j)  = -bR(i,j) / r(i)
                bz(i,j)  = -bz(i,j) / r(i)
            enddo
        enddo
       
        allocate ( xp(nw), yp(nw), temp(nw), ss(nw), fpolRZ(nw,nh), yp_c(nw) ) 

        !   curv1 initialises the spline (fitpack.f)

!        call kurv1 ( nw, fluxGrid, fpol, spl1, spln, 3, xp, yp, temp, ss, 0.0, iErr )
        call curv1 ( nw, fluxGrid, fpol, spl1, spln, 3, yp_c, temp, sigma, iErr )
 
        do i=1,nw
            do j=1,nh

                !   curv2 evaluates the spline (fitpack.f)
               
                t   =  ( psizr(i,j) - simag ) / ( sibry - simag )
                !call kurv2 ( t, xs, ys, nw, fluxGrid, fpol, xp, yp, ss, 0.0 )
                !fPolRZ(i,j) = ys
                fPolRZ(i,j) = curv2 ( psizr(i,j), nw, fluxGrid, fpol, &
                    yp_c, sigma )
                bPhi(i,j)   = fpolRZ(i,j) / r(i)

            end do
        end do
      
!        !   Test the fitpack interpolation
!
!        do i=1,nw
!            
!            t   =  ( fluxGrid(i) - simag ) / ( sibry - simag )
!!            call kurv2 ( t, xs, ys, nw, fluxGrid, fpol, xp, yp, ss, 0.0 )
!            write (*,*) i, t, fluxGrid(i), fpol(i), xs, ys, &
!                curv2 ( fluxGrid(i), nw, fluxGrid, fpol, yp_c, 0.0 )
!         
!        end do

        bMag    = sqrt ( bR**2 + bPhi**2 + bz**2 )

        !!   Plotting

        !if ( present ( plot ) ) then
        !    if ( plot ) then 
        !        
        !        call setFil ( 'eqdsk.eps' ) 
        !        call setPag ( 'DA4P' )
        !        call metaFl ( 'EPS' )
        !        call disIni ()
        !        call graf ( 0.0, 3.0, 0.0, 0.5, -2.0, 2.0, -2.0, 0.5 ) 
        !        nLevs   = 101 
        !        levStep    = ( maxVal ( psizr ) - minVal ( psizr ) ) / nLevs
        !        do i=0,nLevs-1
        !            lev    = i * levStep - nLevs/2*levStep
        !            call contur ( r, nw, z, nh, psizr, lev ) 
        !        end do
        !        call endGrf ()  
        !        call disFin ()

        !    end if 
        !end if

    end subroutine read_geqdsk

end module eqdsk 
