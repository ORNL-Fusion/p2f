module gc_integrate

contains
    function dlg_gc_velocity ( pos, vPer, vPar, NFO )
        use gc_terms
        use interp
        use eqdsk
        implicit none

        real :: dlg_gc_velocity(3)
        real, intent(IN) :: pos(3), vPer, vPar
        real :: grad_R, grad_phi, grad_z, &
            curv_R, curv_phi, curv_z, &
            unitb_R, unitb_phi, unitb_z
        real :: surf2, bMagHere, bHere(3), &
            vgc_R, vgc_phi, vgc_z
        integer, optional, intent(in) :: NFO

        grad_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_R, nw, zp_bGrad_R, sigma )
        grad_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_phi, nw, zp_bGrad_phi, sigma )
        grad_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_z, nw, zp_bGrad_z, sigma )

        curv_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_R, nw, zp_bCurv_R, sigma )
        curv_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_phi, nw, zp_bCurv_phi, sigma )
        curv_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_z, nw, zp_bCurv_z, sigma )

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        unitb_R = bHere(1) / bMagHere
        unitb_phi   = bHere(2) / bMagHere
        unitb_z = bHere(3) / bMagHere

        if ( present ( NFO ) ) then
 
            vgc_R   = vPar * unitb_R 
            vgc_phi   = vPar * unitb_phi 
            vgc_z   = vPar * unitb_z 

        else
            
            vgc_R   = vPar * unitb_R + vPer**2 * grad_R + vPar**2 * curv_R 
            vgc_phi   = vPar * unitb_phi + vPer**2 * grad_phi + vPar**2 * curv_phi
            vgc_z   = vPar * unitb_z + vPer**2 * grad_z + vPar**2 * curv_z

        endif

        dlg_gc_velocity = (/ vgc_R, vgc_phi, vgc_z /)

    end function dlg_gc_velocity

    function dlg_vPer ( pos, u )
        use constants
        implicit none

        real :: bHere(3),bMagHere
        real, intent(IN) :: pos(3), u
        real :: dlg_vPer

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        dlg_vPer   = sqrt ( 2.0 * u * bMagHere / mi )

    end function dlg_vPer

    function dlg_vPar ( pos, u )
        use constants
        use gc_terms
        use interp
        use eqdsk
        implicit none
        
        real :: dlg_vPar, bDotGradB_here, surf2
        real, intent(IN) :: pos(3), u

        bDotGradB_here  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bDotGradB, nw, zp_bDotGradB, sigma )

        dlg_vPar    = -u / mi * bDotGradB_here 

    end function dlg_vPar

    function dlg_interpB ( pos, bMagHere )
        use eqdsk
        use interp
        implicit none
        
        real :: bR_here, bPhi_here, bz_here
        real, intent(IN) :: pos(3)
        real :: dlg_interpB(3)
        real, optional, intent(OUT) :: bMagHere
        real :: surf2

        bR_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bR, nw, zp_bR, sigma )
        bPhi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bPhi, nw, zp_bPhi, sigma )
        bz_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bz, nw, zp_bz, sigma )

        if ( present (bMagHere) ) &
            bMagHere    = sqrt ( bR_here**2 + bPhi_here**2 + bz_here**2 )

        dlg_interpB(1)  = bR_here 
        dlg_interpB(2)  = bPhi_here
        dlg_interpB(3)  = bz_here 
    
    end function dlg_interpB

    subroutine SetNormFac()
        use init_mpi
        use read_namelist, nGV=>nGridPtsGaussianVSpace, nGX=>nGridPtsGaussianXSpace
        use rzvv_grid
        use constants

        implicit none

        real :: GetNormFac
        real(kind=dbl) :: f_vv_update(vPer_nBins,vPar_nBins)
        real(kind=dbl) :: dV(vPer_nBins,vPar_nBins), dV_, this_dV, this_f
        integer :: i,j,m,n
        real :: RCenter, ZCenter

        ! Calculate the normalization factor of a particle 

        normFac_ = 0.0

        if(gParticle)then
    
            v_sigma = particleSize * c
 
            !   Calculate normFac for full particle size,
            !   No bessel function required here since the
            !   vPer offset = 0, i.e., besselI(0,0,)=1.0
  
            f_vv_update = 0.0 
            do i=1,vPer_nBins
                do j=1,vPar_nBins
                
                    f_vv_update(i,j) =  exp ( &
                            - ( vPer_binCenters(i)**2  &
                            + vPar_binCenters(j)**2 ) &
                            / ( 2d0 * v_sigma**2 ) &
                            ) * 2d0 * pi  
                            
                    dV(i,j) = vPer_binCenters(i) * &
                        4.0 * pi**2 * R_binSize * z_binSize * &
                        vPer_binSize * vPar_binSize  
                end do
            end do
            
            normFac_ = sum ( f_vv_update * dV )  

        else if(gParticle4D)then
    
            x_sigma = particleSizeX
            v_sigma = particleSize * c
 
            !do m=max(r_nBins/2-nGX,1),min(r_nBins/2+nGX,r_nBins)
            !do n=max(z_nBins/2-nGX,1),min(z_nBins/2+nGX,z_nBins)
            !do i=max(vPer_nBins/2-nGV,1),min(vPer_nBins/2+nGV,vPer_nBins)
            !do j=max(vPar_nBins/2-nGV,1),min(vPar_nBins/2+nGV,vPar_nBins)
            do m=1,r_nBins
            do n=1,z_nBins
            do i=1,vPer_nBins
            do j=1,vPar_nBins
               
                RCenter = (R_max-R_min)/2 + R_min
                ZCenter = (Z_max-Z_min)/2 + Z_min
     
                this_f =  &
                    exp ( - ( vPer_binCenters(i)**2  + vPar_binCenters(j)**2 ) / ( 2d0 * v_sigma**2 ) ) &
                  * exp ( - ( (RCenter-r_binCenters(m))**2  + (ZCenter-z_binCenters(n))**2 ) / ( 2d0 * x_sigma**2 ) ) 
                        
                this_dV = r_binCenters(m) * vPer_binCenters(i) * &
                    4.0 * pi**2 * R_binSize * z_binSize * &
                    vPer_binSize * vPar_binSize  

                normFac_ = normFac_ + this_f*this_dV

            enddo
            enddo
            enddo
            enddo

        endif

    end subroutine SetNormFac

    subroutine gc_orbit ( start_R, start_z, start_vPer, &
         start_vPar, weight, plot )
        use eqdsk
        use gc_terms
        use constants

        use interp
        use rzvv_grid
        use init_mpi
        use read_namelist, nGV=>nGridPtsGaussianVSpace, nGX=>nGridPtsGaussianXSpace
        implicit none
       
        real, intent(IN) :: start_R, start_z,&
            start_vPer, start_vPar, &
            weight
        real :: weightMod
        integer(kind=LONG) :: stepCnt, var_dt, i, j, ii, ss
        real, allocatable, dimension(:) :: rTrack, zTrack, &
            vPerTrack, vParTrack, distance, dtArray, &
            R_index, z_index, &
            vPer_index, vPar_index, &
            rTrack_nfo, zTrack_nfo, vPerTrack_nfo, &
            vParTrack_nfo
        real :: tau, dTau, pos(3), &
            u, bMagHere, bHere(3), vPer, vPar, &
            vgc(3), vPer_n, vPar_n
        real(kind=dbl) :: dt, dtMin, dtMax
        logical :: stillIn, firstOrbit
        logical, optional, intent(IN) :: plot
        real :: k1_vPar, k1_vgc(3), k2_vPar, k2_vgc(3), &
            k3_vPar, k3_vgc(3), k4_vPar, k4_vgc(3), lastStep, &
            psi_here, surf2, vgc_n(3), k1_vgc_n(3), k2_vgc_n(3), &
            k3_vgc_n(3), k4_vgc_n(3), pos_(3), k1_vPar_n, &
            k2_vPar_n, k3_vPar_n, k4_vPar_n
        
        integer, allocatable, dimension(:) :: vPerL, vPerR, &
            vParL, vParR, R_L, R_R, &
            Z_L, Z_R

        real(kind=dbl) :: f_vv_update(vPer_nBins,vPar_nBins)
	    real(kind=dbl), allocatable :: f_rzvv_update(:,:,:,:)
        real(kind=dbl) :: dV(vPer_nBins,vPar_nBins), dV_, this_dV, this_f

        real(kind=dbl) :: normFac, bArg
        real(kind=dbl) :: bessi, bF, result_D, expTerm_D, expTermX

        !   if you change this you need to alter
        !   the initialisation of phi down further
        integer, parameter :: nPhi = 21         
        
        real :: phi(nPhi), phiBase(nPhi)
        real :: phiRange, phi0, vX(nPhi), vY(nPhi), &
            vZ(nPhi), vX0, vY0, vz0

        real :: theta, start_theta, theta_old, &
            theta_diff, theta_diff_old
        logical :: nearStart, skip_dtUpdate
        real :: EStep, energy, RCenter, ZCenter
        integer :: stride, m, n, i_, j_, m_, n_
       
        !   Initialize variables
      
        allocate(rTrack(MaxSteps))
        allocate(zTrack(MaxSteps), &
            vPerTrack(MaxSteps), vParTrack(MaxSteps), distance(MaxSteps), dtArray(MaxSteps), &
            R_index(MaxSteps), z_index(MaxSteps), &
            vPer_index(MaxSteps), vPar_index(MaxSteps), &
            rTrack_nfo(MaxSteps), zTrack_nfo(MaxSteps), vPerTrack_nfo(MaxSteps), &
            vParTrack_nfo(MaxSteps))
        allocate(vPerL(MaxSteps), vPerR(MaxSteps), &
            vParL(MaxSteps), vParR(MaxSteps), R_L(MaxSteps), R_R(MaxSteps), &
            Z_L(MaxSteps), Z_R(MaxSteps))

        if(gParticle4D)then
            allocate(f_rzvv_update(2*nGX+1,2*nGX+1,2*nGV+1,2*nGV+1))
        endif

        stepCnt = 0
        firstOrbit  = .true.
        stillIn = .true.
        dTau    = 0.0
        tau = 0.0
        var_dt  = 1
        skip_dtUpdate   = .false.
        weightMod   = weight
        
        energy  = 0.5 * mi * ( start_vPer**2 + start_vPar**2 ) / e_ * 1d-3!  [keV]

        if ( energy <= energyThreshold ) then 
            nP_badEnergy = nP_badEnergy + 1
            weightMod = 0.0
        endif

        if ( weight > weightLimit ) then
            nP_badWeight = nP_badWeight + 1
            weightMod = 0.0
        endif

        if (.not.DistributeAlongOrbit) then

            NFO = .false.
            firstOrbit = .false.
            rTrack(1) = start_R
            zTrack(1) = start_z
            vPerTrack(1) = start_vPer
            vParTrack(1) = start_vPar
            tau    = 1.0
            dtArray(1)  = 1.0
            stepCnt = 1
            go to 11 ! yes, i know, shutUp ;-)

        end if

10      continue

        pos(1)  = start_R
        pos(2)  = 0.0
        pos(3)  = start_z

        start_theta = aTan2 ( start_z - zmAxis, start_R -rmAxis ) * 180.0 / pi

        vPer = start_vPer
        vPar = start_vPar
        
        if ( NFO ) then 

            pos_    = pos
            vPar_n  = vPar
            vPer_n = vPer

        endif

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )

        u   = mi * vPer**2 / ( 2.0 * bMagHere ) 

        dtMin   = 1.0e-9
        if ( var_dt == 1 ) then
            dtMax   = 5.0e-7
        else 
            dtMin   = 0.0005e-7
            dtMax   = 0.05e-7
        end if
        
        dt  = dtMin

       do 
  
            ! Include finite orbits

12          continue
 
            if ( dt < dtMin ) dt = dtMin
            if ( dt > dtMax ) dt = dtMax
 
            vPer   = dlg_vPer ( pos, u ) 
            vgc = dlg_gc_velocity ( pos, vPer, vPar )
            k1_vPar   = dt * dlg_vPar ( pos, u ) 
            k1_vgc  = dt * vgc
            
            vPer   = dlg_vPer ( pos + k1_vgc / 2.0, u ) 
            vgc = dlg_gc_velocity ( pos + k1_vgc / 2.0, vPer, vPar + k1_vPar / 2.0 )
            k2_vPar   = dt * dlg_vPar ( pos + k1_vgc / 2.0, u ) 
            k2_vgc  = dt * vgc
           
            vPer   = dlg_vPer ( pos + k2_vgc / 2.0, u ) 
            vgc = dlg_gc_velocity ( pos + k2_vgc / 2.0, vPer, vPar + k2_vPar / 2.0 )
            k3_vPar   = dt * dlg_vPar ( pos + k2_vgc / 2.0, u ) 
            k3_vgc  = dt * vgc
            
            vPer   = dlg_vPer ( pos + k3_vgc, u ) 
            vgc = dlg_gc_velocity ( pos + k3_vgc, vPer, vPar + k3_vPar )
            k4_vPar   = dt * dlg_vPar ( pos + k3_vgc, u ) 
            k4_vgc  = dt * vgc
           
            vPar    = vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar &
                + k4_vPar ) / 6.0
           
            pos   = pos + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0

            if ( NFO ) then
           
            ! No finite orbits 
 
            vPer_n   = dlg_vPer ( pos_, u ) 
            vgc_n = dlg_gc_velocity ( pos_, vPer_n, vPar_n, NFO = 1 )
            k1_vPar_n   = dt * dlg_vPar ( pos_, u ) 
            k1_vgc_n  = dt * vgc_n
            
            vPer_n   = dlg_vPer ( pos_ + k1_vgc_n / 2.0, u ) 
            vgc_n = dlg_gc_velocity ( pos_ + k1_vgc_n / 2.0, vPer_n, vPar_n + k1_vPar_n / 2.0, NFO = 1 )
            k2_vPar_n   = dt * dlg_vPar ( pos_ + k1_vgc_n / 2.0, u ) 
            k2_vgc_n  = dt * vgc_n
           
            vPer_n   = dlg_vPer ( pos_ + k2_vgc_n / 2.0, u ) 
            vgc_n = dlg_gc_velocity ( pos_ + k2_vgc_n / 2.0, vPer_n, vPar_n + k2_vPar_n / 2.0, NFO = 1 )
            k3_vPar_n   = dt * dlg_vPar ( pos_ + k2_vgc_n / 2.0, u ) 
            k3_vgc_n  = dt * vgc_n
            
            vPer_n   = dlg_vPer ( pos_ + k3_vgc_n, u ) 
            vgc_n = dlg_gc_velocity ( pos_ + k3_vgc_n, vPer_n, vPar_n + k3_vPar_n, NFO = 1 )
            k4_vPar_n   = dt * dlg_vPar ( pos_ + k3_vgc_n, u ) 
            k4_vgc_n  = dt * vgc_n
           
            vPar_n    = vPar_n + ( k1_vPar_n + 2.0 * k2_vPar_n + 2.0 * k3_vPar_n &
                + k4_vPar_n ) / 6.0
           
            pos_   = pos_ + ( k1_vgc_n + 2.0 * k2_vgc_n + 2.0 * k3_vgc_n + k4_vgc_n ) / 6.0

            endif
            
            !   Interpolate to get psi at new pos
            !   and if psi(pos) is outside the last 
            !   closed flux surface then dump the particle
            !   The 0.98 is a factor suggested by EFJ to keep
            !   away boundary of the flux grid.

            psi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
                psizr, nw, zp_psi, sigma )

            if ( .not. ascending_flux ) then
                if ( psi_here < sibry * 0.999 ) stillIn = .false.
            else
                if ( psi_here > sibry * 0.999 ) stillIn = .false.
            endif

            if ( (.not. stillIn) .or. (.not. firstOrbit) .or. stepCnt >= maxSteps ) then

                if ( stepCnt+1 >= maxSteps ) then
                    if ( var_dt == 1 ) then 
                        
                        !   The algorithm seems to trace the particles well,
                        !   including finding the orbit end points. Howevere,
                        !   there seem to be two exceptions that I have found so
                        !   far, these being particles whose orbits track very
                        !   near the plasma boundary and those that are on the
                        !   trapped passing boundary. Reducing the step size
                        !   does not seem to fix the trajectories near the
                        !   plasma boundary while the trapped/passing boundary
                        !   orbits are rare.
                           
                        stepCnt = 0
                        var_dt = 0
                        go to 10 
                    else
                        nP_bad  = nP_bad + 1
                        if(FirstOrbit) nP_TookMaxStepsBeforeBounce = nP_TookMaxStepsBeforeBounce + 1
                    end if
                end if
                
                if ( .not. stillIn ) nP_wall = nP_wall + 1

                exit

            end if

            stepCnt = stepCnt + 1

            rTrack(stepCnt) = pos(1)
            zTrack(stepCnt) = pos(3)
            vPerTrack(stepCnt) = vPer
            vParTrack(stepCnt)  = vPar

            if ( NFO ) then 
                rTrack_nfo(stepCnt) = pos_(1)
                zTrack_nfo(stepCnt) = pos_(3)
                vPerTrack_nfo(stepCnt) = vPer_n
                vParTrack_nfo(stepCnt)  = vPar_n
            endif

            tau = tau + dt
            dtArray(stepCnt)    = dt
          
            if ( NFO ) then  
        
                distance(stepCnt)   = sqrt ( ( start_R - pos_(1) )**2 &
                    + ( start_z - pos_(3) )**2 )
            else

                 distance(stepCnt)   = sqrt ( ( start_R - pos(1) )**2 &
                    + ( start_z - pos(3) )**2 )

            endif
            
            !   Detect end of orbit

            theta_old   = theta
            theta   = aTan2 ( zTrack(stepCnt) - zmAxis, rTrack(stepCnt) -rmAxis ) * 180.0 / pi &
                - start_theta
            if ( theta < 0 ) theta = theta + 360.0
            
            if ( stepCnt > 20 ) then

                if ( NFO ) then
       
                    ! NFO end of orbit is small distance plus same sign
                    ! of vPar as starting vPar.
                     
                    if ( distance(stepCnt-1) < 0.0006 &
                        .and. vPar_n / start_vPar > 0 ) firstOrbit = .false.
 
                else

                    theta_diff_old  = theta_diff
                    theta_diff  = theta_old - theta
                    nearStart = .false.
                    if ( theta < 2.5 .or. theta > 357.5 ) nearStart = .true.
                    if ( vParTrack(stepCnt) / start_vPar > 0 .and. &
                        start_vPar / abs ( start_vPar ) * theta_diff > 0 .and. &
                        sin ( theta * pi / 180.0 ) < 0.045 .and. nearStart ) &
                            firstOrbit = .false.  

                endif

            end if

            if ( firstOrbit .and. stepCnt > 1 ) then 

                !   Adjust time step 

                if ( NFO ) then 

                    lastStep    = sqrt ( ( rTrack_nfo(stepCnt) - rTrack_nfo(stepCnt-1) )**2 &
                        + ( zTrack_nfo(stepCnt) - zTrack_nfo(stepCnt-1) )**2 )
                else

                    lastStep    = sqrt ( ( rTrack(stepCnt) - rTrack(stepCnt-1) )**2 &
                    + ( zTrack(stepCnt) - zTrack(stepCnt-1) )**2 )

                endif

                if ( lastStep .eq. lastStep .and. lastStep > 0 .and. &
                        ( .not. skip_dtUpdate ) ) &
                    dt  = traceStepLength / lastStep * dt

            end if 

        end do

11      continue

        !   Introduce a stride over which to select points out
        !   of a track since we need a higher resolution track for
        !   tracing than we do for bounce averaging. Remember to 
        !   shorten the tau by the same factor.

        stride  = int ( stepCnt / usePtsPerTrack )
        if ( stride < 1 ) stride = 1
        if (.not.DistributeAlongOrbit) stride = 1

#if DEBUG_LEVEL>0 
        if(mpi_pId==0) write(*,*) "Stride : ", stride
#endif
        
        tau = tau / stride
        dtArray = dtArray / stride

        !   Calculate f_rzvv indices and add this orbit
        !   to the f_rzvv particle count. At this point the
        !   f_rzvv function will only be the number of particles
        !   but since its a regular grid we can divide by the
        !   volume element after counting all the particles ;-)

        if ( NFO ) then

            !   Overwrite tracks with NFO versions

            rTrack  = rTrack_nfo
            zTrack  = ztrack_nfo
            vPerTrack  = vPerTrack_nfo
            vParTrack   = vParTrack_nfo

        endif

        if ( stillIn ) then 
        
            R_index = ( rTrack - r_min ) / R_range * R_nBins + 1
            z_index = ( zTrack - z_min ) / z_range * z_nBins + 1
            vPer_index = vPerTrack / vPer_range * vPer_nBins + 1
            vPar_index  = ( vParTrack + vPar_range ) / ( 2.0 * vPar_range ) &
                * vPar_nBins + 1

            if (gParticle) then

                if(vPer_nBins<64.or.vPar_nBins<128)then
                    write(*,*) 'ERROR : When using gParticle = .true. you need a fine v-space grid'
                    write(*,*) 'Please set ...'
                    write(*,*) 'vPer_nBins = 64'
                    write(*,*) 'vPar_nBins = 128'
                    stop
                endif

               
                !   Loop only over a few cells around the particle location, 
                !   otherwise this is prohibitavely slow. Its already a factor
                !   of 10 slower than using the delta functions above :-(
    
                vPerL = int(vPer_index)-nGV
                where (vPerL < 1) vPerL = 1
                vPerR  = int(vPer_index)+nGV
                where (vPerR > vPer_nBins) vPerR = vPer_nBins
     
                vParL = int(vPar_index)-nGV
                where (vParL < 1) vParL = 1
                vParR  = int(vPar_index)+nGV
                where (vParR > vPar_nBins) vParR = vPar_nBins
     
                R_L = int(R_index)!-nGX
                where (R_L < 1) R_L = 1
                where (R_L > R_nBins) R_L = R_nBins

                Z_L = int(z_index)!-nGX
                where (Z_L < 1) Z_L = 1
                where (Z_L > z_nBins) Z_L = z_nBins


                !   Setup the base phi grid in case the 
                !   uPer offset is too large for the bessel
                !   function.

                phiBase = (/ (i*1.0,i=-10,10,1) /) / 10.0
 
                do ii=1,stepCnt,stride          

                normFac = normFac_ * rTrack(ii) ! not sure why this Jacobian was here?
   
                f_vv_update = 0.0

                do i=vPerL(ii),vPerR(ii)
                    do j=vParL(ii),vParR(ii)

                        !   If we end up with Infinity type nonsense poping up
                        !   in the distribution function I would suggest
                        !   altering this step so it does not consider particles
                        !   located outside the chosen vPer/vPar range.
 
                        bArg    =  vPerTrack(ii) * vPer_binCenters(i) / v_sigma**2

                        if ( bArg < 650 ) then

                            bF = bessI ( 0, real(bArg,dbl) )
                       
                            expTerm_D   = exp ( real ( &
                                    - ( vPerTrack(ii)**2 + ( vParTrack(ii) - vPar_binCenters(j) )**2  &
                                    + vPer_binCenters(i)**2 ) &
                                    / ( 2.0 * v_sigma**2 ) &
                                    , dbl ) ) 

                            result_D    = expTerm_D * bF * 2d0 * real(pi,dbl)
                            f_vv_update(i,j)  = result_D

                            !   This is a catch for when the modified bessel
                            !   function returns +Infinity

                            if ( bF > 1d300 ) then
                              
                                f_vv_update(i,j)  = 0.0 
 
                            end if

                        else    ! do the gyro angle integral the old fashioned way

                            !   but use an adaptive range for the phi
                            !   integral depending on the value of vPerTrack.
                            !   Oh yes, this is brilliant!

                            phiRange    = 2.5 * v_sigma / vPerTrack(ii)     
                            phi = phiBase * phiRange
                            phi0    = 0.0

                            vX  = vPer_binCenters(i) * cos ( phi )
                            vY  = vPer_binCenters(i) * sin ( phi )
                            vZ  = vPar_binCenters(j)

                            vX0 = vPerTrack(ii) * cos ( phi0 )
                            vY0 = vPerTrack(ii) * sin ( phi0 )
                            vZ0 = vParTrack(ii)

                            f_vv_update(i,j)   = sum ( &
                                exp ( -1.0 * &
                                      ( ( vX - vX0 )**2 + ( vY - vY0 )**2 + ( vZ - vZ0 )**2 ) & 
                                            / ( 2.0 * v_sigma**2 ) ) &
                                ) * abs ( phi(1) - phi(2) )

                        end if
 
                    end do
                end do
               
                !   Normalise and add

                if ( normFac == 0 .or. tau == 0 ) then
                
                    write (*,*) 'DLG: ERROR, normFac == 0 or tau == 0 '
                    write (*,*) normFac, normFac_, tau, v_sigma**2

                else

                    f_vv_update   = f_vv_update / normFac * weightMod

                    f_rzvv(R_L(ii),Z_L(ii),vPerL(ii):vPerR(ii),vParL(ii):vParR(ii))  = &
                        f_rzvv(R_L(ii),Z_L(ii),vPerL(ii):vPerR(ii),vParL(ii):vParR(ii)) &
                        + dtArray(ii) / tau &
                        * f_vv_update(vPerL(ii):vPerR(ii),vParL(ii):vParR(ii)) 

                endif 

                end do 

            else if(gParticle4D) then

                if(vPer_nBins<64.or.vPar_nBins<128)then
                    write(*,*) 'ERROR : When using gParticle = .true. you need a fine v-space grid'
                    write(*,*) 'Please set ...'
                    write(*,*) 'vPer_nBins = 64'
                    write(*,*) 'vPar_nBins = 128'
                    stop
                endif

                !   Loop only over a few cells around the particle location, 
                !   otherwise this is prohibitavely slow. Its already a factor
                !   of 10 slower than using the delta functions above :-(
    
                vPerL = int(vPer_index)-nGV
                where (vPerL < 1) vPerL = 1
                vPerR  = int(vPer_index)+nGV
                where (vPerR > vPer_nBins) vPerR = vPer_nBins
     
                vParL = int(vPar_index)-nGV
                where (vParL < 1) vParL = 1
                vParR  = int(vPar_index)+nGV
                where (vParR > vPar_nBins) vParR = vPar_nBins
     
                R_L = int(R_index)-nGX
                where (R_L < 1) R_L = 1
                R_R  = int(R_index)+nGX
                where (R_R > R_nBins) R_R = R_nBins

                Z_L = int(z_index)-nGX
                where (Z_L < 1) Z_L = 1
                Z_R  = int(z_index)+nGX
                where (Z_R > z_nBins) Z_R = z_nBins

                !   Calculate normFac for full particle size,
                !   No bessel function required here since the
                !   vPer offset = 0, i.e., besselI(0,0,)=1.0
  
                f_rzvv_update = 0.0 
             
                !   Setup the base phi grid in case the 
                !   uPer offset is too large for the bessel
                !   function.

                phiBase = (/ (i*1.0,i=-10,10,1) /) / 10.0
 
                do ii=1,stepCnt,stride          

                normFac = normFac_ * rTrack(ii) ! not sure why this Jacobian was here?
   
                f_rzvv_update = 0.0

                do m=R_L(ii),R_R(ii)
                do n=Z_L(ii),Z_R(ii)
                do i=vPerL(ii),vPerR(ii)
                do j=vParL(ii),vParR(ii)

                    m_ = m-R_L(1)+1
                    n_ = n-Z_L(1)+1
                    i_ = i-vPerL(1)+1
                    j_ = j-vParL(1)+1

                    !   If we end up with Infinity type nonsense poping up
                    !   in the distribution function I would suggest
                    !   altering this step so it does not consider particles
                    !   located outside the chosen vPer/vPar range.

                    bArg = vPerTrack(ii) * vPer_binCenters(i) / v_sigma**2

                    if ( bArg < 650 ) then

                        bF = bessI ( 0, real(bArg,dbl) )
                  
                        expTerm_D = exp ( real ( &
                                - ( vPerTrack(ii)**2 + ( vParTrack(ii) - vPar_binCenters(j) )**2  &
                                + vPer_binCenters(i)**2 ) &
                                / ( 2.0 * v_sigma**2 ) &
                                , dbl ) ) 

                        expTermX = exp ( - ( (rTrack(ii)-r_binCenters(m))**2  + (zTrack(ii)-z_binCenters(n))**2 ) &
                                / ( 2d0 * x_sigma**2 ) ) 

                        result_D = expTermX * expTerm_D * bF * 2d0 * real(pi,dbl)
                        f_rzvv_update(m_,n_,i_,j_) = result_D

                        !   This is a catch for when the modified bessel
                        !   function returns +Infinity

                        if ( bF > 1d300 ) then
                         
                            write(*,*) 'Bad'
 
                            f_rzvv_update(m_,n_,i_,j_)  = 0.0 
 
                        end if

                    else    ! do the gyro angle integral the old fashioned way

                        write(*,*) 'Ouch'

                        !   but use an adaptive range for the phi
                        !   integral depending on the value of vPerTrack.
                        !   Oh yes, this is brilliant!

                        phiRange    = 2.5 * v_sigma / vPerTrack(ii)     
                        phi = phiBase * phiRange
                        phi0    = 0.0

                        vX  = vPer_binCenters(i) * cos ( phi )
                        vY  = vPer_binCenters(i) * sin ( phi )
                        vZ  = vPar_binCenters(j)

                        vX0 = vPerTrack(ii) * cos ( phi0 )
                        vY0 = vPerTrack(ii) * sin ( phi0 )
                        vZ0 = vParTrack(ii)

                        expTermX = exp ( - ( rTrack(ii)**2  + zTrack(ii)**2 ) / ( 2d0 * x_sigma**2 ) ) 

                        f_rzvv_update(m_,n_,i_,j_)   = expTermX * sum ( &
                            exp ( -1.0 * &
                                  ( ( vX - vX0 )**2 + ( vY - vY0 )**2 + ( vZ - vZ0 )**2 ) & 
                                        / ( 2.0 * v_sigma**2 ) ) &
                            ) * abs ( phi(1) - phi(2) )

                    end if
 
                enddo
                enddo
                enddo
                enddo
               
                !   Normalise and add

                if ( normFac == 0 .or. tau == 0 ) then
                
                    write (*,*) 'DLG: ERROR, normFac == 0 or tau == 0 '
                    write (*,*) normFac, normFac_, tau, v_sigma**2
                    stop
                else

                    f_rzvv_update = f_rzvv_update / normFac * weightMod

                    f_rzvv(R_L(ii):R_R(ii),Z_L(ii):Z_R(ii),vPerL(ii):vPerR(ii),vParL(ii):vParR(ii))  = &
                        f_rzvv(R_L(ii):R_R(ii),Z_L(ii):Z_R(ii),vPerL(ii):vPerR(ii),vParL(ii):vParR(ii)) &
                        + dtArray(ii) / tau &
                        * f_rzvv_update(1:m_,1:n_,1:i_,1:j_) 

                endif 

                end do 

            else

                do ss=1,stepCnt
                    
                    if ( int(R_index(ss)) <= R_nBins .and. int(R_index(ss)) >= 1 .and. &
                        int(z_index(ss)) <= z_nBins .and. int(z_index(ss)) >= 1 .and. &
                        int(vPer_index(ss)) <= vPer_nBins .and. int(vPer_index(ss)) >= 1 .and. &
                        int(vPar_index(ss)) <= vPar_nBins .and. int(vPar_index(ss)) >= 1 ) then 

                        dV_ = vPer_binCenters(int(vPer_index(ss))) * rTrack(ss) * &
                            4.0 * pi**2 * R_binSize * z_binSize * &
                            vPer_binSize * vPar_binSize  
                   
                        f_rzvv(int(R_index(ss)),int(z_index(ss)), &
                               int(vPer_index(ss)),int(vPar_index(ss))) &
                           = f_rzvv(int(R_index(ss)),int(z_index(ss)), &
                               int(vPer_index(ss)),int(vPar_index(ss))) + dtArray(ss) / tau * weightMod / dV_
                    
                    else
#if DEBUG_LEVEL>1  
                        write(*,*) 'Off vPer Grid: ', vPer_index(ss), vPer_nBins
                        write(*,*) 'Off vPar Grid: ', vPar_index(ss), vPar_nBins
                        write(*,*) 'Off R Grid: ', R_index(ss), R_nBins
                        write(*,*) 'Off z Grid: ', z_index(ss), z_nBins
                        write(*,*) 'R: ', rTrack(ss), r_min, r_max, r_range 
#endif
                        nP_off_vGrid    = nP_off_vGrid + dtArray(ss) / tau

                    end if
                
                end do

            end if    
        
        end if

    end subroutine gc_orbit 

end module gc_integrate
