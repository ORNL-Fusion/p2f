pro plot_tdep

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 
	if strCmp ( getEnv ( 'HOSTNAME' ), 'benten.gat.com' ) then scratchDir = '/u/greendl/scratch' 

	eqdsk_fileName   = 'data/g129x129_1051206002.01120.cmod'
	;eqdsk_fileName	= 'data/eqdsk.122993'
	;eqdsk_fileName	= '~/data/eqdsk/g122080.03100'

	amu	= 1d0

	e_   = 1.60217646d-19
	q   = 1d0 * e_
	k   = 1.3806504d-23
	mi  = amu * 1.67262158d-27
	c   = 3.0d8

	print, '*** Using ', eqdsk_fileName
	eqdsk   = readGEQDSK ( eqdsk_fileName )

fileList	= file_search ( 'data/p2f_*' )
for ff=0,n_elements(fileList)-1 do begin

	fileName	= fileList[ff] 
	print, fileName

	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )

	; Read netcdf data into variables

	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv
	ncdf_varGet, cdfId, 'R_binCenters', R_binCenters
	ncdf_varGet, cdfId, 'z_binCenters', z_binCenters
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters
	ncdf_varGet, cdfId, 'R_binEdges', R_binEdges
	ncdf_varGet, cdfId, 'z_binEdges', z_binEdges
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges
	ncdf_varGet, cdfId, 'nP', nP
	ncdf_varGet, cdfId, 'R_binSize', R_binSize
	ncdf_varGet, cdfId, 'z_binSize', z_binSize
	ncdf_varGet, cdfId, 'vPerp_binSize', vPerp_binSize
	ncdf_varGet, cdfId, 'vPar_binSize', vPar_binSize
	ncdf_varGet, cdfId, 'density', density_file

	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'R_nBins' ), name, R_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'z_nBins' ), name, z_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPerp_nBins' ), name, vPerp_nBins
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPar_nBins' ), name, vPar_nBins

	ncdf_close, cdfId

	R_nBins	= n_elements ( R_binCenters )
	z_nBins	= n_elements ( z_binCenters )	
	vPerp_nBins	= n_elements ( vPerp_binCenters )
	vPar_nBins	= n_elements ( vPar_binCenters )	

;	Create a density map

	vPerp_3D	= rebin ( vPerp_binCenters, vPerp_nBins, R_nBins, z_nBins )
	vPerp_3D	= transpose ( vPerp_3D, [1,2,0] )
	vPar_4D	= rebin ( vPar_binCenters, vPar_nBins, vPerp_nBins, z_nBins, R_nBins )
	vPar_4D	= transpose ( vPar_4D, [3,2,1,0] )

	R_2D	= rebin ( R_binCenters, R_nBins, z_nBins )

	density	= total ( total ( f_rzvv, 4 ) * vPerp_3D , 3 ) * $
		vPerp_binSize * vPar_binSize * 2.0 * !pi

;	Create nP map

	nParticle_map	= density * z_binSize * R_binSize * R_2D * 2.0 * !pi

;	Create wPerp map ( perp energy PER PARTICLE )

	wPerp	= total ( 0.5 * mi * vPerp_3D^2 $
						* total ( f_rzvv, 4 ) * vPerp_3D , 3 ) $
		* vPerp_binSize * vPar_binSize * 2.0 * !pi

	wPar	= total ( total ( 0.5 * mi * vPar_4D^2 * f_rzvv, 4 ) * vPerp_3D , 3 ) $
		* vPerp_binSize * vPar_binSize * 2.0 * !pi


	wPerp	= wPerp / e_ * 1d-3 / density
	wPar	= wPar / e_ * 1d-3 / density

	iiBad	= where ( density eq 0, iiBadCnt )
	if iiBadCnt gt 0 then wPerp[iiBad] = 0
	if iiBadCnt gt 0 then wPar[iiBad] = 0


;	Create flux surface averaged quantities

;--------------------------------------------
;	Create a density profile

	density_all	= 0.0
	wperp_all	= 0.0
	wpar_all	= 0.0

	R_all	= 0.0
	z_all	= 0.0

	maskIn	= intArr ( size (density) )

	for i = 0, R_nBins-1 do begin
		for j = 0, z_nBins - 1 do begin

	    q1_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q1 )
	    q2_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q2 )
	    q3_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q3 )
	    q4_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q4 )

	    if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then begin

			;if (size(f_vv_all,/dim))[0] eq 0 then f_vv_all = reform(f_rzvv[i,j,*,*]) $
			;	else f_vv_all = [ [[ f_vv_all ]], [[ reform ( f_rzvv[i,j,*,*] ) ]] ]
			
			density_all	= [ density_all, density[i,j] ]
			wperp_all	= [ wperp_all, wperp[i,j] ]
			wpar_all	= [ wpar_all, wpar[i,j] ]

			R_all	= [ R_all, R_binCenters[i] ]
			z_all	= [ z_all, z_binCenters[j] ]

		endif

	    endfor 
	
	endfor

	density_all	= density_all[1:*]
	wperp_all	= wperp_all[1:*]
	wpar_all	= wpar_all[1:*]

	R_all	= R_all[1:*]
	z_all	= Z_all[1:*]
 
	psi_all = interpolate ( eqdsk.psizr, ( R_all - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_all - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psiRange	= abs ( eqdsk.siMag - eqdsk.siBry )
	psiNorm_all	= ( psi_all - eqdsk.siMag ) / psiRange			

	rho_all	= sqrt ( psiNorm_all )	

	R_binCenters2D	= rebin ( R_binCenters, R_nBins, z_nBins )
	z_binCenters2D	= transpose ( rebin ( z_binCenters, z_nBins, R_nBins ) )	
	psi_2D = interpolate ( eqdsk.psizr, ( R_binCenters2D - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_binCenters2D - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

	psi_2D	= ( psi_2D - eqdsk.siMag ) / psiRange

	;	Bin by rho coord.

	rho_nBins	= R_nBins/1.5
	rho_binEdges	= fIndGen ( rho_nBins+1 ) / rho_nBins 
	dRho	= abs(rho_binEdges[1]-rho_binEdges[2])
	rho_binCenters	= rho_binEdges[1:*] - dRho/2.0

	density_rho	= fltArr ( n_elements ( rho_binCenters ) )
	wperp_rho	= fltArr ( n_elements ( rho_binCenters ) )
	wpar_rho	= fltArr ( n_elements ( rho_binCenters ) )

	;f_vv_rho	= fltArr ( n_elements ( rho_binCenters ), vPerp_nBins, vPar_nBins )

	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rho_all ge rho_binCenters[i] - dRho $
							AND rho_all lt rho_binCenters[i] + dRho $
							AND psiNorm_all le 1, cnt)
		if cnt gt 0 then begin
			density_rho[i]	= total ( density_all[ iiDx ] ) / cnt
			wperp_rho[i]	= total ( wperp_all[ iiDx ] ) / cnt
			wpar_rho[i]	= total ( wpar_all[ iiDx ] ) / cnt

			;if cnt eq 1 then $
			;	f_vv_rho[i,*,*]	= f_vv_all[ *, *, iiDx ] $
			;else $
			;	f_vv_rho[i,*,*]	= total ( f_vv_all[ *, *, iiDx ], 3 ) / cnt
			;print, cnt
		endif
		
	endfor

	if (size(wPerpTime,/type) eq 0 ) then begin
		wPerpTime	= wPerp_rho
		wParTime	= wPar_rho
		densityTime	= density_rho
	endif else begin
		wPerpTime	= [ [ wPerpTime ], [ wPerp_rho ] ]
		wParTime	= [ [ wParTime ], [ wPar_rho ] ]
		densityTime	= [ [ densityTime ], [ density_rho ] ]
	endelse

endfor

time	= fIndGen ( n_elements ( wPerpTime[0,*] ) ) * 0.1
iSurface, wPerpTime, rho_binCenters, time
iSurface, wParTime, rho_binCenters, time
iSurface, densityTime, rho_binCenters, time

stop
end

