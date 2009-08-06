pro plot_p2f, $
	COMPARE = compare, $
	LINEAR = linear, $
	CQL3D = cql3d, $
	D3D = d3d, $
	AMU = amu, $
	FILENAME = fileName, $
	TEMPERATURE = temperature, $
	JACOBIAN = jacobian

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 
	if strCmp ( getEnv ( 'HOSTNAME' ), 'benten.gat.com' ) then scratchDir = '/u/greendl/scratch' 

	;fileName2	= scratchDir + '/p2f/d3d_D_Heidbrink_delta/data/fdis.dav.nc'
	fileName2	= scratchDir + '/p2f/cmod_H_minority/2.4MW/t00_delta/data/fdis.dav.nc'
	;eqdsk_fileName   = 'data/g129x129_1051206002.01120.cmod'
	;eqdsk_fileName	= 'data/eqdsk.122993'
	eqdsk_fileName	= 'data/eqdsk'

	if keyword_set(d3d) then eqdsk_fileName	= '~/data/eqdsk/g122080.03100'
	cql3dFileName	= 'data/cql3dData.sav'

	if keyword_set ( amu ) then amu = amu else amu	= 1d0

	if keyWord_set ( cql3d ) then restore, cql3dFileName

	e_   = 1.60217646d-19
	q   = 1d0 * e_
	k   = 1.3806504d-23
	mi  = amu * 1.67262158d-27
	c   = 3.0d8

	print, 'Using ', eqdsk_fileName
	eqdsk   = readGEQDSK ( eqdsk_fileName )

if keyword_set ( compare ) then begin

	cdfId	= ncdf_open ( fileName2, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv2
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges2
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges2
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters2
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters2
	nCdf_varGet, cdfId, 'nP', nP2

	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'R_nBins' ), name, R_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'z_nBins' ), name, z_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPerp_nBins' ), name, vPerp_nBins2
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'vPar_nBins' ), name, vPar_nBins2

	ncdf_close, cdfId

	dVPer2	= abs ( vPerp_binEdges2[0] - vPerp_binEdges2[1] )
	dVPar2	= abs ( vPar_binEdges2[0] - vPar_binEdges2[1] )

endif
	
	; Read netcdf data into variables

	if keyword_set(fileName ) then fileName = fileName else fileName = 'data/fdis.dav.nc'

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

if  keyword_set ( compare ) then begin

	R_nBins2	= n_elements ( R_binCenters2 )
	z_nBins2	= n_elements ( z_binCenters2 )	
	vPerp_nBins2	= n_elements ( vPerp_binCenters2 )
	vPar_nBins2	= n_elements ( vPar_binCenters2 )	

endif

;	Create a density map
	print, 'Create density map'

	vPerp_3D	= rebin ( vPerp_binCenters, vPerp_nBins, R_nBins, z_nBins )
	vPerp_3D	= transpose ( vPerp_3D, [1,2,0] )
	vPar_4D	= rebin ( vPar_binCenters, vPar_nBins, vPerp_nBins, z_nBins, R_nBins )
	vPar_4D	= transpose ( vPar_4D, [3,2,1,0] )

	R_2D	= rebin ( R_binCenters, R_nBins, z_nBins )

	density	= total ( total ( f_rzvv, 4 ) * vPerp_3D , 3 ) * $
		vPerp_binSize * vPar_binSize * 2.0 * !pi

;	Create nP map
	print, 'Create nP map'

	nParticle_map	= density * z_binSize * R_binSize * R_2D * 2.0 * !pi

;	Create wPerp map ( perp energy PER PARTICLE )
	print, 'Create wPerp map'

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

	print, 'Create flux surface average f'
	for i = 0, R_nBins-1 do begin
		for j = 0, z_nBins - 1 do begin

	    q1_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q1 )
	    q2_  = where ( ( R_binCenters[i] - eqdsk.rbbbs gt 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q2 )
	    q3_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs gt 0 ), q3 )
	    q4_  = where ( ( R_binCenters[i] - eqdsk.rbbbs le 0 ) and ( z_binCenters[j] - eqdsk.zbbbs le 0 ), q4 )

	    if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then begin

			if (size(f_vv_all,/dim))[0] eq 0 then f_vv_all = reform(f_rzvv[i,j,*,*]) $
				else f_vv_all = [ [[ f_vv_all ]], [[ reform ( f_rzvv[i,j,*,*] ) ]] ]
			
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

	f_vv_rho	= fltArr ( n_elements ( rho_binCenters ), vPerp_nBins, vPar_nBins )

	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rho_all ge rho_binCenters[i] - dRho $
							AND rho_all lt rho_binCenters[i] + dRho $
							AND psiNorm_all le 1, cnt)
		if cnt gt 0 then begin
			density_rho[i]	= total ( density_all[ iiDx ] ) / cnt
			wperp_rho[i]	= total ( wperp_all[ iiDx ] ) / cnt
			wpar_rho[i]	= total ( wpar_all[ iiDx ] ) / cnt

			if cnt eq 1 then $
				f_vv_rho[i,*,*]	= f_vv_all[ *, *, iiDx ] $
			else $
				f_vv_rho[i,*,*]	= total ( f_vv_all[ *, *, iiDx ], 3 ) / cnt
			;print, cnt
		endif
		
	endfor

	;	Write a netCDF file for f_vv_rho

	print, 'Write netCDF file f_vv_rho.nc'
	nc_id	= nCdf_create ( 'data/f_rho_vv.nc', /clobber )
	nCdf_control, nc_id, /fill
	
	vPer_nBins_id	= nCdf_dimDef ( nc_id, 'vPer_nBins', vPerp_nBins )
	vPar_nBins_id	= nCdf_dimDef ( nc_id, 'vPar_nBins', vPar_nBins )
	rho_nBins_id	= nCdf_dimDef ( nc_id, 'rho_nBins', rho_nBins )
	
	f_rho_vv_id	= nCdf_varDef ( nc_id, 'f_vv_rho', [ rho_nBins_id, vPer_nBins_id, vPar_nBins_id ], /float )
	
	nCdf_control, nc_id, /enDef
	nCdf_varPut, nc_id, f_rho_vv_id, f_vv_rho
	nCdf_close, nc_id

	;	Create the analytical profile AORSA uses	

	n0	= 2.2e18
	nLim	= 2.0e17
	alpha	= 1.5
	beta_	= 1.0

	;!n0	= 1.522e18
	;!nLim	= 0.933e17
	;!alpha	= 5.0
	;!beta_	= 3.0

	;n0	= 3.6e18
	;nLim	= 10.8e17
	;alpha	= 0.6
	;beta_	= 1.4

	print, 'Fit to get aorsa profile parameters'	
	AA	= [nLim,n0,alpha,beta_]
	iiKeep	= where ( psi_2D lt 0.9 )
	densityFit	= curveFit ( sqrt(psi_2D[iiKeep]), density[iiKeep], density[iiKeep]*0+1, AA, sigmaFit, $
		   function_name = 'profile_fit', /noDeriv, status = status ) 

	density_rho_aorsaFIT	= aa[0]+ ( aa[1] - aa[0] )*(1d0-rho_binCenters^aa[3])^aa[2]
	density_rho_aorsa	= nLim + ( n0 - nLim )*(1d0-rho_binCenters^beta_)^alpha

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/profiles.eps'
	print, 'Write data/profiles.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=10, ysize = 16,xoffset=.1, yoffset=.1, /encapsul
	!p.charSize = 1.0
	!p.multi	= [0,1,2]
	plot, sqrt ( psi_2D ), density/1e19, $
		psym = 4, $
		xtitle = 'rho', $
		ytitle = 'density [m!U-3!N] x10!U19!N', $
		xRange = [0, 1], $
		color = 0, $
		yRange = [0,0.3], $
		yStyle = 1
	loadct, 12, /silent
	oPlot, rho_binCenters, density_rho/1e19, $
		color = 8*16-1, $
		thick = 2.0
	oPlot, rho_binCenters, density_rho_aorsa/1e19,$
		thick = 2.0, $
		color = 1*16-1

	if status eq 0 then begin
	print, 'Density fit successfull'
	oPlot, rho_binCenters, density_rho_aorsaFIT/1e19,$
		thick = 2.0, $
		color = 12*16-1
	xyOuts, 0.1, 0.98, 'nLim x19: '+string(aa[0]/1e19,for='(f4.2)') $
			+'  n0 x19: '+string(aa[1]/1e19,for='(f5.2)') $
			+'  a: '+string(aa[2],for='(f4.2)') $
			+'  b: '+string(aa[3],for='(f4.2)'), $
		   /norm, $
		   color = 12*16-1, $
		   charSize = 1.0
    endif else print, 'Density fit UNSUCCESSFUL'

;		plots, [0.7,0.7], [1.8,1.8], $
;			   color = 12*16-1, psym = 4 
;	   plots, [0.6,0.7], [1.6,1.6], $
;			   color = 10*16-1, thick=2
;	   xyOuts, 0.75, 1.8, 'orbit-rf', color = 0, charSize = 1
;	   xyOuts, 0.75, 1.6, 'cql3d', color = 1*16-1, charSize = 1
	
	plot, sqrt ( psi_2D ), wperp, $
		psym = 4, $
		xtitle = 'rho', $
		ytitle = 'wPerp/wPar [keV] per particle', $
		xRange = [0, 1], $
		color = 0, $
		yRange = [0,20], $
		yStyle = 1, /noData
	;oPlot, sqrt ( psi_2D ), wpar, $
	;		psym = 4, $
	;		color = 12*16-1
	oPlot, rho_binCenters, wperp_rho, $
		color = 8*16-1, $
		thick = 2.0
	oPlot, rho_binCenters, wpar_rho, $
		color = 12*16-1, $
		thick = 2.0

	if keyword_set ( cql3d ) then begin
		oPlot, rho_cql3d[1:*], wPerp_cql3d[1:*], $
			   color = 1*16-1, $
			   thick = 2
		oPlot, rho_cql3d[1:*], wPar_cql3d[1:*], $
			   color = 10*16-1, $
			   thick = 2
	   plots, [0.5,0.5], [9,9], $
			   color = 0, psym = 4 
	   plots, [0.4,0.5], [8,8], $
			   color = 1*16-1, thick=2
	   plots, [0.5,0.5], [7,7], $
			   color = 12*16-1, psym = 4 
	   plots, [0.4,0.5], [6,6], $
			   color = 10*16-1, thick=2
	   xyOuts, 0.55, 9, 'orbit-rf wPerp', color = 0, charSize = 1
	   xyOuts, 0.55, 8, 'cql3d wPerp', color = 1*16-1, charSize = 1
	   xyOuts, 0.55, 7, 'orbit-rf wPar', color = 12*16-1, charSize = 1
	   xyOuts, 0.55, 6, 'cql3d wPar', color = 10*16-1, charSize = 1


	
	endif
	device, /close_file

;
;--------------------------------------------


;--------------------------------------------
;	Create a temperature profile by fitting 2D gaussian
;	functions at each pt in space

	if keyword_set ( temperature ) then begin

	tempProfile	= fltArr ( R_nBins, z_nBins )
	psiTemp2D	= fltArr ( R_nBins, z_nBins )
	rhoTemp2D	= fltArr ( R_nBins, z_nBins )

	print, 'Create temperature profile'
	for i = 0, R_nBins-1 do begin
		;print, i
		for j = 0, z_nBins-1 do begin
		
			psiTemp2D[i,j] = interpolate ( eqdsk.psizr, $
				( R_binCenters[i] - eqdsk.rleft ) / eqdsk.rdim * eqdsk.nW, $
    			( z_binCenters[j] - min ( eqdsk.z ) ) / eqdsk.zdim * eqdsk.nH )

			rhoTemp2D[i,j]	= sqrt ( ( psiTemp2D[i,j] - eqdsk.siMag ) / psiRange )	

			tmp	= [reverse(reform(f_rzvv[i,j,*,*]),1),reform(f_rzvv[i,j,*,*])]
			if mean ( tmp ) gt 0 then begin
				fit	= gauss2dFit ( tmp, A )
				tempProfile[i,j]	= mi * ( $
					( A[2] * vPerp_binsize )^2 + ( A[3] * vPar_binSize )^2 $
						) / ( 2.0 * 1e3 * e_ )  
			endif

		endfor
	endfor
	
	temp_rho	= fltArr ( n_elements ( rho_binCenters ) )

	for i = 0, n_elements ( rho_binCenters ) - 1 do begin

		iiDx	= where ( rhoTemp2D[*] ge rho_binCenters[i] - dRho $
			AND rhoTemp2D[*] lt rho_binCenters[i] + dRho, cnt)
		if cnt gt 0 then begin
			temp_rho[i]	= total ( (tempProfile[*])[ iiDx ] )
			temp_rho[i]	= temp_rho[i] / cnt
			;print, cnt
		endif
		
	endfor

	temp_rho	= temp_rho * 1e3


	;	Create the analytical profile AORSA uses	

	t0	= 3.0e3 
	tLim	= 0.052e3
	alpha	= 1.3 
	beta_	= 1.9

	t0	= 24.0e3 
	tLim	= 2.0e3
	alpha	= 1.0 
	beta_	= 2.8

	;	Fit to get aorsa profile parameters

	print, 'Fit to get aorsa profile parameters'	
	maxTmp = 30.0*1e3;[keV]
	iiKeep	= where ( tempProfile*1e3 lt maxTmp and rhoTemp2d lt 1.0 )
	AA	= [tLim,t0,alpha,beta_]
	tempFit	= curveFit ( rhoTemp2D[iiKeep], tempProfile[iiKeep]*1e3, tempProfile[iiKeep]*0+1, AA, sigmaFit, $
		   function_name = 'profile_fit', /noDeriv, status = status ) 

	temp_rho_aorsa	= tLim + ( t0 - tLim )*(1d0-rho_binCenters^beta_)^alpha
	temp_rho_aorsaFIT	= aa[0]+ ( aa[1] - aa[0] )*(1d0-rho_binCenters^aa[3])^aa[2]

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/temp.eps'
	print, 'Wrtie data/temp.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=10, ysize = 10,xoffset=.1, yoffset=.1, /encapsul
	!p.charSize = 1.0
	!p.multi	= 0 
	
	plot, rhoTemp2D, tempProfile, $
		yTitle = 'temperature [keV]', $
		xTitle = 'rho', $
		thick = 1, $
		color = 0, $
		psym = 4, $
		xRange = [0,1.0], $
		xStyle = 1, $
		yStyle = 1, $
		yRange = [0, maxTmp/1e3]

	oPlot, rho_binCenters, temp_rho_aorsa/1e3, $
		psym = -4, $
	 	color = 8*16-1, $
		thick = 2

	if status eq 0 then begin
	print, 'Temperature fit successful'
	oPlot,rho_binCenters, temp_rho_aorsaFIT/1e3, $
		color = 12*16-1, $
		thick = 2
	xyOuts, 0.1, 0.95, 'tLim: '+string(aa[0]/1e3,for='(f4.2)') $
			+'  t0: '+string(aa[1]/1e3,for='(f5.2)') $
			+'  a: '+string(aa[2],for='(f4.2)') $
			+'  b: '+string(aa[3],for='(f4.2)'), $
		   /norm, $
		   color = 12*16-1, $
		   charSize = 1.0
   	endif else print, 'Temperature fit UNSUCCESSFUL'

   	device, /close_file

	endif 
;
;--------------------------------------------


	dVPer	= abs ( vPerp_binEdges[0] - vPerp_binEdges[1] )
	dVPar	= abs ( vPar_binEdges[0] - vPar_binEdges[1] )
	dR	= abs ( R_binEdges[0] - R_binEdges[1] )
	dz	= abs ( z_binEdges[0] - z_binEdges[1] )

	
	; Plot

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/f_rzvv.eps'
	print, 'Write data/f_rzvv.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=25, ysize=37.5,xoffset=.1, yoffset=.1, /encapsul

	nLevs	= 10 
	levScale	= 1d-7
	levels	= 10.0^fIndGen(nLevs)*levScale
	colors	= reverse ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 )

	if keyword_set ( linear ) then begin

		nLevs	= 10
		levels = fIndGen ( nLevs ) / nLevs *1.0d-1
		colors	= reverse ( bytScl ( levels, top = 253 ) + 1 )
		levels[nLevs-1]	= 1000

	endif
	
	loadct, 0, /silent

	!p.multi = [0,R_nBins,z_nBins]
	
	for j=0,z_nBins-1 do begin
		for i=0,R_nBins-1 do begin
	
			f_vv_smooth	= reform(f_rzvv[i,j,*,*])

			contourVar	= transpose ( f_vv_smooth )
			if keyword_set ( jacobian ) then $
				contourVar  = transpose ( f_vv_smooth * rebin(vPerp_binCenters,vPerp_nBins, vPar_nBins)$
					/median(vPerp_binCenters) )

			contour, contourVar, $
				vPar_binCenters / 3.0e6, vPerp_binCenters / 3.0e6, $
				levels = levels, $
				c_colors = colors, $
				color = 0, $
				charSize = 0.01, $
				xRange = [min(vPar_binCenters),max(vPar_binCenters)] / 3.0e6,$
				yRange = [0.0,max(vPerp_binCenters)] / 3.0e6, $
				title = 'R: '+string ( r_binCenters[i], for='(f5.2)' ) $
					+ ' z: '+string ( z_binCenters[j], for='(f5.2)' ), $
				xTicks = 1, $
				yTicks = 1, $
			 	xStyle = 9, $
				yStyle = 9, $
				thick = 0.5, /fill

			;print, i, j
		endfor	
	endfor	

	!p.position = 0
	!p.charSize = 2	
	loadct, 0, /silent
	xyOuts, 0.98, 0.975, 'fileName: '+fileName, color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.95, 'nP: '+string(nP,format='(e7.1)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.925, 'R_nBins: '+string(R_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.9, 'z_nBins: '+string(z_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.875, 'vPerp_nBins: '+string(vPerp_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0
	xyOuts, 0.98, 0.85, 'vPar_nBins: '+string(vPar_nBins,format='(i3.3)'), color = 0, /norm, align = 1.0

	device, /close_file

if keyword_set ( compare ) then begin

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/compare.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=16, ysize=8,xoffset=.1, yoffset=.1, /encapsul

	!p.multi = [0,3,2]
	!p.charsize = 1.0
	loadct, 12, /silent

	pt	= where ( density eq max ( density ) )
	i	= (array_indices(density,pt))[0]
	j	= (array_indices(density,pt))[1]

	densityHere	= density[i,j]
	tempHere	= tempProfile[i,j]	

	f_vv	= reform(f_rzvv[i,j,*,*])
	f_vv2	= reform(f_rzvv2[i,j,*,*])

	vPer2D	= rebin(vPerp_binCenters,vPerp_nBins,vPar_nBins)
	vPer2D2	= rebin(vPerp_binCenters2,vPerp_nBins2,vPar_nBins2)

	;	Create analytical maxwellian


	T_keV    = tempHere ; [keV]  
	T	= T_keV * 1d3 * e_ / k

	f_vv_a	= fltArr ( size ( f_vv, /dim ) )
	
	for i = 0, n_elements ( vPerp_binCenters ) - 1 do begin
		for j = 0, n_elements ( vPar_binCenters ) - 1 do begin
			v	= sqrt ( (vPerp_binCenters[i]^2) $
					+ (vPar_binCenters[j]^2) )
			meanV	= 0.0
			f_vv_a[i,j]	= $
				 ( mi / ( 2.0 * !pi * k * T ))^1.5 * $
					exp ( -  mi * (v-meanV)^2 / (2.0*k*T) )
		endfor
	endfor
	
	f_vv_a	= densityHere * f_vv_a $
		/ total ( f_vv_a * rebin(vPerp_binCenters,vPerp_nBins,vPar_nBins) * dVpar * dVPer * 2.0 * !Pi )

	dfdupar	= dlg_pderiv ( f_vv, 2, dVpar )
	dfdupar2	= dlg_pderiv ( f_vv2, 2, dVpar2 )
	dfdupar_a	= dlg_pderiv ( f_vv_a, 2, dVpar )

	dfduper	= dlg_pderiv ( f_vv, 1, dVper )
	dfduper2	= dlg_pderiv ( f_vv2, 1, dVper2 )
	dfduper_a	= dlg_pderiv ( f_vv_a, 1, dVper )

	;	Interpoate to match grids
	
	nGridPer	= (vPerp_binCenters2 - min ( vPerp_binCenters )) $
		/ ( max ( vPerp_binCenters ) - min ( vPerp_binCenters ) ) * (vPerp_nBins-1)

	nGridPar	= (vPar_binCenters2 - min ( vPar_binCenters )) $
		/ ( max ( vPar_binCenters ) - min ( vPar_binCenters ) ) * (vPar_nBins-1)

	f_vv_a	= interpolate ( f_vv_a, nGridPer, nGridPar, /grid )
	f_vv	= interpolate ( f_vv, nGridPer, nGridPar, /grid )
	dfduper_a	= interpolate ( dfduper_a, nGridPer, nGridPar, /grid )
	dfdupar_a	= interpolate ( dfdupar_a, nGridPer, nGridPar, /grid )
	dfduper	= interpolate ( dfduper, nGridPer, nGridPar, /grid )
	dfdupar	= interpolate ( dfdupar, nGridPer, nGridPar, /grid )

	nLevs	= 10
	levels = fIndGen ( nLevs ) / nLevs * max(f_vv)	
	colors	= reverse ( bytScl ( levels, top = 253 ) + 1 )
	loadct, 1, /silent
	contour, transpose ( f_vv2 ), $
				vPar_binCenters2*1d-6, vPerp_binCenters2*1d-6, $
				levels = levels, $
				c_colors = colors, $
				color = 0, $
				charSize = 1.0, /fill, $
				xTitle = 'vPar [x10!U6!N m/s]', $
				yTitle = 'vPer [x10!U6!N m/s]', $
				xRange = [-1.1,1.1], xStyle = 1

	contour, transpose ( f_vv2 ), $
				vPar_binCenters2*1d-6, vPerp_binCenters2*1d-6, $
				levels = levels, $
				c_colors = (colors*1.2)<254, $
				color = 0, $
				charSize = 1.0, /over
	plots, [-1.1,1.1],$
			[vPerp_binCenters2[1],vPerp_binCenters2[1]]*1d-6, $
			color = 0, $
			lineStyle = 2
	plots, [vPar_binCenters2[vPar_nBins2/2],vPar_binCenters2[vPar_nBins2/2]]*1d-6,$
			[vPerp_binEdges2[0],max(vPerp_binEdges2)]*1d-6, $
			color = 0, $
			lineStyle = 2
	
	loadct, 12 , /silent
	plot, vPar_binCenters2*1d-6, f_vv2[1,*],$
			color = 0 , $
			xTitle = 'vPar [x10!U6!N m/s]', $
			yTitle = 'f [s!U3!N/m!U6!N]', /noData
	oplot, vPar_binCenters2*1d-6, f_vv_a[1,*], color = 0, thick = 2.0
	loadct, 3, /silent
	oplot, vPar_binCenters2*1d-6, f_vv[1,*], color = 150, thick = 1.0
	loadct, 1, /silent
	oplot, vPar_binCenters2*1d-6, f_vv2[1,*], color = 150, thick = 1.0

	plot, vPar_binCenters2*1d-6, dfdupar2[1,*],$
		   	color = 0,$
			xTitle = 'vPar [x10!U6!N m/s]', $
			yTitle = 'dfdvPar', /noData
	oplot, vPar_binCenters2*1d-6, dfdupar_a[1,*], color = 0, thick = 2.0
	loadct, 3, /silent
	oplot, vPar_binCenters2*1d-6, dfdupar[1,*], color = 150, thick = 1.0
	loadct, 1, /silent
	oplot, vPar_binCenters2*1d-6, dfdupar2[1,*], color = 150, thick = 1.0

	loadct, 3, /silent
	contour, transpose ( f_vv ), $
				vPar_binCenters2*1d-6, vPerp_binCenters2*1d-6, $
				levels = levels, $
				c_colors = colors, $
				color = 0, $
				charSize = 1.0, /fill, $
				xTitle = 'vPar [x10!U6!N m/s]', $
				yTitle = 'vPer [x10!U6!N m/s]', $
				xRange = [-1.1,1.1], xStyle = 1


	contour, transpose ( f_vv ), $
				vPar_binCenters2*1d-6, vPerp_binCenters2*1d-6, $
				levels = levels, $
				c_colors = (colors*1.2)<254, $
				color = 0, $
				charSize = 1.0, /over
	plots, [-1.1,1.1],$
			[vPerp_binCenters2[1],vPerp_binCenters2[1]]*1d-6, $
			color = 0, $
			lineStyle = 2
	plots, [vPar_binCenters2[vPar_nBins2/2],vPar_binCenters2[vPar_nBins2/2]]*1d-6,$
			[vPerp_binEdges2[0],max(vPerp_binEdges2)]*1d-6, $
			color = 0, $
			lineStyle = 2
	
	loadct, 12, /silent
	plot, vPerp_binCenters2*1d-6, f_vv2[*,vPar_nBins2/2],$
		   	color = 0,$
			xTitle = 'vPer [x10!U6!N m/s]', $
			yTitle = 'f [s!U3!N/m!U6!N]', /noData
	oplot, vPerp_binCenters2*1d-6, f_vv_a[*,vPar_nBins2/2], color = 0, thick = 2.0
	loadct, 3, /silent
	oplot, vPerp_binCenters2*1d-6, f_vv[*,vPar_nBins2/2], color = 150, thick = 1.0
	loadct, 1, /silent
	oplot, vPerp_binCenters2*1d-6, f_vv2[*,vPar_nBins2/2], color = 150, thick = 1.0

	plot, vPerp_binCenters2*1d-6, dfduper2[*,vPar_nBins2/2], $
			color = 0,$
			xTitle = 'vPer [x10!U6!N m/s]', $
			yTitle = 'dfdvPer', /noData
	oplot, vPerp_binCenters2*1d-6, dfduper_a[*,vPar_nBins2/2], color = 0, thick = 2.0
	loadct, 3, /silent
	oplot, vPerp_binCenters2*1d-6, dfduper[*,vPar_nBins2/2], color = 150, thick = 1.0
	loadct, 1, /silent
	oplot, vPerp_binCenters2*1d-6, dfduper2[*,vPar_nBins2/2], color = 150, thick = 1.0

endif
device, /close_file
set_plot, 'X'

;;	Contour f_rho_vv
;
;for i=0,n_elements(rho_binCenters)-1 do begin
;
;	contour, transpose ( f_vv_rho[i,*,*] ), $
;		levels = 10.0^fIndgen(12)*1d-2, $
;		path_xy = path, path_info = info
;	for j=0,n_elements(info.offset)-1 do begin
;		iSurface, path[0,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1]*0+j, $
;					path[0,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
;					path[1,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
;					/overPlot, $
;					transparency = j*10
;	endfor
;
;endfor
	nRPlot	= 7
	!p.multi = [0,nRPlot,1]
	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'data/f_Rvv_'+file_baseName(fileName)+'.eps'
	print, 'Write data/f_Rvv_ file'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=42, ysize = 7,xoffset=.1, yoffset=.1, /encapsul

	nLevs	= 10 
	levScale	= 1d-7
	levels	= 10.0^fIndGen(nLevs)*levScale
	colors	= reverse ( bytScl ( fIndGen(nLevs), top = 253 ) + 1 ) 
	nRpts	= n_elements(f_rzvv[*,0,0,0])

	cnt = 0
	for j=nRpts/2-1,nRPts-1,(nRpts/2)/nRPlot>1 do begin
			if cnt le nRPlot-1 then begin	
			!p.position	= [cnt*1d0/nRPlot+0.03,0.2,(cnt+1)*1d0/nRPlot,0.95]
			cnt	= cnt+1
			jII	= 	n_elements(f_rzvv[0,*,0,0])/2
			f_vv_smooth	= reform(f_rzvv[j,jII,*,*])
	
			if j eq 0 then ycharSize = 1.0 else ycharSize = 0.01

			contour, transpose ( f_vv_smooth )	, $
				vPar_binCenters / 1.0e6, vPerp_binCenters / 1.0e6, $
				levels = levels, $
				c_colors = colors, $
				color = 0, $
				charSize = 1, $
				ycharSize = 2, $
				xCharSize = 2, $
				;xRange = [min(vPar_binCenters),max(vPar_binCenters)] / 1.0e6,$
				xRange = [-3,3],$
				yRange = [0.0,max(vPerp_binCenters)] / 1.0e6, $
			 	xStyle = 9, $
				yStyle = 9, $
				xTicks = 4, $
				yTitle = 'vPer [m/s] x10!U6!N', $
				xTitle = 'vPar [m/s] x10!U6!N', $
				c_label	= levels * 0 + 1, $
				c_charSize = 0.5, $
				thick = 1.0, $
				xTickFormat = '(i2)'

				;xyouts, 2.0, 5.0, 'R:'+string ( r_binCenters[j], for='(f5.2)' ),$
				;	   align = 1, $
				;	   charSize = 1.5, $
				;	   charThick = 3.0
			    ;xyouts, 2.0, 4.0, 'rho:'+string ( sqrt(psi_2d[j,jII]), for='(f5.2)' ),$
				;	   align = 1, $
				;	   charSize = 1.5, $
				;	   charThick = 3.0
	
			endif
	endfor	

	device, /close_file

old_dev = !D.name
set_plot, 'ps'
outfname	= 'data/wPerp2D.eps'
print, 'Write data/wPerp2D.eps'
device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
	xsize=10, ysize = 10,xoffset=.1, yoffset=.1, /encapsul
	
loadct, 3, /silen
;!p.multi=[0,2,2]
!p.position = 0
!p.multi = 0
!p.charSize = 1.0
levels	= (fIndGen(100)+0.1)/100*80
colors	= reverse (  bytScl( levels, top = 253 ) + 1 )
contour, wPerp, R_binCenters, z_binCenters, $
	levels = levels, $
	c_colors = colors, $
	color = 0, $
	yRange = [min(eqdsk.zbbbs),max(eqdsk.zbbbs)], $
	yStyle = 1, $
	/fill, $
	title = 'wPerp', $
	xTitle = 'R [m]', $
	yTitle = 'z [m]'
plots, eqdsk.rbbbs, eqdsk.zbbbs, $
	color = 0, $
	thick = 2

;	Resonant surfaces

	freq	= 80.0e6 * 2.0 * !pi
	harm	= [1.0]

	B	= freq * mi / ( q * harm )

	for i=0,n_elements(harm)-1 do begin

		for z=0,n_elements(eqdsk.bMag[0,*])-1 do begin
			
				rIndex	= where ( abs(eqdsk.bMag[*,z]-B[i]) eq min(abs(eqdsk.bMag[*,z]-B[i])), iiCnt )
				if (eqdsk.mask)[rIndex,z] gt 0 then $
					plots, eqdsk.r[rIndex],eqdsk.z[z],$
						psym = 3,$
					   	color = 0
		endfor	
	endfor


device, /close_file

outfname	= 'data/density2D.eps'
print, 'Write data/density2D.eps'
device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
	xsize=10, ysize = 10,xoffset=.1, yoffset=.1, /encapsul
	
levels	= (fIndGen(100)+1)/100*6e18
colors	= reverse (  bytScl( levels, top = 253 ) + 1 )
contour, density, R_binCenters, z_binCenters, $
	levels = levels, $
	c_colors = colors, $
	color = 0, $
	yRange = [min(eqdsk.zbbbs),max(eqdsk.zbbbs)], $
	yStyle = 1, $
	/fill, $
	title = 'density', $
	xtitle = 'R [m]', $
	yTitle = 'z [m]'
plots, eqdsk.rbbbs, eqdsk.zbbbs, $
	color = 0, $
	thick = 2

	for i=0,n_elements(harm)-1 do begin

		for z=0,n_elements(eqdsk.bMag[0,*])-1 do begin
			
				rIndex	= where ( abs(eqdsk.bMag[*,z]-B[i]) eq min(abs(eqdsk.bMag[*,z]-B[i])), iiCnt )
				if (eqdsk.mask)[rIndex,z] gt 0 then $
					plots, eqdsk.r[rIndex],eqdsk.z[z],$
						psym = 4,$
					   	color = 0, $
						symSize = 1
		endfor	
	endfor


device, /close_file
for i=1,128 do begin
		free_lun, i
endfor
end

