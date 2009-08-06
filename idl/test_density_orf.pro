pro test_density_orf

	;   Constants
	
	k   = 1.3806504e-23
	e_   = 1.60217646e-19
	amu	= 1.0
	mi  = amu * 1.67262158e-27
	q   = 1.0 * e_
	c   = 3.0e8

	struct	= { psi : 0.0, $
				theta : 0.0, $
				R : 0.0, $
				z : 0.0, $
				E : 0.0, $
				lambda : 0.0, $
				weight : 0.0};, $
				;rho : 0.0 }

	;fileName	= '~/data/particleLists/fdis_axis'
	;fileName	= '~/data/particleLists/fdis_cmod_t0_2kev'
	;fileName	= '~/data/particleLists/fdis_cmod_highResol'
	;fileName	= '~/data/particleLists/fdis_2keV_CMOD_dlg.dav'
	;fileName	= '~/data/particleLists/fdis_122993_highRes'
	;fileName	= '~/data/particleLists/d3d_heidbrink/fdis_122993_t400_heidbrink'
	;fileName	= '~/data/particleLists/fdis_cmod_highsmap'
	;fileName	= '~/data/particleLists/cmod_600kW/fdis_cmod_t1msec_2kev'
	;fileName	= '~/data/particleLists/cmod_600kW/fdis_cmod_1msec_new'
	;fileName	= '~/data/particleLists/cmod_600kW/fdis_p'
	;fileName	= '~/data/particleLists/cmod_600kW/fdis_nop'
	;fileName	= '~/data/particleLists/cmod_600kW/fdis_cmod_2msec'
	;fileName	= '~/data/particleLists/d3d_heidbrink/fdis_122993_t14msec'
	fileName	= '~/data/particleLists/cmod_600kW/fdis_t1msec_rf'
	;fileName	= '~/data/particleLists/d3d_heidbrink/fdis_122993_t14ms_may13'

	eqDskFileName	= '~/code/sMap/run/cmod/g129x129_1051206002.01120.cmod'
	;eqDskFileName	= '~/data/eqdsk/eqdsk.122993'

	mapFileName	= "maptoeqs"

	eqdsk	= readGeqdsk ( eqDskFileName )

	nP	= file_lines ( fileName ) - 1
	;nP	= 100000
	data	= replicate ( struct, nP )

	openR, unit, fileName, /get_lun
		skip_lun, unit, 1, /lines
		readf, unit, data
	free_lun, unit

	data.R	= data.R / 100.0
	data.z	= data.z / 100.0

	weight	= data.weight
	
	;	Use mapToEqs file to convert rho,theta into
	;	R,z coords instead of Orbit-RF native mapping

	plot_mapToEqs, (data.psi), data.theta, rOut, zOut, eqdsk, mapFileName, rz2pt=0
	;plot_mapToEqs, psi_dlg, theta_dlg, rOut, zOut, eqdsk, mapFileName, $
	;	/rz2pt, rIn	= data.r, zIn = data.z

	;R_nBins = 22.0
	R_nBins = 64.0
	R_min   = min ( eqdsk.rbbbs )
	R_max   = eqdsk.rdim + eqdsk.rleft
	R_range = ( R_max - R_min )
	R_binSize   = R_range / R_nBins
	R_binEdges  = fIndGen ( R_nBins + 1 ) * R_binSize + R_min
	R_binCenters    = R_binEdges[0:R_nBins-1] + R_binSize / 2.0
	
	;z_nBins = 11.0
	z_nBins = 64.0
	z_min   = min ( eqdsk.zbbbs ) 
	z_max   = max ( eqdsk.zbbbs ) 
	z_range = ( z_max - z_min )
	z_binSize   = z_range / z_nBins
	z_binEdges  = fIndGen ( z_nBins + 1 ) * z_binSize + z_min
	z_binCenters    = z_binEdges[0:z_nBins-1] + z_binSize / 2.0
	
	r_index =  ( data.R - R_min ) / R_range * R_nBins
	z_index = ( data.z - z_min ) / z_range * z_nBins

	r_indexA =  ( rOut - R_min ) / R_range * R_nBins
	z_indexA = ( zOut - z_min ) / z_range * z_nBins

	f_Rz  = fltArr ( R_nBins, z_nBins )
	f_RzA  = fltArr ( R_nBins, z_nBins )

	f_Rz_  = fltArr ( R_nBins, z_nBins )
	f_RzA_  = fltArr ( R_nBins, z_nBins )

	nP_mhc	= 0L
	nP_dlg	= 0L

	for i = 0L, nP - 1 do begin
	
	    if R_index[i] le R_nBins-1 AND R_index[i] ge 0 AND $
	        z_index[i] le z_nBins-1 AND z_index[i] ge 0 then begin

	        f_Rz [ R_index[i], z_index[i] ] += weight[i]
			nP_mhc	+= 1L
			
		endif

	    if R_indexA[i] le R_nBins-1 AND R_indexA[i] ge 0 AND $
	        z_indexA[i] le z_nBins-1 AND z_indexA[i] ge 0 then begin 
	        
			f_RzA [ R_indexA[i], z_indexA[i] ] += weight[i]
			nP_dlg	+= 1L

		endif


	    if R_index[i] le R_nBins-1 AND R_index[i] ge 0 AND $
	        z_index[i] le z_nBins-1 AND z_index[i] ge 0 then $
		    f_Rz_ [ R_index[i], z_index[i] ] += 1.0 

	    if R_indexA[i] le R_nBins-1 AND R_indexA[i] ge 0 AND $
	        z_indexA[i] le z_nBins-1 AND z_indexA[i] ge 0 then $
	        f_RzA_ [ R_indexA[i], z_indexA[i] ] += 1.0
	
	endfor

	;	Apply jacobian

	f_Rz	= f_Rz / rebin ( R_binCenters, R_nBins, z_nBins ) * 2.0*!pi*1e3
	f_RzA	= f_RzA / rebin ( R_binCenters, R_nBins, z_nBins ) * 2.0*!pi*1e3

	f_Rz_	= f_Rz_ / rebin ( R_binCenters, R_nBins, z_nBins )
	f_RzA_	= f_RzA_ / rebin ( R_binCenters, R_nBins, z_nBins )

	!p.multi=[0,1,2]
	!p.charSize = 1.8
	window, 0, xSize = 1000, ySize = 600
	zii	= z_nBins / 2 
	plot, R_binCenters, f_Rz[*,zii], title = 'With weighting', $
		xTitle = 'R[m]', yTitle='Dens', $
		psym = -4, $
		yRange = [0,max(f_Rz)]
	oplot, R_binCenters, f_Rz[*,zii-1]
	oplot, R_binCenters, f_Rz[*,zii+1]
	
	oplot, R_binCenters, f_RzA[*,zii], thick = 2.0
	oplot, R_binCenters, f_RzA[*,zii-1], thick = 2.0
	oplot, R_binCenters, f_RzA[*,zii+1], thick = 2.0

	oPlot, [eqdsk.rmaxis, eqdsk.rmaxis ], [0.0, max(f_Rz)]

	plot, R_binCenters, f_Rz_[*,zii], title = 'Without weighting', $
		xTitle = 'R[m]', yTitle='Dens [Arb. units]';, yRange = [0,1000]
	oplot, R_binCenters, f_Rz_[*,zii-1]
	oplot, R_binCenters, f_Rz_[*,zii+1]

	oplot, R_binCenters, f_RzA_[*,zii], thick = 2.0
	oplot, R_binCenters, f_RzA_[*,zii-1], thick = 2.0
	oplot, R_binCenters, f_RzA_[*,zii+1], thick = 2.0

;	Contour densities in 2D

	window, 4
	loadct, 39
	!p.background = 255
	!p.multi = [0,2,2]
	device, decomposed = 0
	levels	= (fIndGen(10)+1)*5e18
	colors	= bytScl ( levels, top = 251 ) + 2
	contour, f_Rz, R_binCenters, z_binCenters, $
		levels = levels, $
		c_colors = colors, $
		color = 0, $
		charSize = 2.0, $
		/cell_fill

	levels	= (fIndGen(10)+1)*2e2
	colors	= bytScl ( levels, top = 251 ) + 2
	contour, f_Rz_, R_binCenters, z_binCenters, $
		levels = levels, $
		c_colors = colors, $
		color = 0, $
		charSize = 2.0, $
		/cell_fill


;	Write a modified particle list to use in p2f
;	based on the new R,z coords, not the ORF provided
;	coords.

	v	= sqrt ( 2.0 * data.E * e_ * 1e3 / mi )
	vPar	= data.lambda * v
	vPer	= sqrt ( v^2 - vPar^2 )	

	fName	= '/home/dg6/data/particleLists/dlg_mod_list.nc'
	nc_id	= nCdf_create ( fName, /clobber )
	nCdf_control, nc_id, /fill
	
	nP_id	= nCdf_dimDef ( nc_id, 'nP', nP )	
	R_nBins_id	= nCdf_dimDef ( nc_id, 'R_nBins', R_nBins )	
	z_nBins_id	= nCdf_dimDef ( nc_id, 'z_nBins', z_nBins )	
	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )
	
	vPer_id = nCdf_varDef ( nc_id, 'vPer', nP_id, /float )
	vPar_id = nCdf_varDef ( nc_id, 'vPar', nP_id, /float )
	R_id = nCdf_varDef ( nc_id, 'R', nP_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', nP_id, /float )
	weight_id = nCdf_varDef ( nc_id, 'weight', nP_id, /float )
	density_id = nCdf_varDef ( nc_id, 'density', [R_nBins_id, z_nBins_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, vPer_id, vPer
	nCdf_varPut, nc_id, vPar_id, vPar
	nCdf_varPut, nc_id, R_id, ROut
	nCdf_varPut, nc_id, z_id, zOut
	nCdf_varPut, nc_id, weight_id, data.weight 
	nCdf_varPut, nc_id, density_id, f_Rz;A 

	nCdf_close, nc_id
	
	stop

end
