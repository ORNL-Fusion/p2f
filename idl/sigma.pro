pro sigma

	cdfId	= ncdf_open ( 'f_rzvv.nc', /noWrite )
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
;	ncdf_varGet, cdfId, 'R_binSize', R_binSize
;	ncdf_varGet, cdfId, 'z_binSize', z_binSize
;	ncdf_varGet, cdfId, 'vPerp_binSize', vPerp_binSize
;	ncdf_varGet, cdfId, 'vPar_binSize', vPar_binSize

	ncdf_close, cdfId

	R_nBins	= n_elements ( R_binCenters )
	z_nBins	= n_elements ( z_binCenters )	
	vPerp_nBins	= n_elements ( vPerp_binCenters )
	vPar_nBins	= n_elements ( vPar_binCenters )	

	vPerp_binSize	= vPerp_binCenters[1]-vPerp_binCenters[0]
	vPar_binSize	= vPar_binCenters[1]-vPar_binCenters[0]
	R_binSize	= R_binCenters[1]-R_binCenters[0]
	z_binSize	= abs ( z_binCenters[1]-z_binCenters[0] )

	w	= 30.0e6 * 2.0 * !pi
	kPar_	= 100.0 
	kPerp_	= 100.0 
	c	= 3.0e8
	m	= 2.0 * 1.67262e-27
	e	= 1.602176e-19
	eNorm	= 100.0e3
	n	= 1.0e19
	
	;	Normalisation variables, inherited from Harvey and CQL3D
	;	Stuff that, follow stix eqn. 10.38 to start.

	mu	= m * c^2 / ( 2.0 * e * eNorm )
	vc	= c / sqrt ( mu )

	;	Pick some test location

	rPt	= 12
	zPt	= 25

	f_vv	= reform ( f_rzvv[rPt,zPt,*,*] )


	nf_f	= n / vc^3 / ( total ( f_vv ) * $
		vPerp_binSize * vPar_binSize )

	;f_vv	= f_vv * nf_f

	vPerp_2D	= rebin ( vPerp_binCenters, vPerp_nBins, vPar_nBins )
	vPar_2D		= transpose ( rebin ( vPar_binCenters, vPar_nBins, vPerp_nBins ) )

	nf_v	= vc / max ( sqrt ( vPerp_2D^2 + vPar_2D^2 ) )

	uPerp_binCenters	= vPerp_binCenters
	uPar_binCenters	= vPar_binCenters
	uPerp_binSize	= vPerp_binSize 
	uPar_binSize	= vPar_binSize 
	uPerp_2D	= vPerp_2D 
	uPar_2D	= vPar_2D 

		;	Test maxwellian

	vSigma	= 7e5 
	vA	= max ( f_vv )	
	f_vv_	= vA * exp (-( uPerp_2D^2/(2.0*vSigma^2) + uPar_2D^2/(2.0*vSigma^2) ) )

	;f_vv	= f_vv_

	;	Calculate derivatives

	dfduPerp	= dlg_pDeriv ( f_vv, 1, uPerp_binSize )
	dfduPar		= dlg_pDeriv ( f_vv, 2, uPar_binSize )

	dfduPerp_	= dlg_pDeriv ( f_vv_, 1, uPerp_binSize )
	dfduPar_	= dlg_pDeriv ( f_vv_, 2, uPar_binSize )

	;	Calculate susceptability as a function of kPerp and kPar
	;	for some frequency @ rPt, zPt

	q	= e
	B	= 1.5
	omega	= q * B / ( m )

	kPerp	= (fIndGen ( 100 )+1.0) * 10.0
	kPar	= (fIndGen ( 20 )+1.0) * 50.0

	kPerp_nBins	= n_elements ( kPerp )
	kPar_nBins	= n_elements ( kPar )

	chi	= fltArr ( kPerp_nBins, kPar_nBins )
	chi_	= fltArr ( kPerp_nBins, kPar_nBins )

	lMax	= 10

	Jm1	= fltArr ( kPerp_nBins, lMax*2+1, vPerp_nBins, vPar_nBins )
	J_	= fltArr ( kPerp_nBins, lMax*2+1, vPerp_nBins, vPar_nBins )
	Jp1	= fltArr ( kPerp_nBins, lMax*2+1, vPerp_nBins, vPar_nBins )
	dJ	= fltArr ( kPerp_nBins, lMax*2+1, vPerp_nBins, vPar_nBins )
	T11	= fltArr ( kPerp_nBins, lMax*2+1, vPerp_nBins, vPar_nBins )

	;	Calculate besel functions outside of the kPar loop

	print, 'Calculating besel functions ... '
	for i=0,kPerp_nBins-1 do begin

		x	= kPerp[i] * uPerp_2D / omega

		for l=-lMax,lMax do begin
			
			Jm1[i,l+lMax,*,*]	= beselJ ( x, abs ( l-1 ) )
			J_[i,l+lMax,*,*]	= beselJ ( x, abs ( l ) ) 
			Jp1[i,l+lMax,*,*]	= beselJ ( x, abs ( l+1 ) )
		
			if l lt 0 then begin
				Jm1[i,l+lMax,*,*]	= (-1)^abs ( l-1 ) * Jm1[i,l+lMax,*,*]
				J_[i,l+lMax,*,*]	= (-1)^abs ( l ) * J_[i,l+lMax,*,*]
				Jp1[i,l+lMax,*,*]	= (-1)^abs ( l+1 ) * Jp1[i,l+lMax,*,*]
			endif

			if abs(l) ge 1 then begin
				dJ[i,l+lMax,*,*]	= Jm1[i,l+lMax,*,*] - J_[i,l+lMax,*,*] * l / x
			endif else begin
				dJ[i,l+lMax,*,*]	= l / x * J_[i,l+lMax,*,*] - Jp1[i,l+lMax,*,*]
			endelse
	
			T11[i,l+lMax,*,*]	= uPerp_2D^2 * l^2 * J_[i,l+lMax,*,*]^2 / x^2

		endfor
	endfor

	;	Calculate U outside of the kPerp loop

	U	= fltArr ( kPar_nBins, vPerp_nBins, vPar_nBins )
	U_	= fltArr ( kPar_nBins, vPerp_nBins, vPar_nBins )

	for j=0,kPar_nBins-1 do begin
	
		U[j,*,*]	= dfduPerp + kPar[j] / w * ( uPerp_2D * dfduPar - uPar_2D * dfduPerp )
		U_[j,*,*]	= dfduPerp_ + kPar[j] / w * ( uPerp_2D * dfduPar_ - uPar_2D * dfduPerp_ )

	endfor

	for i=0,kPerp_nBins-1 do begin
	
		for j=0,kPar_nBins-1 do begin
			
			for l=-lmax,lMax do begin

				int_dVPerp	= fltArr ( vPar_nBins )
				for n=0,vPar_nBins-1 do begin	
					int_dVPerp[n]	= int_tabulated ( uPerp_binCenters, U[j,*,n] * T11[i,l+lMax,*,n] )
				endfor

				int_dVPar	= int_tabulated ( uPar_binCenters, $
					int_dVPerp / ( w - kPar[j] * uPar_binCenters - l * omega ) ) 
			
				chi[i,j]	= chi[i,j] + int_dVPar

				int_dVPerp	= fltArr ( vPar_nBins )
				for n=0,vPar_nBins-1 do begin	
					int_dVPerp[n]	= int_tabulated ( uPerp_binCenters, U_[j,*,n] * T11[i,l+lMax,*,n] )
				endfor

				int_dVPar	= int_tabulated ( uPar_binCenters, $
					int_dVPerp / ( w - kPar[j] * uPar_binCenters - l * omega ) ) 
			
				chi_[i,j]	= chi_[i,j] + int_dVPar

			endfor
			print, l, kPerp[i], kPar[j], chi[i,j], min ( abs ( w - kPar[j] * uPar_binCenters - l * omega ) )
		endfor
	endfor

	;	Plot some info

	!p.multi = [0,3,4]
	window, xSize = 1000, ySize = 1200
	@daves_paths
	set_plot, 'X'
	device, decomposed = 0
	!p.background = 0
	!p.charSize	= 2.0
	loadct, 13, file = colortb_path
	levels	= 2.0 ^ fIndGen ( 10 ) * 1e2
	dLevels	= 2.0 ^ fIndGen ( 10 ) * 1e-4
	dLevels	= [ -reverse ( dLevels ), dLevels ]
	colors	= bytScl ( fIndGen ( 10 ), top = 254/2.0 ) + 254/2.0
	dColors	= bytScl ( fIndGen ( 20 ), top = 253 ) + 1

	contour, transpose ( f_vv ), uPar_binCenters, uPerp_binCenters, $
		levels = levels, c_colors = colors, color = 255, charSize = 2.0
	contour, transpose ( dfduPerp ), uPar_binCenters, uPerp_binCenters, $
		color = 255, levels = dLevels, c_colors = dColors, charSize = 2.0
	contour, transpose ( dfduPar ), uPar_binCenters, uPerp_binCenters, $
		color = 255, levels = dLevels, c_colors = dColors, charSize = 2.0
	contour, transpose ( uPerp_2D * dfduPar ), uPar_binCenters, uPerp_binCenters, $
		color = 255, levels = dLevels*1e5, c_colors = dColors, charSize = 2.0
	contour, transpose ( uPar_2D * dfduPerp ), uPar_binCenters, uPerp_binCenters, $
		color = 255, levels = dLevels*1e5, c_colors = dColors, charSize = 2.0
	contour, transpose ( uPar_2D * dfduPerp - uPerp_2D * dfduPar ), uPar_binCenters, uPerp_binCenters, $
		color = 255, levels = dLevels*1e5, c_colors = dColors, charSize = 2.0

	plot, uPar_binCenters, uPar_binCenters * total ( dfduPerp, 1 ) * uPerp_binSize, color = 255
	oplot, uPar_binCenters, uPar_binCenters * total ( dfduPerp_, 1 ) * uPerp_binSize, color = 250, thick = 2.0

	plot, uPar_binCenters, dlg_pDeriv ( total ( uPerp_2D * f_vv, 1 ), 1, uPar_binSize ) * uPerp_binSize, color = 255
	oplot, uPar_binCenters, dlg_pDeriv ( total ( uPerp_2D * f_vv_, 1 ), 1, uPar_binSize ) * uPerp_binSize,$
		color = 250, thick = 2.0 
	;oPlot, uPar_binCenters, total ( uPerp_2D * dfduPar_binSize, 1 ) * uPerp_binSize, color = 50

	plot, uPar_binCenters, total ( uPar_2D * dfduPerp - uPerp_2D * dfduPar, 1 ) * uPerp_binSize, color = 255
	;oPlot, uPar_binCenters, uPar_binCenters * total ( dfduPerp_binSize, 1 ) * uPerp_binSize $
	;	- dlg_pDeriv ( total ( uPerp_2D * f_vv, 1 ), 1, uPar_binSize ) * uPerp_binSize, color = 50
	oplot, uPar_binCenters, total ( uPar_2D * dfduPerp_ - uPerp_2D * dfduPar_, 1 ) * uPerp_binSize,$
		color = 250, thick = 2.0

	xyOuts, 2.0e6*nf_v,1.5e-7*nf_v, $
		string( total (total ( uPar_2D * dfduPerp - uPerp_2D * dfduPar, 1 ) * uPerp_binSize ),$
		for = '(e7.1)' ),$
		color = 255, /data, charSize = 1.0
	xyOuts, 2.0e6*nf_v,1.0e-7*nf_v, $
		string( total (total ( uPar_2D * dfduPerp_ - uPerp_2D * dfduPar_, 1 ) * uPerp_binSize ),$
		for = '(e7.1)' ),$
		color = 250, /data, charSize = 1.0, charThick = 2.0
	
	plot, uPar_binCenters, total ( dfduPerp, 1 ) * uPerp_binSize, color = 255
	oplot, uPar_binCenters, kPar_ / w * uPar_binCenters * total ( dfduPerp, 1 ) * uPerp_binSize,$
		color = 255
	oplot, uPar_binCenters,$
		kPar_ / w * dlg_pDeriv ( total ( uPerp_2D * f_vv, 1 ), 1, uPar_binSize ) * uPerp_binSize,$
		color = 255
	oplot, uPar_binCenters, kPar_ / w * $
		( uPar_binCenters * total ( dfduPerp, 1 ) * uPerp_binSize $
			- dlg_pDeriv ( total ( uPerp_2D * f_vv, 1 ), 1, uPar_binSize ) * uPerp_binSize ),$
		color = 255

	!p.multi = 0

	window, 3
	
	dLevels	= 2.0 ^ fIndGen ( 10 ) * 1e-3
	dLevels	= [ -reverse ( dLevels ), dLevels ]
	dColors	= bytScl ( fIndGen ( 20 ), top = 253 ) + 1


	contour, chi, c_colors=dColors, levels=dLevels, /fill
	iiPos	= where(dLevels gt 0 )
	dColors[iiPos]	= (dColors[iiPos]+20)<254
	iiNeg	= where(dLevels lt 0 )
	dColors[iiNeg]	= (dColors[iiNeg]-20)>1
	contour, chi_, c_colors=((dColors*1.2)<254)>1, levels=dLevels, /over

stop
end
