pro test_hist, subSet = subSet, lowEnergy = lowEnergy

	fileName	= 'fdis0-newtry'

	n	= file_lines ( fileName )

	data = { 	R : 0.0, $
				z : 0.0, $
				E : 0.0, $
				g : 0.0, $
				w : 0.0 }

	data	= replicate ( data, n )

	openr, lun, fileName, /get_lun
	readf, lun, data
	free_lun, lun

	data.R	= data.R * 1e-2
	data.z	= data.z * 1e-2

	ePar	= data.E * data.g
	ePer	= sqrt ( data.E^2 - ePar^2 )
	eMag	= data.E

;	extract a subsect above lowEnergy

	if keyword_set ( lowEnergy ) then begin

		iiSubSet	= where ( eMag gt lowEnergy, iiCntE )
		if iiCntE gt 0 then begin

			ePar	= ePar[iiSubSet]
			ePer	= ePer[iiSubSet]
			eMag	= eMag[iiSubSet]
	
		endif

	endif

;	extract a subset of the data if /subSet keyword is present

	iiSubSet	= where ( data.R gt 1.5 and data.R lt 1.8 $
			and data.z gt -0.05 and data.z lt 0.05, iiCnt )

	if iiCnt gt 0 and keyword_set ( subSet ) then begin
		ePar	= ePar[iiSubSet]
		ePer	= ePer[iiSubSet]
		eMag	= eMag[iiSubSet]
	endif


;	plot histograms at various resolutions

	!p.multi = [0,1,3]
	!p.charSize = 2
	set_plot, 'X'
	device, decomposed = 0
	!p.background = 255
	loadct, 12, /sil

	nBins	= [50,100,200,500]

	for i=0,n_elements(nBins)-1 do begin
		hist	= histogram( eMag, nBins = nBins[i], loca = locs, max = 100, min = 0)
		if i eq 0 then plot, locs, hist*nBins[i], xRange =[-10, 100], color=0 $
		else oplot, locs, hist*nBins[i], color = (4*i)*16-1
	endfor

	for i=0,n_elements(nBins)-1 do begin
		hist	= histogram( ePer, nBins = nBins[i], loca = locs, max = 100, min = 0 )
		if i eq 0 then plot, locs, hist*nBins[i], xRange = [-10,100], color=0 $
		else oPlot, locs, hist*nBins[i], color = (4*i)*16-1
	endfor

	for i=0,n_elements(nBins)-1 do begin
		hist	= histogram( ePar, nBins = nBins[i], loca = locs, max = 100, min = -100 )
		if i eq 0 then plot, locs, hist*nBins[i], xRange = [-100,100], color=0 $
		else oPlot, locs, hist*nBins[i], color = (4*i)*16-1
	endfor

	stop

end
