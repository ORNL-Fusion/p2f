pro p2f_read_particle_list, fName = _fName

	if keyword_set(_fName) then fName = _fName else fName = 'plist.cdf'
	cdfId = ncdf_open(fName)

		ncdf_varget, cdfId, 'R', r
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 'vPer', vPer 
		ncdf_varget, cdfId, 'vPar', vPar 
		ncdf_varget, cdfId, 'weight', weight 

	nCdf_close,	cdfId 

    eqdskFName = 'eqdsk'
    g=readgeqdsk(eqdskFName,/noTor)

    MinR = min(g.rbbbs) 
    MaxR = max(g.rbbbs) 
    MinZ = min(g.zbbbs) 
    MaxZ = max(g.zbbbs) 

    nR = 30
    nZ = 30

    HistR = fIndGen(nR)/(nR-1)*(MaxR-MinR)+MinR
    HistZ = fIndGen(nZ)/(nZ-1)*(MaxZ-MinZ)+MinZ

    dr = HistR[1]-HistR[0]
    dz = HistZ[1]-HistZ[0]

    h = hist_2d(r,z,bin1=dr,bin2=dz,min1=MinR,max1=MaxR,min2=MinZ,max2=MaxZ) 
    Hist = fltArr(nR,nZ)
    Hist[0:n_elements(h[*,0])-1,0:n_elements(h[0,*])-1] = h
    c=contour(Hist,HistR,HistZ,/fill,rgb_table=51,layout=[2,1,1],aspect_ratio=1.0)
    i=image(Hist,margin=0.05, axis_style=2, rgb_table=51, layout=[2,1,2], /current)
stop
end
