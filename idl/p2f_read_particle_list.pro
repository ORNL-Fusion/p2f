pro p2f_read_particle_list

    fName = 'plist.cdf'
	cdfId = ncdf_open(fName)

		ncdf_varget, cdfId, 'R', r
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 'vPer', vPer 
		ncdf_varget, cdfId, 'vPar', vPar 
		ncdf_varget, cdfId, 'weight', weight 

	nCdf_close,	cdfId 

    dr = 0.025
    dz = 0.05
    h = hist_2d(r,z,bin1=dr,bin2=dz) 
    c=contour(h,/fill,rgb_table=51)

stop
end
