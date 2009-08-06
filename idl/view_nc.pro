pro view_nc 

	cdfId	= ncdf_open ( 'phillips_nstx3.5.2.nc', /noWrite )
	glob	= ncdf_inquire ( cdfId )

	;	Print some info about the netcdf file

	print
	print, 'Dimensions: ', glob.ndims
	
	for i=0,glob.ndims-1 do begin

		ncdf_dimInq, cdfId, i, name, size_
		
		if i eq glob.recdim then $
			print, name, size_, '(Unlimited dim)' $
		else $
			print, name, size_
	endfor

	print
	print, 'Variables:'
	for i=0,glob.nvars-1 do begin

		info	= ncdf_varInq ( cdfId, i )
        fmtStr = '(A," (",A," ) Dimension Ids = [ ", 10(I0," "),$)'
		print, format= fmtStr, info.name, info.datatype, info.dim[*]
		print, ']'

		for j=0,info.natts-1 do begin
		
			attname = ncdf_attname ( cdfId, i, j )
			ncdf_attget, cdfId, i, attname, attvalue
			print, 'Attribute ', attname, '=', string ( attvalue )			
		
		endfor
	
	endfor

stop
end
