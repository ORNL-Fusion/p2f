function p2f_read_nml, nmlFile

	openr, lun, nmlFile, /get_lun
	runFileArray = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
		runFileArray = [runFileArray,line]
	endwhile
	free_lun, lun

	for f=0,n_elements(runFileArray)-1 do begin
		if(strMatch(runFileArray[f],'*eqdsk_fileName*'))then eqdsk_fName=(strSplit(runFileArray[f],'"',/extract))[1]
	endfor

	cfg = create_struct ( name='kjCfg', $
			'eqdsk_fName', eqdsk_fName )

	return, cfg
end 


end
