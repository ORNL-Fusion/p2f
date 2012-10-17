pro p2f_convert_nubeam

	runIdent = '69247B08'
	;runIdent = '69247M09'

	fileList = file_search(runIdent+'_debug_nbi_ptcl_state_cpu*.cdf')

	complete_nbi_r = !null
	complete_nbi_z = !null
	complete_nbi_vPer = !null
	complete_nbi_vPar = !null
	complete_nbi_w = !null

	for f=0,n_elements(fileList)-1 do begin
		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'rmjionay', R ; cm
			ncdf_varget, cdfId, 'xzionay', z ; cm
			ncdf_varget, cdfId, 'vay', vAbs ; cm/s
			ncdf_varget, cdfId, 'xksidy', vPar_v
			ncdf_varget, cdfId, 'wghtay', w

		nCdf_close,	cdfId 

		fileName = runIdent+'_debug_nbi_fld_state.cdf'
		cdfId = ncdf_open(fileName)

			ncdf_varget, cdfId, 'nlfprod', nlfprod 

		nCdf_close,	cdfId 

		for i=0,n_elements(nlfprod)-1 do begin
				if nlfprod[i] ne 0 then begin
						print, 'Not beam ion species, must be fusion species?'
						stop
				endif
		endfor


		fileName = runIdent+'_ps_ts1_state.cdf'
		cdfId = ncdf_open(fileName)

			ncdf_varget, cdfId, 'qatom_ALL', qatom_ALL 
			ncdf_varget, cdfId, 'm_ALL', m_ALL 

		nCdf_close,	cdfId 

		ps_mp  = 1.6726e-27
		ps_xe  = 1.6022e-19

		_Z = qatom_all/ps_xe
		amu = m_all/ps_mp

		fileName = runIdent+'_ps_ts2_state.cdf'
		cdfId = ncdf_open(fileName)

			ncdf_varget, cdfId, 'snbi_to_all', snbi_to_all 
			;ncdf_varget, cdfId, 'sfus_to_all', sfus_to_all 

		nCdf_close,	cdfId 

		nbiZ = _Z[snbi_to_all]
		nbiamu = amu[snbi_to_all]
		idx = 0
		nbi_r = R[*,idx]*1e-2
		nbi_z = z[*,idx]*1e-2
		nbivAbs = vAbs[*,idx]*1e-2
		nbivPar_v = vPar_v[*,idx]
		nbi_w = w[*,idx]
		nbi_vPar = nbivAbs*nbivPar_v
		nbi_vPer = sqrt(nbivAbs^2-nbi_vPar^2)

		iiGood = where(nbi_w gt 0.0,iiGoodCnt)
		print, 'Number of non-zero weight particles: ', iiGoodCnt

		nbi_r = nbi_r[iiGood]
		nbi_z = nbi_z[iiGood]
		nbi_vPer = nbi_vPer[iiGood]
		nbi_vPar = nbi_vPar[iiGood]
		nbi_w = nbi_w[iiGood]

		complete_nbi_r    = [complete_nbi_r, nbi_r] 
		complete_nbi_z    = [complete_nbi_z, nbi_z] 
		complete_nbi_vPer = [complete_nbi_vPer, nbi_vPer] 
		complete_nbi_vPar = [complete_nbi_vPar, nbi_vPar] 
		complete_nbi_w    = [complete_nbi_w, nbi_w] 

	endfor

	nP = n_elements(complete_nbi_r)

	print, 'Total nP: ', nP

	; Create sMC compatible particle list

	outputFileName='nubeam_smc.nc'

	nc_id = nCdf_create (outputFileName, /clobber )
		nCdf_control, nc_id, /fill
		
		nP_id = nCdf_dimDef ( nc_id, 'nP', nP )
		scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

		r_id = nCdf_varDef ( nc_id, 'R', nP_id, /float )
		z_id = nCdf_varDef ( nc_id, 'z', nP_id, /float )
		vPer_id = nCdf_varDef ( nc_id, 'vPer', nP_id, /float )
		vPar_id = nCdf_varDef ( nc_id, 'vPar', nP_id, /float )
		w_id = nCdf_varDef ( nc_id, 'weight', nP_id, /float )

		nCdf_control, nc_id, /enDef

		nCdf_varPut, nc_id, r_id, complete_nbi_r
		nCdf_varPut, nc_id, z_id, complete_nbi_z
		nCdf_varPut, nc_id, vPer_id, complete_nbi_vPer
		nCdf_varPut, nc_id, vPar_id, complete_nbi_vPar
		nCdf_varPut, nc_id, w_id, complete_nbi_w

	nCdf_close, nc_id

	stop	

end 
