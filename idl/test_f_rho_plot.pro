pro test_f_rho_plot

restore, 'tmp.sav'

fName	= 'data/f_rho_vv.nc'
nc_id	= nCdf_create ( fName, /clobber )
nCdf_control, nc_id, /fill

nRho_id	= nCdf_dimDef ( nc_id, 'nRho', n_elements(f_vv_rho[*,0,0]) )	
nvPer_id	= nCdf_dimDef ( nc_id, 'nvPer', n_elements(f_vv_rho[0,*,0]) )	
nvPar_id	= nCdf_dimDef ( nc_id, 'nvPar', n_elements(f_vv_rho[0,0,*]) )	

f_rho_vv_id = nCdf_varDef ( nc_id, 'f_rho_vv', [nRho_id, nvPer_id, nvPar_id], /float )
rho_id = nCdf_varDef ( nc_id, 'rho', nRho_id, /float )

nCdf_control, nc_id, /enDef

nCdf_varPut, nc_id, f_rho_vv_id, (f_vv_rho)
nCdf_varPut, nc_id, rho_id, fIndGen( n_elements(f_vv_rho[*,0,0]) ) / (n_elements(f_vv_rho[*,0,0])-1)

nCdf_close, nc_id
stop	
;	Contour f_rho_vv

for i=0,n_elements(f_vv_rho[*,0,0])-1 do begin

	contour, transpose ( f_vv_rho[i,*,*] ), $
		levels = 2.0^fIndgen(100)*1d-2, $
		path_xy = path, path_info = info
	if max ( f_vv_rho[i,*,*] ) gt 0 then begin
		for j=0,n_elements(info.offset)-1 do begin
			iPlot, path[0,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
						path[1,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1], $
						path[1,(info.offset)[j]:(info.offset)[j]+(info.n)[j]-1]*0+i, $
						/overPlot;, transparency = (100-j*10)>0, thick = 3 
		endfor
	endif
endfor

end
