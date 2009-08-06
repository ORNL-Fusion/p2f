pro test_particle_size

amu	= 2d0
mi	= amu * 1.672615d-27
e_	= 1.602176d-19
c	= 3.0d+08
vParRange	= 3.0d-2 * c
vPerRange	= 3.0d-2 * c
sigma	= 2.0 * 1d-04 * c
nPar	= 128 
nPer	= 128 
nPhi	= 21  
density	= 1.0e19

suggestedSigma	= sqrt ( vPerRange^2/700.0 ) 

particle_keV	= sigma^2 * 0.5 * mi / e_ / 1d3
sugg_particle_keV	= suggestedSigma^2 * 0.5 * mi / e_ / 1d3

maxRange_keV	= ( vParRange^2 + vPerRange^2 ) * 0.5 * mi / e_ / 1d3

print, 'sigma: ', sigma / c *1d4
print, 'suggested sigma: ', suggestedSigma / c * 1d4
print, 'particle size in keV: ', particle_keV
print, 'suggested particle size in keV: ', sugg_particle_keV
print, 'max energy range [keV]: ', maxRange_keV 
print, 'per energy range [keV]: ', ( vPerRange^2 ) * 0.5 * mi / e_ / 1d3
print, 'par energy range [keV]: ', ( vParRange^2 ) * 0.5 * mi / e_ / 1d3
stop

;	Build in cylindrical 

vPar	= transpose ( rebin ( (dIndGen ( nPar + 1 )- nPar / 2d0) / (nPar/2d0) $
	* vParRange, nPar+1, nPer, nPhi ), [1,0,2] )
vPer	= rebin ( dIndGen ( nPer ) / nPer * vPerRange, nPer, nPar+1, nPhi ) 
phi		= transpose ( rebin ( dIndGen ( nPhi ) / (nPhi-1) * (2d0*!dpi), nPhi, nPer, nPar+1 ), [1,2,0] )

vX_	= vPer * cos ( phi )
vY_	= vPer * sin ( phi )
vZ_	= vPar

uPar	= 0d0*c
uPer	= 0d0*c
uPhi	= 0d0
uPer2	= 0.004*c

uX_	= uPer * cos ( uPhi )
uY_	= uPer * sin ( uPhi )
uX_2	= uPer2 * cos ( uPhi )
uY_2	= uPer2 * sin ( uPhi )

uZ_	= uPar

uX	= 0d0*c
uY	= 0.025d0*c
uZ	= 0.000d0*c
uY2	= 0.025d0*c

f_rpz	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX_ - uX )^2 / ( 2.0 * sigma^2 ) + $
						( vY_ - uY )^2 / ( 2.0 * sigma^2 ) + $
						( vZ_ - uZ )^2 / ( 2.0 * sigma^2 ) ) ) 
f_rpz2	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX_ - uX )^2 / ( 2.0 * sigma^2 ) + $
						( vY_ - uY2 )^2 / ( 2.0 * sigma^2 ) + $
						( vZ_ - uZ )^2 / ( 2.0 * sigma^2 ) ) ) 


dPhi	= abs ( phi[0,0,0] - phi[0,0,1] )
dvPar	= abs ( vPar[0,0,0] - vPar[0,1,0] )
dvPer	= abs ( vPer[0,0,0] - vPer[1,0,0] )

restore, '~/data/particleLists/fdis_25keV_D_4_flat_D3D.sav'

f_rpz_int	= fltArr ( size ( f_rpz, /dim ) )
f_rpz_int2	= fltArr ( size ( f_rpz, /dim ) )
f_rpz_int3	= fltArr ( size ( f_rpz, /dim ) )
f_vv_int	= fltArr ( nPer, nPar+1) 

for i = 0, 100-1 do begin

	uPar_p	=0.0; p.vPar[i]
	uPer_p	= p.vPerp[i]*2 
	uPhi_p	= 0.0
	uPhi_p2	= !pi/2d0

	print, i, uper_p, upar_p

	parii	= ( uPar_p - min ( vpar[0,*,0] ) ) / (vParRange * 2.0) * nPar
	perii	= uPer_p / vPerRange * nPer
	
	if parii lt npar and perii lt nper and parii ge 0 then $
		++f_rpz_int3[perii,parii,0]	

	uXX	= uPer_p * cos ( uPhi_p )
	uYY	= uPer_p * sin ( uPhi_p )

	uXX2	= uPer_p * cos ( uPhi_p2 )
	uYY2	= uPer_p * sin ( uPhi_p2 )
	
	uZZ	= uPar_p
	
	f_rpz_tmp	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX_ - uXX )^2 / ( 2.0 * sigma^2 ) + $
						( vY_ - uYY )^2 / ( 2.0 * sigma^2 ) + $
						( vZ_ - uZZ )^2 / ( 2.0 * sigma^2 ) ) ) 
	f_rpz_tmp2	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX_ - uXX2 )^2 / ( 2.0 * sigma^2 ) + $
						( vY_ - uYY2 )^2 / ( 2.0 * sigma^2 ) + $
						( vZ_ - uZZ )^2 / ( 2.0 * sigma^2 ) ) ) 
	
	f_vv_tmp	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( ( -uPer_p^2 - (uPar_p-vPar[*,*,0])^2-vPer[*,*,0]^2 ) / ( 2.0 * sigma^2 ) ) $
						 * 2.0 * !pi * beselI (  (uPer_p * vPer[*,*,0]/sigma^2)<700, 0.0, /double )

	;	Try adaptive phi integral range
	
	phiR	= 2.5 * sigma / uPer_p
	phi_	= ( fIndGen ( nPhi ) - nPhi / 2 ) / ( nPhi / 2 ) * phiR
	phi_	= transpose ( rebin ( phi_, nPhi, nPer, nPar + 1 ), [1,2,0] )
	vX__	= vPer * cos ( phi_ )
	vY__	= vPer * sin ( phi_ )

	f_rpz_tmp3	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX__ - uXX )^2 / ( 2.0 * sigma^2 ) + $
						( vY__ - uYY )^2 / ( 2.0 * sigma^2 ) + $
						( vZ_ - uZZ )^2 / ( 2.0 * sigma^2 ) ) ) 
	f_rpz_tmp3	= total ( f_rpz_tmp3, 3 ) * abs ( phi_[0,0,0] - phi_[0,0,1] )

	dummy = max ( f_vv_tmp[*,npar/2], iiMax )
	print, uPer_p * vPer[iiMax,npar/2,0]/sigma^2
	print, phiR 	
	f_rpz_int2	= f_rpz_int2 + f_rpz_tmp2 
	f_vv_int	= f_vv_int + f_vv_tmp 
	f_rpz_int	= f_rpz_int + f_rpz_tmp 
	plot, f_vv_tmp[*,npar/2], psym = -4
	oPlot, f_rpz_tmp3[*,npar/2]
	stop

endfor


;	 Build a vPer/vPar function

;f_vv	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
;	* exp ( ( -uPer^2 - (uPar-vPar[*,*,0])^2-vPer[*,*,0]^2 ) / ( 2.0 * sigma^2 ) ) $
;						 * beselI (  uPer * vPer[*,*,0]/sigma^2, 0.0, /double )
;
;f_vv2	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
;	* exp ( ( -uPer2^2 - (uPar-vPar[*,*,0])^2-vPer[*,*,0]^2 ) / ( 2.0 * sigma^2 ) ) $
;						 * beselI ( uPer2 * vPer[*,*,0]/sigma^2, 0.0, /double )


;f_rpz	= f_rpz * dvPer * dvPar * dPhi * vPer

;	Build in cartesian

nX	= 128 
nY	= 128 
nZ	= 128 

xRange = 0.01*c
yRange = 0.01*c
zRange = 0.01*c

vX	= transpose ( rebin ( (dIndGen ( nX + 1 )- nX / 2d0) / (nX/2d0) * xRange, nX+1, nY+1, nZ+1 ), [0,1,2] )
vY	= transpose ( rebin ( (dIndGen ( nY + 1 )- nY / 2d0) / (nY/2d0) * yRange, nY+1, nX+1, nZ+1 ), [1,0,2] )
vZ	= transpose ( rebin ( (dIndGen ( nZ + 1 )- nZ / 2d0) / (nZ/2d0) * zRange, nZ+1, nX+1, nY+1 ), [2,1,0] )

f_xyz	= density / ( ( 2.0 * !pi ) ^ 3.0/2.0 * sigma^3 ) $
	* exp ( -1.0 * ( ( vX - uX )^2 / ( 2d0 * sigma^2 ) + $
						( vY - uY )^2 / ( 2d0 * sigma^2 ) + $
						( vZ - uZ )^2 / ( 2d0 * sigma^2 ) ) )

dX	= abs ( vX[0,0,0] - vX[1,0,0] )
dY	= abs ( vY[0,0,0] - vY[0,1,0] )
dZ	= abs ( vZ[0,0,0] - vZ[0,0,1] )

;f_xyz	= f_xyz * dx * dz * dy

!p.multi = [0,3,2]
!p.charsize = 2.0
window, 0, xSize = 900, ySize = 600
contour, transpose ( f_rpz[*,*,0] ), vPar(0,*,0), vPer(*,0,0), $
	levels = 10.0^fIndGen(12)*1e-14
plot, vper[*,0,0],total(f_rpz[*,nPer/2.0,*],3),psym=-4
oplot, vper[*,0,0],total(f_rpz2[*,nPer/2.0,*],3),psym=-4

contour, transpose ( f_vv_int[*,*] ), vPar(0,*,0), vPer(*,0,0), $
	levels = 10.0^fIndGen(18)*1e-10
plot, vper[*,0,0],(total(f_rpz_int,3))[*,nPer/2.0] * dPhi
oplot, vper[*,0,0],(total(f_rpz_int2,3))[*,nPer/2.0] * dPhi
oplot, vper[*,0,0],f_vv_int[*,nPer/2.0], thick = 2.0
;plot, vper[*,0,0],total(f_vv_int,2)

a1	= total ( (total(f_rpz_int3,3))[1:*,nPer/2.0] / vPer[1:*,0,0] )
a2	= total ( (total(f_rpz_int2,3))[*,nPer/2.0] ) 
oplot, vper[*,0,0],(total(f_rpz_int3,3))[*,nPer/2.0] / vPer[*,0,0] / a1 * a2

contour, total(f_xyz,3), vX[*,0,0], vY[0,*,0], levels = 10.0^fIndGen(12)*1e-14
plot, vX[*,0,0], f_xyz[*,nY/2,nZ/2], psym=-4

;	Try fancy integration

tmp1d	= fltArr ( nPer )
for i = 0, nPer - 1 do begin
	tmp2d	= fltArr ( nPer, nPar+1 )
	for j = 0, nPar +1 - 1 do begin

		tmp2d[i,j]	= int_tabulated ( reform(phi[0,0,*]), reform(f_rpz[i,j,*]) )

	endfor
	tmp1d[i]	= int_tabulated ( reform(vPar[0,*,0]), reform(tmp2d[i,*]) )
endfor
int_rpz	= int_tabulated ( reform(vPer[*,0,0]), tmp1d * vPer )


tmp1d	= fltArr ( nPer )
for i = 0, nPer - 1 do begin
	tmp2d	= fltArr ( nPer, nPar+1 )
	for j = 0, nPar +1 - 1 do begin

		tmp2d[i,j]	= int_tabulated ( reform(phi[0,0,*]), reform(f_rpz2[i,j,*]) )

	endfor
	tmp1d[i]	= int_tabulated ( reform(vPar[0,*,0]), reform(tmp2d[i,*]) )
endfor
int_rpz2	= int_tabulated ( reform(vPer[*,0,0]), tmp1d * vPer )


tmp1d	= fltArr ( nX+1 )
for i = 0, nX +1- 1 do begin
	tmp2d	= fltArr ( nX+1, nY+1 )
	for j = 0, nY +1- 1 do begin

		tmp2d[i,j]	= int_tabulated ( reform(vZ[0,0,*]), reform(f_xyz[i,j,*]) )

	endfor
	tmp1d[i]	= int_tabulated ( reform(vY[0,*,0]), reform(tmp2d[i,*]) )
endfor
int_xyz	= int_tabulated ( reform(vX[*,0,0]), tmp1d )

print, int_xyz, int_rpz, int_rpz2
print, total(f_xyz*dX*dY*dZ),total(f_rpz*dvPer*dvPar*dPhi*vPer)
print, abs(int_xyz-int_rpz)/int_xyz*100

stop

end
