pro plot_maptoeqs, psiIn, thetaIn, rOut, zOut, eqdsk, fileName, $
	rz2pt = rz2pt, $ ; keyword for going from rOut/zOut -> psiIn/thetaIn
	rIn	= rIn, $
	zIn = zIn

	openr, unit, fileName, /get_lun

	readf, unit, nPsi, nTheta
	readf, unit, xmac, rrc
	readf, unit, n,m,l

	xFit_p	= fltArr ( nPsi )
	readf, unit, xFit_p

	qValCh_p	= fltArr ( nPsi )
	readf, unit, qValCh_p
	
	xFit_g	= fltArr ( nPsi )
	readf, unit, xFit_g

	;psiVal	= fltArr ( nPsi )
	;readf, unit, psiVal 

	xFit_psiVal	= fltArr ( nPsi )
	readf, unit, xFit_psiVal

	x	= fltArr ( nPsi, nTheta+5 )
	readf, unit, x

	z	= fltArr ( nPsi, nTheta+5 )
	readf, unit, z

	grPsSq	= fltArr ( nPsi, nTheta+5 )
	readf, unit, grPsSq 

	xJacob	= fltArr ( nPsi, nTheta+5 )
	readf, unit, xJacob 

	readf, unit, bNorm, xNorm

	free_lun, unit

	print, nPsi, nTheta
	print, xmac, rrc
	print, n, m, l
	print, bNorm, xNorm

	!p.multi = [0,2,3]
	!p.charSize = 2.0
	plot, xFit_p, title='xFit_p'
	plot, qvalCh_p, title='qValCh_p'
	plot, xFit_g, title='xFit_g'
	;plot, psiVal, title='psiVal'
	plot, xFit_psiVal, title='xFit_psiVal'

	;psiRange	= eqdsk.siBry - eqdsk.siMag
	;psiMin	= eqdsk.siMag

	psiHalf	= 1.0 / nPsi / 2.0
	psiMin	= psiHalf
	psiRange	= 1.0 - 2.0*psiHalf

	;psiMin = 0.0
	;psiRange = 1.0

	psiVal	= fIndGen ( nPsi ) * ( psiHalf * 2.0 ) + psiHalf
	theta	= fIndGen ( nTheta+5 ) / (nTheta) * 2.0 * !pi

	;psiI	= ( psiIn - (psiMin+psiHalf) ) / (psiRange+psiHalf*2.0) * nPsi

if keyword_set ( rz2pt ) then begin

;	Search for closest R/z point in the 2D interpolation
;	meshes we have and use that.
	
	psiIn	= fltArr ( n_elements ( rIn ) )
	thetaIn	= fltArr ( n_elements ( rIn ) )

	for i = 0L, n_elements ( rIn ) - 1 do begin

		if ( i mod 1000 ) eq 0 then print, i

		rDiff	= rIn[i]/xmac - x
		zDiff	= zIn[i]/rrc - z
	
		distance	= sqrt ( rDiff^2 + zDiff^2 )
		
		iiUse	= where ( distance eq min ( distance ), iiUseCnt )
		if iiUseCnt ge 1 then begin
	
			arrayIndex	= array_indices ( x, iiUse[0] )
			;psiIn[i]	= xFit_psiVal[arrayIndex[0]]
			psiIn[i]	= 0.5 / nPsi + arrayIndex[0] / (nPsi-1) * ( 1.0 - 1.0 / nPsi )
			thetaIn[i]	= theta[arrayIndex[1]]
		
		endif else begin

			stop
		
		endelse

	endfor

endif 

;	Interpolate from given psi, theta coords to 
;	R,z coords

	psiI	= ( psiIn - (psiMin) ) / (psiRange) * nPsi
	thetaJ	= thetaIn / ( 2.0 * !pi ) * (nTheta)

	rOut	= interpolate ( x, psiI, thetaJ )
	zOut	= interpolate ( z, psiI, thetaJ )
	rzJacob	= interpolate ( xJacob, psiI, thetaJ )

	;iiNearCore	= where ( psiI lt 0, iiNCCnt )
	;if iiNCCnt gt 0 then begin

	;	for i = 0L, iiNCCnt - 1 do begin

	;		coeffs	= poly_fit ( psiVal[0:10], x[0:10,thetaJ[iiNearCore[i]]], 5 )
	;		rOut[iiNearCore[i]]	= coeffs[0] + $
	;								coeffs[1] * psiIn[iiNearCore[i]] + $
	;								coeffs[2] * psiIn[iiNearCore[i]]^2 + $
	;								coeffs[3] * psiIn[iiNearCore[i]]^3 + $
	;								coeffs[4] * psiIn[iiNearCore[i]]^4 + $
	;								coeffs[5] * psiIn[iiNearCore[i]]^5
	;
	;		coeffs	= poly_fit ( psiVal[0:10], z[0:10,thetaJ[iiNearCore[i]]], 5 )
	;		zOut[iiNearCore[i]]	= coeffs[0] + $
	;								coeffs[1] * psiIn[iiNearCore[i]] + $
	;								coeffs[2] * psiIn[iiNearCore[i]]^2 + $
	;								coeffs[3] * psiIn[iiNearCore[i]]^3 + $
	;								coeffs[4] * psiIn[iiNearCore[i]]^4 + $
	;								coeffs[5] * psiIn[iiNearCore[i]]^5

	;	endFor

	;endIf

	rOut	= rOut * xmac
	zOut	= zOut * rrc

end
