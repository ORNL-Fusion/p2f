function interpB,  bStruct, rPos, zPos

    bRHere  = interpolate ( bStruct.bR, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    bPhiHere  = interpolate ( bStruct.bPhi, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    bzHere  = interpolate ( bStruct.bz, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    
    bMag    = sqrt ( bRHere^2 + bPhiHere^2 + bzHere^2 )

    bOut    = { bR : bRHere, $
                bPhi :bPhiHere, $
                bz : bzHere, $
                bMag : bMag }

    return, bOut

end 

function readGeqdsk, fileName, $
	plot_ = plot_, $
	pAngle = pAngle, $
	half = half, $
	reWrite = reWrite

;   Read in data from g-eqdsk file

openr, lun, fileName, /get_lun

case_    = strArr ( 6 )

f1  = '(6a8,3i4)'
f2  = '(5e16.9)'
f3  = '(2i5)'

readf, lun, format = f1, case_, idum, nW, nH
readf, lun, format = f2, rdim, zdim, rcentr, rleft, zmid
readf, lun, format = f2, rmaxis, zmaxis, simag, sibry, bcentr
readf, lun, format = f2, current, simag, xdum, rmaxis, xdum
readf, lun, format = f2, zmaxis, xdum, sibry, xdum, xdum

fpol    = fltArr ( nW )
pres    = fltArr ( nW )
ffprim  = fltArr ( nW )
pprime  = fltArr ( nW )
psizr   = fltArr ( nW, nH )
qpsi    = fltArr ( nW )

readf, lun, format = f2, fpol
readf, lun, format = f2, pres
readf, lun, format = f2, ffprim 
readf, lun, format = f2, pprime 
readf, lun, format = f2, psizr 
readf, lun, format = f2, qpsi

readf, lun, format = f3, nbbbs, limitr

bbbs    = fltArr ( 2, nbbbs )
lim    = fltArr ( 2, limitr )

readf, lun, format = f2, bbbs
readf, lun, format = f2, lim

rbbbs   = bbbs[0,*]
zbbbs    = bbbs[1,*]

rlim    = lim[0,*]
zlim    = lim[1,*]

if keyword_set ( half ) then begin

	psizr	= psizr / 2.0
	fpol	= fpol / 2.0
	simag	= simag / 2.0
	sibry	= sibry / 2.0

endif

if keyword_set ( reWrite ) then begin

	openw, lun, fileName+'.dlgMod', /get_lun
	
	printf, lun, format = f1, case_, idum, nW, nH
	printf, lun, format = f2, rdim, zdim, rcentr, rleft, zmid
	printf, lun, format = f2, rmaxis, zmaxis, simag, sibry, bcentr
	printf, lun, format = f2, current, simag, xdum, rmaxis, xdum
	printf, lun, format = f2, zmaxis, xdum, sibry, xdum, xdum
	printf, lun, format = f2, fpol
	printf, lun, format = f2, pres
	printf, lun, format = f2, ffprim 
	printf, lun, format = f2, pprime 
	printf, lun, format = f2, psizr 
	printf, lun, format = f2, qpsi
	printf, lun, format = f3, nbbbs, limitr
	printf, lun, format = f2, bbbs
	printf, lun, format = f2, lim

endif

;   Calculate other desired quantities

rStep   = rdim / ( nW - 1 )
zStep   = zdim / ( nH - 1 )
fStep   = -( simag - sibry ) / ( nW - 1 )

R   = fIndGen ( nW ) * rStep + rleft
z   = fIndGen ( nH ) * zStep + zmid - zdim / 2.0

fluxGrid    = fIndGen ( nW ) * fStep + simag

;	Adjust to an up down symmetric field

;psizr	= ( psizr + reverse ( psizr, 1 ) ) / 2.0

;	Remember psi = -R * A

bR  = -1.0 * dlg_pDeriv ( psizr, 2, zStep ) / rebin ( R, nW, nH )
bz  = dlg_pDeriv ( psizr, 1, rStep ) / rebin ( R, nW, nH )
print, max(bR), min(bR)
print, max(bz), min(bz)
fPol_spline = spl_init ( fluxGrid, fPol )

fPolRZ  = reform ( spl_interp ( fluxGrid, fPol, fPol_spline, psizr[*] ), nW, nH )
bPhi    = fPolRZ / rebin ( R, nW, nH )
APhi    = -psizr / rebin ( R, nW, nH )

bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

buR = bR / bMag
buPhi   = bPhi / bMag
buZ = bZ / bMag

;   Find points inside the boundary

mask   = intArr ( size ( psizr, /dim ) )

for i = 0, nW - 1 do begin
    for j = 0, nH - 1 do begin 

        q1_  = n_elements ( where ( ( R[i] - rbbbs gt 0 ) and ( z[j] - zbbbs gt 0 ), q1 ) )
        q2_  = n_elements ( where ( ( R[i] - rbbbs gt 0 ) and ( z[j] - zbbbs le 0 ), q2 ) )
        q3_  = n_elements ( where ( ( R[i] - rbbbs le 0 ) and ( z[j] - zbbbs gt 0 ), q3 ) )
        q4_  = n_elements ( where ( ( R[i] - rbbbs le 0 ) and ( z[j] - zbbbs le 0 ), q4 ) )

        if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then $
            mask[i,j]  = 1
               
    endfor
endfor
 
iiInside    = where ( mask gt 0 )
iiOutside   = where ( mask eq 0 )

if keyword_set ( pAngle ) then begin
    
    ;   Calculate a poloidal angle coordinate which labels
    ;   the position along a a flux surface. This will be on 
    ;   the rz grid so we can create grad Chi.
    
    ;   Trace the poloidal field lines at each r (minor) coord.
 
    bInterpS    = { bR : bR, $
                    rleft : rleft, $
                    rdim : rdim, $
                    nW : nW, $
                    z : z, $
                    zdim : zdim, $
                    nH : nH, $
                    bPhi : bPhi, $
                    bz : bz }   

    rMinor   = R - rmaxis
    rMinorRight = max ( rbbbs ) - rmaxis
    iiPositiveRMinor    = where ( rMinor gt 0 and rMinor le rMinorRight, iirMinorCnt )
    
    dl  = 0.001
   
    window, 0, xSize = 800, ySize = 800
    device, decomposed = 0
    loadct, 39, /silent
    !p.background = 255
    plot, [0,0], [0,0], /noData, $
        xRange = [0.0,3.0], yRange = [-1.4,1.3], xStyle = 9, yStyle = 9, color = 0
   
    rChi_all    = fltArr ( 12, iirMinorCnt )
    zChi_all    = fltArr ( 12, iirMinorCnt )

    rPos_all    = 0.0
    zPos_all    = 0.0
    chi_all     = 0.0
    gradChi_R_all   = 0.0
    gradChi_Z_all   = 0.0
    lengthP = 0.0
    lengthP_fluxGrid_R  = rMinor[iiPositiveRMinor] + rmaxis
    lengthP_fluxGrid_z  = lengthP_fluxGrid_R * 0.0 + zmaxis
    lengthP_fluxGrid    = interpolate ( psizr, ( lengthP_fluxGrid_R - rleft ) / rdim * (nW-1.0), $
        ( lengthP_fluxGrid_z - min ( z ) ) / zdim * (nH-1.0), cubic = -0.5 )


    for i = 0, iirMinorCnt - 1 do begin
    
        rStart   = rMinor[iiPositiveRMinor[i]] + rmaxis 
        zStart  = zmaxis
       
        ;   Try RK4 

        rPos    = rStart
        zPos    = zStart
    
        rArray  = rStart
        zArray  = zStart

        lArray  = 0.0
 
        stepCnt = 0
        thetaOld    = 2.0 * !Pi
        dPhi    = -2 * !pi / 100.0
        keepRunning = 1
        while keepRunning do begin

            bHere   = interpB ( bInterpS, rPos, zPos )

            K1_R  = dPhi * bHere.bR / bHere.bMag
            K1_z    = dPhi * bHere.bz / bHere.bMag
 
            bHere   = interpB ( bInterpS, rPos + K1_R / 2.0, zPos + K1_z / 2.0 )

            K2_R    = dPhi * bHere.bR / bHere.bMag
            K2_z    = dPhi * bHere.bz / bHere.bMag 
    
            bHere   = interpB ( bInterpS, rPos + K2_R / 2.0, zPos + K2_z / 2.0 )

            K3_R    = dPhi * bHere.bR / bHere.bMag
            K3_z    = dPhi * bHere.bz / bHere.bMag

            bHere   = interpB ( bInterpS, rPos + K3_R, zPos + K3_z )

            K4_R    = dPhi * bHere.bR / bHere.bMag
            K4_z    = dPhi * bHere.bz / bHere.bMag

            rPos    = rPos + ( K1_R + 2.0 * K2_R + 2.0 * K3_R + K4_R ) / 6.0
            zPos    = zPos + ( K1_z + 2.0 * K2_z + 2.0 * K3_z + K4_z ) / 6.0
 
            if stepCnt gt 0 then thetaOld    = theta
            theta   = aTan ( zPos - zmaxis, rPos - rmaxis )
            if theta lt 0 then theta = theta + 2.0 * !pi 
            
            rArray  = [ rArray, rPos ]
            zArray  = [ zArray, zPos ]
            
            dlO =  sqrt ( ( rStart - rArray[stepCnt] )^2 + ( zStart - zArray[stepCnt] )^2 )
            dl  = sqrt ( ( rPos - rArray[stepCnt] )^2 + ( zPos - zArray[stepCnt] )^2 )
            lArray  = [ lArray, lArray[stepCnt] + dl ]
 
            ++ stepCnt 
            if stepCnt gt 3 then oPlot, rArray, zArray, color = 0
            if ( theta - thetaOld ) gt 0 and stepCnt gt 10 then keepRunning = 0
        
		endWhile

        chiArray    = lArray / max ( lArray ) * ( theta + 2.0 * !pi )
        bHere   = interpB ( bInterpS, rArray, zArray )
        
        lengthP = [ lengthP, max ( lArray ) * 2.0 * !pi / ( theta + 2.0 * !pi ) ]

        dChi_dl = ( theta + 2.0 * !pi ) / max ( lArray )      
        gradChi_R = dChi_dl * bHere.bR / sqrt ( bHere.bR^2 + bHere.bz^2 )
        gradChi_Z   = dChi_dl * bHere.bz / sqrt ( bHere.bR^2 + bHere.bz^2 )

        rPos_all    = [ rPos_all, rArray ]
        zPos_all    = [ zPos_all, zArray ]
        chi_all = [ chi_all, chiArray ]
        gradChi_R_all = [ gradChi_R_all, gradChi_R ]
        gradChi_Z_all = [ gradChi_Z_all, gradChi_Z ]


        chi = 0 
        rChi    = rArray[0]
        zChi    = zArray[0]
        for j = 1, n_elements ( chiArray ) - 1 do begin
           if abs ( fix ( chiArray[j] * !radeg / 30.0 ) $
					- fix ( chiArray[j-1] * !radeg / 30.0 ) ) gt 0 then begin 

                chi = [ chi, chiArray[j] ]
                rChi    = [ rChi, rArray[j] ]
                zChi    = [ zChi, zArray[j] ]

            endif
        endfor
        
        ;   Interpolate the r,z location of a regular chi spacing

        rChi_   = interpol ( rArray, chiArray, fIndGen ( 12 ) * 30.0 * !dtor )
        zChi_   = interpol ( zArray, chiArray, fIndGen ( 12 ) * 30.0 * !dtor )

        rChi_all[*,i]   = rChi_        
        zChi_all[*,i]   = zChi_

        oPlot, rChi, zChi, psym = 5, color = 14 * 16 - 1
        for j = 1, n_elements ( chi ) - 1 do begin
            plots, [rChi[j],rChi[j]], [zChi[j],zChi[j]], $
				psym = 4, color = 4 * 16 - 1, symSize = abs ( j * 30.0 - chi[j] * !radeg ) * 2.0
        endfor
    endfor

rPos_all    = rPos_all[1:*]
zPos_all    = zPos_all[1:*]
chi_all    = chi_all[1:*]

triangulate, rPos_all, zPos_all, triangles, b

rz_chi  = triGrid ( rPos_all, zPos_all, chi_all, triangles, $
    [ rStep, zStep ], [ min ( R ), min ( z ), max ( R ), max ( z ) ] ); , exptrapolate = b )

gradChi_R_all   = gradChi_R_all[1:*]
gradChi_Z_all   = gradChi_Z_all[1:*]
lengthP = lengthP[1:*]

lengthP_spline = spl_init ( fluxGrid, fPol )
lengthPolRZ  = reform ( spl_interp ( lengthP_fluxGrid, lengthP, lengthP_spline, psizr[*] ), nW, nH )

for i = 0, 11 do begin

    oPlot, rChi_all[i,*], zChi_all[i,*], color = 0

endfor
!p.background = 0

endif
if keyword_set ( plot_ ) then begin

    device, decomposed = 0
    !p.multi = [ 0, 2, 2 ]
    
    contour, psizr, R, z
    oplot, rbbbs, zbbbs
    oplot, rlim, zlim
    
    contour, bPhi * mask, R, z, levels = fIndGen ( 21 ) / 5 - 2
    
    veloVect, bR * mask , bz * mask , R, z
    
    !p.multi = 0

endif

;   Create data structure

if keyword_set ( pAngle ) then begin

	eqdsk   = { case_ : case_, $
	            nW : nW, $
	            nH : nH, $
	            rDim : rdim, $
	            zDim : zdim, $
	            rLeft : rleft, $
	            zMid : zmid, $
	            rMAxis : rmaxis, $
	            zMAxis : zmaxis, $
	            siMag : siMag, $
	            siBry : sibry, $
	            rCentr : rcentr, $
	            bCentr : bcentr, $
	            current : current, $
	            fPol : fpol, $
	            pres : pres, $
	            ffPrim : ffprim, $
	            pPrime : pprime, $
	            psizr : psizr, $
	            qPsi : qpsi, $
	            nbbbs : nbbbs, $
	            limitr : limitr, $
	            rbbbs : rbbbs, $
	            zbbbs : zbbbs, $
	            rLim : rlim, $
	            zLim : zlim, $
	            bR : bR, $
	            bz : bz, $
	            bPhi : bPhi, $
	            bMag : bMag, $
	            mask : mask, $
	            iiInside : iiInside, $
	            iiOutside : iiOutside, $
	            fPolRZ : fPolRZ, $
	            chiRZ : rz_chi, $
	            lengthPolRZ : lengthPolRZ, $
	            rStep : rStep, $
	            zStep : zStep, $
	            fStep : fStep, $
	            R : R, $
	            z : z, $
	            fluxGrid : fluxGrid, $
	            lengthP : lengthP, $
	            APhi : APhi, $
	            buR : buR, $
	            buPhi : buPhi, $
	            buZ : buZ }

endif else begin

	eqdsk   = { case_ : case_, $
	            nW : nW, $
	            nH : nH, $
	            rDim : rdim, $
	            zDim : zdim, $
	            rLeft : rleft, $
	            zMid : zmid, $
	            rMAxis : rmaxis, $
	            zMAxis : zmaxis, $
	            siMag : siMag, $
	            siBry : sibry, $
	            rCentr : rcentr, $
	            bCentr : bcentr, $
	            current : current, $
	            fPol : fpol, $
	            pres : pres, $
	            ffPrim : ffprim, $
	            pPrime : pprime, $
	            psizr : psizr, $
	            qPsi : qpsi, $
	            nbbbs : nbbbs, $
	            limitr : limitr, $
	            rbbbs : rbbbs, $
	            zbbbs : zbbbs, $
	            rLim : rlim, $
	            zLim : zlim, $
	            bR : bR, $
	            bz : bz, $
	            bPhi : bPhi, $
	            bMag : bMag, $
	            mask : mask, $
	            iiInside : iiInside, $
	            iiOutside : iiOutside, $
	            fPolRZ : fPolRZ, $
	            rStep : rStep, $
	            zStep : zStep, $
	            fStep : fStep, $
	            R : R, $
	            z : z, $
	            fluxGrid : fluxGrid, $
	            APhi : APhi, $
	            buR : buR, $
	            buPhi : buPhi, $
	            buZ : buZ }

endelse

	
return, eqdsk

end
