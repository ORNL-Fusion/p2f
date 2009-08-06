function dlg_pDeriv, array, dir, dS 

    nX  = n_elements ( array[*,0] )
    nY  = n_elements ( array[0,*] )

    pDeriv  = fltArr ( size ( array, /dim ) )

    if dir eq 2 then begin
        
        for i = 0, nX - 1 do begin
            for j = 0, nY - 1 do begin
                
                if ( j ge 1 ) and ( j le nY-2 ) then $
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        array[i,j+1] - array[i,j-1] ) 

				if ( j ge 2 ) and ( j le nY-3 ) then $
					pDeriv[i,j]	= 1.0 / ( 12.0 * dS ) * $
						( $
						array[i,j-2] - 8.0 * array[i,j-1] $
						+ 8.0 * array[i,j+1] - array[i,j+2] )
        
                if j eq 0 then $
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        - 3.0 * array[i,j] + 4.0 * array[i,j+1] $
                        - array[i,j+2] )

				if j eq 0 then $
					pDeriv[i,j]	= 1.0 / ( 12.0 * dS ) * $
						( $
						-25.0 * array[i,j] + 48.0 * array[i,j+1] $
						-36.0 * array[i,j+2] + 16.0 * array[i,j+3] $
						-3.0 * array[i,j+4] )

                if j eq nY-1 then $ 
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        3.0 * array[i,j] - 4.0 * array[i,j-1] $
                        + array[i,j-2] )

            endfor
        endfor
 
    endif else if dir eq 1 then begin
    
        for i = 0, nX-1 do begin
            for j = 0, nY - 1 do begin
        
                if i ge 1 and i le nX-2 then $
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        array[i+1,j] - array[i-1,j] ) 

				if i ge 2 and i le nX-3 then $
					pDeriv[i,j]	= 1.0 / ( 12.0 * dS ) * $
						( $
						array[i-2,j] - 8.0 * array[i-1,j] $
						+ 8.0 * array[i+1,j] - array[i+2,j] )
        
                if i eq 0 then $
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        - 3.0 * array[i,j] + 4.0 * array[i+1,j] $
                        - array[i+2,j] )

				if i eq 0 then $
					pDeriv[i,j]	= 1.0 / ( 12.0 * dS ) * $
						( $
						-25.0 * array[i,j] + 48.0 * array[i+1,j] $
						-36.0 * array[i+2,j] + 16.0 * array[i+3,j] $
						-3.0 * array[i+4,j] )

                if i eq nX-1 then $
                    pDeriv[i,j] = 1.0 / ( 2.0 * dS ) * $
                        ( $
                        3.0 * array[i,j] - 4.0 * array[i-1,j] $
                        + array[i-2,j] )

            endfor
        endfor
   
    endif
 
    return, pDeriv   

end

