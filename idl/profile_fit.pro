pro profile_fit, x, a, f 

	;	a[0] 	= tlim
	;	a[1]	= t0
	;	a[2]	= alpha
	;	a[3]	= beta
	;	x	= rho

	f = a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]

end
