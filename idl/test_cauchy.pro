FUNCTION f, x  
   RETURN, 1/(5 * x^3 + 6)  
END

pro test_cauchy

ans = imsl_intFcn ( 'f', -1, 5, 0, /cauchy )
PM, 'Computed Answer:', ans
exact = ALOG(125/631.)/18  
PM, 'Exact - Computed:', exact - ans
stop
end
