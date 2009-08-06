pro plot_p2f_tdep

fileList	= file_search ( 'data/p2f_*' )
for ff=0,n_elements(fileList)-1 do begin

	fileName	= fileList[ff] 
	plot_p2f, fileName = fileName, amu = 1

endfor	

end
