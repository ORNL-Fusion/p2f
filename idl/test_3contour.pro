PRO TEST_3CONTOUR, x, y, z, xy, xz, yz, xc, yc, zc 
;   x,  y,  z  : coordinates of points to plot in 3-D scatterplot 
;   xy, xz, yz : functions to be contoured in each of the planes 
;   xc, yx, zc : coordinates of grids for xy, xz, yz 
IF (N_PARAMS() EQ 0L) THEN BEGIN 
   np = 1000L 
   x = RANDOMU(dseed, np) 
   y = RANDOMU(dseed, np) 
   z = RANDOMU(dseed, np) 
   nx = 20L 
   ny = 25L 
   nz = 30L 
   xc = FINDGEN(nx)/(nx-1L) 
   yc = FINDGEN(ny)/(ny-1L) 
   zc = FINDGEN(nz)/(nz-1L) 
   xy = DIST(nx, ny)   
   xz = DIST(nx, nz)   
   yz = DIST(ny, nz)   
ENDIF   
CONTOUR, xy, xc, yc, PATH_XY = xy_path, /PATH_DATA_COORDS, PATH_INFO = xy_info  ;Compute x-y contours 
CONTOUR, xz, xc, zc, PATH_XY = xz_path, /PATH_DATA_COORDS, PATH_INFO = xz_info  ;Compute x-z contours 
CONTOUR, yz, yc, zc, PATH_XY = yz_path, /PATH_DATA_COORDS, PATH_INFO = yz_info  ;Compute y-z contours 
xmin =  0.0 
xmax =  1.0 
ymin =  0.0 
ymax =  1.0 
zmin =  0.0 
zmax =  1.0 
PLOT_3DBOX, x, y, z, PSYM = 3, TITLE = 'Test Plot', $                           ;Plot 3D scatterplot 
   XTITLE = 'X', XRANGE = [xmin, xmax], $ 
   YTITLE = 'Y', YRANGE = [ymin, ymax], $ 
   ZTITLE = 'Z', ZRANGE = [zmin, zmax] 
FOR i = 0L, N_ELEMENTS(xy_info) - 1L DO BEGIN 
   ii = xy_info(i).offset + [LINDGEN(xy_info(i).n), 0L] 
   PLOTS, xy_path(0L,ii), $                                                     ;Draw an x-y contour 
          xy_path(1L,ii), $ 
          REPLICATE(zmin, xy_info(i).n+1L), $ 
          /T3D 
ENDFOR 
FOR i = 0L, N_ELEMENTS(xz_info) - 1L DO BEGIN 
   ii = xz_info(i).offset + [LINDGEN(xz_info(i).n), 0L] 
   PLOTS, xz_path(0L,ii), $                                                     ;Draw an x-z contour 
          REPLICATE(ymax, xz_info(i).n+1L), $ 
          xz_path(1L,ii), $ 
          /T3D 
ENDFOR 
FOR i = 0L, N_ELEMENTS(yz_info) - 1L DO BEGIN 
   ii = yz_info(i).offset + [LINDGEN(yz_info(i).n), 0L] 
   PLOTS, REPLICATE(xmax, yz_info(i).n+1L), $                                   ;Draw an y-z contour 
          yz_path(0L,ii), $ 
          yz_path(1L,ii), $ 
          /T3D 
ENDFOR 
RETURN 
END
