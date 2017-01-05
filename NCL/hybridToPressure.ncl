;-------------------------------------------------------------------------
; Function being put to use to transform from hybrid coords to pressure 
; coords
;-------------------------------------------------------------------------
; TODO: Close all the files
; TODO: Only load parts of nc files that are needed. Most important
; TODO:    with T, because it is huge. 
; TODO: Make generic for any variable or passed file. 
; TODO: Make generic to any passed variable daily or hourly

; Set static directory path variables 
variable = "HEIGHT"
scenario = "2000Base/" 
yearString="200001-201012" ; must match the files under scenario dir 
dataDir  = "/fischer-scratch/sbrey/outputFromYellowstone/"

; get hyam and hybm | hyai and hybi
; NOTE: These variables are static. They do not change as a function of time
; NOTE: or scenario
coordData = dataDir+"cesm122_fmozsoa_f09f09_2000_fires_00.cam.h0.2000-01.nc"
meta = addfile(coordData, "r")

; TODO: If statement on which coordinates to use. Depends on variable.
hyam = meta->hyam ; pertains to lev
hybm = meta->hybm ; pertains to lev

; for variables like HEIGHT that are defined along 'ilev' use?: 
hyai = meta->hyai ; pertains to ilev
hybi = meta->hybi ; petains to ilev 

;reference sea level specified in model, agurment requires [mb]
P0mb = 0.01*meta->P0

; close fcesm 
delete(meta)

; Now get the surface pressure 
PSFileBase = dataDir + "AirQualityData/" + scenario
PSFile = PSFileBase + "cesm122_fmozsoa_f09f09_2000_fires_00.PS.daily."+yearString+".nc"
fPS = addfile(PSFile,"r")
ps = fPS->PS
delete(fPS)

; Surface temperature is always needed as well
variableFileBase = dataDir + "AirQualityData/" + scenario + "cesm122_fmozsoa_f09f09_2000_fires_00."
variableFile = variableFileBase + "T" + ".daily." + yearString + ".nc"
fccT = addfile(variableFile,"r")
T = fccT->T
T_dims = getvardims(T)
srf_index = T_dims(1) ; lev is the second dimension. 

; T is on lev not ilev (HEIGHT), 25 is the closest level to ground. 
; we want all times and points 
T_surface = T(:,srf_index,:,:)  

delete(fcct)
delete(T)

; Now get the variable to be interpolated to a new grid. 
variableFileBase = dataDir + "AirQualityData/" + scenario + "cesm122_fmozsoa_f09f09_2000_fires_00."
variableFile = variableFileBase + variable + ".daily." + yearString + ".nc"
fccm = addfile(variableFile,"r")

; take of the descriptive variables needed for saving
HEIGHT = fccm->HEIGHT
time = fccm->time
levs  = fccm->ilev ; ilev | lev
lat  = fccm->lat
lon  = fccm->lon

delete(fccm)

; Ground level geopotential needs to be included? I do not know how yet.
; This is something
; I want to try. When I did not add ground level geopotential Mountains
; were too low. AFTER regridding I see ~5000m Z500 heights over the ocean
; but values like 400m over the hemaliya 
usgsFile= dataDir + "USGS-gtopo30_0.9x1.25_remap_c051027.nc"
usgs = addfile(usgsFile,"r")
phis = (usgs->PHIS) ;/ 9.81
delete(usgs)

;HEIGHT_plus = new((/4016,27,192,288/), float, 0) ; also slow because height is huge 

; Add this to all levels all days 
;nDays = dimsizes(time)
;nLev  = dimsizes(levs)
;d = 0
;do while(d.le.nDays)
;	l = 0
;	do while(l.le.nLev)
;		HEIGHT_plus(d,l,:,:) = HEIGHT(d,l,:,:) + groundGeopotential
;		l = l + 1 
;	end do
;	d = d + 1
;end do


;new pressure levels we want output on
lev_p = (/ 200,250,300,400,500,600,700,850,925,1000 /)
;lev_p = (/ 200,500 /)

; http://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml
;HEIGHT_P = vinth2p(HEIGHT, hyai, hybi, lev_p, ps, 1, P0mb, 1, False)

;https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml
HEIGHT_P = vinth2p_ecmwf(HEIGHT, hyai, hybi, lev_p, ps, 1, P0mb, 1, True, -1, T_surface, phis)

;===================================================================
; Sanity check plot
; https://www.ncl.ucar.edu/Applications/Scripts/conwomap_2.ncl 
;===================================================================
wks  = gsn_open_wks("png","X")   ; send graphics to PNG file

res                     = True          ; plot mods desired

res@cnFillOn            = True          ; turn on color
res@cnFillPalette       = "gui_default" ; set color map
 
res@vpWidthF            = 0.5           ; change aspect ratio of plot
res@vpHeightF           = 0.3

res@gsnMaximize         = True          ; maximize plot size

topo = phis/9.81

;x = HEIGHT_P(0,4,:,:) + topo; plot the 1000th day
g = 9.81
topo = doubletofloat(phis/g)
x = HEIGHT_P(1000,9,:,:) ;+ topo;
plot = gsn_csm_contour(wks,x,res)      


; comes with dims
;[time | 4016] x [lev_p | 2] x [lat | 192] x [lon | 288]

;===================================================================
; https://www.ncl.ucar.edu/Applications/method_2.shtml
; save output. plot sea level pressure and 500mb pressure to make sure this 
; makes sense. 
;===================================================================
dims  = dimsizes(HEIGHT_P)                 ; get dimension sizes  
ntime = dims(0)
nlev = 10                                            
nlat = dims(2) 
nlon = dims(3)     

diro = "/fischer-scratch/sbrey/" ; Output directory
filo = "test.nc"                 ; Output file, will be same as old name 
                                 ;only with _P added to indicate that this 
								 ;data is on pressure coords

fout  = addfile (diro + filo, "c")  ; open output file

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "NCL Efficient Approach to netCDF Creation"  
fAtt@source_file   = "original-file.nc"
fAtt@Conventions   = "None"   
fAtt@creation_date = systemfunc ("date")        
fileattdef( fout, fAtt )            ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================
dimNames = (/"time", "lat", "lon", "lev_p"/)  
dimSizes = (/ ntime, nlat, nlon, nlev /) 
dimUnlim = (/ True , False, False, False/)   
filedimdef(fout, dimNames, dimSizes, dimUnlim)

;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining 
; each variable's dimension name(s) and type. 
; The following could be replaced with explicit, user defined dimension 
; names different from those associated with the variable in memory. 
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/)) 
;===================================================================
filevardef(fout, "time" ,typeof(time),getvardims(time)) 
filevardef(fout, "lev_p"  ,typeof(lev_p),"lev_p")                           
filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
filevardef(fout, "HEIGHT_P"    ,typeof(HEIGHT_P)  ,getvardims(HEIGHT_P))    

;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================
filevarattdef(fout,"HEIGHT_P",HEIGHT_P)             ; copy HEIGHT_P attributes
filevarattdef(fout,"time" ,time)                    ; copy time attributes
filevarattdef(fout,"lev_p"  ,lev_p)                 ; copy lev attributes
filevarattdef(fout,"lat"  ,lat)                     ; copy lat attributes
filevarattdef(fout,"lon"  ,lon)                     ; copy lon attributes

;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
setfileoption(fout,"DefineMode", False)  
