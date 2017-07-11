;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                             ;;
;; make_geo_wind.ncl                                                           ;;
;;                                                                             ;;
;; AUTHOR                                                                      ;;
;; -> Steven Brey, sjbrey@rams.colostate.edu                                   ;;
;; -> HECTOR LIVES                                                             ;;
;; DESCRIPTION                                                                 ;;
;; -> 	This script is used to calculated geostrophic winds using geopotential ;;
;;       height fields. The NCL function z2geouv() is used to estimate the u_g ;;
;;       and v_g                                                               ;;
;;									       ;;	
;;                                                                             ;; 
;;  									       ;;
;;    									       ;;	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadscript("substring.ncl")

; TODO: Take the if statements that follow and make a function that can be used
; TODO: in all of my NCL scripts. I do not like ncl but there sure are nice built
; TODO: in functions.  

begin ; begin main script 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define variables and paths that are specific to this script. ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;variable   = "Z3" 
;scenario = "2000Firev1" 
fmod = "_fires_01."

nl       = strlen(scenario)          
year     = substring(scenario, 0, 3)
scenRCP  = substring(scenario, 4, nl)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set string variables required for loading files              ; 
; May the coding gods have mercy on me for this repeating mess ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dataDir    = "/pierce-scratch/mariavm/"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set string variables required for loading files              ; 
; May the coding gods have mercy on me for this repeating mess ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; base scenarios
if(scenRCP.eq."Firev1") then

	rcpScenario = ""
	scenario    = scenario + "/" 
	midYear     = year + rcpScenario 

end if

if(scenRCP.eq."NoFirev1") then

	rcpScenario = ""
	scenario    = scenario + "/" 
	midYear     = year + rcpScenario 

end if

; RCP45 scenarios	
if(scenRCP.eq."RCP45Firev1") then

	rcpScenario = "rcp45"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 

if(scenRCP.eq."RCP45NoFirev1") then

	rcpScenario = "rcp45"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 

if(scenRCP.eq."RCP45FireHIv1") then

	rcpScenario = "rcp45"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 

; RCP88 scenarios
if(scenRCP.eq."RCP85Firev1") then

	rcpScenario = "rcp85"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 

if(scenRCP.eq."RCP85NoFirev1") then

	rcpScenario = "rcp85"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 

if(scenRCP.eq."RCP85FireHIv1") then

	rcpScenario = "rcp85"
	scenario    = scenario + "/"
	midYear     = year + "_" + rcpScenario 

end if 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; yearString depnds on only year
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if(year.eq."2000") then

	yearString = "200001-201012" 

end if

if(year.eq."2050") then

	yearString = "204001-205012"

end if

if(year.eq."2100") then

	yearString = "209001-209912"

end if



print("------------------------------------------------------------")
print("User written variables are set.")
print("------------------------------------------------------------")

print("------------------------------------------------------------")
print("Loading geopotential heights................................")
print("------------------------------------------------------------")


variableFileBase = dataDir + scenario
variableFileMid  = variableFileBase + "cesm122_fmozsoa_f09f09_" + midYear
variableFile = variableFileMid + fmod +"Z3_P"+".daily." + yearString + ".nc"

Z_nc = addfile(variableFile, "r") ; NEED CREATED FILES TO BE MOVED HERE
print("Working on loading Z3_P")
print("file="+variableFile)

; take the descriptive variables needed for saving
Z      = Z_nc->Z3 
time   = Z_nc->time
lev_p  = Z_nc->plevel ; ilev | lev, determined above in if statement
lat    = Z_nc->lat
lon    = Z_nc->lon

delete(Z_nc)

print("------------------------------------------------------------")
print("Using z2geouv() to estimate geostophic winds.")
print("------------------------------------------------------------")
; http://www.ncl.ucar.edu/Document/Functions/Built-in/z2geouv.shtml
; Assume z is a global grid with dimension sizes (ntim,klev,nlat,mlon). 
; There are no missing values and the data are periodic in longitude. Then:

uv = z2geouv(Z, lat, lon, 1) ; 1 indicates Z cyclic in longitude


; create u component variable 
u = uv(0, :, :, :, :)

u!0    = "time"
u!1    = "plevel"
u!2    = "lat"
u!3    = "lon"

u&time   =  time
u&plevel =  lev_p
u&lat    =  lat
u&lon    =  lon

; create v component variable
v = uv(1, :, :, :, :)

v!0    = "time"
v!1    = "plevel"
v!2    = "lat"
v!3    = "lon"

v&time   =  time
v&plevel =  lev_p
v&lat    =  lat
v&lon    =  lon


;===================================================================
; Sanity check plot for this variable 
; https://www.ncl.ucar.edu/Applications/Scripts/conwomap_2.ncl 
;===================================================================
print("plotting geo winds at the surface")
wks  = gsn_open_wks("png", "1000_mb_Ug")   ; send graphics to PNG file

res                     = True          ; plot mods desired
res@tiMainString        = "test plot of u wind field"
res@cnFillOn            = True          ; turn on color
res@cnFillPalette       = "gui_default" ; set color map
res@vpWidthF            = 0.5           ; change aspect ratio of plot
res@vpHeightF           = 0.3
res@gsnMaximize         = True          ; maximize plot size

; always plot the 0th day for sanitck check plot. 
xn = u(1000, 9, 150:180, 150:180) ;
plot = gsn_csm_contour(wks, xn, res)


print("------------------------------------------------------------")
print("Saving u_g and v_g as nc file with nicely labelled attribute")
print("------------------------------------------------------------")
;===================================================================
; http://www.ncl.ucar.edu/Applications/netcdf4.shtml
; save output 
;===================================================================

print("Working on saving u component as netcdf")


; Make the nc file save name the same name as input but append _P
; cesm122_fmozsoa_f09f09_2050_rcp45_fires_00.HEIGHT.daily.204001-205012.nc
fon =  variableFileBase + "cesm122_fmozsoa_f09f09_" + midYear + fmod + "ug_P.daily." + yearString + ".nc"

if (fileexists(fon)) then
	system("rm "+ fon) ; remove it because if it exists this wont write
end if
fo = addfile(fon, "c")

setfileoption("nc", "DefineMode", True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "u component of geostrophic winds created with Z3 regridded to pressure coords."  
fAtt@source_file   = "Z3 regridded to pressure coords."
fAtt@Conventions   = "Geostrophic winds calculated using z2geouv()"   
fAtt@creation_date = systemfunc("date")        
fileattdef( fo, fAtt )            ; copy file attributes   


; predefine dimensionality
dimNames = (/"time", "plevel", "lat", "lon"/)
dimSizes = (/dimsizes(time), dimsizes(lev_p), dimsizes(lat), dimsizes(lon)/)
dimUnlim = (/ False , False, False, False/)
filedimdef(fo, dimNames, dimSizes, dimUnlim)

; predefine the the dimensionality of the variables to be written out
filevardef(fo, "time", typeof(time), getvardims(time))
filevarattdef(fo,"time", time)                   ; copy time attributes
fo->time   = (/time/)

filevardef(fo, "plevel", typeof(lev_p), getvardims(lev_p) )
filevarattdef(fo,"plevel", lev_p)                     ; copy lev attributes
fo->plevel    = (/lev_p/)

filevardef(fo, "lat", typeof(lat), getvardims(lat))
filevarattdef(fo,"lat", lat)                     ; copy lat attributes
fo->lat    = (/lat/)

filevardef(fo, "lon", typeof(lon), getvardims(lon))
filevarattdef(fo,"lon", lon)                     ; copy lon attributes
fo->lon    = (/lon/)

filevardef(fo, "ug", typeof(u), getvardims(u))
filevarattdef(fo,"ug", u)                         ; copy U attributes
fo->ug      = (/u/)

delete(fo)

print("------------------------------------------------------------")
print("Done writing ug")
print("------------------------------------------------------------")


print("------------------------------------------------------------")
print("Working on saving v component as netcdf")
print("------------------------------------------------------------")

; Make the nc file save name the same name as input but append _P
; cesm122_fmozsoa_f09f09_2050_rcp45_fires_00.HEIGHT.daily.204001-205012.nc
fon =  variableFileBase + "cesm122_fmozsoa_f09f09_" + midYear + fmod + "vg_P.daily." + yearString + ".nc"

if (fileexists(fon)) then
	system("rm "+ fon) ; remove it because if it exists this wont write
end if
fo = addfile(fon, "c")

setfileoption("nc", "DefineMode", True)

;===================================================================
; create global attributes of the file
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = "v component of geostrophic winds created with Z3 regridded to pressure coords."  
fAtt@source_file   = "Z3 regridded to pressure coords."
fAtt@Conventions   = "Geostrophic winds calculated using z2geouv()"   
fAtt@creation_date = systemfunc("date")        
fileattdef( fo, fAtt )            ; copy file attributes   


; predefine dimensionality
dimNames = (/"time", "plevel", "lat", "lon"/)
dimSizes = (/dimsizes(time), dimsizes(lev_p), dimsizes(lat), dimsizes(lon)/)
dimUnlim = (/ False , False, False, False/)
filedimdef(fo, dimNames, dimSizes, dimUnlim)

; predefine the the dimensionality of the variables to be written out
filevardef(fo, "time", typeof(time), getvardims(time))
filevarattdef(fo,"time", time)                   ; copy time attributes
fo->time   = (/time/)

filevardef(fo, "plevel", typeof(lev_p), getvardims(lev_p) )
filevarattdef(fo,"plevel", lev_p)                     ; copy lev attributes
fo->plevel    = (/lev_p/)

filevardef(fo, "lat", typeof(lat), getvardims(lat))
filevarattdef(fo,"lat", lat)                     ; copy lat attributes
fo->lat    = (/lat/)

filevardef(fo, "lon", typeof(lon), getvardims(lon))
filevarattdef(fo,"lon", lon)                     ; copy lon attributes
fo->lon    = (/lon/)

filevardef(fo, "vg", typeof(v), getvardims(v))
filevarattdef(fo,"vg", v)                         ; copy V attributes
fo->vg      = (/v/)

delete(fo)

end   ; end main script
