;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                             ;;
;; hybridToPressure.ncl                                                        ;;
;;                                                                             ;;
;; AUTHOR                                                                      ;;
;; -> Steven Brey, sjbrey@rams.colostate.edu                                   ;;
;;                                                                             ;;
;; DESCRIPTION                                                                 ;;
;; -> 	Function being put to use to transform from hybrid coords to pressure  ;;
;; 		coords								                                   ;;
;;                                                                             ;;	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; TODO: Only load parts of nc files that are needed. Most important
;       with T, because it is huge. 
; TODO: Make generic enough to pass any variable or passed file. 
; TODO: Make generic to any passed variable daily or hourly.

begin ; begin main script 

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Define variables and paths that are specific to this script. ;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	; Set static directory path variables 
	; TODO: PASS THESE AS ARGUMENTS
	variable   = "HEIGHT" 
	lev_p      = (/ 200,250,300,400,500,600,700,850,925,1000 /) ; make grid on these 
	plot_lev   = 4 ; level that will be plotted on sanity plot 
	onLev      = False  ; Set to False if variable vertical dim is on ilev, 
                        ; True if on lev. for
	scenario   = "2000Base/" ; 2000Base, 2050RCP45, 2050RCP85, 2100RCP45, 2100RCP85
	midYear    = "2000" ; 2000, 2050, 2100
	yearString = "200001-201012" ; 200001-201012, 204001-205012, 209001-209912
	dataDir    = "/fischer-scratch/sbrey/outputFromYellowstone/"

    ;===================================================================
	; Begin loading required model info 
	;===================================================================
    
	; NOTE: These variables are static. They do not change as a function of date
	; NOTE: or scenario
	coordData = dataDir+"cesm122_fmozsoa_f09f09_2000_fires_00.cam.h0.2000-01.nc"
	meta = addfile(coordData, "r")

	; get hyam and hybm | hyai and hybi
	if(onLev) then
		a_coef = meta->hyam ; pertains to lev
		b_coef = meta->hybm ; pertains to lev
		print("Using lev for vertical coordinates on " + variable)
		levStr = "lev"
	else
		; for variables like HEIGHT that are defined along 'ilev'
		a_coef = meta->hyai ; pertains to ilev
		b_coef = meta->hybi ; petains to ilev 
		print("Using ilev for vertical coordinate on " + variable)
		levStr = "ilev"
	end if

	;reference sea level specified in model, agurment requires [mb]
	P0mb = 0.01*meta->P0

	; close fcesm 
	delete(meta)

	; Now get the surface pressure. This is always required.  
	PSFileBase = dataDir + "AirQualityData/" + scenario
	PSFile = PSFileBase + "cesm122_fmozsoa_f09f09_"+midYear+"_fires_00.PS.daily."+yearString+".nc"
	print("Working on loading sea level pressure file")
	fPS = addfile(PSFile,"r")
	ps = fPS->PS
	delete(fPS)

	; Surface temperature is always needed as well
	variableFileBase = dataDir + "AirQualityData/" + scenario + "cesm122_fmozsoa_f09f09_"+midYear+"_fires_00."
	variableFile = variableFileBase + "T" + ".daily." + yearString + ".nc"
	print("Working on loading temperature. Ground level temperature needed.")
	fccT = addfile(variableFile,"r")
	T = fccT->T
	T_lev = fccT->lev
	nTLev = dimsizes(T_lev)
	srf_index = nTLev-1 ; lev is the second dimension. 

	; T is on lev not ilev (HEIGHT), 25 is the closest level to ground. 
	; we want all times and points 
	T_surface = T(:,srf_index,:,:)  

	delete(fccT)
	delete(T)

	; Now get the variable to be interpolated to a new grid. 
	variableFileBase = dataDir + "AirQualityData/" + scenario + "cesm122_fmozsoa_f09f09_2000_fires_00."
	variableFile = variableFileBase + variable + ".daily." + yearString + ".nc"
	fccm = addfile(variableFile,"r")
	print("Working on loading " + variable)
	print("file="+variableFile)

	; take the descriptive variables needed for saving
	VAR    = fccm->$variable$ 
	time   = fccm->time
	levs   = fccm->$levStr$ ; ilev | lev, determined above in if statement
	lat    = fccm->lat
	lon    = fccm->lon

	delete(fccm)

	; Ground level geopotential is needed for regridding function. 
	print("Working on loading ground level groupotential PHIS")
	usgsFile= dataDir + "USGS-gtopo30_0.9x1.25_remap_c051027.nc"
	usgs = addfile(usgsFile,"r")
	phis = (usgs->PHIS) 
	delete(usgs)

	; calculate ground level in meters from PHIS
	topo = doubletofloat(phis/9.81) ; [m]
	
	; Add topo to HEIGHT in order to make true geopotential
	; Figure out what this means for other variables to be tansformed

	; last arg is what dimensions match
	if(variable == "HEIGHT") then
		print("Working on adding ground level geopotential to each level of height")
		geopotential = VAR + conform(VAR, topo, (/2,3/))  
		copy_VarCoords(VAR, geopotential)
		VAR_ready = geopotential
		delete(geopotential)
	else 
		print("no addition of topo need for " + variable)
		VAR_ready = VAR
	end if 
	delete(VAR)

	;===================================================================
	; Transform from Hybrid vertical levels to pressure levels	
	;https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml
	print("Regridding "+ variable)
	VAR_P = vinth2p_ecmwf(VAR_ready, a_coef, b_coef, lev_p, ps, 1, P0mb, 1, True, -1, T_surface, phis)

	;===================================================================
	; Sanity check plot for this variable 
	; https://www.ncl.ucar.edu/Applications/Scripts/conwomap_2.ncl 
	;===================================================================
	;print("plotting "+variable)
	;wks  = gsn_open_wks("png",variable+lev_p(plot_lev)+"_Topo")   ; send graphics to PNG file

	;res                     = True          ; plot mods desired
	;res@tiMainString        = Variable+"on"+lev_p(plot_lev)+"pressure surface"
	;res@cnFillOn            = True          ; turn on color
	;res@cnFillPalette       = "gui_default" ; set color map
	;res@vpWidthF            = 0.5           ; change aspect ratio of plot
	;res@vpHeightF           = 0.3
	;res@gsnMaximize         = True          ; maximize plot size
	
	; always plot the 0th day for sanitck check plot. 
	;x = VAR_P(0,plot_lev,:,:) + topo; add topo if needed in analysis phases
	;plot = gsn_csm_contour(wks, VAR_P(0,4,:,:), res)

	;===================================================================
	; http://www.ncl.ucar.edu/Applications/netcdf4.shtml
	; save output 
	;===================================================================

	; get dim sizes
	;ntime = dimsizes(time)
	;nlev = dimsizes(lev_p)                                            
	;nlat = dimsizes(lat)
	;nlon = dimsizes(lon) 
	print("Working on saving output as netcdf")

	setfileoption("nc", "FileStructure", "Advanced")
	; Make the nc file save name the same name as input but append _P
	fon =  scenario + "cesm122_fmozsoa_f09f09_2000_fires_00." + variable + ".daily." + yearString + "_P.nc"
	system("rm"+ fon) ; remove it because if it exists this wont write
	fo = addfile(fon, "c")

	;===================================================================
	; create global attributes of the file
	;===================================================================
	fAtt               = True            ; assign file attributes
	fAtt@title         = variable + " regridded to pressure coords"  
	fAtt@source_file   = variableFile
	fAtt@Conventions   = "hybrid coords to pressure coords using vinth2p_ecmwf"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fo, fAtt )            ; copy file attributes    

	VAR_P!1="plevel"
	VAR_P&plevel = lev_p
	fo->$variable$ = VAR_P
	
	; close the file connection
	delete(fo)
	;===================================================================
	; explicitly declare file definition mode. Improve efficiency.
	;===================================================================
 	;setfileoption(fo,"DefineMode",True)



	;===================================================================
	; predefine the coordinate variables and their dimensionality
	; Note: to get an UNLIMITED record dimension, we set the dimensionality
	; to -1 (or the actual size) and set the dimension name to True.
	;===================================================================
 	;dimNames = (/"time", "lev_p", "lat", "lon"/)
 	;dimSizes = (/ 1 ,   nlev,  nlat,  nlon /)
 	;dimUnlim = (/ True , False, False, False/)
 	;filedimdef(fo, dimNames, dimSizes, dimUnlim)

	; predefine the the dimensionality of the variables to be written out
	;===================================================================
	; Here we are using NCL functions to facilitate defining
	; each variable's dimension name(s) and type.
	;===================================================================
	;filevardef(fo, "time", typeof(time), getvardims(time))
	;filevarattdef(fo,"time", time)                   ; copy time attributes
	;fo->time   = (/time/)

	;lev = levs(0:9)
	;lev(0:9) = (/lev_p/)
	;filevardef(fo, "lev", typeof(lev), getvardims(lev) )
	;filevarattdef(fo,"lev", lev)                     ; copy lev attributes
	;fo->lev    = (/lev/)

	;filevardef(fo, "lat", typeof(lat), getvardims(lat))
	;filevarattdef(fo,"lat", lat)                     ; copy lat attributes
	;fo->lat    = (/lat/)

	;filevardef(fo, "lon", typeof(lon), getvardims(lon))
	;filevarattdef(fo,"lon", lon)                     ; copy lon attributes
	;fo->lon    = (/lon/)

	;filevardef(fo, variable, typeof(VAR_P), getvardims(VAR_P))
	;filevarattdef(fo,variable, VAR_P)                         ; copy T attributes
	;filevarchunkdef(fo, variable, chunkSizes)
	;filevarcompressleveldef(fo, variable, 2)
	;fo->variable     = (/VAR_P/)
	;print(fo)
	;printVarSummary(t)
	;exit

    


end   ; end main script
