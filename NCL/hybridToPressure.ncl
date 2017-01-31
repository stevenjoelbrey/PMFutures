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
loadscript("substring.ncl")

; TODO: Only load parts of nc files that are needed. Most important
;       with T, because it is huge. 
; TODO: Make generic enough to pass any variable or passed file. 
; TODO: Make generic to any passed variable daily or hourly.
; TODO: Add units to written netcdf information. 

begin ; begin main script 


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Define variables and paths that are specific to this script. ;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	; Set static directory path variables 
	; TODO: PASS THESE AS ARGUMENTS
	variable   = "T" 

	; make grid on these
	lev_p      = (/ 200,250,300,400,500,600,700,850,925,1000 /)  
	plot_lev   = 4 ; level that will be plotted on sanity plot 
	onLev      = True  ; Set to False if variable vertical dim is on ilev, 
		                ; True if on lev. 

	; scenario options
	; 2000Base, 2050RCP45, 2050RCP85, 2100RCP45, 2100RCP85
	scenario = "2050RCP45" 
	nl       = strlen(scenario)          
	year     = substring(scenario, 0, 3)
	scenRCP  = substring(scenario, 4, nl)

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Set string variables required for loading files              ; 
	; May the coding gods have mercy on me for this repeating mess ;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	dataDir    = "/fischer-scratch/sbrey/outputFromYellowstone/"
	if(scenRCP.eq."Base") then
	
		rcpScenario = ""; rcp45 | rcp85 | "" (blank) for base scenario
		scenario    = scenario + "/" ; 2000Base, 2050RCP45, 2050RCP85, 2100RCP45, 2100RCP85
		midYear     = year + rcpScenario ; 2000, 2050_rcp85, 2100_rcp45

	end if

	if(scenRCP.eq."RCP45") then

		rcpScenario = "rcp45"
		scenario    = scenario + "/"
		midYear     = year + "_" + rcpScenario 

	end if 
	
	if(scenRCP.eq."RCP85") then
	
		rcpScenario = "rcp485"
		scenario    = scenario + "/" 
		midYear     = year + "_" + rcpScenario 

	end if

	; yearString depnds on only year
	if(year.eq."2000") then
	
		yearString = "200001-201012" 

	end if

	if(year.eq."2050") then

		yearString = "204001-205012"

	end if

	if(year.eq."2100") then

		yearString = "209001-209912"

	end if


	; In the olden days you had to use your brain to make sure these all worked
	; together when setting manually. 
	;rcpScenario= ""; rcp45 | rcp85 | "" (blank) for base scenario
	;scenario   = "2000Base/" ; 2000Base, 2050RCP45, 2050RCP85, 2100RCP45, 2100RCP85
	;midYear    = "2000"+rcpScenario ; 2000, 2050_rcp85, 2100_rcp45
	;yearString = "200001-201012" ; 200001-201012, 204001-205012, 209001-209912
	;year       = "2000"
	;dataDir    = "/fischer-scratch/sbrey/outputFromYellowstone/"
	print("------------------------------------------------------------")
	print("User written variables are set.")
	print("------------------------------------------------------------")
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

	print("------------------------------------------------------------")
	print("Air Quality model data loaded.")
	print("------------------------------------------------------------")

	; Now get the surface pressure. This is always required for vinth2p_ecmwf()  
	PSFileBase = dataDir + "AirQualityData/" + scenario
	PSFile = PSFileBase + "cesm122_fmozsoa_f09f09_" + midYear + "_fires_00.PS.daily."+yearString+".nc"
	print("Working on loading sea level pressure file")
	fPS = addfile(PSFile,"r")
	ps = fPS->PS
	delete(fPS)

	print("------------------------------------------------------------")
	print("Surface pressure loaded into the workspace.")
	print("------------------------------------------------------------")


	print("------------------------------------------------------------")
	print("Working on loading T into the workspace. Really big file....")
	print("------------------------------------------------------------")
	; Surface temperature is always needed as well
	; e.g.cesm122_fmozsoa_f09f09_2050_rcp45_fires_00.T.daily.204001-205012.nc
	variableFileBase = dataDir + "AirQualityData/" + scenario + "cesm122_fmozsoa_f09f09_"+midYear+"_fires_00."
	variableFile = variableFileBase + "T" + ".daily." + yearString + ".nc"
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

	print("------------------------------------------------------------")
	print("Working on loading " + variable +  " into the workspace.....")
	print("------------------------------------------------------------")

	; Now get the variable to be interpolated to a new grid. 
	variableFileBase = dataDir + "AirQualityData/" + scenario
	variableFileMid  = variableFileBase + "cesm122_fmozsoa_f09f09_" + midYear
	variableFile = variableFileMid + "_fires_00."+variable+".daily."+yearString+".nc"

	fccm = addfile(variableFile, "r")
	print("Working on loading " + variable)
	print("file="+variableFile)

	; take the descriptive variables needed for saving
	VAR    = fccm->$variable$ 
	time   = fccm->time
	levs   = fccm->$levStr$ ; ilev | lev, determined above in if statement
	lat    = fccm->lat
	lon    = fccm->lon

	delete(fccm)


	; last arg is what dimensions match
	if(variable.eq."HEIGHT") then		

		print("------------------------------------------------------------")
		print("Wierd units of HEIGHT require regridding.")
		print("------------------------------------------------------------")

		; Ground level geopotential is needed for regridding function. 
		print("Working on loading ground level groupotential PHIS")
		usgsFile= dataDir + "USGS-gtopo30_0.9x1.25_remap_c051027.nc"
		usgs = addfile(usgsFile,"r")
		phis = (usgs->PHIS) 
		delete(usgs)

		; calculate ground level in meters from PHIS. We convert this to m
		; because HEIGHT (geopotential) is saved in units of m not m^2/s^s.
		; When HEIGHT is passed we want to get from units of m above ground
		; level to units of sea level, so we add topo height in meters to
		; every level. That is only for when the var to be transformed is
		; HEIGHT. 
		topo = doubletofloat(phis/9.81) ; [m]

		; Add topo to HEIGHT in order to make true geopotential
		; Figure out what this means for other variables to be tansformed

		print("Working on adding ground level geopotential to each level of height")
		; https://www.ncl.ucar.edu/Document/Functions/Built-in/conform.shtml
		geopotential = VAR + conform(VAR, topo, (/2,3/))  
		copy_VarCoords(VAR, geopotential)
		VAR_ready = geopotential
		VAR_ready@units = "m"
		VAR_ready@long_name = "geopotential above sea level in meters"
		delete(geopotential)
		varflg = -1

	else 
		print("No addition of topo need for: " + variable)
		VAR_ready = VAR
		varflg = 0 ; for interpreting all non geopotential and temperature layers. 
	end if 


	;===================================================================
	; Transform from Hybrid vertical levels to pressure levels	
	;https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml
	print("------------------------------------------------------------")
	print("Regridding: " + variable + " to pressure coords.............")
	print("------------------------------------------------------------")
	VAR_P = vinth2p_ecmwf(VAR_ready, a_coef, b_coef, lev_p, ps, 1, P0mb, 1, True, varflg, T_surface, phis)

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
	; cesm122_fmozsoa_f09f09_2050_rcp45_fires_00.HEIGHT.daily.204001-205012.nc
	fon =  "../NCL_output/" + scenario + "cesm122_fmozsoa_f09f09_" + midYear + "_fires_00." + variable + "_P.daily." + yearString + ".nc"

	if (fileexists(fon)) then
		system("rm "+ fon) ; remove it because if it exists this wont write
	end if
	fo = addfile(fon, "c")

	;===================================================================
	; create global attributes of the file
	;===================================================================
	fAtt               = True            ; assign file attributes
	fAtt@title         = variable + " regridded to pressure coords"  
	fAtt@source_file   = variableFile
	fAtt@Conventions   = "hybrid coords to pressure coords using vinth2p_ecmwf"   
	fAtt@creation_date = systemfunc("date")        
	fileattdef( fo, fAtt )            ; copy file attributes    


	print("------------------------------------------------------------")
	print("Wring " + variable + " on pressure coords to nc file........")
	print("------------------------------------------------------------")

	; label the pressure level dimension
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

