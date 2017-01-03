; get hyam and hybm
coordData = ("/fischer-scratch/sbrey/outputFromYellowstone/cesm122_fmozsoa_f09f09_2000_fires_00.cam.h0.2000-01.nc")

fcesm = addfile(coordData, "r")
hyam = fcesm->hyam ; pertains to lev
hybm = fcesm->hybm ; pertains to lev

; for variables like HEIGHT that are defined along 'ilev' use?: 
; hyai = fcesn->hyai
; hybi = fcesm->hybi

;reference sea level specified in model, agurment requires [mb]
P0mb = 0.01*fcesm->P0

; Now get the surface pressure 
fPS = addfile("/fischer-scratch/sbrey/outputFromYellowstone/AirQualityData/2000Base/cesm122_fmozsoa_f09f09_2000_fires_00.PS.daily.200001-201012.nc","r")

ps = fPS->PS

; Now get the HEIGHT field to be interpolated to a new grid. 
fccm = addfile("/fischer-scratch/sbrey/outputFromYellowstone/AirQualityData/2000Base/cesm122_fmozsoa_f09f09_2000_fires_00.HEIGHT.daily.200001-201012.nc","r")

;new pressure levels we want output on
lev_p = (/ 200,250, 300,400,500,600,700,850,925,1000 /)

HEIGHT_P = vinth2p(fccm->HEIGHT, hyam, hybm, lev_p, ps,1,P0mb,1, True)

; save output. plot sea level pressure and 500mb pressure to make sure this 
; makes sense. 

ntim  = dimsizes(time)                 ; get dimension sizes  
klev  = dimsizes(lev)                                               
nlat  = dimsizes(lat)  
nlon  = dimsizes(lon)      
