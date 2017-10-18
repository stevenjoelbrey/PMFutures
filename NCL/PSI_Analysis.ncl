;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                             ;;
;; PSI_Analysis.ncl                                                            ;;
;;                                                                             ;;
;; AUTHOR                                                                      ;;
;; -> Cory Baggett, cfb128@psu.edu                                             ;;
;;                                                                             ;;
;; DESCRIPTION                                                                 ;;
;; -> 									       ;;
;;                                                                             ;;	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load NCL function scripts. ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Load constants.ncl and energetics.ncl. ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "/charney/s0/cfb128/NCL_FUNCTIONS/constants.ncl"
load "/charney/s0/cfb128/NCL_FUNCTIONS/energetics.ncl"

begin   ; begin main script

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Define variables and paths that are specific to this script. ;
   ; -> The user must chose which storm to analyze by defining    ;
   ;    storm_num.                                                ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   nlat         = 73      ; Values that the output data is interpolated to
   nlon         = 144   
   ti           = 0       ; IC time step
   tf           = 19      ; Final time step
   levout       = LEV23   ; ouput pressure levels
   level        = 250

   figdate        = "20150218"
   cnMaxLevelValF = 2.                 ; maximum contour value to be plotted
   plot_type      = "pdf"
   nplot          = 5

   mdl_run_name = "bcc_csm1_1_P2mP1"
   clm_name     = "bcc_csm1_1_no"
   obs_name     = "erai.DJFM"

   mdl_path     = "/acc/s1/cfb128/GFDL/test/" + mdl_run_name + "/ensemble0/atmos_daily.nc"
   clm_path     = "/acc/s1/cfb128/GFDL/test/" + clm_name + "/ensemble0/atmos_daily.nc"
   obs_path     = "/acc/s1/cfb128/GFDL/data/" + obs_name + ".nc"
   outpath      = "/acc/s1/cfb128/GFDL/data/output/PSI_ANOM_" + mdl_run_name + ".nc"

   figpath      = "/acc/s1/cfb128/FIGURES/" + figdate + "/"
   figname      = "PSI_Anom_vs_" + mdl_run_name

   mdl_file     = addfile(mdl_path, "r")
   clm_file     = addfile(clm_path, "r")
   obs_file     = addfile(obs_path, "r")

   print("Loading data.")
   ucomp_mdl        = mdl_file->ucomp(ti:tf,:,:,:)
   vcomp_mdl        = mdl_file->vcomp(ti:tf,:,:,:)
   temp_mdl         = mdl_file->temp(ti:tf,:,:,:)
   height_mdl       = mdl_file->height(ti:tf,:,:,:)
   ps_mdl           = mdl_file->ps(ti:tf,:,:)

   ucomp_clm        = clm_file->ucomp(ti:tf,:,:,:)
   vcomp_clm        = clm_file->vcomp(ti:tf,:,:,:)
   temp_clm         = clm_file->temp(ti:tf,:,:,:)
   height_clm       = clm_file->height(ti:tf,:,:,:)
   ps_clm           = clm_file->ps(ti:tf,:,:)

   print("Converting data to a fixed grid and horizontally interpolating.")
   ucompF_mdl  = g2fsh_Wrap(ucomp_mdl,  (/nlat,nlon/))
   vcompF_mdl  = g2fsh_Wrap(vcomp_mdl,  (/nlat,nlon/))
   tempF_mdl   = g2fsh_Wrap(temp_mdl,   (/nlat,nlon/))
   heightF_mdl = g2fsh_Wrap(height_mdl, (/nlat,nlon/))
   psF_mdl     = g2fsh_Wrap(ps_mdl,     (/nlat,nlon/))

   ucompF_clm  = g2fsh_Wrap(ucomp_clm,  (/nlat,nlon/))
   vcompF_clm  = g2fsh_Wrap(vcomp_clm,  (/nlat,nlon/))
   tempF_clm   = g2fsh_Wrap(temp_clm,   (/nlat,nlon/))
   heightF_clm = g2fsh_Wrap(height_clm, (/nlat,nlon/))
   psF_clm     = g2fsh_Wrap(ps_clm,     (/nlat,nlon/))

   hyB         = mdl_file->pfull / 1000.
   nhyB        = dimsizes(hyB)
   hyA         = hyB
   hyA         = 0.
   tempF_bot   = tempF_mdl(:,nhyB-1,:,:)
   ZF_bot      = G * heightF_mdl(:,nhyB-1,:,:) 

   print("Vertically interpolating to ECMWF pressure levels.")
   ucompFp_mdl  = vinth2p_ecmwf(ucompF_mdl,  hyA, hyB, levout, psF_mdl, 1, 1000, 1, True,  0, tempF_bot, ZF_bot)
   vcompFp_mdl  = vinth2p_ecmwf(vcompF_mdl,  hyA, hyB, levout, psF_mdl, 1, 1000, 1, True,  0, tempF_bot, ZF_bot)

   hyB         = clm_file->pfull / 1000.
   nhyB        = dimsizes(hyB)
   hyA         = hyB
   hyA         = 0.
   tempF_bot   = tempF_clm(:,nhyB-1,:,:)
   ZF_bot      = G * heightF_clm(:,nhyB-1,:,:)

   ucompFp_clm  = vinth2p_ecmwf(ucompF_clm,  hyA, hyB, levout, psF_clm, 1, 1000, 1, True,  0, tempF_bot, ZF_bot)
   vcompFp_clm  = vinth2p_ecmwf(vcompF_clm,  hyA, hyB, levout, psF_clm, 1, 1000, 1, True,  0, tempF_bot, ZF_bot)

   levelist  = obs_file->levelist
   latitude  = obs_file->latitude
   longitude = obs_file->longitude

   ucompFp_mdl!1         = "levelist"
   ucompFp_mdl&levelist  = levelist
   ucompFp_mdl!2         = "latitude"
   ucompFp_mdl&latitude  = latitude
   ucompFp_mdl!3         = "longitude"
   ucompFp_mdl&longitude = longitude

   vcompFp_mdl!1         = "levelist"
   vcompFp_mdl&levelist  = levelist
   vcompFp_mdl!2         = "latitude"
   vcompFp_mdl&latitude  = latitude
   vcompFp_mdl!3         = "longitude"
   vcompFp_mdl&longitude = longitude

   ucompFp_clm!1         = "levelist"
   ucompFp_clm&levelist  = levelist
   ucompFp_clm!2         = "latitude"
   ucompFp_clm&latitude  = latitude
   ucompFp_clm!3         = "longitude"
   ucompFp_clm&longitude = longitude

   vcompFp_clm!1         = "levelist"
   vcompFp_clm&levelist  = levelist
   vcompFp_clm!2         = "latitude"
   vcompFp_clm&latitude  = latitude
   vcompFp_clm!3         = "longitude"
   vcompFp_clm&longitude = longitude

   sfvp_mdl = uv2sfvpF(ucompFp_mdl(:,{level},:,:),vcompFp_mdl(:,{level},:,:))
   psi_mdl  = sfvp_mdl(0,:,:,:)

   sfvp_clm = uv2sfvpF(ucompFp_clm(:,{level},:,:),vcompFp_clm(:,{level},:,:)) 
   psi_clm  = sfvp_clm(0,:,:,:)

   psi_mdl_anom = psi_mdl - psi_clm

   psi_mdl_anom!0         = "time"
   psi_mdl_anom&time      = ispan(ti,tf,1)
   psi_mdl_anom!1         = "latitude"
   psi_mdl_anom&latitude  = latitude
   psi_mdl_anom!2         = "longitude"
   psi_mdl_anom&longitude = longitude
   psi_mdl_anom@level     = level + " hPa"
   psi_mdl_anom@units     = "m**2 s**-1"

   printVarSummary(psi_mdl_anom)

   print("Saving " + outpath + ".")
   system ("rm " + outpath)
   outfile               = addfile(outpath, "c")
   outfile->psi_mdl_anom = psi_mdl_anom

   psi_mdl_anom_temp = psi_mdl_anom
   psi_mdl_anom      = psi_mdl_anom / 1000000.

   copy_VarMeta(psi_mdl_anom_temp, psi_mdl_anom)

   tad        = new((/nplot,nlat,nlon/),float)
   tad(0,:,:) = dim_avg_n_Wrap(psi_mdl_anom(0:2,:,:),0)
   tad(1,:,:) = dim_avg_n_Wrap(psi_mdl_anom(3:7,:,:),0)
   tad(2,:,:) = dim_avg_n_Wrap(psi_mdl_anom(8:12,:,:),0)
   tad(3,:,:) = dim_avg_n_Wrap(psi_mdl_anom(13:17,:,:),0)
   tad(4,:,:) = dim_avg_n_Wrap(psi_mdl_anom(18:19,:,:),0)

   wks   = gsn_open_wks (plot_type, figpath + figname)              ; open workstation
   gsn_define_colormap(wks,"cmp_flux")                              ; choose colormap
   plot = new(nplot,graphic)                                        ; create an array to hold the individual panels of the panel plot

   res                      = True                                  ; plot mods desired
   res@gsnDraw              = False                                 ; don't draw
   res@gsnFrame             = False                                 ; don't advance frame
   res@cnLevelSelectionMode = "ManualLevels"                        ; manually select levels
   res@cnLineLabelsOn       = False                                 ; turn on line labels
   res@cnFillOn             = True                                  ; turn on color fill
   res@cnLinesOn            = False                                 ; turn off countour lines
   res@gsnSpreadColors      = True                                  ; use full range of map
   res@gsnLeftString        = ""
   res@gsnRightString       = ""                                    ; define right subtitle
   res@mpGridAndLimbOn      = False                                 ; turn off lat/lon lines
   res@mpCenterLonF         = 180.0                                 ; Center the map on longitude 180E
   res@mpMinLatF            = -30.0                                 ; Minimum Latitude
   res@tmYLMode   = "Explicit"                                      ; explicit Y-axis tick marks
   res@tmYLValues = ((/-90,-60,-30,0,30,60,90/))                    ; selects left y-axis tick marks
   res@tmYRValues = ((/-90,-60,-30,0,30,60,90/))                    ; selects right y-axis tick marks
   res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)     ; labels y-axis tick marks
   res@tmXBMode   = "Explicit"                                      ; explicit x-axis tick marks
   res@tmXBValues = ((/0,60,120,180,240,300,360/))                  ; selects bottom x-axis tick marks
   res@tmXBLabels = (/"0","60E","120E","180E","240E","300E","0"/)   ; labels x-axis tick marks
   res@tmXTValues = ((/0,60,120,180,240,300,260/))                  ; selects top x-axis tick marks
   res@tmXBLabelFontHeightF = 0.02                                 ; bottom x-axis tick mark label font height
   res@tmYLLabelFontHeightF = 0.02                                 ; left y-axis tick mark label font height
   res@tiYAxisString        = "Latitude (degrees)"                  ; y-axis label
   res@tiXAxisString        = "Longitude (degrees)"                 ; x-axis label
   res@tiMainFontHeightF    = 0.02                                 ; title font height
   res@tiYAxisFontHeightF   = 0.02                                 ; y-axis font height
   res@tiXAxisFontHeightF   = 0.02                                 ; x-axis font height
   res@gsnLeftStringFontHeightF = 0.02
   res@tiYAxisOffsetXF      = -0.045                              ; offset y-axis title
   res@tiXAxisOffsetYF      = 0.085                               ; offset x-axis title
   res@tmXBLabelsOn         = False

   res@lbLabelBarOn         = False                                 ; turn off individual contour bars
   res@cnMaxLevelValF       = cnMaxLevelValF
   res@cnMinLevelValF       = -cnMaxLevelValF
   res@cnLevelSpacingF      = 2. * cnMaxLevelValF / 20

   res@tiXAxisOn                = False

   lags = (/"~F22~a~F21~ lags 0 to +2", \
            "~F22~b~F21~ lags +3 to +7",  \
            "~F22~c~F21~ lags +8 to +12",  \
            "~F22~d~F21~ lags +13 to +17",  \
            "~F22~e~F21~ lags +18 to +19"/)

   txres                       = True
   txres@txFontHeightF         = 0.02
   txres@txPerimOn             = True
   txres@txBackgroundFillColor = "White"

   loncoord = (/311,314,313,311,308/)

   do p = 0,nplot-1
     if (p .eq. nplot - 1) then
        res@tmXBLabelsOn             = True
        res@tiXAxisOn                = True
        res@lbLabelBarOn             = True
        res@lbBoxLinesOn             = True
        res@lbTitlePosition          = "Bottom"
        res@lbTitleFontHeightF       = 0.02                         ; contour title font height
        res@lbLabelFontHeightF       = 0.02
        res@lbTitleOffsetF           = 0.05                          ; offsets the label title from the contour labels
        res@lbLabelStride            = 5
        res@lbTitleString            = "10~S~6~N~ m~S~2~N~ s~S~-1~N~"
        res@pmLabelBarOrthogonalPosF = 0.35
     end if
     plot(p)               = gsn_csm_contour_map(wks,tad(p,:,:),res)
     plottxt               = gsn_add_text(wks,plot(p),lags(p),loncoord(p),-18,txres)
   end do

   resP                            = True                          ; panel mods desired
   resP@gsnPanelYWhiteSpacePercent = 2                             ; add 5% white space in the y-direction
   resP@gsnPanelBottom             = 0.12                          ; shrink panel plot by setting bottom edge of plot
   resP@gsnMaximize                = True

   gsn_panel(wks,plot,(/nplot,1/),resP)

end   ; end main script

