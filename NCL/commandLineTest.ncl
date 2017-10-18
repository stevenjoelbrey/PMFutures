;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                             ;;
;;commandLineTest.ncl                                                          ;;
;;                                                                             ;;
;; AUTHOR                                                                      ;;
;; -> Steven Brey, sjbrey@rams.colostate.edu                                   ;;
;;                                                                             ;;
;; DESCRIPTION                                                                 ;;
;; -> 	This script is used to see if I can run scripts from command line	   ;;
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
	rcpScenario= "rcp45"; can be "" (blank) for base scenario
	scenario   = "2050RCP45/" ; 2000Base, 2050RCP45, 2050RCP85, 2100RCP45, 2100RCP85
	midYear    = "2050_"+rcpScenario ; 2000, 2050_rcp85, 2100_rcp45
	yearString = "204001-205012" ; 200001-201012, 204001-205012, 209001-209912
	year       = "2050"
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

	print("MUUUUUUUUUUUUUUUUURRRRRPPPPPPHHHHH")
	print("WHY DID YOU LEAVE ME MURPH")
	print("Interstellar is extemely avarege.")

end   ; end main script
