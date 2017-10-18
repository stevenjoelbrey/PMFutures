;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                             ;;
;; calculateRelativeVorticity.ncl                                              ;;
;;                                                                             ;;
;; AUTHOR                                                                      ;;
;; -> Steven Brey, sjbrey@rams.colostate.edu                                   ;;
;;                                                                  	       ;;
;; DESCRIPTION                                                                 ;;
;; -> Calculates MERRA2 vorticity from u and v using procedure uv2vrdvf()      ;;
;;                           						       ;;	
;;                                                                             ;;	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


MERRA2Path = "/barnes-scratch/REANALYSES/MERRA2/inst6_3d_ana_Np"
outputDir = "/fischer-scratch/sbrey/MERRA2_vorticity"


