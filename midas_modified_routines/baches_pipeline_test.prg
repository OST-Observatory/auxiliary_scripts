! $Id
!
! AUTHOR: C. Guirao
!
! NAME: baches_pipeline.prg
!   Standard reduction MIDAS procedure for BACHES echelle spectra
!
! INPUT: 
!   CCD image in FITS format with science spectrum (e.g. sirius.fits)
!
! OUTPUT: 
!   Bi-dimensional image (e.g. sirius_wm.fit) alredy reduced. You can plot
!   each row with the command "PLOT/ECH sirius_wm.fit <row_number>"
!  
!   One-dimensional image (e.g. sirius_wmr.fit) alredy reduced
!   
! CONVENTION:
!   b with bias subtraction
!   d with dark-current subtraction
!   p with flat-field division
!   w calibrated in wavelength
!   r rebined
!   m orders merged
!   f calibrated in flux with a standard star
!
! USAGE: PIPELINE/BACHES frame.fit
!
!*******************************************************************************
! ECHO/FULL
!
! INPUT PARAMS:
! p1 REQUIRED: FITS image of spectrum to be reduced
!
!
! Edit here to change defaults:
!
 DEFINE/PARA p1 ?	  IMA "Enter FITS image to be calibrated: "
 DEFINE/PARA p2 1000      NUM
 DEFINE/PARA p3 3,3       NUM
 DEFINE/PARA p4 7,2,3     NUM
!
! Convert paramaters into variables:
!
 DEFINE/LOCAL TIME/D/1/1        1.
 DEFINE/LOCA cmd/c/1/20		"pipeline"
 DEFINE/LOCA fitsframe/c/1/40   {p1}
 DEFINE/LOCA original/c/1/40
 DEFINE/LOCA frame/c/1/40
 DEFINE/LOCA exposure/i/1/1	0

 DEFINE/LOCA bias/c/1/40         "master_bias.bdf"
 DEFINE/LOCA dark_table/c/1/40   "master_dark.tbl"
 DEFINE/LOCA flat/c/1/40         "master_flat.bdf"
 DEFINE/LOCA session/c/1/40      "baches"
 DEFINE/LOCA flux/c/1/40         "master_flux.bdf"
 DEFINE/LOCA printer/c/1/40      "postscript"
 DEFINE/LOCA session_file1/c/1/40
 DEFINE/LOCA session_file2/c/1/40
 DEFINE/LOCA tmpname/c/1/40
!
 DEFINE/LOCA reply/c/1/1        y
 DEFINE/LOCA fcentral/r/1/1	0
 DEFINE/LOCA fcentrum/i/1/1	0
 DEFINE/LOCA stat_flux/r/1/1    0.0
 DEFINE/LOCA basename/c/1/40    "" ? +lower_levels
 DEFINE/LOCA ext/c/1/10         "" ? +lower_levels
 DEFINE/LOCA next/i/1/1		1
!
! MAIN, main:
!
! Ask before any remove
!
 WRITE/OUT "{cmd}: Version: baches_pipeline"
 WRITE/OUT "{cmd}: Bias={bias}"
 WRITE/OUT "{cmd}: Dark table={dark_table}"
 WRITE/OUT "{cmd}: Flat={flat}"
 WRITE/OUT "{cmd}: Echelle session={session}"
 WRITE/OUT "{cmd}: Flux={flux}"
 WRITE/OUT 
!
! Set context ECHELLE, if not yet already
! and initia session
!
 SHOW/CONT echelle > NULL
 IF  OUTPUTI(5) .eq. 0 THEN
   SET/CONTEXT echelle >NULL
 ENDIF
!
! Gets basename and extension of input file
!
 FILE/BACHES {fitsframe}
! ECHO/FULL
 original = basename
 frame = original
!
! Convert FITS format into MIDAS BDF format
!
 IF ext .ne. ".bdf" THEN
   COMPUTE/BACHES {fitsframe} {frame}
 ENDIF
 frame = frame + "_"
 $rm -f {frame}.bdf
 $ln {original}.bdf {frame}.bdf 
!
! Get the exposure time
!
exposure = m$value({frame},O_TIME(7)) 
WRITE/OUT "{cmd}: exposure time for {frame}: {exposure}"
!
! Bias subtraction
!
 IF m$exist(bias) .eq. 1 THEN 
   COMPUTE/IMA {frame}b = {frame} - {bias}
   frame = frame + "b"
   WRITE/OUT "{cmd}: creating {frame} with {bias} subtracted"
! ELSE
!   WRITE/OUT "{cmd}: {bias} not found. No bias subtracted"
 ENDIF
!
! Dark subtraction
!
 IF m$exist(dark_table) .eq. 1 THEN 
   DARK/BACHES {frame} {frame}d
   frame = frame + "d"
   WRITE/OUT "{cmd}: creating {frame} with dark_current subtracted"
! ELSE
!   WRITE/OUT "{cmd}: {dark_table} not found. No dark_current subtracted"
 ENDIF
!
! Flat-fielding (pixel to pixel variations)
!
 IF m$exist(flat) .eq. 1 THEN
   WRITE/KEYW reply/c/1/1 y
   INQUIRE/KEYW reply "{cmd}: Do you want to correct flux with the flatfield [yn] (y)?"
   IF reply(:1) .eq. "y" THEN
     COMPUTE/IMA {frame}f = {frame} / {flat} >NULL
     frame = frame + "f"
     WRITE/OUT "{cmd}: creating {frame} with {flat} divided"
   ENDIF
 ENDIF
!
! Wavelength calibration 
!
 session_file1 = session + "ORDE.tbl"
 session_file2 = session + "LINE.tbl"

 IF m$exist(session_file1) .eq. 1 THEN
   IF m$exist(session_file2) .eq. 1 THEN
     ! ECHO/FULL
     INITIA/ECHE {session}
     SET/ECHELLE EXTMTD=OPTIMAL
!     SET/ECHELLE MRGMTD=OPTIMAL
!     SET/ECHELLE MRGMTD=NOAPPEND
!     SET/ECHELLE MRGORD=1,26
     REDUCE/ECHE {frame} {frame}wrm
     $cp middummr.bdf  {frame}wr.bdf
     frame = frame + "wr"
     OUTDISK/FITS {frame} {frame}.fit >NULL
     WRITE/OUT "{cmd}: {frame} with all orders separated"
     WRITE/OUT "{cmd}: converting {frame}.bdf to {frame}.fit"
     frame = frame + "m"
    
     ! Flat-field rough flattening
     IF RESPOPT(1:1) .EQ. "Y" .AND. RESPMTD(1:1) .EQ. "F" THEN
       WRITE/OUT "{cmd}: {frame} blazed corrected"
       COMPUTE/IMA  {frame} = (50 * {frame}) / {BLAZE}m
     ENDIF

     OUTDISK/FITS {frame} {frame}.fit >NULL
     WRITE/OUT "{cmd}: {frame} with all orders merged"
   ELSE
     WRITE/OUT "{cmd}: ERROR {session} incomplete:"
     WRITE/OUT "{cmd}: {session_file2} does not exist."
     RETURN/EXIT
   ENDIF
 ELSE
   WRITE/OUT "{cmd}: ERROR {session} incomplete:"
   WRITE/OUT "{cmd}: {session_file1} does not exist."
   RETURN/EXIT
 ENDIF

!
! Converting result into FITS
!
 WRITE/OUT "{cmd}: converting {frame}.bdf to {frame}.fit"
 OUTDISK/FITS {frame} {frame}.fit >NULL
!
 CREATE/GRAPH
 SET/GRAPH default=
 SET/GRAPH XAXIS=4050,7300
! ECHO/FULL
 PLOT {frame}
! $sleep 5
! PLOT/AXES 6450,6650 0,20000
! OVERPLOT {frame}
 ECHO/OFF
