! $Id: baches_calib.prg,v 1.5 2012/11/25 18:33:45 midasmgr Exp $
! $Name:  $
!
! AUTHOR: C. Guirao
!
! NAME: baches_calib.prg
!   Semi-automatic wavelength procedure for echelle spectra BACHES.
!   The calibration is based in the ECHELLE package of MIDAS which assumes
!   raw-images oriented so that orders appears along the rows, withe low orders 
!   up and wavelength increases from left to right, i.e. red is up and left.
!   START and STEP descriptors of the image must be positive.
!
! INPUT:
!   - 2-dimensional FITS/BDF image with the spectrum of flat field lamp.  
!   - 2-dimensional FITS/BDF image with the spectrum of a calibration lamp.
!   - MIDAS table of the calibration lamp in FITS format. Default "thar.fit"
!   - Number of orders to be calibrated. Default 24
!   - Slit width for order identification. Default 20
!   - Slit width for extraction. Default 20
!   - Final tolerance on RMS error. Default 0.4
!   - Final degree of polynomial defining dispersion. Default 4
!
! OUTPUT:
!   - ECHELLE session files for BACHES. Result of the "SAVE/ECHE baches"
!     generates: bachesLINE.tbl, bachesORDE.tbl and bachesback.tbl
!     
!
! SUFIX CONVENTION:
!   b with bias subtraction
!   d with dark-current subtraction
!   f with flat-field division
!   1 averaged to one-dimensional image
!   w calibrated in wavelength
!   r rebined in wavelength
!   s calibrated in flux with a standard star
!   m all orders merged 
!
! USAGE: BACHES/CALI <ff> <lamp> [table] [orders] [withi] [slit] [offset] [extract_method] [tol] [poly]
!
! EXAMPLES:
!      BACHES/CALI ff.fit th               (0 orders, 0 width, 0.2 tol, 3 poly)
!      BACHES/CALI ff.fit th ? 25 ? 5 0.5 4  (25 orders, 5 width, 0.5 tol, 5 poly)
!      BACHES/CALI ff.fit th neon ? 5 5 0.5 4 (neon table, 5 orders, 5 width, 0.5 tol, 4 poly)
!      BACHES/CALI ff.fit th neon ? 10 5 0.5 4 (neon table, 5 orders, 10 width for order identification, 5 width for extraction, 0.5 tol, 4 poly)
!
!*******************************************************************************
! ECHO/FULL
!
! INPUT PARAMS:
! p1 REQUIRED: FITS or BDF 2-dimension of image of spectrum of a flat field lamp
! p2 REQUIRED: FITS or BDF 2-dimension of image of spectrum of a calibration lamp
! p3 OPTIONAL: Orders to be calibrated (default 0 -> auto)
! p4 OPTIONAL: Slit width in pixels for order identification (default 0 -> auto)
! p5 OPTIONAL: Slit width in pixels for extraction  (default 15 -> BACHES)
! p6 OPTIONAL: Offset in pixels of the center of the slit (default 0)
! p7 OPTIONAL: Final tolerance on RMS error
! p8 OPTIONAL: Final degree of polynomial defining dispersion
! p9 OPTIONAL: Threshold above background for wavelength calibration
!
! Edit here to change defaults:
!
 DEFINE/PARA p1 ? I "Enter 2-dimensional image with the spectrum of a flat field lamp? "
 DEFINE/PARA p2 ? I "Enter 2-dimensional image with the spectrum of a calibration lamp? "
 DEFINE/PARA p3 24   NUM
 DEFINE/PARA p4 20   NUM
 DEFINE/PARA p5 20  NUM
 DEFINE/PARA p6 0   NUM
 DEFINE/PARA p7 1   NUM
 DEFINE/PARA p8 4   NUM 
! DEFINE/PARA only permits 8 parameter max.
! DEFINE/PARA p9 10  NUM
!
! Convert paramaters into variables:
!
 DEFINE/LOCA cmd/c/1/20         "baches_calib"
 DEFINE/LOCA ff_frame/c/1/40    {p1}
 DEFINE/LOCA lamp_frame/c/1/40  {p2}
 DEFINE/LOCA calib_table/c/1/8  "thar.fit"
 DEFINE/LOCA orders/i/1/1       {p3}
 DEFINE/LOCA widthi/i/1/1       {p4}
 DEFINE/LOCA width/i/1/1        {p5}
 DEFINE/LOCA slit_offset/r/1/1  {p6}
 DEFINE/LOCA tolerance/r/1/1    {p7}
 DEFINE/LOCA poly/i/1/1         {p8}
 DEFINE/LOCA threshold/i/1/1    10
 DEFINE/LOCA std_frame/c/1/40   "std_frame.fit"
 DEFINE/LOCA std_table/c/1/40   "std_table.tfit"

!
! Other local definitions
!
 DEFINE/LOCA ff/c/1/40
 DEFINE/LOCA lamp/c/1/40
 DEFINE/LOCA basename/c/1/40    "" ? +lower_levels
 DEFINE/LOCA ext/c/1/10         "" ? +lower_levels
 DEFINE/LOCA session/c/1/40	"baches"
 DEFINE/LOCA session_file/c/1/40
 DEFINE/LOCA tmpname/c/1/40
 DEFINE/LOCA total_first/i/1/1   0        
 DEFINE/LOCA total_last/i/1/1    0    
 DEFINE/LOCA iteration/i/1/1     0    

!
! MAIN, main:
!
 WRITE/OUT
 WRITE/OUT "PARAMETERS FOR THIS CALIBRATION:"
 WRITE/OUT "==============================="
 WRITE/OUT "Flat field  = {ff_frame}"
 WRITE/OUT "Calibration lamp = {lamp_frame}"
 WRITE/OUT "Calibration table = {calib_table}"
 WRITE/OUT "Num. of orders  = {orders}"
 WRITE/OUT "Slit order width = {widthi}"
 WRITE/OUT "Slit extraction width = {width}"
 WRITE/OUT "Offset from center slit = {slit_offset}"
 WRITE/OUT "Tolerance on RMS = {tolerance}"
 WRITE/OUT "Polynomial degree  = {poly}"
 WRITE/OUT "Threshold = {threshold}"

 IF m$exist(std_frame) .eq. 1 THEN
   INDISK/FIT {std_frame} std_frame >NULL
   WRITE/OUT "Standard photometric star = {std_frame}"
 ELSE
   WRITE/OUT "Spectrum of a standard photometric star = NOT FOUND"
 ENDIF

 IF m$exist(std_table) .eq. 1 THEN
   INDISK/FIT {std_table} std_table >NULL
   WRITE/OUT "Flux table for standard photometric star = {std_table}"
 ELSE
   WRITE/OUT "Flux table for the standard photometric star = NOT FOUND"
 ENDIF
 WRITE/OUT

!
! Ask before continuing
!
WRITE/KEYW reply/c/1/1 y
INQUIRE/KEYW reply "{cmd}: Do you want to continue [yn] (y)?"
IF reply(:1) .eq. "y" THEN
  GOTO continue
ELSE
  RETURN/EXIT
ENDIF

continue:
!
! Set context ECHELLE, if not yet already
!
 SHOW/CONT echelle >NULL
 IF  OUTPUTI(5) .eq. 0 THEN
   SET/CONTEXT echelle >NULL
 ENDIF
 INITIA/ECHE 
 RESET/DIS >NULL

!
! Gets basename and extension of input file
!
 FILE/BACHES {ff_frame}
 ff = basename
!
! Convert FITS format into MIDAS BDF format
!
 IF ext .ne. ".bdf" THEN
   COMPUTE/BACHES {ff_frame} {ff}
 ENDIF

!
! Check if an order identification exist from a previous session
!
 session_file = session + "ORDE.tbl"
 tmpname = session + "ORDE"
 IF m$exist(session_file) .eq. 1 THEN
   WRITE/OUT "{cmd}: An order identification already exist from a previous session"
   WRITE/KEYW reply/c/1/1 y
   INQUIRE/KEYW reply "{cmd}: Do you want to use it [yn] (y)?"
   IF reply(:1) .eq. "y" THEN
     INITIA/ECHE {session}
     GOTO continue_with_line_identification
   ENDIF
 ENDIF

!
! Start creating a display
!
! CREATE/DISP 0 590,390,0,25
! CREATE/DISP 0 1090,730,0,25
 CREATE/DISP 0 1417,949,0,25
! CREATE/DISP 0 1635,1095,0,25

! CREATE/DISP 0 2180,1460,0,0
 LOAD/LUT heat

 CUTS/IMAG {ff} =sigma
 DISPLA/ECHE {ff}
 SET/ECHE WIDTHI={widthi}
 SET/ECHE NBORDI={orders}

! The file order.tbl is created but never removed
! If corrupted, DEFINE/HOUG will never work. So remove it:
 $rm -f $order.tbl

! Next command can also be used to reduced the scan area from rows 1 to 1000
! (bottom-up) in an image of 1024 rows, and with the purpose of eliminating
! a truncated red order. A truncated order may cause hanging in DEFINE/HOUGH
! SCAN/ECHE {ff} 1,1000

! Next command could be used to identifiy with the cursor and a graphical
! diplay the area of scanning:
! SCAN/ECHE {ff} cursor

! DEFINE/HOUG could use a low threshold in an attempt to idenitify weak orders
! DEFINE/HOUG {ff} ? ? ? 500 2,5

! Parameters 1,5 requires the usage of line functions
 DEFINE/HOUG {ff} ? ? ? ? 2,5
 SAVE/ECHE {session} 
!
! Show current ECHELLE parameters before asking to continue
!
 SHOW/ECHE w
 WRITE/KEYW reply/c/1/1 y
 WRITE/OUT
 WRITE/OUT "{cmd}: Order identification finished"
 INQUIRE/KEYW reply "{cmd}: Do you want to continue [yn] (y)?"
 IF reply(:1) .eq. "y" THEN
   GOTO continue_with_line_identification
 ELSE
   RETURN/EXIT
 ENDIF

continue_with_line_identification:

!
! Check if an line identification exist from a previous session
!
 session_file = session + "LINE.tbl"
 tmpname = session + "LINE"
 IF m$exist(session_file) .eq. 1 THEN
   WRITE/OUT "{cmd}: A line identification already exits from a previous session"
   WRITE/KEYW reply/c/1/1 y
   INQUIRE/KEYW reply "{cmd}: Do you want to use it [yn] (y)?"
   IF reply(:1) .eq. "y" THEN
     INITIA/ECHE {session}
     GOTO continue_with_flat_calibration
   ENDIF
 ENDIF

!
! Create two displays: 1 for reference lamp, 0 for input lamp
!
 CLEAR/CHAN OVER
! CREATE/DISP 0 590,390,0,25
! CREATE/DISP 0 1090,730,0,25
 CREATE/DISP 0 1417,949,0,25
! CREATE/DISP 0 1635,1095,0,25
 LOAD/LUT heat
 CREATE/GRAP 0 590,440,590,25
 CREATE/DISP 1 590,390,590,25
 LOAD/LUT heat

 INDISK/FITS MID_HOME:contrib/baches/demo/thar_ref.fit thar_ref.bdf
 INDISK/FITS MID_HOME:contrib/baches/demo/thar_ref_ident.fit thar_ref_ident.tbl
 CUTS/IMAG thar_ref =sigma
 DISPLA/ECHE thar_ref
 CLEAR/CHAN OVER
 LOAD/TABLE thar_ref_ident :X :Y :IDENT 0 6
!
 ASSIGN/DISP 0

!
! Gets basename and extension of calibration lamp frame
!
 FILE/BACHES {lamp_frame}
 lamp = basename
!
! Convert FITS format into MIDAS BDF format
!
 IF ext .ne. ".bdf" THEN
   COMPUTE/BACHES {lamp_frame} {lamp}
 ENDIF

 CUTS/IMAG {lamp} =sigma
! CUTS/IMAG {lamp} 100.0,1420.4
 DISPLA/ECHE {lamp}
!
 SET/ECHELLE WLC={lamp}
 SET/ECHELLE SLIT={width}
 SET/ECHELLE OFFSET={slit_offset}
 SET/ECHELLE WLCOPT=2D
 EXTRAC/ECHELLE {wlc} extwlc
 SET/ECHELLE THRES2={threshold}
 SEARCH/ECHELLE extwlc 
!
 LOAD/SEAR
 WRITE/OUT
 WRITE/OUT "{cmd}: line identification finished"
 WRITE/OUT "{cmd}: if not enough lines, try changing threshold"

 WRITE/KEYW reply/c/1/1 y
 INQUIRE/KEYW reply "{cmd}: Do you want to continue [yn] (y)?"
 IF reply(:1) .eq. "y" THEN
   GOTO continue_next
 ELSE
   WRITE/OUT "{cmd}: SET/ECHE THRES2={threshold}"
   WRITE/OUT "{cmd}: SEARCH/ECHE extwlc"
   WRITE/OUT "{cmd}: LOAD/SEAR"
   RETURN/EXIT
 ENDIF
continue_next:
 SAVE/ECHE {session}
 CLEAR/CHAN OVER

!
! Convert FITS format into MIDAS BDF format
!
 INDISK/FITS MID_HOME:contrib/baches/demo/thar.fit thar.tbl
 LOAD/LUT heat
 SET/ECHE LINCAT = thar.tbl
!
! Start with tol=2 and dc=3 (poly)
! WIDTH2 represents the analysis window width of the calibration lamp;
! equivalent to the FHWM of the emission lines (+2 pixels)
! For BACHES WIDTH2=5, For FLECHAS WIDTH2=15
! wlciter(1) = 1.5
! wlciter(2) = 0.2
 SET/ECHELLE WIDTH2 = 5
 SET/ECHE TOL = 2
 SET/ECHE DC = 3

!
! Change size of the square cursor
!
 SET/CURSOR 0 RECTANGLE 1,1,5,12

!
! First identification using reference file
!
 WRITE/OUT
 WRITE/OUT "{cmd}: Executing IDENTI/ECHE"
 IDENTI/ECHE 
 SET/ECHE GUESS = {session}
 SAVE/ECHE {session}
 SET/ECHE WLCMTD = guess
 DELETE/DISPL 1
!
! Set the tolerance 1.0 higher than requested
!
 tol = {tolerance} + 1
! dc = {poly} + 1

begin_loop1:
! dc = {poly} + 1
 iteration = {iteration} + 1
 $rm -f IDENTI.res
 IDENTI/ECHE >>IDENTI.res
 SHOW/ECHE w
 $tail -13 IDENTI.res
 $grep TOTAL IDENTI.res | awk '{print $7}' | write/key total_last
 IF total_last .gt. total_first THEN
   SAVE/ECHE {session}
   total_first = total_last
   GOTO begin_loop1
 ENDIF

begin_loop2:
 dc = {poly}
 iteration = {iteration} + 1
 $rm -f IDENTI.res
 IDENTI/ECHE >>IDENTI.res
 SHOW/ECHE w
 $tail -13 IDENTI.res
 $grep TOTAL IDENTI.res | awk '{print $7}' | write/key total_last
 IF total_last .gt. total_first THEN
   SAVE/ECHE {session}
   total_first = total_last
   GOTO begin_loop2
 ENDIF

! IF tol .gt. tolerance then
!   tol = {tol} - 0.2
!   total_first = 0
!   GOTO begin_loop1
!   GOTO begin_loop2
! The method ROBUST works great in Cygwin but not on other Linux
! ELSE
   tol = {tolerance}
   dc = {poly}
   SET/ECHE WLCREG = ROBUST
   IDENTI/ECHE 
   SAVE/ECHE {session}
! ENDIF
 WRITE/OUT "{cmd}: Identification completed after {iteration} iterations."

continue_with_flat_calibration:
!
! Background subtraction from flat field
! Set FFOPT=NO if you do not want flat-field division
! SET/ECHELLE FLAT={ff} SAMPLE=0.15 SLIT=5.0 BKGMTD=SPLINE BKGSMO=50.
 SET/ECHELLE BKGRAD=6,2 BKGDEG=3 BKGSMO=500. BKGMTD=POLY 
 SET/ECHELLE BKGVIS=NO
 SET/ECHELLE FLAT={ff} 
 FLAT/ECHELLE 
 SAVE/ECHE {session}
 SET/ECHELLE DELTA = 3 
 SET/ECHELLE MRGMTD  = AVERAGE
 MERGE/ECHE {BLAZE} {BLAZE}m
! CREA/GRAPH
! SET/GRAPH default=
! PLOT  {BLAZE}m

!
! Flat correction (pixel-to-pixel) is a very risky task as 
! the division of images may cause very weird results.
 WRITE/KEYW reply/c/1/1 y
 INQUIRE/KEYW reply "{cmd}: Do you want flat correction [yn] (y)?"
 IF reply(:1) .eq. "y" THEN
   SET/ECHELLE FFOPT = YES
 ELSE
   SET/ECHELLE FFOPT = NO
 ENDIF
! Default: no flat correction
! SET/ECHELLE FFOPT = NO

flux_calibration:

 WRITE/KEYW reply/c/1/1 y
 INQUIRE/KEYW reply "{cmd}: Do you want to correct flux with a standard star [yn] (y)?"
 RESET/DISPLAY
 IF reply(:1) .eq. "y" THEN
   SET/ECHELLE RESPOPT = YES
 ELSE
   SET/ECHELLE RESPOPT = NO
   GOTO ask_blackbody
 ENDIF
!
! Flux calibration if we have both the raw standard star and its table
! recognized with the filenames {std_frame} and {std_table}
!
 SET/ECHELLE RESPMTD = STD

 IF m$exist(std_frame) .eq. 1 THEN
   INDISK/FIT {std_frame} std_frame
   WRITE/OUT "{cmd}: converting {std_frame} to std_frame.bdf"
   SET/ECHELLE STD=std_frame
 ELSE
   WRITE/OUT "{cmd}: {std_frame} not found"
   GOTO ask_blackbody
 ENDIF

 IF m$exist(std_table) .eq. 1 THEN
   INDISK/FIT {std_table} std_table
   WRITE/OUT "{cmd}: converting {std_table} to std_table.tbl"
   SET/ECHELLE FLUXTAB=std_table
   GOTO response
 ELSE
   WRITE/OUT "{cmd}: {std_table} not found"
   GOTO ask_blackbody
 ENDIF

ask_blackbody:
 WRITE/KEYW reply/c/1/1 y
 INQUIRE/KEYW reply "{cmd}: Do you want to correct flux with {ff} as a blackbody [yn] (y)?"
 IF reply(:1) .eq. "y" THEN
   SET/ECHELLE RESPOPT = YES
   SET/ECHELLE RESPMTD = STD
   SET/ECHELLE STD={ff}
   BLACKBODY/BACHES 3200
   SET/GRAPH default= 
   SET/ECHELLE FLUXTAB=blackbody_3200
   GOTO response
 ELSE
   WRITE/KEYW reply/c/1/1 y
   INQUIRE/KEYW reply "{cmd}: Do you want to remove the blaze (flatten) [yn] (y)?"
   IF reply(:1) .eq. "y" THEN
     SET/ECHELLE RESPOPT = YES
     SET/ECHELLE RESPMTD = FLAT
     GOTO no_response
   ELSE
     SET/ECHELLE RESPOPT = NO
     GOTO no_response
   ENDIF
 ENDIF

response:
 RESPONSE/ECHELLE
 CLEAR/CHAN OVER
! LOAD response  cuts=0.,55000.  scale=-2,16

no_response:
 SAVE/ECHE {session}
! MERGE/ECHELLE &a name 3.0 AVERAGE

 WRITE/KEYW reply/c/1/1 y
 INQUIRE/KEYW reply "{cmd}: Do you want to calculate R for {lamp_frame} [yn] (y)?"
 IF reply(:1) .eq. "y" THEN
   GOTO calculateR
 ELSE
   GOTO nocalculateR
 ENDIF

!
! calculateR/nocalculateR
!
calculateR:
!
! Gets basename and extension of calibration lamp frame
!
 FILE/BACHES {lamp_frame}
 lamp = basename
!
! Convert FITS format into MIDAS BDF format
!
 IF ext .ne. ".bdf" THEN
   COMPUTE/BACHES {lamp_frame} {lamp}
 ENDIF

 WRITE/OUT "REDUCE/ECHELLE {lamp} {lamp}_wrm"
 REDUCE/ECHELLE {lamp} {lamp}_spectrum
! @@ baches_pipeline {lamp}
 WRITE/OUT "RESOLV/BACHES {lamp}_spectrum"
 RESOLV/BACHES {lamp}_spectrum
!
! Plot the four horses:
!
! SET/GRAPH YAXIS=0.,400.
! SET/GRAPH YAXIS=6570.,6600.
! PLOT {lamp}_wrm

!
 ECHO/OFF
 RETURN/EXIT

nocalculateR:
! RESET/DISP
 ECHO/OFF
 RETURN/EXIT
