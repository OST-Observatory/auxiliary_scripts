!
! USAGE: @@ dados_stack.prg tostack_1.fits tostack_2.fits ... tostack_8.fits
! Used to stack 2d images (and save result as .fit)
!


!!
!  clean:
!!

$ rm -f *.bdf *.tbl *.cat *.ascii *.plt *.KEY


!!
!  read input
!!

DEFINE/PAR P1 + f "Image to stack"
DEFINE/PAR P2 + f "Image to stack"
DEFINE/PAR P3 + f "Image to stack"
DEFINE/PAR P4 + f "Image to stack"
DEFINE/PAR P5 + f "Image to stack"
DEFINE/PAR P6 + f "Image to stack"
DEFINE/PAR P7 + f "Image to stack"
DEFINE/PAR P8 + f "Image to stack"

DEFINE/LOCA basename/c/1/80   ""
DEFINE/LOCA indx/i/1/1         0
DEFINE/LOCA name/c/1/80        {P1}
indx = m$index(name,".")-1
basename = name(1:{indx})

SET/CONT    long

CREATE/ICAT scienceimages null

IF P1(1:1) .ne. "+" THEN
    INDISK/FITS 'P1' stack_1
    COPY/DKEY stack_1 O_TIME texp
    WRITE/DESCR stack_1 WEIGHT/d/1/10 {texp(7)}
    READ/DESCR stack_1 WEIGHT
    ADD/ICAT scienceimages stack_1.bdf

    WRITE/KEYW  revnorm/R/1/1 1
ENDIF
IF P2(1:1) .ne. "+" THEN
    INDISK/FITS 'P2' stack_2
    COPY/DKEY stack_2 O_TIME texp 
    WRITE/DESCR stack_2 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_2.bdf
    WRITE/KEYW  revnorm/R/1/1 2
ENDIF
IF P3(1:1) .ne. "+" THEN
    INDISK/FITS 'P3' stack_3
    COPY/DKEY stack_3 O_TIME texp 
    WRITE/DESCR stack_3 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_3.bdf
    WRITE/KEYW  revnorm/R/1/1 3
ENDIF
IF P4(1:1) .ne. "+" THEN
    INDISK/FITS 'P4' stack_4
    COPY/DKEY stack_4 O_TIME texp 
    WRITE/DESCR stack_4 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_4.bdf
    WRITE/KEYW  revnorm/R/1/1 4
ENDIF
IF P5(1:1) .ne. "+" THEN
    INDISK/FITS 'P5' stack_5
    COPY/DKEY stack_5 O_TIME texp 
    WRITE/DESCR stack_5 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_5.bdf
    WRITE/KEYW  revnorm/R/1/1 5
ENDIF
IF P6(1:1) .ne. "+" THEN
    INDISK/FITS 'P6' stack_6
    COPY/DKEY stack_6 O_TIME texp 
    WRITE/DESCR stack_6 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_6.bdf
    WRITE/KEYW  revnorm/R/1/1 6
ENDIF
IF P7(1:1) .ne. "+" THEN
    INDISK/FITS 'P7' stack_7
    COPY/DKEY stack_7 O_TIME texp 
    WRITE/DESCR stack_7 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_7.bdf
    WRITE/KEYW  revnorm/R/1/1 7
ENDIF
IF P8(1:1) .ne. "+" THEN
    INDISK/FITS 'P8' stack_8
    COPY/DKEY stack_8 O_TIME texp 
    WRITE/DESCR stack_8 WEIGHT/d/1/10 {texp(7)}
    ADD/ICAT scienceimages stack_8.bdf
    WRITE/KEYW  revnorm/R/1/1 8
ENDIF


!!
!  stack images
!!

!COMBINE/LONG scienceimages.cat stacked_0 AVERAGE
COMBINE/LONG scienceimages.cat stacked_0 MEDIAN
! this uses texp as weight for stacking
!AVERAGE/WEIGHTS stacked_0 = scienceimages.cat

! this can be used to get cat weigths from med, sig in window
!CREATE/DISPLAY ? 920,430,1000,1080
!LOAD/IMAG   stack_1 cuts=f scale=full
!COMPUTE/WEIGHTS scienceimages.cat

!AVERAGE/IMAGE stacked_0 = scienceimages.cat
! this does not make sense
!COMPUTE/IMAGE stacked = stacked_0*{revnorm}
COMPUTE/IMAGE stacked = stacked_0

OUTDISK/FITS stacked {basename}_stacked.fit

!!
!  clean
!!

$ rm -f *.bdf *.tbl *.cat *.KEY *.plt *.ascii

ECHO/OFF
RETURN/EXIT
