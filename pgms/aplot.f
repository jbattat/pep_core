C       This file contains both "PLOT" and "AIPLOT"
C       initial implementation 4/27/87 by SEB@crnlns.bitnet
C
C       Modified 11/6/87 to decrese duplicate lines and
C       to limit the number of vectors in a "stroke" to 100 lines
C       Modified Jan 15, 1988 to add %%BoundingBox. seb
C       Modified 1992 for non-VAX Fortran, OFFSET, bug fixes.
      SUBROUTINE PLOT(X,Y,IPENIN)
C
C       This is the Adobe PostScript device dependent interface routine.
C       It is still one level above the hardware dependent routine.
C       This is where the current floating point plotting offsets
C       scale factors, and other standard information is kept.
C
C       PLOT is called with relative coordinates in units of inches.
C       It calls AIPLOT with absolute coordinates in units of 1/300 in.
C
      IMPLICIT INTEGER(A-W)
      IMPLICIT REAL(X-Z)
      PARAMETER (SYMSHFT=128,ROMSHFT=129)
      REAL*4 RIX,RIY,RPWID,XW,YW,SC
      CHARACTER*8 ATIME
Csys      CHARACTER*9 ADATE
      CHARACTER*8 ADATE
Cend
      CHARACTER*20 DATE_WORK
      CHARACTER*80 FILE_NAME
      CHARACTER*(*) UQNAME
      CHARACTER*(*) HWSTR
      REAL*4 HWSIZE,HWANG
      INCLUDE '../peputil/plttyp.inc'
C       Common block to keep track of info:
      INCLUDE 'pltopnfo.inc'
C       MUNITS: Currently selected ploter devices
C       LPOPEN: An open has been done (and data written to files, maybe)
C       LTRACE: .true. = type subr. name and args
C       MARGIN: 1 = shrink plot to fit 1-inch margins
C     DATA LPOPEN/.FALSE./
      COMMON/AdobePS0/LQUNIT,NPAGES
      COMMON/AdobePS1/LQNAME,APEN
      CHARACTER*80 LQNAME
C dash patterns: solid . _ = ._ .= _= .._ ..= ._= ..._ ...= ..=_ _...=
C                          1   5   10   15   20   25   30
C                          |   |    |    |    |    |    |
      CHARACTER*32 PENPAT/'20 15 5 15 5 15 5 15 55 15 20 15'/
      INTEGER PENSTRT(14),PENEND(14),P1,P2
      DATA PENSTRT/3, 7, 1, 22,  1, 17, 22,  1, 12, 17, 1,  7, 12, 1/
      DATA PENEND/ 3,10, 5, 26, 10, 26, 32, 15, 26, 32, 20, 26,32,26/
      CHARACTER*1 APEN
C output default unit and file name
C LQUNIT should be 53, LQNAME should be 'PLOT.PS', APEN should be '0'
C Default pen color: 0 = black, 1 = white = erase
C Default pen width
      DATA IPWID/1/
C
      IF(.NOT.LPOPEN) THEN
        ASSIGN 88882 TO I88888
        GOTO 88888
        ENDIF
88882 CONTINUE
      IPEN=IPENIN
      IF(LTRACE) WRITE(6,*) '%PLOT:  X= ',X,', Y= ',Y,', pen= ',ipen
      IPENCTRL=IABS(IPENIN)/10
      IF(IPENCTRL.LE.2) IPEN=MOD(IPEN,10)
      IF(IPEN.LT.0 .AND. IPEN.NE.-999) IPENCTRL=2
      IPENABS=IABS(IPEN)
      IF(IPENABS.EQ.1) IPEN=LASTPEN*IPEN
      IF(IPENABS.EQ.2.OR.IPENABS.EQ.3) LASTPEN=IPENABS
C
C SAVE CURRENT POSITION (with offsets, if any)
      IF(IPENCTRL.EQ.1) THEN
        X=(X-XOFF)/XFACT
        Y=(Y-YOFF)/YFACT
        ENDIF
      X1=X
      Y1=Y
C
C IPEN=999 to close file
      IF(IPEN.EQ.999)THEN
        CALL APOUT
        WRITE(UNIT=LQUNIT,FMT='(1X,A)') 'showpage'
        NPAGES=NPAGES+1
        WRITE(UNIT=LQUNIT,FMT='(A)') '%%Trailer'
        WRITE(UNIT=LQUNIT,FMT='(1X,A)') 'grestore'
        WRITE(UNIT=LQUNIT,FMT='(A,I4)')'%%Pages:',NPAGES
C       1.3 add "eof"
        WRITE(UNIT=LQUNIT,FMT='(A)') char(4)
        CLOSE(UNIT=LQUNIT)
        LPOPEN=.FALSE.
        RETURN
        ENDIF
C
C IPEN=-999 to go to next page
      IF(IPEN.EQ.-999)THEN
C Finish this page's output
         CALL APOUT
C display current plot and do a form-feed
         WRITE(UNIT=LQUNIT,FMT='(A)') ' showpage'
C note what page this is for possible post processing
         NPAGES=NPAGES+1
         WRITE(UNIT=LQUNIT,FMT='(A,I4,A)') '%%Page: ',NPAGES, ' ?'
C reinitialize plotting parameters
         XLT = -10.9
         YLT = 0.25
C       cvt from default of 1/72 inch to 1/300 inch
         SC = 72.0/300.0
         IF(MARGIN.EQ.1) THEN
            XLT = -9.9
            YLT = 1.0
            SC = 72.0/380.0
         ENDIF
         XLT = XLT * 72.0 / SC
         YLT = YLT * 72.0 / SC
C Move origin to other edge of paper (rotated 90 deg)
         IXLT = XLT
         IYLT = YLT
         WRITE(UNIT=LQUNIT,FMT=505) SC,SC, RPWID,APEN,IXLT,IYLT
  505    FORMAT(
     1  2(' ',F4.3),' scale'/
     1  ' ',F4.1,' setlinewidth 1 setlinecap 1 setlinejoin'/
     1  ' [] 0 setdash ',A1,' setgray'/
     1  ' -90 rotate ',I6,1X,I6,' translate')
         CALL AIPLOT(0.,0.,-3)
C Reset everything to 0,0
         X0=0.
         Y0=0.
         X1=0.
         Y1=0.
         X2=0.
         Y2=0.
         LASTPEN=3
         RETURN
      ENDIF
C IPEN=0 reset current location to (X,Y)
      IF(IPENIN.EQ.0) THEN
        X0=X2-X
        Y0=Y2-Y
        RETURN
        ENDIF
C Plot user's data
C Convert to f.p. absolute coords
      X2=X+X0
      Y2=Y+Y0
C Scale to user's size
      X3=X2*XYSCL
      Y3=Y2*XYSCL
C Convert to integer screen positions
      RIX=X3*XS
      RIY=Y3*YS
      IF(LTRACE)THEN
        WRITE(6,*) '%PLOT: (X+X0)*XYSCL*XS = (',
     1  X,' +',X0,' )*',XYSCL,' *',XS,' =',RIX
        WRITE(6,*) '%PLOT: (Y+Y0)*XYSCL*YS = (',
     1  Y,' +',Y0,' )*',XYSCL,' *',YS,' =',RIY
        ENDIF
C
C IPEN=1,2,3 = plot vector
C call hardware dependent routine
      CALL AIPLOT(RIX,RIY,IPEN)
      IF(IPENCTRL.EQ.2)THEN   !RESET ORIGIN
C set origin to current physical location
        X0 =X2
        Y0 =Y2
C set current position to 0,0
        X1 =0.
        Y1 =0.
        CALL APOUT
        ENDIF
      RETURN
C----------------------------------------------------------------------
      ENTRY PLTOPN
C       Select another AdobePostScript output unit: ain't none
      IF(LTRACE) WRITE(6,*) '%PLTOPN'
      RETURN
C----------------------------------------------------------------------
      ENTRY PLTFIL(UQNAME)
C       Set name of AdobePostScript plot file
      IF(LTRACE) WRITE(6,*) '%PLTFIL: UQNAME=',UQNAME
      LQNAME=UQNAME
      RETURN
C----------------------------------------------------------------------
      ENTRY NEWPEN(IPNUM)
      IF(LTRACE) WRITE(6,*) '%NEWPEN: IPNUM=',IPNUM
C       ensure all pending points are out,
C       since this writes directly to the file
      CALL APOUT
C       Set line pattern
      P1=PENSTRT(IPNUM+1)
      P2=PENEND(IPNUM+1)
      WRITE(UNIT=LQUNIT,FMT='(A)')
     . ' ['//PENPAT(P1:P2)//'] 0 setdash'
      RETURN
C----------------------------------------------------------------------
      ENTRY PENWID(IPNUM)
      IF(LTRACE) WRITE(6,*) 'PENWID: IPNUM=',IPNUM
C       ensure all pending points are out,
C       since this writes directly to the file
      CALL APOUT
C       Set line width
      IPWID=ABS(IPNUM)
      RPWID = IPWID
      IF(IPNUM.LT.0)THEN
        APEN='1'
      ELSE
        APEN='0'
        ENDIF
      WRITE(UNIT=LQUNIT,FMT='(1X,F4.1,A)')
     1RPWID,' setlinewidth '//APEN//' setgray'
      RETURN
C----------------------------------------------------------------------
      ENTRY PENUP
      IF(LTRACE) WRITE(6,*) '%PENUP'
C Raise pen for subsequent calls with IPEN=1
      LASTPEN=3
      RETURN
C----------------------------------------------------------------------
      ENTRY PENDN
      IF(LTRACE) WRITE(6,*) '%PENDN'
C Lower pen for subsequent calls with IPEN=1
      LASTPEN=2
      RETURN
C----------------------------------------------------------------------
      ENTRY PLOTS    !INITIALIZE SELECTED UNITS
      IF(LTRACE) WRITE(6,*) '%PLOTS'
      ASSIGN 88881 TO I88888
88888 CONTINUE
C
C       AdobePostScript units are 1/72 of an inch (aprox 1 "point")
C       XS = 72.0 * 1.005
C       YS = 72.0 * 0.996 ! make 1 inch be so on Apple LaserWriterPlus
      XS = 300.0 * 1.005
      YS = 300.0 * 0.996 ! make 1 inch be so on Apple LaserWriterPlus
      X0=0.
      Y0=0.
      APEN='0'
      LASTPEN=3           ! Start with pen up
C User's SCALE FACTOR (TO BE MODIFIED BY FACTOR)
      XYSCL=1.0
C Initialize offsets and scale factors
      XOFF=0.0
      XFACT=1.0
      YOFF=0.0
      YFACT=1.0
C Initialize other things as required
C according to Adobe's PostScript Language Reference Manual
C appendix D.4, the Apple LaserWriter's page is 10.92 by 8.00
      XLT = -10.9
      YLT = 0.25
C       cvt from default of 1/72 inch to 1/300 inch
      SC = 72.0/300.0
      IF(MARGIN.EQ.1) THEN
         XLT = -9.9
         YLT = 1.0
         SC = 72.0/380.0
      ENDIF
      XLT = XLT * 72.0 / SC
      YLT = YLT * 72.0 / SC
C       Move origin to other edge of paper (rotated 90 deg)
      IXLT = XLT
      IYLT = YLT
C initialize pen width
      RPWID = ipwid
C Open output file
      CALL ASGN_PLT_FILE(LQUNIT,LQNAME)
      OPEN(UNIT=LQUNIT,STATUS='NEW',FILE=LQNAME)
Csys       CARRIAGECONTROL='LIST',
      INQUIRE(UNIT=LQUNIT,NAME=FILE_NAME)
      NPAGES=0
C Initialize laser printer:
C start with Adobe standard comments
C Trim file name.
      LF = LEN(FILE_NAME)
      DO WHILE (LF.GT.1 .AND.
     1 (FILE_NAME(LF:LF).EQ.' '.OR.FILE_NAME(LF:LF).EQ.CHAR(0)))
        LF=LF-1
        ENDDO
  498 CONTINUE
C       Get current date & time
Csys      CALL DATE(ADATE)
Csys      CALL TIME(ATIME)
      CALL TODAY(DATE_WORK)
      ADATE=DATE_WORK(1:8)
      ATIME=DATE_WORK(13:20)
Cend
      WRITE(UNIT=LQUNIT,FMT=500) FILE_NAME(:LF), ADATE,ATIME
  500 FORMAT(
     1 '%!PS-Adobe-1.0'/
     1 '%%Creator: APLOT(v1.7x)'/
     1 '%%Title: ',A/
     1 '%%CreationDate: ',A,1X,A/
     1 '%%Pages: (atend)'/
     1 '%%BoundingBox: 000 000 612 792'/
     1 '%%EndComments')
      WRITE(UNIT=LQUNIT,FMT=504) SC,SC,RPWID,APEN,IXLT,IYLT
C The following definitions are needed for the commands to work
C mostly abbreviations to reduce plot file size.
  504 FORMAT (
     1 ' gsave '/
     1 ' /ps {stroke} def'/
     1 ' /pu {moveto} def'/
     1 ' /pd {lineto} def'/
C       dt = make a tiny dot at provided x,y (not used yet)
C       1' /dt {moveto currentpoint lineto stroke} def'/
C       ny = draw line using previous x value : prefix by new y
     1 ' /ny {currentpoint pop exch lineto} def'/
C       nx = draw line using previous y value : prefix by new x
     1 ' /nx {currentpoint exch pop lineto} def'/
     . ' /fm {exch dup .7 mul exch 0 exch 0 exch 0 0 6 array astore'/
     . 'exch 6 array rotate 6 array concatmatrix} bind def'/
     . ' /rf {fm /Courier exch selectfont} def'/
     . ' /sf {fm /Symbol exch selectfont} def'/
     . '%%Page: 0 ?'/
     1 2(' ',f4.3),' scale'/
     1 ' ',f4.1,' setlinewidth 1 setlinecap 1 setlinejoin'/
     1 ' [] 0 setdash ',A1,' setgray'/
C
C       The laser printer will display things in "portrait" orientation
C       while the upper levels of this code assume that it will be
C       in "landscape" orientation.
C       so transform the coordinates so they'll come out "right"
C
     1 ' -90 rotate ',i6,' ',i6,' translate')
      CALL AIPLOT(0.,0.,-3)
      LPOPEN=.TRUE.
      GOTO I88888,(88881,88882)
88881 CONTINUE
      RETURN
C----------------------------------------------------------------------
      ENTRY FACTOR(YXSCL)
      IF(LTRACE) WRITE(6,*) '%FACTOR: YXSCL=',YXSCL
C       Set scale factor for selected devices
      XYSCL=YXSCL
      RETURN
C----------------------------------------------------------------------
      ENTRY OFFSET(ZXOFF,ZXFACT,ZYOFF,ZYFACT)
      IF(LTRACE) WRITE(6,*) '%OFFSET: X OFF,FACT=',ZXOFF,ZXFACT,
     1 ' Y OFF,FACT=',ZYOFF,ZYFACT
C Set offsets and scale factors
      XOFF=ZXOFF
      XFACT=ZXFACT
      YOFF=ZYOFF
      YFACT=ZYFACT
C For safety...
      IF(XFACT.EQ.0.0) XFACT=1.0
      IF(YFACT.EQ.0.0) YFACT=1.0
      RETURN
C----------------------------------------------------------------------
      ENTRY ERASE
C ERASE SELECTED DISPLAYS
      IF(LTRACE) WRITE(6,*) '%ERASE'
      X0=0.
      Y0=0.
      X1=0.
      Y1=0.
      X2=0.
      Y2=0.
      CALL APOUT
      RETURN
C----------------------------------------------------------------------
      ENTRY WHERE(XW,YW)
C       Return current relative position on selected device
      XW=X1
      YW=Y1
      IF(LTRACE) WRITE(6,*) '%WHERE: XW=',XW,', YW=',YW
      RETURN
C----------------------------------------------------------------------
      ENTRY AORIGN
C       Reset coords to origin, but don't draw anything
      IF(LTRACE) WRITE(6,*) '%AORIGN'
      CALL AIPLOT(0.,0.,-3)
      X0=0.
      Y0=0.
      X1=0.
      Y1=0.
      X2=0.
      Y2=0.
      CALL AIPLOT(0.,0.,3)
      RETURN
C--------------------------------------------------------------------
      ENTRY SYMBHW(HWSTR,HWLEN,HWSIZE,HWANG)
      IF(LTRACE) WRITE(6,*) '%SYMBHW: L=',HWLEN,' STR=',HWSTR(1:HWLEN),
     . '  SIZE=S*XYSCL*YS=',HWSIZE,'*',XYSCL,'*',YS,'  ANGLE=',HWANG
      CSIZE=HWSIZE*XYSCL*YS*1.54+0.5
      CANG=AMOD(HWANG,360.)+360.5
      IF(CANG.GE.360) CANG=CANG-360
      CALL AISYMB(HWSTR,HWLEN,CSIZE,CANG)
      PLTLEN=0
      DO IPTR = 1,HWLEN
         IF(HWSTR(IPTR:IPTR).NE.CHAR(SYMSHFT) .AND.
     .    HWSTR(IPTR:IPTR).NE.CHAR(ROMSHFT)) PLTLEN=PLTLEN+1
      END DO 
      X1=X1+HWSIZE*XYSCL*PLTLEN*COS(HWANG*1.7453292E-2)*PROPOR
      Y1=Y1+HWSIZE*XYSCL*PLTLEN*SIN(HWANG*1.7453292E-2)
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE ASGN_PLT_FILE(LQUNIT,LQNAME)
      CHARACTER*(*) LQNAME
      DATA NITUNT/52/

      NITUNT=NITUNT+1
      LQUNIT=NITUNT
      WRITE(LQNAME,100) LQUNIT
  100 FORMAT('plot',i2,'.ps')
      RETURN
      END
C--------------------------------------------------------------------
      SUBROUTINE AIPLOT(RIX,RIY,ICPEN)
C       DEVICE DEPENDENT AdobePostScript PLOTTING ROUTINE
      IMPLICIT INTEGER(A-Z)
      PARAMETER (MAXVEC=20)
      PARAMETER (SYMSHFT=128,ROMSHFT=129)
      REAL*4 RIX,RIY,RMIX,RMIY,RMIXP,RMIYP
      CHARACTER*255 LINE2WRITE,PLINE2WRITE
      CHARACTER*70 OLINE(0:MAXVEC)
      CHARACTER*23 FVEC
      CHARACTER*15 AX,AY,APX,APY
      CHARACTER*(*) HWSTR
      INTEGER NVPS                    !# vec in current stroke
      INTEGER NDUP                    !# duplicate lines in a row
      INTEGER NSTORE(0:MAXVEC)
      LOGICAL LUAGAIN,NOPS
      INCLUDE 'pltopnfo.inc'
      COMMON/AdobePS0/ LQUNIT,NPAGES
      COMMON/AdobePS1/ LQNAME,APEN
      CHARACTER*80 LQNAME
      CHARACTER*1 APEN
      DATA LASTP/3/,NVEC/0/,NVPS/0/,NDUP/0/
      CHARACTER*2 CURFONT,HWFONT
      DATA CURSIZE/-999/,CURANG/-999/,CURFONT/'  '/
      IF(LTRACE)
     1  WRITE(6,*)  '%AIPLOT: X= ',RIX,', Y= ',RIY,', pen= ',ICPEN
C       IPEN = 2, PEN DOWN
C       IPEN = 3, PEN UP
C       IF NEGATIVE, RE-ORIGIN: ASSUME HANDLED BY PLOT
C       this routine expects absolute values in units of 1/300 in.
      IPEN=IABS(ICPEN)
      IF (IPEN.GT.3 .OR. IPEN.LT.2) RETURN
      NOWPEN=IPEN
C       ENSURE X and y POSITIVE
      RMIX=AMAX1(0.0,RIX)
      RMIY=AMAX1(0.0,RIY)
C don't redraw a line to the current position
      IF(RMIX.EQ.RMIXP .AND. RMIY.EQ.RMIYP .AND. NOWPEN.EQ.LASTP) RETURN
C TO STORE-VECTOR
C try to keep the vector list from wrapping when typing the file
      IF(NVEC.GT.0)THEN
        LLN = 0
        DO NC=1,NVEC
          LLN = LLN+NSTORE(NC)
          ENDDO
        IF(LLN.GT.72)THEN
C          WRITE-BUFFER
          ASSIGN 99992 TO I99993
          GO TO 99993
          ENDIF
        ENDIF
99992 CONTINUE
C A "stroke" must start with a "moveto" setting the "currentpoint"
C so always start a new stroke when we see a new "moveto"
C (not optimal, but simplifies other optimizations)
      IF(LASTP.EQ.2 .AND. NOWPEN.EQ.3)THEN
         IF(NSTORE(NVEC).GE.67) THEN
            NVEC=NVEC+1
            NSTORE(NVEC)=0
         ENDIF
         IF(NSTORE(NVEC).LE.0) THEN
            NSTORE(NVEC)=1
            OLINE(NVEC) = ' '
         ENDIF
         OLINE(NVEC) = OLINE(NVEC)(:NSTORE(NVEC))//' ps'
         NSTORE(NVEC) = NSTORE(NVEC) + 3
         NVPS  = 0
C        WRITE-BUFFER
         ASSIGN 99991 TO I99993
         GO TO 99993
      ENDIF
99991 CONTINUE
C Ascii coordinates
      MIX = RMIX + 0.5
      IF(MIX.LT.1)THEN
        AX ='0'
        NWX = 1
      ELSE
        WRITE(UNIT=AX,FMT='(1X,I10)') MIX
        NWX = 11
        CALL APLOT_UNPAD(AX,NWX)
        ENDIF
      MIY = RMIY + 0.5
      IF(MIY.LT.1)THEN
        AY ='0'
        NWY = 1
      ELSE
        WRITE(UNIT=AY,FMT='(1X,I10)') MIY
        NWY = 11
        CALL APLOT_UNPAD(AY,NWY)
        ENDIF
      IF(LTRACE)
     1  WRITE(6,*) '%AIPLOT: (ascii) X= ',ax(:nwx),', Y= ',ay(:nwy)
C NVEC points to the plot command previously entered
C (or maybe to 0 - before the buffer)
C Pen down=2, pen up=3
      IF(NOWPEN.EQ. 3)THEN
C Don't gen a new command if we are just moving around,
C just overwrite previous positioning command.
        LUAGAIN=.FALSE.
        IF(NVEC.GE.1) THEN
          IF(NSTORE(NVEC).GT.2)
     1      LUAGAIN=(OLINE(NVEC)(NSTORE(NVEC)-2:NSTORE(NVEC)).EQ.' pu')
          ENDIF
        IF(.NOT.LUAGAIN) NVEC=NVEC+1
        OLINE(NVEC) = ' '//AX(:NWX)//' '//AY(:NWY)//' pu'
        NSTORE(NVEC) = NWX+NWY+5
      ELSE
C Always enter a "pen down" command as a new entry
        NVEC=NVEC+1
C Optimize data transmission somewhat:
C if doing a "pd" ("lineto"),
C then only send changed part of coordinate.
        IF(AX.EQ.APX)THEN
C AY .ne. APY or wouldn't get here
          OLINE(NVEC)  = ' '//AY(:NWY)//' ny'
          NSTORE(NVEC) = NWY+4
        ELSEIF(AY.EQ.APY)THEN
C AX .ne. APX or wouldn't get here
          OLINE(NVEC)  = ' '//AX(:NWX)//' nx'
          NSTORE(NVEC) = NWX+4
        ELSE
C AX AND AY both new
          OLINE(NVEC) = ' '//AX(:NWX)//' '//AY(:NWY)//' pd'
          NSTORE(NVEC) = NWX+NWY +5
          ENDIF
        ENDIF
      APX = AX
      APY = AY
      RMIXP = RMIX
      RMIYP = RMIY
      NVPS = NVPS + 1
C Save pen state for comparisons next time
      LASTP=NOWPEN
C Start a new series of vectors after re-origining.
C Try to keep from overflowing printer buffer too.
      IF(ICPEN.LT.0 .OR. NVPS.GT.100) GOTO 99997
      RETURN
C==================================================================
      ENTRY APOUT
      IF(LTRACE) WRITE(6,*) '%APOUT'
C
C        TO WRITE-STROKE
99997 CONTINUE
C Guarantee that the previous vector list is properly terminated
C (but code loses track if there are multiple "pen up" commands)
      IF(NSTORE(NVEC).GE.67) THEN
         NVEC=NVEC+1
         NSTORE(NVEC)=0
      ENDIF
      IF(NSTORE(NVEC).LE.0)THEN
         NSTORE(NVEC)=1
         OLINE(NVEC) = ' '
      ENDIF
      OLINE(NVEC) = OLINE(NVEC)(:NSTORE(NVEC))//' ps'
      NSTORE(NVEC) = NSTORE(NVEC) + 3
      NVPS = 0
C        WRITE-BUFFER
      ASSIGN 99990 TO I99993
      GO TO 99993
99990 CONTINUE
C and that the next one (IF ANY) will be properly initiated
C Ascii coordinates
      MIX = IFIX(RMIXP + 0.5)
      WRITE(UNIT=ax,FMT='(1X,i10)') MIX
      NWX = 11
      CALL APLOT_UNPAD(AX,NWX)
      MIY = IFIX(RMIYP + 0.5)
      WRITE(UNIT=AY,FMT='(1X,i10)') MIY
      NWY = 11
      CALL APLOT_UNPAD(AY,NWY)
C starting a new vector list
      NVEC=1
      OLINE(NVEC) = ' '//AX(:NWX)//' '//AY(:NWY)//' pu'
      NSTORE(NVEC) = NWX+NWY+5
      RETURN
C==================================================================
      ENTRY AISYMB(HWSTR,HWLEN,HWSIZE,HWANG)
      IF(LTRACE) WRITE(6,*) '%AISYMB:  STR=',HWSTR(1:HWLEN),
     . '  SIZE=',HWSIZE,' ANGLE=',HWANG
      IF(HWSTR(1:1).EQ.CHAR(SYMSHFT)) THEN
         HWFONT='sf'
      ELSE
         HWFONT='rf'
      ENDIF
      IF(HWSIZE.NE.CURSIZE .OR. HWANG.NE.CURANG .OR.
     . HWFONT.NE.CURFONT) THEN
         WRITE(UNIT=ax,FMT='(1X,i10)') HWSIZE
         NWX = 11
         CALL APLOT_UNPAD(AX,NWX)
         WRITE(UNIT=AY,FMT='(1X,i10)') HWANG
         NWY = 11
         CALL APLOT_UNPAD(AY,NWY)
         FVEC=' '//AX(:NWX)//' '//AY(:NWY)//' '
         FLEN=NWX+NWY+3
         IF(NSTORE(NVEC).GT.1) NVEC=NVEC+1
         OLINE(NVEC) = FVEC(:FLEN)//HWFONT
         NSTORE(NVEC) = FLEN+2
         CURSIZE=HWSIZE
         CURANG=HWANG
         CURFONT=HWFONT
      ENDIF

C try to keep the vector list from wrapping when typing the file
      IF(NVEC.GT.0)THEN
         LLN = 0
         DO NC=1,NVEC
            LLN = LLN+NSTORE(NC)
         ENDDO
         IF(LLN.GT.72)THEN
C  WRITE-BUFFER
            ASSIGN 99988 TO I99993
            GO TO 99993
         ENDIF
      ENDIF
99988 CONTINUE

      IF(NSTORE(NVEC)+HWLEN.GT.50 .OR. NVEC.EQ.0) THEN
         NVEC=NVEC+1
         NSTORE(NVEC)=0
      ENDIF
      OLPTR=NSTORE(NVEC)+1
      OLINE(NVEC)(OLPTR:OLPTR)='('
      OLPTR=OLPTR+1
      DO IPTR=1,HWLEN
         IF(HWSTR(IPTR:IPTR).EQ.CHAR(SYMSHFT) .OR. 
     .    HWSTR(IPTR:IPTR).EQ.CHAR(ROMSHFT)) THEN
            IF(IPTR.EQ.HWLEN) GOTO 188
            IF(HWSTR(IPTR:IPTR).EQ.CHAR(SYMSHFT)) THEN
               HWFONT='sf'
            ELSE
               HWFONT='rf'
            ENDIF
            IF(HWFONT.NE.CURFONT) THEN
               OLINE(NVEC)(OLPTR:OLPTR+4)=')show'
               OLPTR=OLPTR+5
               NSTORE(NVEC)=OLPTR-1

C try to keep the vector list from wrapping when typing the file
               LLN = 0
               DO NC=1,NVEC
                  LLN = LLN+NSTORE(NC)
               ENDDO
               IF(LLN.GT.72)THEN
C     WRITE-BUFFER
                  ASSIGN 99989 TO I99993
                  GO TO 99993
               ENDIF
99989          CONTINUE

               NVEC=NVEC+1
               OLINE(NVEC)=FVEC(:FLEN)//HWFONT
               OLPTR=1+FLEN+2
               CURFONT=HWFONT
               OLINE(NVEC)(OLPTR:OLPTR)='('
               OLPTR=OLPTR+1
            ENDIF
            GOTO 188
         ENDIF
         IF(OLPTR.GT.65) THEN
            NSTORE(NVEC)=OLPTR-1
            NVEC=NVEC+1
            OLPTR=1
         ENDIF
         IF(HWSTR(IPTR:IPTR).EQ.'(') THEN
C the following line is compiler-dependent
            OLINE(NVEC)(OLPTR:OLPTR+3)='\050'
            OLPTR=OLPTR+4
         ELSE IF(HWSTR(IPTR:IPTR).EQ.')') THEN
C the following line is compiler-dependent
            OLINE(NVEC)(OLPTR:OLPTR+3)='\051'
            OLPTR=OLPTR+4
         ELSE
            OLINE(NVEC)(OLPTR:OLPTR)=HWSTR(IPTR:IPTR)
            OLPTR=OLPTR+1
         ENDIF
  188 END DO
      IF(OLPTR.GT.65) THEN
         NSTORE(NVEC)=OLPTR-1
         NVEC=NVEC+1
         OLPTR=1
      ENDIF
      OLINE(NVEC)(OLPTR:OLPTR+4)=')show'
      NSTORE(NVEC)=OLPTR+4
      RETURN
         
C----------------------------------------
C
C        TO WRITE-BUFFER
99993 CONTINUE
C send buffer to output file
      IF(NVEC.GT.0 .AND. OLINE(1).NE.' ')THEN
C element 0 is just for safety - don't try to plot it
        LLN = 1
        LINE2WRITE=' '
C if this line starts with a 'pu' and the previous line
C didn't end with a 'ps', then emit a 'ps' for good luck
        NOPS = LLNP.GT.1
        IF(NOPS)THEN
          NOPS =
     1    (OLINE(1)(NSTORE(1)-2:NSTORE(1)).EQ.' pu') .AND.
     1    (PLINE2WRITE(LLNP-2:LLNP) .NE. ' ps')
          ENDIF
        IF(NOPS)THEN
          WRITE(UNIT=LQUNIT,FMT='(A)') ' ps'
          NVPS = 0
          ENDIF
        DO I=1,NVEC
          LINE2WRITE =
     1    LINE2WRITE(:LLN) // OLINE(I)(:NSTORE(I))
          LLN = LLN + NSTORE(I)
          ENDDO
C eliminate all extra blanks, ".0" and the like
        CALL APLOT_UNPAD(LINE2WRITE,LLN)
C don't emit any line with only pu and ps in it but no pd
C (fixing yet another optimization bug. sigh.)
        IF(LLN.LT.17)THEN
          IPU = INDEX (LINE2WRITE(:LLN),'pu')
          IPS = INDEX (LINE2WRITE(:LLN),'ps')
          IPD = INDEX (LINE2WRITE(:LLN),'pd')
          IF (IPD.EQ.0)   IPD = INDEX (LINE2WRITE(:LLN),'nx')
          IF (IPD.EQ.0)   IPD = INDEX (LINE2WRITE(:LLN),'ny')
          IF(.NOT.((IPU.NE.0) .AND. (IPS.NE.0) .AND. (IPD.EQ.0)))THEN
            WRITE(UNIT=LQUNIT,FMT='(A)') LINE2WRITE(:LLN)
            PLINE2WRITE = LINE2WRITE
            LLNP = LLN
            ENDIF
        ELSE
C not so minor optimization:
C don't write more than two duplicate lines in a row
C (eg: scatter plots)
          IF(LINE2WRITE .EQ. PLINE2WRITE)THEN
            NDUP = NDUP + 1
          ELSE
            NDUP = 0
            ENDIF
          IF(NDUP.LT.2)THEN
            WRITE(UNIT=LQUNIT,FMT='(A)') LINE2WRITE(:LLN)
            ENDIF
          PLINE2WRITE = LINE2WRITE
          LLNP = LLN
          ENDIF
        LINE2WRITE = ' '
        DO I=0,NVEC
          OLINE(I)=' '
          NSTORE(I)=0
          ENDDO
        NVEC=0
        ENDIF
      GO TO I99993, (99988,99989,99990,99991,99992)
      END
      SUBROUTINE APLOT_UNPAD(STR,LENS)
C Remove extra padding from string - fixes string "in place"
      CHARACTER*(*) STR
      IF(LENS.GT.0)THEN
        LENS = MIN(LEN(STR),LENS)
      ELSE
        LENS = LEN(STR)
        ENDIF
C Remove extra tabs and spaces
C 0. Convert NULLS to spaces
      IPT = INDEX(STR(:LENS),CHAR(0))
      DO WHILE(IPT.NE.0)
        STR(IPT:IPT)=' '
        IPT = INDEX(STR,CHAR(0))
        ENDDO
C 1. Convert tabs to spaces
      IPT = INDEX(STR(:LENS),CHAR(9))
      DO WHILE(IPT.NE.0)
        STR(IPT:IPT)=' '
        IPT = INDEX(STR,CHAR(9))
        ENDDO
C 2. remove leading spaces
      DO WHILE(STR(1:1).EQ.' ' .AND. LENS.GT.1)
        STR(1:LENS)=STR(2:LENS)
        LENS=LENS-1
        ENDDO
      IF(STR(1:1).EQ.' ' .AND. LENS.LE.1)THEN
        LENS = 0
        RETURN
        ENDIF
C 2A. remove trailing spaces
      DO WHILE(STR(LENS:LENS).EQ.' ' .AND. LENS.GT.0)
        LENS=LENS-1
        ENDDO
cC 3. convert multiple spaces to a single space
c      IPS=1
c      DO WHILE(IPS.GT.0 .AND. LENS.GT.0)
c        IPS=INDEX(STR(:LENS),'  ')
c        IF(IPS.GT.0)THEN
c          IPN=IPS
cC Find the next non-space character (Must be one)
c          DO WHILE(STR(IPN:IPN).EQ.' ' .AND. IPN.LT.LENS)
c            IPN=IPN+1
c            ENDDO
cC Concatenate end of string with part including leading space
c          STR(:LENS)=STR(1:IPS)//STR(IPN:LENS)
c          LENS= IPS + (LENS - IPN +1)
c          IF(LENS.LT.0) LENS = 0
c          ENDIF
c        ENDDO
      RETURN
      END
