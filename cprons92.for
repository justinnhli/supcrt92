*** cprons92 - Read a sequential-access sprons.dat file and 
***            generate its direct-access dprons.dat equivalent.
***
*********************************************************************
***
*** Author:     James W. Johnson
***             Earth Sciences Department, L-219
***             Lawrence Livermore National Laboratory
***             Livermore, CA 94550
***             johnson@s05.es.llnl.gov
***
*** Abandoned:  8 November 1991
***
*********************************************************************

      PROGRAM cprons

      CALL opener
      CALL trnsfr

      END

*********************************************************************

*** consts - Constants.

      BLOCK DATA consts

      INTEGER rterm, wterm, saf, daf

      COMMON /io/     rterm, wterm, saf, daf

      SAVE

      DATA rterm, wterm, saf, daf
     1   /   5,     6,    43,  44  /

      END

*********************************************************************

*** opener - Open an existing sequential-access sprons.dat file (saf)
***          and a new direct-access dprons.dat file (daf).

      SUBROUTINE opener

      LOGICAL openf
      CHARACTER*20 fname
      INTEGER rterm, wterm, saf, daf

      COMMON /io/ rterm, wterm, saf, daf

      SAVE


      CALL banner

  1   WRITE(wterm,10)
 10   FORMAT(/,' specify name of INPUT sequential-access'
     1        ,' THERMODYNAMIC DATABASE: ',/)
      READ(rterm,20) fname
 20   FORMAT(a20)
      IF (.NOT. openf(wterm,saf,fname,1,1,1,90)) GO TO 1

  2   WRITE(wterm,30)
 30   FORMAT(/,' specify name of OUTPUT direct-access'
     1        ,' THERMODYNAMIC DATABASE: ',/)
      READ(rterm,20) fname
      IF (.NOT. openf(wterm,daf,fname,2,2,1,90)) GO TO 2


      END 

*******************************************************************

*** banner - Write program banner to the terminal screen.

      SUBROUTINE banner

      INTEGER rterm, wterm, saf, daf

      COMMON /io/ rterm, wterm, saf, daf

      SAVE


      WRITE(wterm,10)
 10   FORMAT(/,5x,' Welcome to CPRONS92 Version 1.0',
     1       /,5x,' Author:    James W. Johnson',
     2       /,5x,' Abandoned: 8 November 1991',/)

      END 

*******************************************************************

*** trnsfr - Transfer data from saf to daf.

      SUBROUTINE trnsfr

      INTEGER rterm, wterm, saf, daf,
     1        rec1m1, rec1m2, rec1m3, rec1m4, rec1g, rec1a

      COMMON /io/ rterm, wterm, saf, daf

      SAVE


      CALL tail(nmin1,nmin2,nmin3,nmin4,ngas,naqs,
     1          rec1m1,rec1m2,rec1m3,rec1m4,rec1g,rec1a)

      WRITE(wterm,*) ' starting min1 transfer '
      CALL mgtran(0,nmin1,rec1m1)

      WRITE(wterm,*) ' starting min2 transfer '
      CALL mgtran(1,nmin2,rec1m2)

      WRITE(wterm,*) ' starting min3 transfer '
      CALL mgtran(2,nmin3,rec1m3)

      WRITE(wterm,*) ' starting min4 transfer '
      CALL mgtran(3,nmin4,rec1m4)

      WRITE(wterm,*) ' starting gas transfer '
      CALL mgtran(0,ngas,rec1g)

      WRITE(wterm,*) ' starting aqueous species transfer '
      CALL aqtran(naqs,rec1a)

      END

******************************************************************
  
*** tail - Read nmin1..4, ngas, naqs, lskip from tail of saf;
***        calculate corresponding rec1m1..4, rec1g, rec1a;  
***        transfer these bookkeeping parameters to head of daf;
***        read-skip to top of min1 block of saf.

      SUBROUTINE tail(nmin1,nmin2,nmin3,nmin4,ngas,naqs,
     1                rec1m1,rec1m2,rec1m3,rec1m4,rec1g,rec1a)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NBACK = 12)
***** NBACK = # lines to backspace before reading nmin[1..4],ngas,naqs.
*****         NOTE!!!  On some installations, NBACK must be 
*****                  incremented to 11

      INTEGER rterm, wterm, saf, daf,
     1        rec1m1, rec1m2, rec1m3, rec1m4, rec1g, rec1a

      COMMON /io/ rterm, wterm, saf, daf

      SAVE


*** go to EOF and backpeddle to read current species totals and lskip

      REWIND(saf)
      DO 31 i = 1,10000
           READ(saf,32,END=3)
 32        FORMAT()
 31        CONTINUE

  3   DO 33 i = 1,NBACK
 33        BACKSPACE(saf)

      READ(saf,*) nmin1
      print *, nmin1
      READ(saf,*) nmin2
      print *, nmin2
      READ(saf,*) nmin3
      print *, nmin3
      READ(saf,*) nmin4
      print *, nmin4
      READ(saf,*) ngas
      print *, ngas
      READ(saf,*) naqs
      print *, naqs
      READ(saf,32)
      READ(saf,32)
      READ(saf,32)
      READ(saf,*) lskip
      print *, lskip

*** return to TOF and skip over comment lines,
*** i.e. move forward to top of min1 block

      REWIND(saf)
      DO 40  l = 1,lskip
           READ(saf,50) 
 50        FORMAT()
 40        CONTINUE

*** calculate rec1m1..4, rec1g, rec1a; 
*** write these parameters to line 2 of daf.

      rec1m1 = 3
      rec1m2 = rec1m1 + 6*nmin1
      rec1m3 = rec1m2 + 7*nmin2
      rec1m4 = rec1m3 + 8*nmin3
      rec1g  = rec1m4 + 9*nmin4
      rec1a  = rec1g  + 6*ngas

      WRITE(daf,60,REC=1) nmin1,  nmin2,  nmin3,  nmin4,  ngas,  naqs
      WRITE(daf,60,REC=2) rec1m1, rec1m2, rec1m3, rec1m4, rec1g, rec1a
 60   FORMAT(6(1x,i4))

      END

******************************************************************
      
*** mgtran - Read (from saf) data for current block of nmg mineral
***          or gas species, each of which contains ntran phase 
***          transitions; write these data to daf, 
***          starting at REC=rec1. 

      SUBROUTINE mgtran(ntran,nmg,rec1)
 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (NTPLUS = 4)

      INTEGER      rterm, wterm, saf, daf, rec1, r
      CHARACTER*9  date
      CHARACTER*12 ref
      CHARACTER*15 abrev
      CHARACTER*20 name
      CHARACTER*30 scform
      CHARACTER*40 ecform

      DOUBLE PRECISION  a(NTPLUS), b(NTPLUS), c(NTPLUS),
     1                  Ttr(NTPLUS), Htr(NTPLUS),
     2                  Vtr(NTPLUS), tslope(NTPLUS)


      COMMON /io/ rterm, wterm, saf, daf

      SAVE


***** skip global 3-line header for min1..4 or gas data block

      DO 1 i = 1,3
           READ(saf,2)
  2        FORMAT()
  1        CONTINUE

***** transfer data for all nmg individual species blocks 

      DO 5 i = 1,nmg

*****      set r = first record for current species block

           r = rec1 + (i-1)*(6 + ntran)
  
*****      read data from saf 

           READ(saf,10) name, scform
 10        FORMAT(1x,a20,a30)
           READ(saf,20) abrev, ecform
 20        FORMAT(1x,a15,5x,a40)
           READ(saf,30) ref, date
 30        FORMAT(1x,a12,8x,a9)
           READ(saf,40)  Gf, Hf, Sref, Vref
 40        FORMAT(4x,2(2x,f12.1),2(2x,f8.3))

           DO 45 j = 1,ntran                 
                READ(saf,50) a(j), b(j),
     1          c(j), Ttr(j), Htr(j), Vtr(j), tslope(j)
 50             FORMAT(4x,3(2x,f12.6),2x,f7.2,2x,f8.1,2(2x,f10.3))
 45             CONTINUE
 
           j = ntran + 1
           READ(saf,60) a(j), b(j), c(j)
 60        FORMAT(4x,3(2x,f12.6))
           READ(saf,70) Tmax
 70        FORMAT(8x,f7.2)

*****      write data to daf *****

           WRITE(daf,10,REC=r)    name, scform
           WRITE(daf,20,REC=r+1)  abrev, ecform
           WRITE(daf,30,REC=r+2)  ref, date 
           WRITE(daf,40,REC=r+3)  Gf, Hf, Sref, Vref

           DO 65 j = 1,ntran                 
                WRITE(daf,50,REC=r+3+j) a(j), b(j),
     1          c(j), Ttr(j), Htr(j), Vtr(j), tslope(j)
 65             CONTINUE
           j = ntran + 1
           WRITE(daf,60,REC=r+4+ntran) a(j), b(j), c(j)
           WRITE(daf,70,REC=r+5+ntran) Tmax
 
  5        CONTINUE

      END

************************************************************

*** aqtran - Read (from saf) data for aqueous species block; 
***          write these data to daf, starting at REC=rec1.

      SUBROUTINE aqtran(naqs,rec1)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER      rterm, wterm, saf, daf, rec1, r
      CHARACTER*9  date
      CHARACTER*12 ref
      CHARACTER*15 abrev
      CHARACTER*20 name
      CHARACTER*30 scform
      CHARACTER*40 ecform

      COMMON /io/ rterm, wterm, saf, daf

      SAVE


***** skip global 3-line header for aqueous species data block

      DO 1 i = 1,3
           READ(saf,2)
  2        FORMAT()
  1        CONTINUE

***** transfer data for all naqs individual aqueous species blocks 

      DO 5 i = 1,naqs

*****      set r = first record for current species block

           r = rec1 + 6*(i-1)
  
*****      read data from saf

           READ(saf,10) name, scform
 10        FORMAT(1x,a20,a30)
           READ(saf,20) abrev, ecform
 20        FORMAT(1x,a15,5x,a40)
           READ(saf,30) ref, date
 30        FORMAT(1x,a12,8x,a9)
           READ(saf,40)  Gf, Hf, Sref
 40        FORMAT(4x,2(2x,f10.0),4x,f8.3)
           READ(saf,50)  a1, a2, a3, a4
 50        FORMAT(4x,4(2x,f8.4,2x))
           READ(saf,60)  c1, c2, omega, charge
 60        FORMAT(4x,3(2x,f8.4,2x),9x,f3.0)

*****      write data to daf

           WRITE(daf,10,REC=r)    name, scform
           WRITE(daf,20,REC=r+1)  abrev, ecform
           WRITE(daf,30,REC=r+2)  ref, date
           WRITE(daf,40,REC=r+3)  Gf, Hf, Sref
           WRITE(daf,50,REC=r+4)  a1, a2, a3, a4
           WRITE(daf,60,REC=r+5)  c1, c2, omega, charge

  5        CONTINUE

      END

***************************************************************
***
*** openf -  Returns .TRUE. and opens the file specified by 
***          fname, fstat, facces, fform, and frecl  
***          if this file exists and is accessible; otherwise, 
***          returns .FALSE. and prints an appropriate error 
***          message to the device specified by iterm.
***
*** Author:     James W. Johnson
***
*** Abandoned:  8 October 1990
***
***************************************************************

      LOGICAL FUNCTION openf(iterm,iunit,fname,istat,
     1                       iacces,iform,irecl)      

      CHARACTER*11  fform(2)
      CHARACTER*10  facces(2)
      CHARACTER*20  fname
      CHARACTER*3   fstat(2)

      SAVE

      DATA fform  / 'FORMATTED  ',  'UNFORMATTED' /
      DATA facces / 'SEQUENTIAL',   'DIRECT    '  /
      DATA fstat  / 'OLD',          'NEW'         /


      openf = .FALSE.
      
      IF ((iacces .LT. 1) .OR. (iacces .GT. 2) .OR.
     1    (iform  .LT. 1) .OR. (iform  .GT. 2) .OR.
     2    (istat  .LT. 1) .OR. (istat  .GT. 2)) GO TO 10

      IF (iacces .EQ. 1) THEN
           OPEN(UNIT=iunit,FILE=fname,ACCESS=facces(iacces),
     1          FORM=fform(iform),STATUS=fstat(istat),ERR=10)
           openf = .TRUE.
           RETURN
      ELSE
           OPEN(UNIT=iunit,FILE=fname,ACCESS=facces(iacces),
     1          FORM=fform(iform),STATUS=fstat(istat),RECL=irecl,
     2          ERR=10)
           openf = .TRUE.
           RETURN
      END IF

 10   WRITE(iterm,20)
 20   FORMAT(/,' nonexistant file or invalid specifications',
     1         ' ... try again',/)
      RETURN

      END
