*** mprons92 - Open a sequential-access sprons.dat file, 
***            modify according to user specifications,
***            then generate the updated file.
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

      PROGRAM mprons

      CALL opener(nmin1,nmin2,nmin3,nmin4,ngas,naqs,lskip)
      CALL mingas(nmin1,'min',0)
      CALL mingas(nmin2,'min',1)
      CALL mingas(nmin3,'min',2)
      CALL mingas(nmin4,'min',3)
      CALL mingas(ngas,'gas',0)
      CALL aqueus(naqs,'aqs',0)
      CALL summry(nmin1,nmin2,nmin3,nmin4,ngas,naqs,lskip)

      END

********************************************************************

*** consts - Constants.

      BLOCK DATA consts

      INTEGER rterm, wterm, safold, safnew

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE

      DATA rterm, wterm, safold, safnew
     1   /   5,     6,     43,     44  /

      END

**********************************************************************

*** opener - Open existing and new sequential-access sprons.dat files;
***          read header information from old file; read/write to 
***          top of min1 data.

      SUBROUTINE opener(nmin1,nmin2,nmin3,nmin4,ngas,naqs,lskip)

      PARAMETER (NBACK = 12)
***** NBACK = # lines to backspace before reading nmin[1..4],ngas,naqs.
*****         NOTE!!!  On some installations, NBACK must be 
*****                  incremented to 11

      LOGICAL openf
      CHARACTER*20 fname
      CHARACTER*80 line
      INTEGER rterm, wterm, safold, safnew

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


      CALL banner

  1   WRITE(wterm,10)
 10   FORMAT(/,' specify name of an EXISTING sequential-access'
     1        ,' THERMODYNAMIC DATABASE:',/)
      READ(rterm,20) fname
 20   FORMAT(a20)
      IF (.NOT. openf(wterm,safold,fname,1,1,1,90)) GO TO 1

  2   WRITE(wterm,30)
 30   FORMAT(/,' specify name of a NEW sequential-access'
     1        ,' THERMODYNAMIC DATABASE:',/)
      READ(rterm,20) fname
      IF (.NOT. openf(wterm,safnew,fname,2,1,1,90)) GO TO 2

*** go to EOF and backpeddle to read current species totals and lskip

      REWIND(safold)
      DO 31 i = 1,20000
           READ(safold,32,END=3)
 32        FORMAT()
 31        CONTINUE

  3   DO 33 i = 1,NBACK 
 33        BACKSPACE(safold)

      READ(safold,*) nmin1
      READ(safold,*) nmin2
      READ(safold,*) nmin3
      READ(safold,*) nmin4
      READ(safold,*) ngas
      READ(safold,*) naqs
      READ(safold,32)
      READ(safold,32)
      READ(safold,32)
      READ(safold,*) lskip

*** return to TOF and transfer comment lines,
*** i.e. move forward to top of min1 block

      REWIND(safold)
      DO 40  l = 1,lskip
           READ(safold,50) line
 50        FORMAT(a80)
           WRITE(safnew,50) line
 40        CONTINUE

      END 

*******************************************************************

*** banner - Write program banner to the terminal screen.

      SUBROUTINE banner

      INTEGER rterm, wterm, safold, safnew

      COMMON /io/ rterm, wterm, safold, safnew
     
      SAVE


      WRITE(wterm,10)
 10   FORMAT(/,5x,' Welcome to MPRONS92 Version 1.0',
     1       /,5x,' Author:    James W. Johnson',
     2       /,5x,' Abandoned: 8 November 1991',/)

      END 

*******************************************************************

*** mingas - Transfer nmg mineral or gas species from the old to
***          new sprons.dat file; if requested, revise the data for,
***          delete, or add the specified species; if necessary, 
***          update nmg and/or realphabetize.

      SUBROUTINE mingas(nmg,stype,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXCHG = 100, MAXMG = 200, NTPLUS = 4)

      LOGICAL match, rename, realph

      CHARACTER*3  stype
      CHARACTER*9  date(MAXMG)
      CHARACTER*12 ref(MAXMG)
      CHARACTER*15 abrev(MAXMG)
      CHARACTER*20 name(MAXMG), alist(MAXCHG), 
     1             dlist(MAXCHG), mlist(MAXCHG)
      CHARACTER*30 scform(MAXMG)
      CHARACTER*40 ecform(MAXMG)

      INTEGER rterm, wterm, safold, safnew, order(MAXMG)

      DOUBLE PRECISION Gf(MAXMG), Hf(MAXMG), Sref(MAXMG), Vref(MAXMG),
     1                 a(MAXMG,NTPLUS), b(MAXMG,NTPLUS), 
     2                 c(MAXMG,NTPLUS), Ttr(MAXMG,NTPLUS), 
     3                 Htr(MAXMG,NTPLUS), Vtr(MAXMG,NTPLUS),
     4                 Tslope(MAXMG,NTPLUS), Tmax(MAXMG)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


      CALL header(stype,ntran)
      CALL getchg(stype,ntran,nadd,ndel,nmod,alist,dlist,mlist)

      ngone  = 0
      ndone  = 0
      realph = .FALSE.

      DO 10 i = 1,nmg
           j = i-ngone
           CALL readmg(j,name,scform,abrev,ecform,ref,date,Gf,Hf,Sref,
     1                 Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran)
           IF (match(ndel-ngone,ndel,name(j),dlist)) THEN
***             delete the current mineral or gas 
                ngone = ngone + 1
                GO TO 10
           END IF

           IF (match(nmod-ndone,nmod,name(j),mlist)) THEN
***             modify data for the current mineral or gas 
                ndone = ndone + 1
                CALL modmg(j,name,scform,abrev,ecform,ref,date,Gf,
     1               Hf,Sref,Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran,
     2               rename)
		IF (rename) THEN
		     realph = .TRUE.
                END IF
           END IF

 10        CONTINUE

      IF (ndel .NE. ngone) THEN
***        one or more species to be deleted was not found
           WRITE(wterm,20) ndel-ngone
 20        FORMAT(/,i5,' of the species to be deleted were not found',/)
      END IF

      IF (nmod .NE. ndone) THEN
***        one or more species to be modified was not found
           WRITE(wterm,30) nmod-ndone
 30        FORMAT(/,i5,' of the species to be modified were not found',
     1            /)
      END IF

      nmg = nmg - ngone

      IF (nadd .GT. 0) THEN
	   realph = .TRUE.
           DO 40 i = 1,nadd
***             obtain data for the new mineral or gas
                nmg = nmg + 1
                name(nmg) = alist(i)
                CALL addmg(nmg,name,scform,abrev,ecform,ref,date,
     1                     Gf,Hf,Sref,Vref,a,b,c,Ttr,Htr,Vtr,Tslope,
     2                     Tmax,ntran)
 40             CONTINUE
      END IF

*** initialize order
      
      DO 50 i = 1,nmg
           order(i) = i
 50        CONTINUE

*** re-alphabetize if necessary 

      IF (realph) CALL alpha(nmg,name,order)

*** transfer data 

      DO 60 i = 1,nmg
           CALL tranmg(order(i),name,scform,abrev,ecform,ref,date,Gf,Hf,
     1                 Sref,Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran)
 60        CONTINUE

      END

*********************************************************************

*** alpha - Use a straight-selection sort to alphabetize the n-length
***         list of names; the index sequence representing their 
***         alphabetical order is returned in order.

      SUBROUTINE alpha(n,name,order)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMG = 200)

      INTEGER order(MAXMG)

      CHARACTER*20 name(MAXMG)

      SAVE


      DO 10 i = n,2,-1
	   iend = 1
           DO 20 j = 2,i
		IF (name(order(iend)) .LT. name(order(j))) iend = j
 20             CONTINUE
	   itemp = order(i)
	   order(i) = order(iend)
	   order(iend) = itemp
 10        CONTINUE

      END 

***********************************************************************

*** modmg - modify data for the mineral of gas refereced by name

      SUBROUTINE modmg(j,name,scform,abrev,ecform,ref,date,Gf,Hf,
     1                 Sref,Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran,
     2                 rename)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMG = 200, NTPLUS = 4)

      LOGICAL modify, rename
      INTEGER rterm, wterm, safold, safnew

      CHARACTER*9  date(MAXMG)
      CHARACTER*12 ref(MAXMG)
      CHARACTER*15 abrev(MAXMG)
      CHARACTER*20 name(MAXMG)
      CHARACTER*30 scform(MAXMG)
      CHARACTER*40 ecform(MAXMG)

      DOUBLE PRECISION Gf(MAXMG), Hf(MAXMG), Sref(MAXMG), Vref(MAXMG),
     1                 a(MAXMG,NTPLUS), b(MAXMG,NTPLUS), 
     2                 c(MAXMG,NTPLUS), Ttr(MAXMG,NTPLUS), 
     3                 Htr(MAXMG,NTPLUS), Vtr(MAXMG,NTPLUS),
     4                 Tslope(MAXMG,NTPLUS), Tmax(MAXMG)
    
      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


***** modify data 

      WRITE(wterm,1) name(j)
  1   FORMAT(/,' modify data for ',a20,/)

***** specify / [update] name

      WRITE(wterm,10) name(j)
 10   FORMAT(/,' name: ',a20,/)
      IF (.NOT. modify()) THEN
	   rename = .FALSE.
      ELSE
	   rename = .TRUE.
           WRITE(wterm,20) 
 20        FORMAT(/,' specify new name: ',/)
           READ(rterm,30) name(j)
 30        FORMAT(a20)
      END IF

***** specify / [update] structural chemical formula

      WRITE(wterm,11) scform(j)
 11   FORMAT(/,' structural chemical formula: ',a30,/)
      IF (modify()) THEN
           WRITE(wterm,21) 
 21        FORMAT(/,' specify new structural chemical formula: ',/)
           READ(rterm,31) scform(j)
 31        FORMAT(a30)
      END IF

***** specify / [update] abbreviation

      WRITE(wterm,12) abrev(j)
 12   FORMAT(/,' abbreviation: ',a15,/)
      IF (modify()) THEN
           WRITE(wterm,22) 
 22        FORMAT(/,' specify new abbreviation: ',/)
           READ(rterm,32) abrev(j)
 32        FORMAT(a15)
      END IF

***** specify / [update] elemental chemical formula

      WRITE(wterm,13) ecform(j)
 13   FORMAT(/,' elemental chemical formula: ',a40,/)
      IF (modify()) THEN
           WRITE(wterm,23) 
 23        FORMAT(/,' specify new elemental chemical formula: ',/)
           READ(rterm,33) ecform(j)
 33        FORMAT(a40)
      END IF

***** specify / [update] reference

      WRITE(wterm,14) ref(j)
 14   FORMAT(/,' reference entry: ',a12,/)
      IF (modify()) THEN
           WRITE(wterm,24) 
 24        FORMAT(/,' specify new reference entry: ',/)
           READ(rterm,34) ref(j)
 34        FORMAT(a12)
      END IF

***** specify / update date

      WRITE(wterm,15) date(j)
 15   FORMAT(/,' date of last update: ',a9,/)
      WRITE(wterm,25) 
 25   FORMAT(/,' current date: ',/)
      READ(rterm,35) date(j)
 35   FORMAT(a9)

***** specify / [update] Standard molal Gf 

      WRITE(wterm,16) Gf(j)
 16   FORMAT(/,' standard molal Gibbs free energy of formation'
     1        ,' (cal/mol)',
     2       /,' at 25 C, 1 bar [NOTE: 999999 = NULL value]:  ',f12.1,/)
      IF (modify()) THEN
           WRITE(wterm,26) 
 26        FORMAT(/,' specify new Gf: ',/)
           READ(rterm,*) Gf(j)
      END IF

***** specify / [update] Standard molal Hf 

      WRITE(wterm,17) Hf(j)
 17   FORMAT(/,' standard molal enthalpy of formation (cal/mol)',
     1       /,' at 25 C, 1 bar (NOTE: 999999 = NULL value):  ',f12.1,/)
      IF (modify()) THEN
           WRITE(wterm,27) 
 27        FORMAT(/,' specify new Hf: ',/)
           READ(rterm,*) Hf(j)
      END IF

***** specify / [update] Standard molal S 

      WRITE(wterm,18) Sref(j)
 18   FORMAT(/,' standard molal entropy (cal/mol/K) at 25 C, 1 bar: '
     1        ,f8.3,/)
      IF (modify()) THEN
           WRITE(wterm,28) 
 28        FORMAT(/,' specify new SPrTr: ',/)
           READ(rterm,*) Sref(j)
      END IF

***** specify / [update] Standard molal V

      WRITE(wterm,19) Vref(j)
 19   FORMAT(/,' standard molal volume (cc/mol) at 25 C, 1 bar: '
     1        ,f8.3,/)
      IF (modify()) THEN
           WRITE(wterm,29) 
 29        FORMAT(/,' specify new VPrTr: ',/)
           READ(rterm,*) Vref(j)
      END IF

***** specify / [update] Maier-Kelly Cp coeffs

      IF (ntran .EQ. 0) THEN
           WRITE(wterm,40) a(j,1)
 40        FORMAT(/,' Maier-Kelly Cp coeff  a(10**0): ',f12.6,/)
           IF (modify()) THEN
                WRITE(wterm,50) 
 50             FORMAT(/,' specify new  a(10**0): ',/)
                READ(rterm,*) a(j,1)
           END IF  
           WRITE(wterm,42) b(j,1)
 42        FORMAT(/,' Maier-Kelly Cp coeff  b(10**3): ',f12.6,/)
           IF (modify()) THEN
                WRITE(wterm,52) 
 52             FORMAT(/,' specify new  b(10**3): ',/)
                READ(rterm,*) b(j,1)
           END IF  
           WRITE(wterm,43) c(j,1)
 43        FORMAT(/,' Maier-Kelly Cp coeff  c(10**-5): ',f12.6,/)
           IF (modify()) THEN
                WRITE(wterm,53) 
 53             FORMAT(/,' specify new  c(10**-5): ',/)
                READ(rterm,*) c(j,1)
           END IF  
           WRITE(wterm,44) Tmax(j)
 44        FORMAT(/,' maximum temperature (K): ',f7.2,/)
           IF (modify()) THEN
                WRITE(wterm,54) 
 54             FORMAT(/,' specify new maximum temperature (K): ',/)
                READ(rterm,*) Tmax(j)
           END IF  
           RETURN
      END IF

      DO 111 k = 1,ntran
           WRITE(wterm,65) k
 65        FORMAT(/,15x,' data for pre-transition ',i1,/)
           WRITE(wterm,40) a(j,k)
           IF (modify()) THEN
                WRITE(wterm,50) 
                READ(rterm,*) a(j,k)
           END IF  
           WRITE(wterm,42) b(j,k)
           IF (modify()) THEN
                WRITE(wterm,52) 
                READ(rterm,*) b(j,k)
           END IF  
           WRITE(wterm,43) c(j,k)
           IF (modify()) THEN
                WRITE(wterm,53) 
                READ(rterm,*) c(j,k)
           END IF  
           WRITE(wterm,45) k, Ttr(j,k)
 45        FORMAT(/,' temperature (K) of transition ',i1,': ',f7.2,/)
           IF (modify()) THEN
                WRITE(wterm,56) 
 56             FORMAT(/,' specify new transition temperature (K): ',/)
                READ(rterm,*) Ttr(j,k)
           END IF  
           WRITE(wterm,300) k, Htr(j,k)
 300       FORMAT(/,' enthalpy of transition ',i1,' (cal/mol): ',f8.1,
     1            /,' [NOTE: 999999 = NULL value]',/)
           IF (modify()) THEN
                WRITE(wterm,301) k
 301            FORMAT(/,' specify new enthalpy of transition ',i1
     1                  ,' (cal/mol): ',/)
                READ(rterm,*) Htr(j,k)
           END IF  
           WRITE(wterm,302) k, Vtr(j,k)
 302       FORMAT(/,' volume of transition ',i1,' (cc/mol): ',f10.3,
     1            /,' [NOTE: 999999 = NULL value]',/)
           IF (modify()) THEN
                WRITE(wterm,303) k
 303            FORMAT(/,' specify new volume of transition ',i1
     1                  ,' (cc/mol): ',/)
                READ(rterm,*) Vtr(j,k)
           END IF  
           WRITE(wterm,304) k, Tslope(j,k)
 304       FORMAT(/,' dP/dT of transition ',i1,' (bars/K): ',f10.3,
     1            /,' [NOTE: 999999 = NULL value]',/)
           IF (modify()) THEN
                WRITE(wterm,305) k
 305            FORMAT(/,' specify new (dP/dT) of transition ',i1
     1                  ,' (bar/K): ',/)
                READ(rterm,*) Tslope(j,k)
           END IF  
 111       CONTINUE

      WRITE(wterm,95) ntran
 95   FORMAT(/,15x,' data for post-transition ',i1,/)
      WRITE(wterm,40) a(j,ntran+1)
      IF (modify()) THEN
           WRITE(wterm,50) 
           READ(rterm,*) a(j,ntran+1)
      END IF  
      WRITE(wterm,42) b(j,ntran+1)
      IF (modify()) THEN
           WRITE(wterm,52) 
           READ(rterm,*) b(j,ntran+1)
      END IF  
      WRITE(wterm,43) c(j,ntran+1)
      IF (modify()) THEN
           WRITE(wterm,53) 
           READ(rterm,*) c(j,ntran+1)
      END IF  
      WRITE(wterm,44) Tmax(j)
      IF (modify()) THEN
           WRITE(wterm,54) 
           READ(rterm,*) Tmax(j)
      END IF  

      END 

*********************************************************************

*** modify - Returns .TRUE. if the user wishes to modify the current
***          value of a parameter; otherwise, returns .FALSE.

      LOGICAL FUNCTION modify()

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      CHARACTER*1   ans
      INTEGER       rterm, wterm, safold, safnew

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


  1   WRITE(wterm,10) 
 10   FORMAT(' Would you like to update this value?  (y/n)',/)

      READ(rterm,15) ans
 15   FORMAT(a1)
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) THEN
           GO TO 1
      ELSE  
           modify = ((ans .EQ. 'y') .OR. (ans .EQ. 'Y')) 
      END IF
      
      RETURN
      END

*********************************************************************

*** addmg - Add the mineral or gas referenced by name, which undergoes
***         ntran phase transitions, to the database.

      SUBROUTINE addmg(j,name,scform,abrev,ecform,ref,date,Gf,Hf,Sref,
     1                 Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
    
      PARAMETER (MAXMG = 200, NTPLUS = 4)

      CHARACTER*9  date(MAXMG)
      CHARACTER*12 ref(MAXMG)
      CHARACTER*15 abrev(MAXMG)
      CHARACTER*20 name(MAXMG)
      CHARACTER*30 scform(MAXMG)
      CHARACTER*40 ecform(MAXMG)

      INTEGER rterm, wterm, safold, safnew

      DOUBLE PRECISION Gf(MAXMG), Hf(MAXMG), Sref(MAXMG), Vref(MAXMG),
     1                 a(MAXMG,NTPLUS), b(MAXMG,NTPLUS), 
     2                 c(MAXMG,NTPLUS), Ttr(MAXMG,NTPLUS), 
     3                 Htr(MAXMG,NTPLUS), Vtr(MAXMG,NTPLUS),
     4                 Tslope(MAXMG,NTPLUS), Tmax(MAXMG)

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


***** prompt for / read requisite data for new species "name"

      WRITE(wterm,10) name(j)
 10   FORMAT(/,' specify requisite data for ',a20,/)

      WRITE(wterm,15) 
 15   FORMAT(/,' abbreviation [e.g., Qtz]:',/)
      READ(rterm,20) abrev(j)
 20   FORMAT(a15)

      WRITE(wterm,16) 
 16   FORMAT(/,' structural chemical formula [e.g., SiO2]:',/)
      READ(rterm,21) scform(j) 
 21   FORMAT(a30)

      WRITE(wterm,25) 
 25   FORMAT(/,' elemental chemical formula [e.g., Si(1)O(2)]:',/)
      READ(rterm,30) ecform(j) 
 30   FORMAT(a40)

      WRITE(wterm,35) 
 35   FORMAT(/,' reference [e.g., ref:5]:',/)
      READ(rterm,40) ref(j)
 40   FORMAT(a12)

      WRITE(wterm,45) 
 45   FORMAT(/,' date [e.g., 15.Jun.90]:',/)
      READ(rterm,50) date(j)
 50   FORMAT(a9)

      WRITE(wterm,46) 
 46   FORMAT(/,' standard molal G(cal/mol), H(cal/mol),'
     1        ,' S(cal/mol/K), V(cc/mol) at 25 C, 1 bar'
     2       /,' [NOTE: separate by commas, blanks, or <cr>;'
     3        ,' 999999 = NULL value for G, H]',/)
      READ(rterm,*) Gf(j), Hf(j), Sref(j), Vref(j)

      IF (ntran .EQ. 0) THEN
           WRITE(wterm,55)
 55        FORMAT(/,' Cp coeffs: Maier-Kelly  a(10**0),'
     1             ,' b(10**3), c(10**-5), and maximum temp (K):',/)
           READ(rterm,*) a(j,1), b(j,1), c(j,1), Tmax(j) 
           RETURN
      END IF

      DO 111 k = 1,ntran
           WRITE(wterm,65) k
 65        FORMAT(/,15x,' pre-transition ',i1,/)
           WRITE(wterm,75)
 75        FORMAT(/,' Cp coeffs: Maier-Kelly  a(10**0),'
     1             ,' b(10**3), c(10**-5), and maximum temp (K):',/)
           READ(rterm,*) a(j,k), b(j,k), c(j,k), Ttr(j,k) 
           WRITE(wterm,85) 
 85        FORMAT(/,' Htran(cal/mol), Vtran(cc/mol), dP/dT tran(bar/K)',
     2       /,' [NOTE: separate by commas, blanks, or <cr>;'
     3        ,' 999999 = NULL value]',/)
           READ(rterm,*) Htr(j,k), Vtr(j,k), Tslope(j,k)
 111       CONTINUE

      WRITE(wterm,95) ntran
 95   FORMAT(/,15x,' post-transition ',i1,/)
      WRITE(wterm,105) 
 105  FORMAT(/,' Cp coeffs: Maier-Kelly  a*(10**0),'
     1        ,' b(10**3), c(10**-5), and maximum temp (K):',/)
      READ(rterm,*) a(j,ntran+1), b(j,ntran+1), c(j,ntran+1), Tmax(j)

      END

***********************************************************************

*** match - Returns .TRUE. and if nlist > 0 AND the species referenced 
***         by name is found within list; otherwise, returns .FALSE.

      LOGICAL FUNCTION match(nleft,nlist,name,list)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXCHG = 100)

      CHARACTER*20  name, list(MAXCHG)

      SAVE


      match = .FALSE.
      
      IF (nleft .EQ. 0 ) RETURN

      DO 10 i = 1,nlist
           IF (name .EQ. list(i)) THEN
                match = .TRUE. 
                RETURN
           END IF
 10        CONTINUE

      RETURN
    
      END

***********************************************************************
 
*** header - Transfer the three header lines of 
***          the current species block.

      SUBROUTINE header(stype,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
    
      CHARACTER*3   stype
      CHARACTER*46  prompt(6)
      CHARACTER*80  line
      INTEGER       rterm, wterm, safold, safnew

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


      DATA prompt /
     1 'minerals that do not undergo phase transitions',
     2 'minerals that undergo one phase transition    ',
     3 'minerals that undergo two phase transitions   ',
     4 'minerals that undergo three phase transitions ',
     5 'gases                                         ',
     6 'aqueous species                               ' /


      READ(safold,10) line
 10   FORMAT(a80)
      WRITE(safnew,10) line

      READ(safold,20) 
 20   FORMAT()
      IF (stype .EQ. 'min') THEN
           index = ntran + 1
      ELSE IF (stype .EQ. 'gas') THEN
           index = 5
      ELSE
           index = 6     
      END IF
      WRITE(safnew,30) prompt(index)
 30   FORMAT(10x,a46)

      READ(safold,10) line
      WRITE(safnew,10) line

      END

**********************************************************************

*** readmg - Read all data for the current mineral or gas. 

      SUBROUTINE readmg(j,name,scform,abrev,ecform,ref,date,Gf,Hf,Sref,
     1                  Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMG = 200, NTPLUS = 4)

      INTEGER rterm, wterm, safold, safnew

      CHARACTER*9  date(MAXMG)
      CHARACTER*12 ref(MAXMG)
      CHARACTER*15 abrev(MAXMG)
      CHARACTER*20 name(MAXMG)
      CHARACTER*30 scform(MAXMG)
      CHARACTER*40 ecform(MAXMG)

      DOUBLE PRECISION Gf(MAXMG), Hf(MAXMG), Sref(MAXMG), Vref(MAXMG),
     1                 a(MAXMG,NTPLUS), b(MAXMG,NTPLUS), 
     2                 c(MAXMG,NTPLUS), Ttr(MAXMG,NTPLUS), 
     3                 Htr(MAXMG,NTPLUS), Vtr(MAXMG,NTPLUS),
     4                 Tslope(MAXMG,NTPLUS), Tmax(MAXMG)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** read data from old sprons.dat file 

      READ(safold,10) name(j), scform(j)
 10   FORMAT(1x,a20,a30)

      READ(safold,20) abrev(j), ecform(j)
 20   FORMAT(1x,a15,5x,a40)

      READ(safold,30) ref(j), date(j)
 30   FORMAT(1x,a12,8x,a9)

      READ(safold,*) Gf(j), Hf(j), Sref(j), Vref(j)
           
      DO 15 k = 1,ntran
           READ(safold,*) a(j,k), b(j,k), c(j,k), 
     1     Ttr(j,k), Htr(j,k), Vtr(j,k), Tslope(j,k)
 15        CONTINUE

      k = ntran + 1

      READ(safold,*) a(j,k), b(j,k), c(j,k)
      READ(safold,*) Tmax(j)

      END

**********************************************************************

*** getchg - Prompt for and read lists of species to be modified,
***          deleted, or added.

      SUBROUTINE getchg(stype,ntran,nadd,ndel,nmod,alist,dlist,mlist)
    
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXCHG = 100)

      CHARACTER*1   ans
      CHARACTER*3   stype
      CHARACTER*20  mlist(MAXCHG), dlist(MAXCHG), alist(MAXCHG)
      CHARACTER*55 prompt(6)
      INTEGER       rterm, wterm, safold, safnew

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE

      DATA prompt /
     1 ' MINERALS that DO NOT undergo phase transitions?  (y/n)',
     2 ' MINERALS that undergo ONE PHASE TRANSITION?  (y/n)    ',
     3 ' MINERALS that undergo TWO PHASE TRANSITIONS?  (y/n)   ',
     4 ' MINERALS that undergo THREE PHASE TRANSITIONS?  (y/n) ',
     5 ' GASES?  (y/n)                                         ',
     6 ' AQUEOUS SPECIES?  (y/n)                               ' /

      
      nadd = 0
      ndel = 0
      nmod = 0

      IF (stype .EQ. 'min') THEN
           index = ntran + 1
      ELSE IF (stype .EQ. 'gas') THEN
           index = 5
      ELSE
           index = 6     
      END IF

  1   WRITE(wterm,10) prompt(index)  
 10   FORMAT(/,' would you like to ADD, DELETE, or MODIFY data for any',
     1       /,a55,/)

      READ(rterm,15) ans
 15   FORMAT(a1)
      IF ((ans .NE. 'y') .AND. (ans .NE. 'Y') .AND.
     1    (ans .NE. 'n') .AND. (ans .NE. 'N')) GO TO 1
     
      IF ((ans .EQ. 'n') .OR. (ans .EQ. 'N')) RETURN

      WRITE(wterm,20) 
 20   FORMAT(/,' Specify number of species to ADD: ',/)
      READ(rterm,*) nadd

      IF (nadd .GT. 0) CALL getlst(nadd,alist)

      WRITE(wterm,30) 
 30   FORMAT(/,' Specify number of species to DELETE: ',/)
      READ(rterm,*) ndel

      IF (ndel .GT. 0) CALL getlst(ndel,dlist)

      WRITE(wterm,40) 
 40   FORMAT(/,' Specify number of species to MODIFY: ',/)
      READ(rterm,*) nmod

      IF (nmod .GT. 0) CALL getlst(nmod,mlist)

      END

***********************************************************************  

*** getlst - Prompt for and read the n-element list of species 

      SUBROUTINE getlst(n,list)

      PARAMETER (MAXCHG = 100)

      CHARACTER*20 list(MAXCHG) 
      INTEGER       rterm, wterm, safold, safnew

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


      WRITE(wterm,10) n
 10   FORMAT(/,' Specify ',i3,' species names; one per line',/)
      
      DO 20 i = 1,n
           READ (rterm,30) list(i)
 30        FORMAT(a20)
 20        CONTINUE 
           
      END

***********************************************************************  

*** tranmg - Transfer (the possibly revised) data for one mineral 
***          or gas to the new sprons.dat file.

      SUBROUTINE tranmg(j,name,scform,abrev,ecform,ref,date,Gf,Hf,
     1                  Sref,Vref,a,b,c,Ttr,Htr,Vtr,Tslope,Tmax,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMG = 200, NTPLUS = 4)

      INTEGER rterm, wterm, safold, safnew

      CHARACTER*9  date(MAXMG)
      CHARACTER*12 ref(MAXMG)
      CHARACTER*15 abrev(MAXMG)
      CHARACTER*20 name(MAXMG)
      CHARACTER*30 scform(MAXMG)
      CHARACTER*40 ecform(MAXMG)

      DOUBLE PRECISION Gf(MAXMG), Hf(MAXMG), Sref(MAXMG), Vref(MAXMG),
     1                 a(MAXMG,NTPLUS), b(MAXMG,NTPLUS), 
     2                 c(MAXMG,NTPLUS), Ttr(MAXMG,NTPLUS), 
     3                 Htr(MAXMG,NTPLUS), Vtr(MAXMG,NTPLUS),
     4                 Tslope(MAXMG,NTPLUS), Tmax(MAXMG)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** write data to new sprons.dat file 

      WRITE(safnew,40)  name(j), scform(j)
 40   FORMAT(1x,a20,a30)
      WRITE(safnew,50)  abrev(j), ecform(j)
 50   FORMAT(1x,a15,5x,a40)
      WRITE(safnew,60)  ref(j), date(j) 
 60   FORMAT(1x,a12,8x,a9)
      WRITE(safnew,70)  Gf(j), Hf(j), Sref(j), Vref(j)
 70   FORMAT(4x,2(2x,f12.1),2(2x,f8.3))
               
      DO 25 k = 1,ntran
           WRITE(safnew,75) a(j,k), b(j,k), c(j,k), Ttr(j,k),
     1                      Htr(j,k), Vtr(j,k), Tslope(j,k)
 75        FORMAT(4x,3(2x,f12.6),2x,f7.2,2x,f8.1,2(2x,f10.3))
 25        CONTINUE

      k = ntran + 1
      WRITE(safnew,80) a(j,k), b(j,k), c(j,k)
 80   FORMAT(4x,3(2x,f12.6))
      WRITE(safnew,90)  Tmax(j)
 90   FORMAT(8x,f7.2)

      END

**********************************************************************

*** aqueus - Transfer naqs aqueous species from the old to new 
***          sprons.dat file; if requested, revise the data for,
***          delete, or add the specified species; if necessary, 
***          update naqs and/or realphabetize.

      SUBROUTINE aqueus(naqs,stype,ntran)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXCHG = 100, MAXAQ = 500)

      LOGICAL       match, rename, realph

      CHARACTER*3   stype
      CHARACTER*9   date(MAXAQ)
      CHARACTER*12  ref(MAXAQ) 
      CHARACTER*15  abrev(MAXAQ)
      CHARACTER*20  name(MAXAQ), alist(MAXCHG), 
     1              dlist(MAXCHG), mlist(MAXCHG)
      CHARACTER*30  scform(MAXAQ)
      CHARACTER*40  ecform(MAXAQ)

      INTEGER       rterm, wterm, safold, safnew, order(MAXAQ)

      DOUBLE PRECISION Gf(MAXAQ), Hf(MAXAQ), Sref(MAXAQ),
     1                 a1(MAXAQ), a2(MAXAQ), a3(MAXAQ), a4(MAXAQ),
     2                 c1(MAXAQ), c2(MAXAQ), omega(MAXAQ), charge(MAXAQ)

      COMMON /io/   rterm, wterm, safold, safnew

      SAVE


      CALL header(stype,ntran)
      CALL getchg(stype,ntran,nadd,ndel,nmod,alist,dlist,mlist)

      ngone = 0
      ndone = 0
      realph = .FALSE.

      DO 10 i = 1,naqs
           j = i-ngone
           CALL readaq(j,name,scform,abrev,ecform,ref,date,
     1                 Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)
           IF (match(ndel-ngone,ndel,name(j),dlist)) THEN
***             delete the current aqueous species
                ngone = ngone + 1
                GO TO 10
           END IF
           IF (match(nmod-ndone,nmod,name(j),mlist)) THEN
***             modify data for the current aqueous species
                ndone = ndone + 1
                CALL modaq(j,name,scform,abrev,ecform,ref,date,
     1                     Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge,
     2                     rename)
                IF (rename) THEN
                     realph = .TRUE.
                END IF
           END IF
 10        CONTINUE

      IF (ndel .NE. ngone) THEN
***        one or more species to be deleted were not found
           WRITE(wterm,20) ndel-ngone
 20        FORMAT(/,i5,' of the species to be deleted were not found',/)
      END IF

      IF (nmod .NE. ndone) THEN
***        one or more species to be modified were not found
           WRITE(wterm,30) nmod-ndone
 30        FORMAT(/,i5,' of the species to be modified were not found',
     1            /)
      END IF

      naqs = naqs - ngone

      IF (nadd .GT. 0) THEN
           realph = .TRUE.
           DO 40 i = 1,nadd
***             obtain data for the new aqueous species
                naqs = naqs + 1
                name(naqs) = alist(i)
                CALL addaq(naqs,name,scform,abrev,ecform,ref,date,
     1                     Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)
 40             CONTINUE
      END IF

*** initialize order

      DO 50 i = 1,naqs
           order(i) = i
 50        CONTINUE

*** re-alphabetize if necessary 

      IF (realph) CALL alpha(naqs,name,order)

*** transfer data 

      DO 60 i = 1,naqs
           CALL tranaq(order(i),name,scform,abrev,ecform,ref,date,
     1                 Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)
 60        CONTINUE

      END

***********************************************************************

*** readaq - Read all data for current aqueous species. 

      SUBROUTINE readaq(j,name,scform,abrev,ecform,ref,date,
     1                  Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXAQ = 500)

      CHARACTER*9   date(MAXAQ)
      CHARACTER*12  ref(MAXAQ) 
      CHARACTER*15  abrev(MAXAQ)
      CHARACTER*20  name(MAXAQ)
      CHARACTER*30  scform(MAXAQ)
      CHARACTER*40  ecform(MAXAQ)

      INTEGER rterm, wterm, safold, safnew

      DOUBLE PRECISION Gf(MAXAQ), Hf(MAXAQ), Sref(MAXAQ),
     1                 a1(MAXAQ), a2(MAXAQ), a3(MAXAQ), a4(MAXAQ),
     2                 c1(MAXAQ), c2(MAXAQ), omega(MAXAQ), charge(MAXAQ)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** read data from old prons 

      READ(safold,10) name(j), scform(j)
 10   FORMAT(1x,a20,a30)
      READ(safold,20) abrev(j), ecform(j)
 20   FORMAT(1x,a15,5x,a40)
      READ(safold,30) ref(j), date(j)
 30   FORMAT(1x,a12,8x,a9)
      READ(safold,*)  Gf(j), Hf(j), Sref(j)
      READ(safold,*)  a1(j), a2(j), a3(j), a4(j)
      READ(safold,*)  c1(j), c2(j), omega(j), charge(j)

      END

***********************************************************************

*** modaq - Modify data for the aqueous species referenced by name.

      SUBROUTINE modaq(j,name,scform,abrev,ecform,ref,date,
     1                 Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge,
     2                 rename)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
    
      PARAMETER (MAXAQ = 500)

      LOGICAL modify, rename 

      CHARACTER*9   date(MAXAQ)
      CHARACTER*12  ref(MAXAQ) 
      CHARACTER*15  abrev(MAXAQ)
      CHARACTER*20  name(MAXAQ)
      CHARACTER*30  scform(MAXAQ)
      CHARACTER*40  ecform(MAXAQ)

      INTEGER rterm, wterm, safold, safnew

      DOUBLE PRECISION Gf(MAXAQ), Hf(MAXAQ), Sref(MAXAQ),
     1                 a1(MAXAQ), a2(MAXAQ), a3(MAXAQ), a4(MAXAQ),
     2                 c1(MAXAQ), c2(MAXAQ), omega(MAXAQ), charge(MAXAQ)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** modify data 

      WRITE(wterm,1) name(j)
  1   FORMAT(/,' modify data for ',a20,/)

***** specify / [update] name

      WRITE(wterm,10) name(j)
 10   FORMAT(/,' name: ',a20,/)
      IF (.NOT. modify()) THEN
	   rename = .FALSE.
      ELSE
	   rename = .TRUE.
           WRITE(wterm,20) 
 20        FORMAT(/,' specify new name: ',/)
           READ(rterm,30) name(j)
 30        FORMAT(a20)
      END IF

***** specify / [update] structural chemical formula

      WRITE(wterm,11) scform(j)
 11   FORMAT(/,' structural chemical formula: ',a30,/)
      IF (modify()) THEN
           WRITE(wterm,21) 
 21        FORMAT(/,' specify new structural chemical formula: ',/)
           READ(rterm,31) scform(j)
 31        FORMAT(a30)
      END IF

***** specify / [update] abbreviation

      WRITE(wterm,12) abrev(j)
 12   FORMAT(/,' abbreviation: ',a15,/)
      IF (modify()) THEN
           WRITE(wterm,22) 
 22        FORMAT(/,' specify new abbreviation: ',/)
           READ(rterm,32) abrev(j)
 32        FORMAT(a15)
      END IF

***** specify / [update] elemental chemical formula

      WRITE(wterm,13) ecform(j)
 13   FORMAT(/,' elemental chemical formula: ',a40,/)
      IF (modify()) THEN
           WRITE(wterm,23) 
 23        FORMAT(/,' specify new elemental chemical formula: ',/)
           READ(rterm,33) ecform(j)
 33        FORMAT(a40)
      END IF

***** specify / [update] reference

      WRITE(wterm,14) ref(j)
 14   FORMAT(/,' reference entry: ',a12,/)
      IF (modify()) THEN
           WRITE(wterm,24) 
 24        FORMAT(/,' specify new reference entry: ',/)
           READ(rterm,34) ref(j)
 34        FORMAT(a12)
      END IF

***** specify / update date

      WRITE(wterm,15) date(j)
 15   FORMAT(/,' date of last update: ',a9,/)
      WRITE(wterm,25) 
 25   FORMAT(/,' current date: ',/)
      READ(rterm,35) date(j)
 35   FORMAT(a9)

***** specify / [update] Standard partial molal Gf 

      WRITE(wterm,16) Gf(j)
 16   FORMAT(/,' standard partial molal Gibbs free energy',
     1       /,' of formation (cal/mol) at 25 C, 1 bar:  ',f10.0,/)
      IF (modify()) THEN
           WRITE(wterm,26) 
 26        FORMAT(/,' specify new Gf: ',/)
           READ(rterm,*) Gf(j)
      END IF

***** specify / [update] Standard partial molal Hf 

      WRITE(wterm,17) Hf(j)
 17   FORMAT(/,' standard partial molal enthalpy of formation',
     1       /,' (cal/mol) at 25 C, 1 bar:                   ',f10.0,/)
      IF (modify()) THEN
           WRITE(wterm,27) 
 27        FORMAT(/,' specify new Hf: ',/)
           READ(rterm,*) Hf(j)
      END IF

***** specify / [update] Standard partial molal S 

      WRITE(wterm,18) Sref(j)
 18   FORMAT(/,' standard partial molal entropy'
     1       /,' (cal/mol/K) at 25 C, 1 bar:   ',f8.3,/)
      IF (modify()) THEN
           WRITE(wterm,28) 
 28        FORMAT(/,' specify new SPrTr: ',/)
           READ(rterm,*) Sref(j)
      END IF

***** specify / [update] equation-of-state parameters

      WRITE(wterm,19) a1(j)
 19   FORMAT(/,' equation-of-state parameter  a1(10**1): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,29) 
 29        FORMAT(/,' specify new value for  a1(10**1): ',/)
           READ(rterm,*) a1(j)
      END IF
      WRITE(wterm,90) a2(j)
 90   FORMAT(/,' equation-of-state parameter  a2(10**-2): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,40) 
 40        FORMAT(/,' specify new value for  a2(10**-2): ',/)
           READ(rterm,*) a2(j)
      END IF
      WRITE(wterm,91) a3(j)
 91   FORMAT(/,' equation-of-state parameter  a3(10**0): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,41) 
 41        FORMAT(/,' specify new value for  a3(10**0): ',/)
           READ(rterm,*) a3(j)
      END IF
      WRITE(wterm,92) a4(j)
 92   FORMAT(/,' equation-of-state parameter  a4(10**-4): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,42) 
 42        FORMAT(/,' specify new value for  a4(10**-4): ',/)
           READ(rterm,*) a4(j)
      END IF
      WRITE(wterm,93) c1(j)
 93   FORMAT(/,' equation-of-state parameter  c1(10**0): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,43)
 43        FORMAT(/,' specify new value for  c1(10**0): ',/)
           READ(rterm,*) c1(j)
      END IF
      WRITE(wterm,94) c2(j)
 94   FORMAT(/,' equation-of-state parameter  c2(10**-4): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,44)
 44        FORMAT(/,' specify new value for  c2(10**-4): ',/)
           READ(rterm,*) c2(j)
      END IF
      WRITE(wterm,95) omega(j)
 95   FORMAT(/,' equation-of-state parameter  omega (cal/mol): ',f8.4,/)
      IF (modify()) THEN
           WRITE(wterm,45)
 45        FORMAT(/,' specify new value for  omega (cal/mol): ',/)
           READ(rterm,*) omega(j)
      END IF
      WRITE(wterm,96) charge(j)
 96   FORMAT(/,' valence state: ',f3.0,/)
      IF (modify()) THEN
           WRITE(wterm,46)
 46        FORMAT(/,' specify new value for valence state: ',/)
           READ(rterm,*) charge(j)
      END IF

      END 

************************************************************************

*** addaq - Add the aqueous species referenced by name 
***         to the new sprons.dat.

      SUBROUTINE addaq(j,name,scform,abrev,ecform,ref,date,
     1                 Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXAQ = 500)

      CHARACTER*9   date(MAXAQ)
      CHARACTER*12  ref(MAXAQ) 
      CHARACTER*15  abrev(MAXAQ)
      CHARACTER*20  name(MAXAQ)
      CHARACTER*30  scform(MAXAQ)
      CHARACTER*40  ecform(MAXAQ)

      INTEGER rterm, wterm, safold, safnew

      DOUBLE PRECISION Gf(MAXAQ), Hf(MAXAQ), Sref(MAXAQ),
     1                 a1(MAXAQ), a2(MAXAQ), a3(MAXAQ), a4(MAXAQ),
     2                 c1(MAXAQ), c2(MAXAQ), omega(MAXAQ), charge(MAXAQ)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** prompt for / read requisite data for new species "name"

      WRITE(wterm,10) name(j)
 10   FORMAT(/,' specify requisite data for ',a20,/)

      WRITE(wterm,15) 
 15   FORMAT(/,' abbreviation [e.g., Mg+2, TOL(aq)]:',/)
      READ(rterm,20) abrev(j)
 20   FORMAT(a15)

      WRITE(wterm,16) 
 16   FORMAT(/,' structural chemical formula [e.g., Mg(+2), C6H5CH3]:',
     1       /)
      READ(rterm,21) scform(j) 
 21   FORMAT(a30)

      WRITE(wterm,25) 
 25   FORMAT(/,' elemental chemical formula [e.g., Mg(1)+(2),'
     1        ,' C(7)H(8)+(0)]:',/)
      READ(rterm,30) ecform(j) 
 30   FORMAT(a40)

      WRITE(wterm,35) 
 35   FORMAT(/,' reference [e.g., ref:5]:',/)
      READ(rterm,40) ref(j)
 40   FORMAT(a12)

      WRITE(wterm,45) 
 45   FORMAT(/,' date [e.g., 15.Jun.90]:',/)
      READ(rterm,50) date(j)
 50   FORMAT(a9)

      WRITE(wterm,55) 
 55   FORMAT(/,' standard partial molal G(cal/mol), H(cal/mol),'
     1        ,' S(cal/mol/K) at 25 C, 1 bar',/)
      READ(rterm,*) Gf(j), Hf(j), Sref(j)

      WRITE(wterm,65)
 65   FORMAT(/,' equation-of-state coeffs:  a1(10**1), a2(10**-2),'
     1        ,' a3(10**0), a4(10**-4)',/)
      READ(rterm,*) a1(j), a2(j), a3(j), a4(j) 

      WRITE(wterm,75)
 75   FORMAT(/,' equation-of-state coeffs:  c1(10**0), c2(10**-4)',/)
      READ(rterm,*) c1(j), c2(j)

      WRITE(wterm,85)
 85   FORMAT(/,' omega(10**-5), charge',/)
      READ(rterm,*) omega(j), charge(j)

      END

************************************************************************

*** tranaq - Transfer data for one aqueous species from
***          the old to new sprons.dat file.

      SUBROUTINE tranaq(j,name,scform,abrev,ecform,ref,date,
     1                  Gf,Hf,Sref,a1,a2,a3,a4,c1,c2,omega,charge)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXAQ = 500)
 
      CHARACTER*9   date(MAXAQ)
      CHARACTER*12  ref(MAXAQ) 
      CHARACTER*15  abrev(MAXAQ)
      CHARACTER*20  name(MAXAQ)
      CHARACTER*30  scform(MAXAQ)
      CHARACTER*40  ecform(MAXAQ)

      INTEGER rterm, wterm, safold, safnew

      DOUBLE PRECISION Gf(MAXAQ), Hf(MAXAQ), Sref(MAXAQ),
     1                 a1(MAXAQ), a2(MAXAQ), a3(MAXAQ), a4(MAXAQ),
     2                 c1(MAXAQ), c2(MAXAQ), omega(MAXAQ), charge(MAXAQ)

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE


***** write data to new sprons.dat 

      WRITE(safnew,40)  name(j), scform(j)
 40   FORMAT(1x,a20,a30)
      WRITE(safnew,50)  abrev(j), ecform(j)
 50   FORMAT(1x,a15,5x,a40)
      WRITE(safnew,60)  ref(j), date(j)
 60   FORMAT(1x,a12,8x,a9)
      WRITE(safnew,70)  Gf(j), Hf(j), Sref(j)
 70   FORMAT(4x,2(2x,f10.0),4x,f8.3)
      WRITE(safnew,80)  a1(j), a2(j), a3(j), a4(j)
 80   FORMAT(4x,4(2x,f8.4,2x))
      WRITE(safnew,90)  c1(j), c2(j), omega(j), charge(j)
 90   FORMAT(4x,3(2x,f8.4,2x),9x,f3.0)

      END

*******************************************************************

*** summry - Write a succinct summary of file contents.

      SUBROUTINE summry(nmin1,nmin2,nmin3,nmin4,ngas,naqs,lskip)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER rterm, wterm, safold, safnew

      COMMON /io/ rterm, wterm, safold, safnew

      SAVE

      
      WRITE(safnew,30) 
 30   FORMAT(/)

      WRITE(safnew,35) nmin1
 35   FORMAT(5x,i4,'  minerals that do not undergo phase transition',
     1             '  (nmin1)')

      WRITE(safnew,45) nmin2
 45   FORMAT(5x,i4,'  minerals that undergo one phase transition',
     1             '     (nmin2)')

      WRITE(safnew,55) nmin3
 55   FORMAT(5x,i4,'  minerals that undergo two phase transitions',
     1             '    (nmin3)')

      WRITE(safnew,65) nmin4
 65   FORMAT(5x,i4,'  minerals that undergo three phase transitions',
     1             '  (nmin4)')

      WRITE(safnew,75) ngas
 75   FORMAT(5x,i4,'  gases                                        ',
     1             '  (ngas)')

      WRITE(safnew,85) naqs
 85   FORMAT(5x,i4,'  aqueous species                              ',
     1             '  (naqs)')

      WRITE(safnew,95) 
 95   FORMAT(5x,'----',2x,'----------------------------')

      WRITE(safnew,105) nmin1 + nmin2 + nmin3 + nmin4 + ngas + naqs
 105  FORMAT(5x,i4,2x,'total',/)

      WRITE(safnew,115) lskip
 115  FORMAT(5x,i4,2x,'comment lines preceed top of min1 block')

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
