*** reac92 - Calculates the standard molal Gibbs free energy, enthalpy,
***          entropy, heat capacity, and volume of the i[th] reaction
***          (specified in common blocks /icon/ and /reac/) among
***          <= MAXMIN minerals, <= MAXAQS aqueous species, <= MAXGAS
***          gases, and H2O using equations and data given by Helgeson
***          et al. (1978), Tanger and Helgeson (1988), Shock and 
***          Helgeson (1988, 1990), Shock et al. (1989, 1991), Johnson
***          and Norton (1991), Johnson et al. (1991), and Sverjensky
***          et al. (1991).
***
***          Computed reaction properties are stored in COMMON blocks
***          /minsp/, /gassp/, /aqsp/, /solvn/, and /fmeq/.
*** 
*******************************************************************
***
*** Author:     James W. Johnson
***             Earth Sciences Department, L-219
***             Lawrence Livermore National Laboratory
***             Livermore, CA 94550
***             johnson@s05.es.llnl.gov
***
*** Abandoned:  8 November 1991
***
*******************************************************************

      SUBROUTINE reac92(i,P,TC,Dw,Vw,betaw,alphaw,daldTw,
     1                  Sw,Cpw,Hw,Gw,Zw,Qw,Yw,Xw,geqn)

*******************************************************************

*** argument units:  (w suffixes denote H2O properties)
***
***         P ............ bars
***         TC ........... degC
***         Dw ........... g/cm**3
***         Vw ........... cm**3/mol
***         betaw, Qw .... bars**(-1)
***         alphaw, Yw ... K**(-1)
***         daldTw, Xw ... K**(-2)
***         Sw, Cpw ...... cal/(mol*K)
***         Hw, Gw ....... cal/mol
***         Zw ........... dimensionless

*********************************************

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXGAS = 10, MAXAQS = 10, MAXRXN = 10)

      INTEGER geqn

      LOGICAL m2reac(MAXRXN)

      CHARACTER*80  rtitle(MAXRXN)
      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS)
      DOUBLE PRECISION  coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                  coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      COMMON /reac1/ rtitle
      COMMON /reac2/ coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1               rec1m, rec1a, rec1g, m2reac

      SAVE


      TK = TC + 273.15d0
      CALL solids(nm(i),P,TK)
      CALL gases(ng(i),TK)
      CALL aqsps(na(i),P,TK,Dw,betaw,alphaw,daldTw,Zw,Qw,Yw,Xw,geqn)
      CALL reactn(i,TK,Vw,Sw,Cpw,Hw,Gw)

      END

********************************************************************

*** Solids - Computes the standard molal thermodynamic properties of
*            nmin minerals at P,T using equations given by
*            Helgeson et al. (1978).


      SUBROUTINE solids(nmin,P,T)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3)

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      INTEGER  ntran(MAXMIN), phaser(MAXMIN), getphr, getCpr, Cpreg

      DOUBLE PRECISION  mwH2O, TtranP(MXTRAN,MAXMIN),
     1                         PtranT(MXTRAN,MAXMIN)

      DOUBLE PRECISION  Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                  Hmin(MAXMIN), Gmin(MAXMIN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), SPrTrm(MAXMIN),
     2                  VPrTrm(MAXMIN), 
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN)

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
 
      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran

      COMMON /minsp/  Vmin, Smin, Cpmin, Hmin, Gmin, phaser

      COMMON /PTtran/ TtranP, PtranT

      SAVE


      DO 10  i = 1,nmin
           phaser(i) = getphr(i,P,T,TtranP)
           Cpreg = getCpr(i,T)
           CALL Vterms(i,P,T,phaser(i),Vmin(i),VdP,PtranT)
           CALL Cptrms('min',i,Cpreg,T,Cpmin(i),CprdT,CprdlT)
           CALL pttrms(i,phaser(i),T,Spttrm,Hpttrm,Gpttrm)
           Smin(i) = SPrTrm(i) + CprdlT + Spttrm
           Hmin(i) = Hfmin(i) + CprdT + VdP + Hpttrm
           Gmin(i) = Gfmin(i) - SPrTrm(i)*(T-Tref) +
     1               CprdT - T*CprdlT + VdP + Gpttrm
           IF ((mname(i) .EQ. 'QUARTZ') .OR.
     1         (mname(i) .EQ. 'COESITE')) THEN
                Hmin(i) = Hmin(i) - VdP
                Gmin(i) = Gmin(i) - VdP
                CALL quartz(mname(i),P,T,Ttran(1,i),
     1                      Vmin(i),Smin(i),Hmin(i),Gmin(i))
           END IF
 10        CONTINUE

      RETURN
      END

************************************************************************

*** getCpr - Returns the effective phase region for temperature 
***          integration of Cpr(T) for mineral imin (i.e., the
***          phase region specified by T at 1 bar).


      INTEGER FUNCTION getCpr(imin,T)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3)

      INTEGER  ntran(MAXMIN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), SPrTrm(MAXMIN),
     2                  VPrTrm(MAXMIN), 
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN)

      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran

      SAVE


      getCpr = 1
      DO 10 i = 1,ntran(imin)
           IF (T .GT. Ttran(i,imin)) getCpr = getCpr + 1
 10        CONTINUE

      RETURN
      END 

***********************************************************************

*** quartz - Revises the standard molal Gibbs free energy (G), enthalpy
***          (H), entropy (S), and volume (V) of quartz or coesite to
***          account for V(T) > 0 using equations (109) through (115), 
***          Helgeson et al. (1978). 


      SUBROUTINE quartz(mname,P,T,TtPr,V,S,H,G)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      CHARACTER*20      mname 
      INTEGER           qphase
      DOUBLE PRECISION  k, mwH2O

      COMMON /refval/ mwH2O, R, Pr, Tr, ZPrTr, YPrTr
      COMMON /qtzcon/ aa, ba, ca, VPtTta, VPrTtb, Stran

      SAVE
      
*** VPrTra = VPrTr(a-quartz) 
*** Vdiff  = VPrTr(a-quartz) - VPrTr(coesite)
*** k      = dPdTtr(a/b-quartz)  

      DATA VPrTra, Vdiff, k / 22.688d0, 2.047d0, 38.5d0 /


***** set qphase = phase region of quartz

      IF ((T .LE. TtPr) .OR. (P .GE. (Pr + k*(T-TtPr)))) THEN
           qphase = 1
      ELSE
           qphase = 2
      END IF

***** set Pstar and Sstar *****

      IF (T .LE. TtPr) THEN
           Pstar = Pr 
           Sstar = 0.0d0
      ELSE
           IF (qphase .EQ. 2) THEN
                Pstar = P
                Sstar = 0.0d0
           ELSE
                Pstar = Pr + k*(T-TtPr)
                Sstar = Stran
           END IF
      END IF

      IF (qphase .EQ. 2) THEN
***** set volume to beta-quartz *****
           V = VPrTtb
      ELSE
***** calculate volume of alpha-quartz per eqn (109) *****
      V = VPrTra + ca*(P-Pr) + (VPtTta - VPrTra - ca*(P-Pr))*(T-Tr) /
     1    (TtPr + (P-Pr)/k - Tr)
      END IF

      IF (mname .EQ. 'COESITE') V = V - Vdiff

***** leading constant for [G,S]Vterm below 
***** is a coversion factor (cal/cm**3/bar) 

      IF (mname .EQ. 'QUARTZ') THEN
        GVterm = 0.23901488d-1 * (VPrTra*(P-Pstar) + VPrTtb*(Pstar-Pr) -
     1        0.5d0*ca*(2.0d0*Pr*(P-Pstar) - (P**2-Pstar**2)) -
     2        ca*k*(T-Tr)*(P-Pstar) + 
     3        k*(ba + aa*ca*k)*(T-Tr)*DLOG((aa + P/k)/(aa + Pstar/k)))
      ELSE
        GVterm = 0.23901488d-1 * ((VPrTra-Vdiff)*(P-Pstar) +
     1        (VPrTtb-Vdiff)*(Pstar-Pr) - 0.5d0*ca*(2.0d0*Pr*(P-Pstar) - 
     2        (P**2-Pstar**2)) - ca*k*(T-Tr)*(P-Pstar) + 
     3        k*(ba + aa*ca*k)*(T-Tr)*DLOG((aa + P/k)/(aa + Pstar/k)))
      END IF

      SVterm = 0.23901488d-1 * (-k*(ba + aa*ca*k)*
     1         DLOG((aa + P/k)/(aa + Pstar/k)) + ca*k*(P-Pstar)) - 
     2         Sstar

      G = G + GVterm
      S = S + SVterm
      H = H + GVterm + T*SVterm

      END

************************************************************************

*** getphr - Returns phase region for mineral imin at P, T; and, as a
***          side effect, TtranP(1..MXTRAN,imin) as f(P).
***
***          getphr = 1 ... TtranP(1,imin) > T  [or imin lacks transn]
***          getphr = 2 ... TtranP(1,imin) < T  [< TtranP(2,imin)]
***          getphr = 3 ... TtranP(2,imin) < T  [< TtranP(3,imin)]
***          getphr = 4 ... TtranP(3,imin) < T


      INTEGER FUNCTION getphr(imin,P,T,TtranP)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3)

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      INTEGER  ntran(MAXMIN)

      DOUBLE PRECISION  TtranP(MXTRAN,MAXMIN), mwH2O

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     2                  SPrTrm(MAXMIN), 
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN), 
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN)
 
      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr

      SAVE


****** phase region 1 ******

      getphr = 1
      IF (ntran(imin) .EQ. 0) RETURN

      IF (dPdTtr(1,imin) .EQ. 0.0d0) THEN
           TtranP(1,imin) = Ttran(1,imin)
      ELSE
           TtranP(1,imin) = Ttran(1,imin) + (P-Pref)/dPdTtr(1,imin)
      END IF
      IF (T .LE. TtranP(1,imin)) RETURN

****** phase region 2 ******

      getphr = 2
      IF (ntran(imin) .EQ. 1)   RETURN

      IF (dPdTtr(2,imin) .EQ. 0.0d0) THEN
           TtranP(2,imin) = Ttran(2,imin)
      ELSE
           TtranP(2,imin) = Ttran(2,imin) + (P-Pref)/dPdTtr(2,imin)
      END IF
      IF (T .LE. TtranP(2,imin)) RETURN

****** phase region 3 ******

      getphr = 3
      IF (ntran(imin) .EQ. 2)   RETURN

      IF (dPdTtr(3,imin) .EQ. 0.0d0) THEN
           TtranP(3,imin) = Ttran(3,imin)
      ELSE
           TtranP(3,imin) = Ttran(3,imin) + (P-Pref)/dPdTtr(3,imin)
      END IF
      IF (T .LE. TtranP(3,imin)) RETURN

****** phase region 4 ******

      getphr = 4
      RETURN

      END

************************************************************************

*** Vterms - Computes Vmin(P,T), Vmin*dP, and (if necesary) PtranT.


      SUBROUTINE Vterms(imin,P,T,phaser,Vmin,VdP,PtranT)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3)

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      INTEGER  phaser, ntran(MAXMIN)

      DOUBLE PRECISION  mwH2O, PtranT(MXTRAN,MAXMIN)

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     2                  SPrTrm(MAXMIN), 
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN),
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN),
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN)
 
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
 
      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran

      SAVE


      Vmin = VPrTrm(imin)
      DO 10 i = 1,phaser-1
 10        Vmin = Vmin + Vtran(i,imin)
      VdP  = Vmin*(P - Pref)*0.23901488d-1

****** return if Pressure integration does not cross 
****** phase transition boundaries

      IF (ntran(imin) .EQ. 0)        RETURN
      IF (dPdTtr(1,imin) .EQ. 0.0d0) RETURN 
      IF (T .LE. Ttran(1,imin))      RETURN  
      IF ((ntran(imin) .EQ. 1) .AND. (phaser .EQ. 2))  RETURN
      IF ((ntran(imin) .EQ. 2) .AND. (phaser .EQ. 3))  RETURN
      IF ((ntran(imin) .EQ. 2) .AND. (phaser .EQ. 2) .AND.
     1   (T .LT. Ttran(2,imin)))     RETURN

****** take account of cross-boundary pressure integration 

      IF ((ntran(imin) .EQ. 1) .OR. 
     1   ((phaser .EQ. 1) .AND.(T .LT. Ttran(2,imin)))) THEN
           PtranT(1,imin) = Pref + (T - Ttran(1,imin))*dPdTtr(1,imin)
           VdP  = 0.23901488d-1 * (
     1            VPrTrm(imin)*(P - Pref) + 
     2            Vtran(1,imin)*(PtranT(1,imin) - Pref))
           RETURN
      END IF

****** ntran(imin) = 2 and T .GE. Ttran(2,imin) ******

      PtranT(2,imin) = Pref + (T - Ttran(2,imin))*dPdTtr(2,imin)

      IF (phaser .EQ. 2) THEN
           VdP  = 0.23901488d-1 * (
     1            (VPrTrm(imin) + Vtran(1,imin))*(P - Pref) + 
     2            Vtran(2,imin)*(PtranT(2,imin) - Pref))
      ELSE
           PtranT(1,imin) = Pref + (T - Ttran(1,imin))*dPdTtr(1,imin)
           VdP  = 0.23901488d-1 * (
     1            VPrTrm(imin)*(P - Pref) + 
     1            Vtran(1,imin)*(PtranT(1,imin) - Pref) + 
     2            Vtran(2,imin)*(PtranT(2,imin) - Pref))
      END IF

      RETURN

      END

************************************************************************

*** Cptrms - Computes the standard molal heat capacity and heat capacity 
***          temperature integrals, evaluated from Tref to T at 1 bar.


      SUBROUTINE Cptrms(phase,i,Cpreg,T,Cpr,CprdT,CprdlT)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXGAS = 10, MAXMIN = 10, MXTRAN = 3, IABC   = 3)

      CHARACTER*3   phase
      CHARACTER*20  mname(MAXMIN), gname(MAXGAS)
      CHARACTER*30  mform(MAXMIN), gform(MAXGAS)

      INTEGER  ntran(MAXMIN), Cpreg

      DOUBLE PRECISION mwH2O

      DOUBLE PRECISION Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     2                 SPrTrm(MAXMIN), 
     3                 MK1(IABC,MAXMIN), MK2(IABC,MAXMIN), 
     4                 MK3(IABC,MAXMIN), MK4(IABC,MAXMIN), 
     5                 Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                 Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                 Tmaxm(MAXMIN)

      DOUBLE PRECISION Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     2                 SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)
 
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran
      COMMON /gnames/ gname, gform 
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg

      SAVE


      IF (phase .EQ. 'gas') THEN 
         Cpr = Cp(T,MKg(1,i),MKg(2,i),MKg(3,i))
         CprdT = CpdT(Tref,T,MKg(1,i),MKg(2,i),MKg(3,i))
         CprdlT = CpdlnT(Tref,T,MKg(1,i),MKg(2,i),MKg(3,i))
         RETURN
      END IF

***** phase = "min" *****

      IF (Cpreg .EQ. 1) THEN
         Cpr = Cp(T,MK1(1,i),MK1(2,i),MK1(3,i))
         CprdT = CpdT(Tref,T,MK1(1,i),MK1(2,i),MK1(3,i))
         CprdlT = CpdlnT(Tref,T,MK1(1,i),MK1(2,i),MK1(3,i))
         RETURN
      END IF

      IF (Cpreg .EQ. 2) THEN
         Cpr = Cp(T,MK2(1,i),MK2(2,i),MK2(3,i))
         CprdT = CpdT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i)) +
     2           CpdT(Ttran(1,i),T,MK2(1,i),MK2(2,i),MK2(3,i))
         CprdlT = CpdlnT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i)) +
     2            CpdlnT(Ttran(1,i),T,MK2(1,i),MK2(2,i),MK2(3,i))
         RETURN
      END IF

      IF (Cpreg .EQ. 3) THEN
           Cpr = Cp(T,MK3(1,i),MK3(2,i),MK3(3,i))
           CprdT = CpdT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i)) +
     2        CpdT(Ttran(1,i),Ttran(2,i),MK2(1,i),MK2(2,i),MK2(3,i)) +
     3        CpdT(Ttran(2,i),T,MK3(1,i),MK3(2,i),MK3(3,i))
           CprdlT = CpdlnT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i))+
     2       CpdlnT(Ttran(1,i),Ttran(2,i),MK2(1,i),MK2(2,i),MK2(3,i))+
     3       CpdlnT(Ttran(2,i),T,MK3(1,i),MK3(2,i),MK3(3,i))
           RETURN
      END IF

***** Cpreg = 4 *****

      Cpr = Cp(T,MK4(1,i),MK4(2,i),MK4(3,i))
      CprdT = CpdT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i)) +
     2        CpdT(Ttran(1,i),Ttran(2,i),MK2(1,i),MK2(2,i),MK2(3,i)) +
     3        CpdT(Ttran(2,i),Ttran(3,i),MK3(1,i),MK3(2,i),MK3(3,i)) +
     4        CpdT(Ttran(3,i),T,MK4(1,i),MK4(2,i),MK4(3,i)) 
      CprdlT = CpdlnT(Tref,Ttran(1,i),MK1(1,i),MK1(2,i),MK1(3,i)) +
     2       CpdlnT(Ttran(1,i),Ttran(2,i),MK2(1,i),MK2(2,i),MK2(3,i))+
     3       CpdlnT(Ttran(2,i),Ttran(3,i),MK3(1,i),MK3(2,i),MK3(3,i))+
     4       CpdlnT(Ttran(3,i),T,MK4(1,i),MK4(2,i),MK4(3,i))
      RETURN

      END

*********************************************************************

*** Cp - Returns the standard molal heat capacity at T.


      DOUBLE PRECISION FUNCTION Cp(T,a,b,c)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      Cp = a + b*T + c/T**2

      RETURN
      END

*********************************************************************

*** CpdT - Returns the integral CpdT evaluated from T1 to T2.


      DOUBLE PRECISION FUNCTION CpdT(T1,T2,a,b,c)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      CpdT = a*(T2 - T1) + b/2.0d0*(T2**2 - T1**2) - 
     2       c*(1.0d0/T2 - 1.0d0/T1)

      RETURN
      END

*********************************************************************

*** CpdlnT - Returns the integral CpdlnT evaluated from T1 to T2.


      DOUBLE PRECISION FUNCTION CpdlnT(T1,T2,a,b,c)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      SAVE

      CpdlnT = a*DLOG(T2/T1) + b*(T2 - T1) - 
     2         c/2.0d0*(1.0d0/T2**2 - 1.0d0/T1**2)

      RETURN
      END

*********************************************************************

*** pttrms - Computes phase transition terms for Smin, Hmin, and Gmin.


      SUBROUTINE pttrms(imin,phaser,T,Spttrm,Hpttrm,Gpttrm)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MXTRAN = 3, IABC = 3)

      CHARACTER*20  mname(MAXMIN)
      CHARACTER*30  mform(MAXMIN)

      INTEGER  ntran(MAXMIN), phtran, phaser

      DOUBLE PRECISION  Gfmin(MAXMIN), Hfmin(MAXMIN), VPrTrm(MAXMIN),
     2                  SPrTrm(MAXMIN), 
     3                  MK1(IABC,MAXMIN), MK2(IABC,MAXMIN), 
     4                  MK3(IABC,MAXMIN), MK4(IABC,MAXMIN), 
     5                  Ttran(MXTRAN,MAXMIN), Htran(MXTRAN,MAXMIN),
     6                  Vtran(MXTRAN,MAXMIN), dPdTtr(MXTRAN,MAXMIN),
     7                  Tmaxm(MAXMIN)

      COMMON /mnames/ mname, mform
      COMMON /minref/ Gfmin, Hfmin, SPrTrm, VPrTrm, MK1, MK2, MK3, MK4, 
     1                Ttran, Htran, Vtran, dPdTtr, Tmaxm, ntran

      SAVE


      Spttrm = 0.0d0
      Hpttrm = 0.0d0
      Gpttrm = 0.0d0
      DO 10 phtran = 1,phaser-1
           Spttrm = Spttrm + Htran(phtran,imin)/Ttran(phtran,imin)
           Hpttrm = Hpttrm + Htran(phtran,imin)
           Gpttrm = Gpttrm + 
     1              Htran(phtran,imin)*(1.0d0 - T/Ttran(phtran,imin))
 10        CONTINUE

      RETURN
      END

**********************************************************************

*** gases - Computes the standard molal thermodynamic properties of
*           ngas gases at P,T using equations given by
*           Helgeson et al. (1978).


      SUBROUTINE gases(ngas,T)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXGAS = 10, IABC = 3, TS1BAR = 99.6324d0)

      LOGICAL  error

      CHARACTER*20  gname(MAXGAS)
      CHARACTER*30  gform(MAXGAS)

      INTEGER specs(10)

      DOUBLE PRECISION  mwH2O

      DOUBLE PRECISION  Vgas(MAXGAS), Sgas(MAXGAS), Cpgas(MAXGAS), 
     2                  Hgas(MAXGAS), Ggas(MAXGAS)

      DOUBLE PRECISION  Gfgas(MAXGAS), Hfgas(MAXGAS), VPrTrg(MAXGAS), 
     2                  SPrTrg(MAXGAS), MKg(IABC,MAXGAS), Tmaxg(MAXGAS)
 
      DOUBLE PRECISION  states(4), props(46)

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /gnames/ gname, gform 
      COMMON /gasref/ Gfgas, Hfgas, SPrTrg, VPrTrg, MKg, Tmaxg
      COMMON /gassp/  Vgas, Sgas, Cpgas, Hgas, Ggas 
 
      SAVE

      DATA specs  / 2,2,2,5,1,0,2,0,4,0 /
      DATA states / 0.0d0, 1.0d0, 0.0d0, 0.0d0 /
 

      TC = T - 273.15d0

      DO 10  i = 1,ngas
 
           IF ((gname(i) .EQ. 'H2O,g') .AND. (TC .GE. TS1BAR)) THEN
***             use Haar et al. (1984) equation of state to
***             compute H2O,g properties at 1 bar, T > Tsat(1 bar) =
***             99.6324 C.  Note that for T < Tsat(1 bar), 
***             thermodynamic properties of metastable H2O,g are
***             calculated using parameters estimated by J. W. Johnson
***             (3/90) that facilitate smooth transition into the 
***             Haar et al. (1984) equation at Tsat.
***
***             Beacuse (1) P = 1 bar, and (2) thermodynamic properties
***             of steam are independent of dielectric properties,
***             specs(8..9) can be safely hardwired, as above.

                states(1) = TC
                CALL H2O92(specs,states,props,error)
            
                Vgas(i)  = VPrTrg(i)
                Sgas(i)  = props(5)
                Hgas(i)  = props(9)
                Ggas(i)  = props(3)
                Cpgas(i) = props(13)
           ELSE 
                Vgas(i) = VPrTrg(i)

                CALL Cptrms('gas',i,1,T,Cpgas(i),CprdT,CprdlT)
 
                Sgas(i) = SPrTrg(i) + CprdlT

                Hgas(i) = Hfgas(i) + CprdT

                Ggas(i) = Gfgas(i) - SPrTrg(i)*(T - Tref) + 
     1                    CprdT - T*CprdlT
           END IF

 10        CONTINUE

      RETURN
      END

************************************************************************

*** aqsps - Computes the standard partial molal thermodynamic properties
***         of naqs aqueous species at P,T using equations given by
***         Tanger and Helgeson (1988), Shock et al. (1991), and 
***         Johnson et al. (1991).


      SUBROUTINE aqsps(naqs,P,T,Dw,betaw,alphaw,daldTw,Z,Q,Y,X,geqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER  (MAXAQS = 10)

      CHARACTER*20  aname(MAXAQS)
      CHARACTER*30  aform(MAXAQS)

      INTEGER geqn

      DOUBLE PRECISION mwH2O

      DOUBLE PRECISION Vaqs(MAXAQS), Saqs(MAXAQS), Cpaqs(MAXAQS),
     2                 Haqs(MAXAQS), Gaqs(MAXAQS),
     7                 VQterm(MAXAQS), SYterm(MAXAQS), CpXtrm(MAXAQS),
     8                 HYterm(MAXAQS), GZterm(MAXAQS)

      DOUBLE PRECISION  Gfaqs(MAXAQS), Hfaqs(MAXAQS), SPrTra(MAXAQS), 
     2                  a(4,MAXAQS), c(2,MAXAQS), wref(MAXAQS), 
     3                  chg(MAXAQS)
 
      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
      COMMON /aqscon/ eta, theta, psi, anion, cation, gref
      COMMON /anames/ aname, aform
      COMMON /aqsref/ Gfaqs, Hfaqs, SPrTra, c, a, wref, chg
      COMMON /aqsp/   Vaqs, Saqs, Cpaqs, Haqs, Gaqs
      COMMON /solvn/  VQterm, SYterm, CpXtrm, HYterm, GZterm

      SAVE

    
      IF (naqs .EQ. 0) THEN
           RETURN
      ELSE
           CALL gfun92(T-273.15d0,P,Dw,betaw,alphaw,daldTw,
     1                 g,dgdP,dgdT,d2gdT2,geqn)
      END IF

      DO 10  j = 1,naqs
****** compute w, dwdP, dwdT, d2wdT2 ******
           CALL omeg92(g,dgdP,dgdT,d2gdT2,wref(j),chg(j),
     1                 w,dwdP,dwdT,d2wdT2,aname(j))

           VQterm(j) = 0.4184004d2 * (-w*Q + (-Z - 1.0d0)*dwdP)
*** the leading constant converts cal/(mol*bar) -> cm3/mol
           Vaqs(j) = 0.4184004d2 * (a(1,j) + 
     1               a(2,j)/(psi+P) + 
     2               a(3,j)/(T-theta) + 
     3               a(4,j)/(psi+P)/(T-theta)) +
     4               VQterm(j)

           SYterm(j) = w*Y - (-Z - 1.0d0)*dwdT - wref(j)*YPrTr
           Saqs(j) = SPrTra(j) + c(1,j)*DLOG(T/Tref) -
     2               c(2,j)/theta* (1.0d0/(T-theta) - 
     3                              1.0d0/(Tref-theta) +
     4                              (1.0d0/theta)*
     5                          DLOG(Tref*(T-theta)/T/(Tref-theta))) +
     6               (a(3,j)*(P-Pref) + 
     7                a(4,j)*DLOG((psi+P)/(psi+Pref))) * 
     8               (1.0d0/(T-theta))**2 + 
     9               SYterm(j)

           CpXtrm(j) = w*T*X + 2.0d0*T*Y*dwdT + T*(Z + 1.0d0)*d2wdT2
           Cpaqs(j)  = c(1,j) + c(2,j)/(T-theta)**2 -
     1                 (2.0d0*T/(T-theta)**3) * (a(3,j)*(P-Pref) + 
     2                 a(4,j)*DLOG((psi+P)/(psi+Pref))) +
     3                 CpXtrm(j)

           HYterm(j) = w*(-Z - 1.0d0) + w*T*Y - T*(-Z - 1.0d0)*dwdT -
     1                 wref(j)*(-ZPrTr - 1.0d0) -
     2                 wref(j)*Tref*YPrTr
           Haqs(j) = Hfaqs(j) + c(1,j)*(T-Tref) -
     1             c(2,j)*(1.0d0/(T-theta) - 1.0d0/(Tref-theta)) +
     2             a(1,j)*(P-Pref) + a(2,j)*DLOG((psi+P)/(psi+Pref)) +
     3            (a(3,j)*(P-Pref) + a(4,j)*DLOG((psi+P)/(psi+Pref))) *
     4            ((2.0d0*T - theta)/(T - theta)**2) + 
     5             HYterm(j)


           GZterm(j) = w*(-Z - 1.0d0) - wref(j)*(-ZPrTr - 1.0d0) +
     1                 wref(j)*YPrTr*(T-Tref)
           Gaqs(j) = Gfaqs(j) - SPrTra(j)*(T-Tref) -
     1             c(1,j)*(T*DLOG(T/Tref)-T+Tref) +
     2             a(1,j)*(P-Pref) + a(2,j)*DLOG((psi+P)/(psi+Pref)) -
     3             c(2,j)* ( (1.0d0/(T-theta) - 1.0d0/(Tref-theta)) *
     4             ((theta-T)/theta) - T/theta**2 * 
     5             DLOG((Tref*(T-theta))/(T*(Tref-theta))) ) +
     6             (1.0d0/(T-theta)) * (a(3,j)*(P-Pref) + 
     7             a(4,j)*DLOG((psi+P)/(psi+Pref))) +
     8             GZterm(j)
           GZterm(j) = w*(-Z - 1.0d0)

  10       CONTINUE

      RETURN
      END

************************************************************************

*** omeg92 - Computes the conventinal Born coefficient (w) of the 
***          current aqueous species, dwdP, dwdP, and dw2dT2 as a 
***          function of g, dgdP, dgdT, d2gdT2, wref, and Z using
***          equations given by Johnson et al. (1991).


      SUBROUTINE omeg92(g,dgdP,dgdT,d2gdT2,wref,Z,
     1                  w,dwdP,dwdT,d2wdT2,aname)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      CHARACTER*20 aname

      COMMON /aqscon/ eta, theta, psi, anion, cation, gref

      SAVE


      IF ((Z .EQ. 0.0d0) .OR. (aname .EQ. 'H+')) THEN
***        neutral aqueous species or H+ 
           w      = wref
           dwdP   = 0.0d0
           dwdT   = 0.0d0
           d2wdT2 = 0.0d0
           RETURN
      ELSE
***        charged aqueous species other than H+
           reref = Z**2 / (wref/eta + Z/(3.082d0 + gref))
           re = reref + DABS(Z) * g
           w  = eta * (Z**2/re - Z/(3.082d0 + g))
           Z3 = DABS(Z**3)/re**2 - Z/(3.082d0 + g)**2
           Z4 = DABS(Z**4)/re**3 - Z/(3.082d0 + g)**3
           dwdP   = -eta * Z3 * dgdP
           dwdT   = -eta * Z3 * dgdT
           d2wdT2 = 2.0d0 * eta * Z4 * dgdT**2 - eta * Z3 * d2gdT2
      END IF

      END

************************************************************************

*** reactn - Computes the standard molal thermodynamic properties
***          of the i[th] reaction.


      SUBROUTINE reactn(i,T,Vw,Sw,Cpw,Hw,Gw)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER (MAXMIN = 10, MAXAQS = 10, MAXGAS = 10, MAXRXN = 10)

      CHARACTER*80  rtitle(MAXRXN)

      LOGICAL  m2reac(MAXRXN)

      INTEGER  nm(MAXRXN), na(MAXRXN), ng(MAXRXN), nw(MAXRXN),
     1         rec1m(MAXRXN,MAXMIN), rec1a(MAXRXN,MAXAQS), 
     2         rec1g(MAXRXN,MAXGAS), phaser(MAXMIN)

      DOUBLE PRECISION mwH2O

      DOUBLE PRECISION coefm(MAXRXN,MAXMIN), coefa(MAXRXN,MAXAQS),
     1                 coefg(MAXRXN,MAXGAS), coefw(MAXRXN)

      DOUBLE PRECISION Vmin(MAXMIN), Smin(MAXMIN), Cpmin(MAXMIN),
     2                 Hmin(MAXMIN), Gmin(MAXMIN),
     3                 Vgas(MAXGAS), Sgas(MAXGAS), Cpgas(MAXGAS),
     4                 Hgas(MAXGAS), Ggas(MAXGAS),
     5                 Vaqs(MAXAQS), Saqs(MAXAQS), Cpaqs(MAXAQS),
     6                 Haqs(MAXAQS), Gaqs(MAXAQS)

      DOUBLE PRECISION VQterm(MAXAQS), SYterm(MAXAQS), CpXtrm(MAXAQS),
     1                 HYterm(MAXAQS), GZterm(MAXAQS), logKr

      COMMON /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
 
      COMMON /reac1/ rtitle
      COMMON /reac2/ coefm, coefa, coefg, coefw, nm, na, ng, nw,
     1               rec1m, rec1a, rec1g, m2reac

      COMMON /minsp/ Vmin, Smin, Cpmin, Hmin, Gmin, phaser
      COMMON /gassp/ Vgas, Sgas, Cpgas, Hgas, Ggas
      COMMON /aqsp/  Vaqs, Saqs, Cpaqs, Haqs, Gaqs
      COMMON /fmeq/  dVr,  dSr,  dCpr,  dHr,  dGr,
     2               logKr, dlogKT, dlogKP
      COMMON /solvn/ VQterm, SYterm, CpXtrm, HYterm, GZterm

      SAVE


***** sum mineral contributions *****

      dVrm  = 0.0d0
      dCprm = 0.0d0
      dSrm  = 0.0d0
      dHrm  = 0.0d0
      dGrm  = 0.0d0
      DO 10 j = 1,nm(i)
           dVrm  = dVrm  + coefm(i,j)*Vmin(j)
           dCprm = dCprm + coefm(i,j)*Cpmin(j)
           dSrm  = dSrm  + coefm(i,j)*Smin(j)
           dHrm  = dHrm  + coefm(i,j)*Hmin(j)
           dGrm  = dGrm  + coefm(i,j)*Gmin(j)
 10        CONTINUE

***** sum gas contributions *****

      dVrg  = 0.0d0
      dCprg = 0.0d0
      dSrg  = 0.0d0
      dHrg  = 0.0d0
      dGrg  = 0.0d0
      DO 20 j = 1,ng(i)
           dVrg  = dVrg  + coefg(i,j)*Vgas(j)
           dCprg = dCprg + coefg(i,j)*Cpgas(j)
           dSrg  = dSrg  + coefg(i,j)*Sgas(j)
           dHrg  = dHrg  + coefg(i,j)*Hgas(j)
           dGrg  = dGrg  + coefg(i,j)*Ggas(j)
 20        CONTINUE

***** sum aqueous species contributions *****

      dVra  = 0.0d0
      dCpra = 0.0d0
      dSra  = 0.0d0
      dHra  = 0.0d0
      dGra  = 0.0d0
      DO 30 j = 1,na(i)
           dVra  = dVra  + coefa(i,j)*Vaqs(j)
           dCpra = dCpra + coefa(i,j)*Cpaqs(j)
           dSra  = dSra  + coefa(i,j)*Saqs(j)
           dHra  = dHra  + coefa(i,j)*Haqs(j)
           dGra  = dGra  + coefa(i,j)*Gaqs(j)
 30        CONTINUE

***** calculate H2O contributions *****

      dVrw  = coefw(i) * Vw
      dSrw  = coefw(i) * Sw
      dCprw = coefw(i) * Cpw
      dHrw  = coefw(i) * Hw
      dGrw  = coefw(i) * Gw

***** calculate reaction properties *****

      dVr  = dVrm  + dVrg  + dVra  + dVrw
      dSr  = dSrm  + dSrg  + dSra  + dSrw
      dCpr = dCprm + dCprg + dCpra + dCprw
      dHr  = dHrm  + dHrg  + dHra  + dHrw
      dGr  = dGrm  + dGrg  + dGra  + dGrw

      logKr  = -dGr / (2.302585d0 * R * T)
      dlogKT =  dHr / (2.302585d0 * R * T**2) 
      dlogKP = -0.23901488d-1 * dVr / (2.302585d0 * R * T)

      RETURN
      END

******************************************************************

*** gfun92 - Computes the g function (Tanger and Helgeson, 1988;
***          Shock et al., 1991) and its partial derivatives
***          (dgdP, dgdT, d2gdT2) at TdegC, Pbars using the 
***          computational algorithm specified by geqn.
***        
***        geqn = 1 ...... use Tanger-Helgeson (1988) equations 
***        geqn = 2 ...... use Shock et al. (1991) equations
***                        without the f(P,T) difference function
***        geqn = 3 ...... use Shock et al. (1991) equations 
***                        with the f(P,T) difference function


      SUBROUTINE gfun92(TdegC,Pbars,Dgcm3,betab,alphaK,daldT,
     1                  g,dgdP,dgdT,d2gdT2,geqn)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      PARAMETER  (TMAX = 1000.0d0,  PMAX = 5000.0d0,  TOL=1.0d-4)
 
      INTEGER  geqn

      SAVE


****** initialize g and derivatives to zero
      g      = 0.0d0
      dgdP   = 0.0d0
      dgdT   = 0.0d0
      d2gdT2 = 0.0d0

      IF ((TdegC .GT. TMAX+TOL) .OR. (Pbars .GT. PMAX+TOL)) RETURN

*     IF (geqn .EQ. 1) THEN
******     use Tanger-Helgeson (1988) equations 
*          CALL gTangr(Pbars,TdegC+273.15d0,Dgcm3,betab,alphaK,daldT,
*    2                 g,dgdP,dgdT,d2gdT2)
*          RETURN
*     END IF


*     IF (geqn .EQ. 2) THEN
******     use Shock et al. (1991) equations 
******     without f(P,T) difference function
*          CALL gShok1(TdegC,Pbars,Dgcm3,betab,alphaK,daldT,
*    2                 g,dgdP,dgdT,d2gdT2)
*          RETURN
*     END IF

      IF (geqn .EQ. 3) THEN
******     use Shock et al. (1991) equations 
******     with f(P,T) difference function
           CALL gShok2(TdegC,Pbars,Dgcm3,betab,alphaK,daldT,
     2                 g,dgdP,dgdT,d2gdT2)
           RETURN
      END IF

      END

*****************************************************************

*** gShok2- Computes g, dgdP, dgdT, and d2gdT2 using equations given 
***         by Shock et al. (1991) 
***
*** units:   T ................. C
***          D ................. g/cm**3
***          beta, dgdP ........ bars**(-1)
***          alpha, dgdT ....... K**(-1)
***          daldT, d2gdT2 ..... K**(-2)

      SUBROUTINE gShok2(T,P,D,beta,alpha,daldT,g,dgdP,dgdT,d2gdT2)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION c(6), cc(3)

      SAVE

      DATA c /  -0.2037662D+01,  0.5747000D-02, -0.6557892D-05,
     1           0.6107361D+01, -0.1074377D-01,  0.1268348D-04 /
 
      DATA cc /  0.3666666D+02, -0.1504956D-9,   0.5017997D-13 /


      IF (D .GE. 1.0d0) RETURN

      a = c(1) + c(2)*T + c(3)*T**2
      b = c(4) + c(5)*T + c(6)*T**2
      g = a*(1.0d0 - D)**b

      dgdD   = - a * b * (1.0d0 - D)**(b - 1.0d0)
      dgdD2  =   a * b * (b - 1.0d0) * (1.0d0 - D)**(b - 2.0d0)

      dadT   =   c(2) + 2.0d0*c(3)*T
      dadTT  =   2.0d0 * c(3)
      dbdT   =   c(5) + 2.0d0*c(6)*T
      dbdTT  =   2.0d0 * c(6) 

      dDdT   = - D * alpha
      dDdP   =   D * beta
      dDdTT  = - D * (daldT - alpha**2)

      Db     = (1.0d0 - D) ** b

      dDbdT  = -b * (1.0d0 - D)**(b-1.0d0) * dDdT +
     1         DLOG(1.0d0 - D) * Db  * dbdT

      dDbdTT = -(b * (1.0d0 - D)**(b-1.0d0) * dDdTT +
     1           (1.0d0 - D)**(b-1.0d0) * dDdT * dbdT + b * dDdT *
     2           (-(b-1.0d0) * (1.0d0 - D)**(b-2.0d0) * dDdT +
     3           DLOG(1.0d0 - D) * (1.0d0 - D)**(b-1.0d0) * dbdT)) +
     4           DLOG(1.0d0 - D) * (1.0d0 - D)**b * dbdTT -
     5           (1.0d0 - D)**b * dbdT * dDdT / (1.0d0 - D) +
     6           DLOG(1.0d0 - D) * dbdT * dDbdT

      dgdP   = dgdD * dDdP
      dgdT   = a*dDbdT + Db*dadT
      d2gdT2 = a*dDbdTT + 2.0d0*dDbdT*dadT + Db*dadTT

      IF ((T .LT. 155.0d0) .OR. (P .GT. 1000.0d0) .OR. 
     1    (T .GT. 355.0d0)) RETURN

      ft     = ((T - 155.0d0)/300.0d0)**4.8 +
     1         cc(1)*((T - 155.0d0)/300.0d0)**16

      dftdT  = 4.8d0/300.0d0*((T - 155.0d0)/300.0d0)**3.8 +
     1        16.0d0/300.0d0*cc(1)*((T - 155.0d0)/300.0d0)**15

      dftdTT = 3.8d0*4.8d0/300.0d0**2*((T - 155.0d0)/300.0d0)**2.8 +
     1        15.0d0*16.0d0/300.0d0**2*cc(1)*((T - 155.0d0)/300.0d0)**14

      fp     = cc(2)*(1000.0d0 - P)**3 + cc(3)*(1000.0d0 - P)**4

      dfpdP  = -3.0d0*cc(2)*(1000.0d0 - P)**2 -
     1          4.0d0*cc(3)*(1000.0d0 - P)**3 

      f      = ft * fp
      dfdP   = ft * dfpdP
      dfdT   = fp * dftdT
      d2fdT2 = fp * dftdTT

      g      = g      - f
      dgdP   = dgdP   - dfdP
      dgdT   = dgdT   - dfdT
      d2gdT2 = d2gdT2 - d2fdT2

      RETURN

      END
