# SUPCRT

A thin wrapper for building SUPCRT. This repository provides a Makefile to simplify the build process.

NOTE: The source code of SUPCRT is listed as abandoned, hence the inclusion of the source code in this repository. I am happy to remove the source code if this violates any copyright or licenses.

The main differences between this repository and the original SUPCRT92 source code are:

* `DPRONS92.DAT` has been removed (it will be generated in the build process)
* the `README.DOC` has been converted to Markdown and included in this document
* all filenames are in lower case

## Build Instructions

To build SUPCRT, make sure you have [`gfortran`](https://`en.wikipedia`.org/wiki/GNU_Fortran). Once that is installed, run:

```sh
make
```

You can also check that the compiled code is correct by checking that the output matches the test cases:

```sh
make test
```

# The SUPCRT92 Software Package (DOS disk)

[This is the original README documentation from 1991-11-08.]

## Introduction

Welcome to SUPCRT92; we hope you will enjoy using the package. In addition to the double-sided/high-density (1.2 mb) 5.25" DOS disk from which you obtained the present file, you should now have in your possession a copy of "SUPCRT92: A Software Package for Calculating the Standard Molal Thermodynamic Properties of Minerals, Gases, Aqueous Species, and Reactions from 1 to 5000 bars and 0 to 1000 C" (Johnson, Oelkers, and Helgeson, 1992). Questions regarding the capabilities, operation, and effective use of the software are referred to this manuscript, which presents a summary of the equations and data encoded in SUPCRT92, detailed user's documentation for the package, and appendices that contain sample interactive sessions with the programs. The following discussion is largely restricted to contents of the distribution disk, the library of test problems, and installation of the package.

Contents of the Distribution Disk
---------------------------------

The distribution disk contains 50 files that together occupy slightly less that one megabyte of the disk:

    cprons92.f  rep92pc.f   skdef2.gxy  sktd2.con
    h2o92D.f    skarn.rxn   skdef2.hxy  sktd2.tab
    mprons92.f  skdef1.cxn  skdef2.kxy  sktp1.con
    qc.rxn      skdef1.dxy  skdef2.sxy  sktp1.tab
    qcgetp.con  skdef1.gxy  skdef2.tab  sktp2.con
    qcgetp.tab  skdef1.hxy  skdef2.vxy  sktp2.tab
    qcgett.con  skdef1.kxy  skdef3.tab  species.rxn
    qcgett.tab  skdef1.sxy  skodd1.con  species.tab
    qcgett.uxy  skdef1.tab  skodd1.tab  sprons92.dat
    reac92D.f   skdef1.vxy  skodd2.con  sup92D.f
    reac92pc.f  skdef2.2xy  skodd2.tab  sup92pc.f
    readme.doc  skdef2.cxy  sktd1.con
    rep92D.f    skdef2.dxy  sktd1.tab

Among these are nine FORTRAN77 source files:

    cprons92.f  reac92D.f   rep92pc.f
    h2o92D.f    reac92pc.f  sup92D.f
    mprons92.f  rep92D.f    sup92pc.f

One thermodynamic database in sequential-access format:

    sprons92.dat

A test library that includes three input reaction (`*.rxn`) files:

    qc.rxn
    skarn.rxn
    species.rxn

Eight input files that specify the conditions over which thermodynamic properties of the reactions are to be calculated (`*.con` files):

    qcgetp.con  sktd1.con
    qcgett.con  sktd2.con
    skodd1.con  sktp1.con
    skodd2.con  sktp2.con

Twelve output files that tabulate the calculated reaction properties (`*.tab` files):

    qcgetp.tab  skdef3.tab  sktd2.tab
    qcgett.tab  skodd1.tab  sktp1.tab
    skdef1.tab  skodd2.tab  sktp2.tab
    skdef2.tab  sktd1.tab   species.tab

And sixteen files that list the calculated properties for `qcgett.tab`, `skdef1.tab`, and `skdef2.tab` in a format suitable for x-y plotting (`*.*xy` files):

    qcgett.uxy  skdef2.2xy
    skdef1.cxy  skdef2.cxy
    skdef1.dxy  skdef2.dxy
    skdef1.gxy  skdef2.gxy
    skdef1.hxy  skdef2.hxy
    skdef1.kxy  skdef2.kxy
    skdef1.sxy  skdef2.sxy
    skdef1.vxy  skdef2.vxy

## The Test Library

The test library consists of twelve problems; the relationship between the foregoing input and output files is summarized in the following table:

               INPUT                           OUTPUT
       ---------------------           ----------------------
       *.con           *.rxn           *.tab           *.*xy
    ____________    ____________    ____________    ____________

       default1        skarn           skdef1          skdef1
       default2        skarn           skdef2          skdef2
       default3        skarn           skdef3          ------
       sktp1           skarn           sktp1           ------
       sktp2           skarn           sktp2           ------
       sktd1           skarn           sktd1           ------
       sktd2           skarn           sktd2           ------
       skodd1          skarn           skodd1          ------
       skodd2          skarn           skodd2          ------
       default2        species         species         ------
       qcgett          qc              qcgett          qcgett
       qcgetp          qc              qcgetp          ------

In this table, default[1..3] refers to default options that are selected during an interactive session with SUPCRT92. 

The test library provides a representative sampling of the types of geochemical calculations that can be carried out using SUPCRT92. Brief explanations are given below for each of the twelve sample problems.

       *.tab                   summary
    ____________    _______________________________________________

       skdef1          Calculations across the default P-T grid of
                       uniform increments in the single-phase regions
                       of fluid H2O. Because the skarn reaction
                       contains minerals that do and do not undergo
                       phase transitions, gases, and charged and
                       neutral aqueous species, this run produces
                       examples of most of the data and explanatory
                       comments that may be written to a *.tab file.

       skdef2          Calculations across the default temperature
                       range of uniform increments along the liquid
                       side of the H2O vaporization boundary.

       skdef3          Calculations across the eight-element grid
                       used by the EQ3/6 software package.

       sktp1           Calculations along evenly-spaced isobars as
                       a function of uniform temperature increments
                       within the single-phase regions of fluid H2O.

       sktp2           Calculations along evenly-spaced isotherms as
                       a function of uniform pressure increments
                       within the single-phase regions of fluid H2O.

       sktd1           Calculations along evenly-spaced isochores as
                       a function of uniform temperature increments
                       within the single-phase regions of fluid H2O.

       sktd2           Calculations along evenly-spaced isotherms as
                       a function of uniform increments in H2O density
                       within the single-phase regions of fluid H2O.

       skodd1          Calculations for a specified set of [H2O
                       density, temperature] coordinates that are
                       separated by nonuniform increments in the
                       single-phase regions of fluid H2O.

       skodd2          Calculations for a specified set of H2O
                       temperatures separated by nonuniform
                       increments along the liquid side of the
                       H2O vaporization boundary.

       species         Calculations across the default temperature
                       range of uniform increments along the liquid
                       side of the H2O vaporization boundary.
                       The species.rxn file contains nine reactions,
                       each specifying one of the species in the
                       preceding skarn reaction. This run
                       illustrates retrieval of the thermodynamic
                       properties of individual species as well as
                       evaluation of multiple reactions.

       qcgett          Calculation of the P-T projection of
                       univariant quartz-coesite equilibrium as
                       a function of specified pressures.

       qcgetp          Calculation of the P-T projection of
                       univariant quartz-coesite equilibrium as
                       a function of specified temperatures.

We encourage you to exercise SUPCRT92 on each of these test problems as a means of ensuring that the code is operating correctly in your local computing environment. Please report any significant discrepancies to James W. Johnson (address below).

## Compiling and Linking MPRONS92, CPRONS92, and SUPCRT92 in Most Environments

The SUPCRT92 software package is written in FORTRAN77 ANSI standard; hence, we anticipate minimal portability problems. In most environments, installation of the software can be achieved with straightforward use of the resident FORTRAN77 compiler and linker. (THOSE READERS LACKING THIS EXPERTISE ARE URGED TO CONSULT THEIR LOCAL EXPERTS.) Briefly, the generic installation procedure is as follows:

1. Compile and link `cprons92.f` to create the executable CPRONS92.

2. Execute CPRONS92 once, specifying `sprons92.dat` as the input sequential-access thermodynamic database and `dprons92.dat` as the output direct-access database.

3. Compile and link `sup92D.f`, `rep92D.f`, `reac92D.f`, and `h2o92D.f` to create the executable SUPCRT92.

4. Execute SUPCRT92 using the default thermodynamic database (`dprons92.dat`) and verify its correct operation by running the library of test problems described above.

5. If you desire to modify the thermodynamic database, compile and link `mprons92.f` to create the executable MPRONS92, which you can then use to revise `sprons92.dat`.

6. Pipe the revised `sprons92.dat` through CPRONS92 and return to step 4, now using this new, revised thermodynamic database instead of `dprons92.dat`.

## Compiling and Linking MPRONS92, CPRONS92, and SUPCRT92 on a DOS-based PC

The procedure for installing the SUPCRT92 package into the DOS-based PC (Personal Computer) environment is more involved than that outlined above. Because accessible memory on such machines is typically less than 640 kilobytes, it is necessary to use alternate source files for three SUPCRT92 modules; specifically, one should use `sup92pc.f`, `rep92pc.f`, and `reac92pc.f` in lieu of their `*D.f` counterparts. As noted by Johnson, Oelkers, and Helgeson (1992), these `*pc.f` source files differ from the corresponding `*D.f` files only in certain array dimensions, which limits the scope of problems that can be run. Otherwise, however, the two sets of files are identical. In addition to using these alternate `*pc.f` source files, a greater level of expertise with the resident FORTRAN77 compiler is required to transform these files into the operational programs. THOSE READERS LACKING THIS EXPERTISE ARE URGED TO CONSULT THEIR LOCAL EXPERTS.

The following procedure was used successfully to compile, link, and execute the SUPCRT92 package on an standard-configuration IBM PC AT equipped with Microsoft's FORTRAN Optimizing Compiler Version 4.00 and the associated Overlay Linker Version 3.55.

1. Rename the source files such that the current default extension ".f" is changed to ".for".

2. Compile and link `cprons92.f` to create the executable CPRONS92:

        FL /Gt /AH cprons92.f

3. Execute CPRONS92 once, specifying `sprons92.dat` as the input sequential-access thermodynamic database and `dprons92.dat` as the output direct-access database.

4. Compile and link `sup92pc.f`, `rep92pc.f`, `reac92pc.f`, and `h2o92D.f` to create the executable SUP92PC:

        FL /c /Gt /AH /Od sup92pc.f
        FL /c /Gt /AH /Od rep92pc.f
        FL /c /Gt /AH /Od reac92pc.f
        FL /c /Gt /AH /Od h2o92D.f

        LINK sup92pc rep92pc reac92pc h2o92D /I

5. Execute SUP92PC using the default thermodynamic database (`dprons92.dat`) and verify its correct operation by running the library of test problems described above.

6. If you desire to modify the thermodynamic database, compile and link `mprons92.f` to create the executable MPRONS92,

        FL /Gt /AH mprons92.f

    which you can then use to revise `sprons92.dat`.

7. Pipe the revised `sprons92.dat` through CPRONS92 and return to step 5, now using this new, revised thermodynamic database instead of `dprons92.dat`.

## Compiling and Linking MPRONS92, CPRONS92, and SUPCRT92 on an Apple Macintosh

THOSE READERS LACKING FAMILIARITY WITH THE MAC ENVIRONMENT AND EXPERIENCE USING FORTRAN COMPILERS IN THIS ENVIRONMENT ARE URGED TO CONSULT THEIR LOCAL EXPERTS BEFORE ATTEMPTING TO INSTALL THE SUPCRT92 PACKAGE ONTO A MAC.

The following procedure was used successfully to compile, link, and execute the SUPCRT92 package on a Mac IIci equipped with Absoft's MacFortran/020 Compiler (Version 2.3).

1. Update the values of the two variables --- rterm and wterm --- that specify the unit numbers associated with reading data from and writing data to the terminal screen. This must be done in source files `mprons92.f`, `cprons92.f`, and `sup92D.f`. In each of these files, rterm and wterm are assigned global values by the appropriate DATA statement within BLOCK DATA consts. Use your local editor to change the values of rterm and wterm from 5 and 6, respectively, to 9 and 9.

2. Compile `cprons92.f` to create the executable CPRONS92 application. (Compiler Options for "ignore character case", "68020/68030 instructions", and "68881/68882 instructions" should be activated.)

3. Execute CPRONS92 once, specifying `sprons92.dat` as the input sequential-access thermodynamic database and `dprons92.dat` as the output direct-access database.

4. Concatenate `sup92D.f`, `rep92D.f`, `reac92D.f`, and `h2o92D.f` into a single file (named, `e.g`., `supcrt92.f`).

5. Compile `supcrt92.f` to create the executable SUPCRT92 application. (Use the same Compiler Options given in step 2.)

6. Augment the default application memory size for SUPCRT92 by clicking on the SUPCRT92 application, selecting "Get Info" from the Multifinder menu bar, and increasing the value given in the "Application Memory Size (K)" field to 1024. (If necessary on your machine, you may be able to set this memory size to a slightly smaller value.)

7. Execute SUPCRT92 using the default thermodynamic database (`dprons92.dat`) and verify its correct operation by running the library of test problems described above.

8. If you desire to modify the thermodynamic database, compile `mprons92.f` to create the executable MPRONS92 application, which you can then use to revise `sprons92.dat`. (Use the same Compiler Options given in step 2.)

9. Pipe the revised `sprons92.dat` through CPRONS92 and return to step 7, now using this new, revised thermodynamic database instead of `dprons92.dat`.

## Storage Requirements of the Executable SUPCRT92 Package

The executable SUPCRT92 package consists of the interactive programs CPRONS92, MPRONS92, and SUPCRT92 together with the sequential- and direct-access thermodynamic databases, `sprons92.dat` and `dprons92.dat`. The exact size of these files will, of course, vary slightly with the local computing environment. The required storage on an Alliant FX/80 and VAX 11/750 are summarized below as examples.

                    Alliant FX/80   VAX 11/750
      module        (kilobytes)      (blocks)
    ------------    -------------   ----------
    CPRONS92            28.672           15
    MPRONS92            86.016           63
    SUPCRT92           471.040          580
    sprons92.dat       155.398          313
    dprons92.dat       270.810          529
                    -------------   ----------
                      1011.936         1500

## Troubleshooting in the Local Computing Environment

As noted above, portability problems are expected to be minimal. Nevertheless, two commonly encountered difficulties warrant mention. The first of these is related to the database codes MPRONS92 and CPRONS92. If these programs terminate abnormally with a run-time error message such as "incomprehensible list input", the problem almost certainly lies with the current value of PARAMETER NBACK in SUBROUTINE tail of `cprons92.f` and SUBROUTINE opener of `mprons92.f`. As indicated by documentation within these routines, the appropriate value of NBACK is 10 on some systems, 11 on others; recently, we discovered that the appropriate value is 12 on at least one system. Hence, if this problem occurs, edit the files to incorporate the indicated adjustment and recompile and link the codes; this should eliminate the difficulty.

The second potential problem relates to reading [coefficient, species] pairs while interactively building a reaction during SUPCRT92 execution. If the program permanently stalls during these reads, edit sup92[D,pc].f, move to LOGICAL FUNCTION parse, activate the commented-out line

    READ(numstr,*) r8num

and deactivate the next five executable statements, which are

    OPEN(UNIT=tempf,FILE='zero.dat')
    WRITE(tempf,*) numstr
    BACKSPACE(tempf)
    READ(tempf,*) r8num
    CLOSE(UNIT=tempf)

Now recompile, link, and execute the program; all should run smoothly.

## Questions Regarding the SUPCRT92 Software Package
-------------------------------------------------

Questions regarding SUPCRT92 software should be addressed to

    James W. Johnson (email address: johnson@s05.es.llnl.gov)
    Earth Sciences Department, L-219
    Lawrence Livermore National Laboratory
    Livermore, CA 94550

Questions regarding the equations and data used in SUPCRT92 may be addressed to Johnson or to

    Harold C. Helgeson
    Department of Geology and Geophysics
    University of California
    Berkeley, CA 94720
