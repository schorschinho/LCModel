c  To Do:
c  Check for '<<<<<<<<<<<<<<<'
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PROGRAM LCMODL
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C LCModel license
C
C    Copyright (c) 1992-2021, Stephen Provencher
C    All rights reserved.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are met:
C
C  1. Redistributions of source code must retain the above copyright notice, this
C   list of conditions and the following disclaimer.
C
C  2. Redistributions in binary form must reproduce the above copyright notice,
C     this list of conditions and the following disclaimer in the documentation
C     and/or other materials provided with the distribution.
C
C  3. Neither the name of the copyright holder nor the names of its
C     contributors may be used to endorse or promote products derived from
C     this software without specific prior written permission.
C
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
C  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
C  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
C  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
C  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
C  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
C  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
C  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
C  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvFor user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     
C  LCMODL.  Automatic estimation of metabolite concentrations from in vivo NMR
C             spectra.
C           The data spectrum is analyzed as a linear combination from a basis
C             set of model in vitro spectra, with automatic corrections for
C             phase, field inhomogeneity, eddy currents, background, and
C             differences in 1/T2.
C
C
C  To use this program, you need the following references:
C    S.W. Provencher (1993) Magn. Reson. Med. 30, 672-679.
C                    (1993--2021) LCModel Users Manual.
C
C  The first reference above is referred to as simply "MRM" in the comment
C    statements here.
C  In general, the comment statements are mainly for my use.  You need the
C    User's Manual and MRM.
C  Comment statements for you are set off by 2 comment lines, the top with a
C    string of "v" characters and the bottom with a string of "^", both with
C    "For user" in the middle (as these comments are).  Comments for you occur
C    only in subprograms MYCONT, MYDATA, and MYBASI.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
C  Input:
C
C    On unit LCONTR (must be standard input):
C       NAMELIST LCMODL.
C
C    On unit LBASIS (file FILBAS) (in subprogram MYBASI):
C       NAMELIST BASIS1
C       NAMELIST BASIS
C       Frequency-domain basis spectra from zero-filling
C                                      (2*NUNFIL=NDATA COMPLEX values).
C
C    On unit LRAW (file FILRAW) (in subprogram MYDATA):
C       NAMELIST NMID
C       Time-domain data (NUNFIL COMPLEX values)
C
C    If DOECC or DOWS or UNSUPR, on unit LH2O (file FILH2O) (in subprogram
c      MYDATA):
C       NAMELIST NMID
C       Time-domain non-H2O-suppressed data
C
C
C  Output:
C
C    If LPRINT>0, on unit LPRINT (file FILPRI):
C       Detailed output for printing.
C
C    If LCOORD>0, on unit LCOORD (file FILCOO):
C       Output for possible subsequent external plot program.
C
C    If LCOraw>0, on unit LCoraw (file FILCOr):
C       Output of RAW data, after ECC, *FCALIB*TRAMP/VOLUME, correction for
C         BRUKER=T or SEQACQ=T; corrections for phase & referencing shift.
C
C    If LPS>0, on unit LPS (file FILPS):
C       .PS file for PostScript printer.
C
C    If LTABLE>0, on unit LTABLE (file FILTAB):
C       short .TABLE file, containing only the 4 tables from 1-page output.
C
C
C  The fit is to the real part of the spectrum obtained with zero filling.
C     The spectra are rearranged; therefore, the analysis also works for
C     PPM>PPMCEN.
c  Old grid positions + NUNFIL.
c  Rearrangement in CFFT_r and DCFFT_r, but not in CFFT or CFFTIN (since
c    BASISF is not rearranged)
c  SEQACC = T not tested.
c  Had to redefine INCDIM, DELPPM, and test for INITIA 3.
C
C  Fortran 77 extensions: COMPLEX*16, INCLUDE, NAMELIST,
C  (Obsolete) Sun version (marked ! followed directly by Sun): Time, date
C  IRIX Version (marked ! followed directly by IRIX): READONLY, TIME, DATE,
C                                                     UNKNOWN, ZEXP
c  ***************************************************************************
C                 markers preceded by a !
c                 =======================
C  Cyg: (formerly Cygwin g77 version) Intel version for MS Windows.  License
c       tests are skipped.  It is assumed that a HASP envelope protects it.
c  sun: g77 Sun cross-compiler (important to only leave this not commented out;
C       its preparation will only comment out Linux preceded by ! and could
c       otherwise produce the statement twice.  The others comment out all
C       lines with ! string, except its own).
c  IRIX: IRIX
c  OSF: DEC Fortran
c  Linux: g77 PC Linux
c  ***************************************************************************
C  Version-dependent parts are marked with C!!!!...
C
C  Extra tests (COMMENTed out):
C    !CT1: Checks (quadratic) convergence with exact data from previous fit
C          when ISTAGE=1.
C    !CT2: As !CT1, but for ISTAGE=2.
C
c  LDUMP(1) = T to dump weights, etc, for phased-array averaging in AVERAGE.
c       (2) = T to dump details of gaps and baseline regions on ppm-axis.
c       (3) = T to dump compacted lines of CONTROL file.
c       (4) = T & SUBBAS=T & NEACH=99 to dump ratios of peak areas to Concs to
c               check normalization and scaling of simulated basis spectra.
c       (5) = T to dump the real part of basis around WSPPM in AREABA.
c
      INCLUDE 'lcmodel.inc'
      external ilen
      logical iok
      data nanalyses_done/0/, nvoxels_done/0/, nvoxels_done_in/0/
      chsubp = 'MAIN'
      version_lcm = '6.3-1N'
      lversion_lcm = ilen(version_lcm)
      VERSIO='LCModel (Version ' // version_lcm(:lversion_lcm) //
     1       ') Copyright: S.W. Provencher.' //
     2       '          Ref.: Magn. Reson. Med. 30:672-679 (1993).'
C     -------------------------------------------------------------------------
C     Get changes to Control Variables.
C     -------------------------------------------------------------------------
      CALL MYCONT ()
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Uncomment the following to initialize license KEYs and other Control
c       Parameters for Frahm.
c
c      include 'frahm.inc'
c
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     -------------------------------------------------------------------------
c     Load ZERO_VOXEL array, which has values of T for zero voxels.
c     Skip this if BASCAL=T.
C     -------------------------------------------------------------------------
      if (.not.bascal) then
         call check_zero_voxels ()
C        ----------------------------------------------------------------------
c        IAVERG >= 1 to call AVERAGE to average series of "voxels" in CSI
c                    format in DATAT (as in phased-array or Philips or Elscint
c                    series of scans) and then set control parameters for
c                    single-voxel analysis.
C        ----------------------------------------------------------------------
         if (iaverg .ge. 1) call average ()
      end if
 
      if (lcsi_sav_1 .eq. 12) then
         if (lcsi_sav_2 .ne. 13) go to 801
         open (12, file=filcsi_sav_1, status='old', err=801)
         open (13, file=filcsi_sav_2, status='old', err=801)
         read (12, 5030, err=801, end=801) nvoxels_done, nanalyses_done,
     1                                     ioffset_current
 5030    format (3i5)
         nvoxels_done_in = nvoxels_done
         ioffset_current_in = ioffset_current
         if (nanalyses_done .gt. 0) read (12, 5040, err=801, end=40)
     1                              degppm, degzer, dgppmn, dgppmx,
     2                              sddegp, sddegz, shifmn, shifmx
 5040       format (1p5e16.6)
      else
         nxoxels_done_in = 0
         nanalyses_done = 0
         ioffset_current_in = 0
      end if
 
 40   IF (Lcsv .GT. 0   .and.   FILcsv .NE. ' ') then
         if (nanalyses_done .le. 0) then
            OPEN (Lcsv, FILE=FILcsv, STATUS='UNKNOWN', err=802)
         else
            OPEN (Lcsv, FILE=FILcsv, STATUS='OLD', err=802)
            do 30 jline = 1, 9999
 5020          format (a)
               read (lcsv, 5020, err=803, end=35) chline
 30         continue
 35         backspace lcsv
            go to 50
         end if
      end if
C     -------------------------------------------------------------------------
C     Main loop for all analyses of all voxels.
c     Dimensions ND* I* have been checked and ordered in MYCONT
c     Round up or down to put i*_center closest to the overall center.
C     -------------------------------------------------------------------------
 50   center_whole = float(ndrows + 1) / 2.
      i1 = (irowst + irowen) / 2
      i2 = (irowst + irowen + 1) / 2
      if (abs(float(i1) - center_whole) .lt.
     1    abs(float(i2) - center_whole)) then
         irow_center = i1
      else
         irow_center = i2
      end if
      center_whole = float(ndcols + 1) / 2.
      i1 = (icolst + icolen) / 2
      i2 = (icolst + icolen + 1) / 2
      if (abs(float(i1) - center_whole) .lt.
     1    abs(float(i2) - center_whole)) then
         icol_center = i1
      else
         icol_center = i2
      end if
      single_voxel = max0(ndslic, ndrows, ndcols) .eq. 1
      noffset = max0(irowen - irow_center, irow_center - irowst,
     1               icolen - icol_center, icol_center - icolst)
C     -------------------------------------------------------------------------
c     If there is only one voxel (NOFFSET=0), then force an analysis by
c        setting ZERO_VOXEL(1)=F (and let LCModel abort if only voxel really
c        is zero).
C     -------------------------------------------------------------------------
      if (noffset .le. 0) zero_voxel(1) = .false.
      voxel1 = .true.
      do 100 ioffset = ioffset_current_in, noffset
         if (.not. voxel1) then
            rewind lraw
            lraw_at_top = .true.
            IF (filh2o .ne. ' ') rewind lh2o
         end if
         if (ioffset .gt. ioffset_current_in) nvoxels_done_in = 0
         ivoxel = 0
         do 110 idslic = 1, ndslic
            do 120 idrow = 1, ndrows
               do 130 idcol = 1, ndcols
                  ivoxel = ivoxel + 1
                  ir = iabs(idrow - irow_center)
                  ic = iabs(idcol - icol_center)
                  iok = (ir .eq. ioffset   .and.   ic .le. ioffset) .or.
     1                  (ir .le. ioffset   .and.   ic .eq. ioffset)
                  skip_voxel = .not.iok   .or.
     1                         idrow .lt. irowst   .or.
     2                         idrow .gt. irowen   .or.
     3                         idcol .lt. icolst   .or.
     4                         idcol .gt. icolen   .or.
     5                         idslic .ne. islice   .or.
     6                         zero_voxel(ivoxel)   .or.
     7                         ivoxel .le. nvoxels_done_in
                  do 140 j = 1, nvoxsk
                     skip_voxel = skip_voxel   .or.
     1                            (idrow .eq. irowsk(j)   .and.
     2                             idcol .eq. icolsk(j))
 140              continue
                  call restore_settings ()
C
C                 -------------------------------------------------------------
c                 Open output files.
C                 -------------------------------------------------------------
                  call open_output ()
c
C                 -------------------------------------------------------------
                  if (.not.skip_voxel) then
c                    ----------------------------------------------------------
C                    Load changes to Control Variables in array CHANGE for
c                       later output.
C                    ----------------------------------------------------------
                     CALL LOADCH ()
c
C                    ----------------------------------------------------------
C                    Initialize global quantities for later use.
C                    DELPPM(JY) = PPM - PPMCEN on the frequency grid.
C                    ----------------------------------------------------------
                     CALL INITIA ()
c
                  end if
C                 -------------------------------------------------------------
C                 Get DATAT = raw time-domain data.
C                 -------------------------------------------------------------
                  if (.not.bascal) CALL DATAIN ()
                  voxel1 = .false.
                  lraw_at_top = .false.
                  if (skip_voxel) go to 130
 
                  if (lcsi_sav_1 .eq. 12) then
                     nanalyses_done = nanalyses_done + 1
                     nvoxels_done = ivoxel
                     ioffset_current = ioffset
                     rewind 12
                     write (12, 5030, err= 143)
     1                  nvoxels_done, nanalyses_done, ioffset_current
                     go to 145
 143                 call errmes (1, 4, chsubp)
                  end if
 145              initialize_solve = .true.
C                 -------------------------------------------------------------
C                 BASIST(JDATA,JMETAB) = time-domain basis vectors (COMPLEX).
C                 -------------------------------------------------------------
                  CALL MYBASI (1)
c
C                 -------------------------------------------------------------
C                 Compute NCOMPO and LCOMPO for combinations of metabolites.
C                 -------------------------------------------------------------
                  CALL COMBIS ()
C
C                 -------------------------------------------------------------
C                 Preliminary analysis to get starting estimates for the phase
C                   corrections and the referencing shift in the
c                   frequency-domain data.
C                 -------------------------------------------------------------
                  CALL STARTV (1)
C                 -------------------------------------------------------------
C                 Repeat Prel with fewer metabolites.  Typically, L20 & L09
c                    are omitted in CHECK_CHLESS.  This is skipped when
c                    NCHLES <= 0 (default in Block Data).
C                 -------------------------------------------------------------
                  call check_chless ()
                  if (omit_chless) then
                     CALL MYBASI (1)
                     CALL COMBIS ()
                     CALL STARTV (2)
                  end if
C
C                 -------------------------------------------------------------
C                 Analyses with Regula Falsi searches for ALPHAB and ALPHAS
C                   with PROB1 between PRMNMX(1,1) and PRMNMX(2,1).
C                 -------------------------------------------------------------
                  IF (DOFULL) THEN
                     CALL MYBASI (2)
                     CALL COMBIS ()
                     CALL TWOREG ()
                  END IF
C
C                 -------------------------------------------------------------
C                 Final output.
C                 -------------------------------------------------------------
                  CALL FINOUT ()
                  if (.not.single_voxel) call update_priors ()
 130           continue
 120        continue
 110     continue
 100  continue
      if (lcsv .gt. 0) close (lcsv)
      if (lcsi_sav_1 .eq. 12) then
         rewind 12
         nanalyses_done = 0
         write (12, 5030) nvoxels_done, nanalyses_done, ioffset_current
         close (12)
         close (13)
      end if
      stop
 801  call errmes(1, -4, chsubp)
 802  call errmes(2, -4, chsubp)
 803  call errmes(3, -4, chsubp)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BLOCK DATA
      INCLUDE 'lcmodel.inc'
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DATA RRANGE/1.E37/
c     -----------------------------------------------------------------------
c     DRANGE = SQRT(big number), because it can be squared.
c     -----------------------------------------------------------------------
      DATA DRANGE/1.D153/
      DATA CCNTRL/.FALSE./
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
c     ------------------------------------------------------------------------
c     The ordering of the Metabolite Names for the Cho & Cr families in
c       CHCOMB below must be that in COMBIS and in the tests for NAMREL in
c       FINOUT.
c     ------------------------------------------------------------------------
      DATA change/mlines*' ', mlines*' '/,
     1     chbcal/' '/, chcali/mmetab*' '/, chcol/'-'/,
     1     CHCOMB/'GPC+PCh', 'GPC+Cho', 'PCh+Cho', 'Cho+GPC+PCh',
     1            'NAA+NAAG', 'Ins+Glyc', 'mI+Glyc',  'Ins+Gly',
     1            'mI+Gly', 'Cr+PCr', 'Cre+PCr',  'Glu+Gln',
     1            'Lip13a+Lip13b', 'MM14+Lip13a+Lip13b+MM12',
     1            'MM09+Lip09',    'MM20+Lip20', MPMET16*' '/,
     1     chcom2/mpmet*' '/, chdate/' '/,
     1     CHEXT2/MMETAB*' '/, chgam/'gamma'/, chgrsh/mgroup_shift*' '/,
     1     CHKEEP/MMETAB_extra*' '/, chless/'L09', 'L20', mmeta2*' '/,
     1     chlsha/mmetab*' '/, chmore/'L13'/,
     1     chnols/'Lac', 'Ala', mmet_ex2*' '/,
     1     chnot1/'Acn',   'Act',   'Bet',   'bHb',   'Car',
     2            'Cit',   '-CrCH2','Eth',   'Fuc',   'GABA',
     3            'Gcn',  'Gcr',    'Glc',   'Gly',   'Glyc',
     4            'Gua',  'Ilc',    'Leu',   'Lip20', 'MM12',
     5            'MM14', 'MM17',   'MM20',  'NAAG',  'Pal',
     6            'Pgc',  'Pyr',    'Scyllo','Suc',   'notTau',
     7            'Thr',  'TMPO',   'Val',   mmet_ex33*' '/,
     1     chnot2/'Lip13c', 'Lip13d', 'Lip13e', mmet_ex3*' '/
     1     CHOMIT/MPMET*' '/
      DATA CHPMET/
     1            'Acn',    'Act',    'Ala',
     2            'Asp',    'Bet',    'bHb',
     3            'Car',    'Cho',    'Cit',
     4            'Cr',     'Cys',    'Eth',
     5            'Fuc',    'GABA',   'Gcn',
     6            'Gcr',    'Glc',    'Gln',
     7            'Glu',    'Glyc',   'GPC',
     8            'Gua',    'ILc',    'Ins',
     9            'Lac',    'Leu',    'NAA',
     T            'NAAG',   'PAl',    'PCh',
     1            'Pgc',    'Pyr',    'Scyllo',
     2            'Suc',    'Tau',    'Thr',
     3            'TMPO',   'Val',    'PCr',
     4            'Cre',    'Gly',    'mI',
     5            '-CrCH2', 'Lip20',  'MM12',
     6            'MM14',   'MM17',   'MM20',
     7            'Lip13a', 'Lip13b', 'Lip13c',
     8            'Lip13d', 'Lip13e', 'MM09',
     9            'Lip09',  MPME55*' '/
c     -------------------------------------------------------------------------
c     If CHRATO is changed below, then NORATO, NRATIO, etc, may have to be
c       adjusted in assignments according to SPTYPE in MYCONT and liver-1.inc
c       and lipid-1.inc.
c     .286+-.126 below corresponds to GSH/Gln=.4+-.3
c     -------------------------------------------------------------------------
 
      DATA chrati/mmetab*' '/,
     1     chrato/
     1            'Lip09/Lip13* = .267 +- .1 +WT= MM12',
     1            'Lip20/Lip13* = .15 +- .07 +WT= MM12',
     1            'MM20/MM09* = 1.5 +- .375 +WT= Lip09*',
     1            'MM12/MM09* = .3 +- .1 +WT= Lip09*',
     1            'MM14/MM09* = .75 +- .45 +WT= Lip09*',
     1            'MM17/MM09* = .375 +- .3 +WT= Lip09*',
     1            '-CrCH2/totCr = .1 +- .25',
     1            'Asp/Big3 = .05 +- .05',
     1            'GABA/Big3 = .03 +- .03',
     1            'Glc/Big3 = .03 +- .03',
     1            'Scyllo/Big3 = .01 +- .01',
     1            'Tau/Big3 = .05 +- .05',
     1            'GSH/Gln+GSH = .286 +- .126',
     1            mmet13*' '/,
     1     chratw/mmetab*' '/,
     1     chrow/'_'/, CHSDSH/'NAA', 'NAAG', 'Cho', MMETA3*' '/,
     1     CHSDT2/MMETAB*' '/, chsim/mmetab*' '/
      data chsimu/
     1          'Lip13a @ 1.28 +- .01 FWHM= .15 < .2 +- .035 AMP= 2.',
     1          'Lip13b @ 1.28 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.',
     1          'Lip13c @ 1.30 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.',
     1          'Lip13d @ 1.26 +- .01 FWHM= .089 < .09 +- .035 AMP= 2.',
     1          'Lip09 @ .89 +- .02 FWHM= .14 < .19 +- .035 AMP= 3.',
     1          'MM09 @ .91 +- .02 FWHM= .14 < .17 +- .015 AMP= 3.',
     1          'Lip20 @ 2.04 +- .005 FWHM=.15 < .2 +- .025 AMP=1.33
     1             @ 2.25 FWHM=.15 AMP=.67  @ 2.8 FWHM=.2 AMP=.87',
     1          'MM20 @ 2.08 +- .005 FWHM=.15 < .18 +- .01 AMP=1.33
     1             @ 2.25 FWHM=.2 AMP=.33  @1.95 FWHM=.15 AMP=.33
     1             @ 3. FWHM=.2 AMP=.4',
     1          'MM12 @ 1.21 +- .01 FWHM= .15 < .2 +- .02 AMP= 2.',
     1          'MM14 @ 1.43 +- .02 FWHM= .17 < .2 +- .02 AMP= 2.',
     1          'MM17 @ 1.67 +- .03 FWHM= .15 < .17 +- .02 AMP= 2.',
     1          '-CrCH2 @ 3.93 +- 0. FWHM= -9. < 0. +- 0. AMP= -2.',
     1          'Gua @ 3.78 +- 0. FWHM=-9.<0. +- 0. AMP=2.',
     1          'Glyc @ 3.55 +- 0. FWHM=-9.<0. +- 0. AMP=2.',
     1          mmet14*' '/,
     1     chslic/'sl'/,
     1     CHUSE1/'NAA', 'Cr', 'GPC', 'Glu', 'Ins', 'Lac', 'Glc',
     1            'Gln', 'NAAG', 'Tau', 'Asp', 'Ala', 'GABA', 'Scyllo',
     1            MMET_ex14*' '/,
     1     FILBAS/' '/,  FILCOO/' '/, filcor/' '/,
     1     filcsi_sav_1/' '/, filcsi_sav_2/' '/, filcsv/' '/,
     1     FILH2O/' '/, FILPS/' '/, FILPRI/' '/, FILRAW/' '/,
     1     FILTAB/' '/, NAMEAC/MMETAB*' '/, NAMREL/'Cr'/,
     1     norato/mmetab*' '/,
     1     OWNER/' '/, ownout/' '/, PGNORM/' '/, savdir/' '/,
     1     sptype/' '/, srch2o/' '/, srcraw/' '/,
     1     SYNUS1/'GPC','PCh',  'GPC','Cho',  'PCh','Cho',
     1            'Ins','mI',   'Glyc','Gly', 'Cr','Cre',
     1            MPMET6*' ', MPMET6*' '/
     1     TITLE/' '/, title_line/2*' '/, wsmet/'Cr'/
      data power/2.d0, 1.d0, 1.5d0/
      DATA iareaw/2/,
     3     iauto/1/, iaverg/0/, icolen/1/, icolsk/mvoxsk*0/, icolst/1/,
     3     idgppm/-1/, IDUMP/MSTAGE*1/,
     3     IETCOU/3/, imethd/0/, INCSID/1/, INCSMX/15/,
     3     ioffset_current_in/0/,
     3     IPAGE2/1/, ipdump/0/, ipowph/6/, ipowrg/0/, irowen/1/,
     3     irowsk/mvoxsk*0/,  irowst/1/, ISDBOL/15/,
     3     ISHIFW/0/, islice/1/,
     3     ISTAGO/0/, iter_dump/-1/,
     3     KEY/MKEY*0/, LBASIS/3/, LCONTR/5/, lcontr_scratch/29/,
     3     LCOORD/0/, lcoraw/0/, lcsi_sav_1/0/, lcsi_sav_2/0/, lcsv/0/,
     3     LDWFFT/0/,  lett/0/, LH2O/2/, LINCHG/2*0/, LINERR/0/,
     3     linerr_mycont/0/,
     3     LINETC/0/, LINTBL/0/, LPRINT/0/, LPS/8/, LRAW/1/,
     3     LTABLE/0/, LWFFT/0/,
     3     MDALPB/12/, mdegp3/41/,
     3     MERMES/9000/, MFNDAL/25/, MINTER/3*1/,
     3     MITER/20, 40, 30/, mnsamp/9/, mpower/1/,
     3     MREPHA/3, 0/, n1hmet/3/, nback/75, 50/, nbas_ccf/10/,
     3     nbckmn/6/, ncalib/0/, nchgam/5/, nchles/0/, nchlin/35, 49/,
     3     NCOMBI/16/, ndcols/1/, NDEGZ/36, 4/,
     3     NDGPPM/21, 3/, ndrows/1/, ndslic/1/,
     3     NEACH/0/, nermes/0/, NERROR/MMERM9*0/,
     3     nextr/2*0/, NEXT2/0/, ngau/mmetab*1/, ngrsh/0/,
     3     ninfl/2*0/, NKEEP/0/, nlin/0/,
     3     nlshap/0/, nnolsh/0/,
     3     nnot1/33/, nnot2/3/, nnorat/0/, NOMIT/0/,
     3     npar/0/, npower/mmetab*0/,  nratio/13/,
     3     nratio_used/0/, NREFPK/1, 3/,
     3     NRF2MN/2/, NSDSH/3/, NSDT2/0/, NSHIFT/8/, NSIDMN/5/,
     3     NSIDMX/11/, nsimul/13/, NSUBTK/2/, ntitle/2/,
     3     NUNFIL/2048/, NUSE1/5/, nvoxsk/0/,
     3     nwsend/50/, nwsst/10/
      DATA absval/.false./, accept_alpbmn/.false./,
     4     accept_step2/.false./,
     4     areaba_orig_basisf/.true./, asymlp/.false./,
     4     badref/.false./, bascal/.false./,
     4     basout/.false./, biglip/.false./,
     4     chksim/.true./, conc3f/.true./,
     4     DOECC/.FALSE./, DOFULL/.TRUE./,
     4     DOREFS/.FALSE., .TRUE./,
     4     dowatr/.true./, dows/.false./,
     4     DOZERO/2*.FALSE., .TRUE./, eccdon/.false./, endpha/.false./,
     4     fixshf/.false./, forecc/.false./,
     4     gauss_rt2/.true./, gshgua/.true./, LANDSC/.TRUE./,
     4     ldump/10*.false./, lcy_skip/my*.false./,
     4     nobase/.false./, nobasi/.false./,
     4     onlyco/.false./, QUICK/.FALSE./,
     4     reflac/.FALSE./, roomt/.false./, scafwh/.false./,
     4     scasim/.true./, sidump/mmetab*.false./,
     4     sitayl/mmetab*.false./, skip_step3/.false./,
     4     solgrd/.true./, smtail/.true./,
     4     SUBBAS/.FALSE./, unsupr/.false./,
     4     useany/.FALSE./, useglc/.FALSE./, usemxb/.true./,
     4     USINFL/.TRUE./, VITRO/.FALSE./, wsdone/.false./,
     4     year4d/.true./
      DATA ALEXT2/MMETAB*0./, ALPBMN/1.7D-3/, ALPBMX/.84D0/,
     5     ALPBPN/.0021/,
     5     ALPBST/.0025/, alphab_dump/-1./, alphas_dump/-1./,
     5     ALPSMN/.1D0/, ALPSMX/2.D3/, ALPSST/10.0001/,
     5     ALSDSH/2*.002, .008, MMETA3*0./, ALSDT2/MMETAB*0./,
     5     area_met_norm/0./,
     5     atth2o/.7/, attmet/1./, black/3*0./, bwtolr/.001/,
     5     conc_expect/mmetab*0./,
     5     CONREL/1./, COSMIN/3*1.E-2/, ddegp3/5./, DDGPMQ/1., 1./,
     5     DDGZMQ/3., 3./, DEEXT2/2./, DEGMAX/12.5, 13./,
     5     DEGPPM/0./, DEGZER/0./, DELTAT/5.E-4/, DESDSH/.004/,
     5     DESDT2/.4/, DFLDMQ/.1/, DGPPMN/-30./, DGPPMX/30./,
     5     dkngam/.35/, DKNTMN/.15, .6/, dkntmn_standard/.075/,
     5     DSHPAT/.05,.1, 2*0./, echot/-1./, exrati/mmetab*-1./,
     5     EXRT2/MMETAB*2./, fcalib/1./, fcsum/1.e-3/, FH2OMX/.25/,
     5     fmain_power/.6/, fother_power/.1/,
     5     frepha/.5/, FSTPMQ/5./, FWHH2O/.4/, FWHMBA/.013/,
     5     fwhmmn/.1/, fwhmmx/.15/,
     5     fwhmsm/-1./, FWHMST/.04/, hifmm/190./,
     5     hwdwat/1., 2./, hzpgam/90., 210./,
     5     HZPPPM/84.47/, HZREF/MREFPK*0., MREFPk*0./,
     5     PAGEHT/27.9/, PAGEWD/21./, PHITOT/2*0./,
     5     PMQST/3*1./, PMQSTL/3*.6/, PNALPB/9./, ppmbas/.1, .2/,
     5     PPMCEN/4.65/, PPMEND/9999./, ppmend_phalip/-1./,
     5     ppmgap/mgap*1.e37, mgap*1.e37/, ppmh2o/4.65/
c     -------------------------------------------------------------------------
c          With PPMMET below:
c     -CrCH2 is omitted if the CH3 peak around 3ppm is not included.
c     MM12 MM14 MM17 MM20 MM09 require MM09 and baseline down to 0.61.
c        However, these are relaxed as much as possible in MYBASI if
c        TE <= TEMM (typically 99.99), with error message if PPMEND > .6, so
c        that they (especially MM20) can fit under NAA peak, etc.
c
c     Lip* are not excluded if PPMEND>0.6, since (especially mobile) lipids
c        can still be present, even at long TE.
c     Changes below must be also made in MYBASI.
c     -------------------------------------------------------------------------
      DATA PPMMET/
     1 2.27,2.17,2.27,2.17,  1.97,1.87,1.97,1.87,  1.57,1.29,1.57,1.29,
     2 2.91,2.59,2.91,2.59,  3.31,3.21,3.31,3.21,  1.26,1.09,1.26,1.09,
     3 3.71,2.35,3.71,2.35,  3.25,3.22,3.22,3.19,  2.71,2.49,2.71,2.49,
     4 3.92,3.89,3.03,3.00,  3.71,3.05,3.71,3.05,  1.29,1.01,1.29,1.01,
     5 1.31,1.09,1.31,1.09,  3.17,2.10,3.17,2.10,  3.81,3.59,3.81,3.59,
     6 3.99,3.09,3.99,3.09,  3.84,2.99,3.84,2.99,  3.91,3.59,2.65,2.01,
     7 3.84,3.59,2.51,2.26,  3.84,3.49,3.84,3.49,  3.25,3.22,3.22,3.19,
     8 3.83,3.73,3.83,3.73,  1.09,0.79,1.09,0.79,  3.59,3.56,3.56,3.51,
     9 1.41,1.21,1.41,1.21,  1.81,0.79,1.01,0.79,  2.63,2.55,2.01,2.01,
     T 2.81,2.45,2.31,2.01,  3.99,2.99,3.99,2.99,  3.25,3.22,3.22,3.19,
     1 1.31,1.01,1.31,1.01,  2.43,2.33,2.43,2.33,  3.40,3.30,3.40,3.30,
     2 2.45,2.35,2.45,2.35,  3.61,3.06,3.61,3.06,  1.41,1.19,1.41,1.19,
     3 3.95,3.85,3.78,3.68,  1.21,0.79,1.21,0.79,  3.92,3.89,3.03,3.00,
     4 3.92,3.89,3.03,3.00,  3.84,3.49,3.84,3.49,  3.59,3.56,3.56,3.51,
     5 3.03,3.00,3.03,3.00,  2.31,1.11,2.31,1.11,  1.25,0.61,1.25,0.61,
     6 1.50,0.61,1.50,0.61,  1.75,0.61,1.75,0.61,  2.31,0.61,2.31,0.61,
     7 1.50,1.11,1.50,1.11,  1.50,1.11,1.50,1.11,  1.50,1.11,1.50,1.11,
     8 1.50,1.11,1.50,1.11,  1.50,1.11,1.50,1.11,  1.20,0.61,1.20,0.61,
     9 1.20,0.61,1.20,0.61,
     T MPME55*0., MPME55*0.,  MPME55*0., MPME55*0./
      DATA PPMPOS/-1.E6, 1.E6/,
     5     PPMREF/4.65, MREFP1*0.,
     5            2.01, 3.03, 3.22, 3.56, MREFP4*0./,
     5     ppmsep/mgap*1.e37/,
     5     PPMSHF/1.E37/, ppmsig/1.e9, -1.e9/, PPMST/-9999./,
     5     ppmst_phalip/7./,
     5     ppm_truncate_max/-999./, ppm_truncate_min/999./,
     5     ppm_water_range/1./, ppm_water_tol/.08/,
     5     PRMNMX/.02,.08, .2,.5, .02,.08/, PTLABL/7.8/,
     5     PTOUTP/7.8/, PTTITL/11./, PTVERS/9./, RALIMN/1.02/,
     5     RALINC/2./, r_areaba/20./, ratipm/4./, RBACKG/6., 2.5/,
     5     rbasmx/2*.25/, RCONVR/3*1.E-3/,
     5     RDALPB/2./, rfwbas/10./, RFWHCC/.5/, RFWHM/1.8/,
     5     rfwhmst_ccf/1.5/, RFWHST/.5/,
     5     RGBBOL/2*0., .999/, rgberr/.999, 2*0./,
     5     RGBLIN/.999, 17*0./, rgbrat/.999, 2*0./,
     5     RHLABL/1.2/, RHOUTP/1.3/, RHTITL/1.5/, RHVERS/2./,
     5     RINCRS/3*1.E-5/, rincsh/.5/, rlesmo/.35/, rlrntz/1./,
     5     RMQDEC/3*.5/, RMQINC/6*4./, RPENMX/1.00001/,
     5     RPMQMN/3*1./, rpowmq/.1/, RRT2MQ/.1/, rsdgp3/1.3/,
     5     rsdsam/2.5, 1.5, 3./, rsdsmq/.125/,
     5     RSHFMQ/.025/, RSTPMN/3*0./,
     5     RSTPMX/3*1./, rt2min/mmetab*0./,
     5     RWFONT/.62/, SDDEGP/20./, SDDEGZ/999./,
     5     sdgrsh/mgroup_shift*0./, SDMSHF/.02/, sdrati/mmetab*0./,
     5     sdshmn/.002/, sdshmx/.004/, SDSMOO/.01, .02, .03, .01/,
     5     SHIFMN/-.3, -.1/, SHIFMX/.3, .3/, siamp/mmet_mgau*0./,
     5     sifwex/mmetab*0./, sifwmn/mmet_mgau*0./, sifwsd/mmetab*0./,
     5     sippm/mmet_mgau*999./, sisdsh/mmetab*-1./,
     5     temm/99.99/, THRLIN/.1/,
     5     wconc/35880./, WDLINE/.06, .01, .005, .01, .04, .005/,
     5     wsppm/3.027/, XLEFT/1.3/,
     5     XRIGHT/1.3/, XSTEP/.2/, XTRPMX/2./, YBOTT/1.3/, YTOP/1.5/
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      subroutine y2k !OSF
C      INCLUDE 'lcmodel.inc' !OSF
C      character chstr*40 !OSF
c
C      CALL IDATE (KNUM(2), KNUM(1), KNUM(3)) !OSF
C      CALL DATE (CHSTR) !OSF
C      if (year4d) then !OSF
C         chdate(1:7) = chstr !OSF
C         if (knum(3) .lt. 99) then !OSF
C            knum(3) = knum(3) + 2000 !OSF
C         end if !OSF
C         write (chstr, 5100) knum(3) !OSF
C 5100    format (i4) !OSF
C         chdate(8:11) = chstr !OSF
C         CALL TIME (CHSTR) !OSF
C         CHDATE(15:19)=CHSTR !OSF
C      else !OSF
C         CHDATE(1:9)=CHSTR !OSF
C         CALL TIME (CHSTR) !OSF
C         CHDATE(13:17)=CHSTR !OSF
C      end if !OSF
C      end !OSF
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      subroutine y2k !IRIX
C      INCLUDE 'lcmodel.inc' !IRIX
C      character chstr*40 !IRIX
c
C      CALL IDATE (KNUM(2), KNUM(1), KNUM(3)) !IRIX
C      CALL DATE (CHSTR) !IRIX
C      if (year4d) then !IRIX
C         chdate(1:7) = chstr !IRIX
C         if (knum(3) .lt. 99) then !IRIX
C            knum(3) = knum(3) + 2000 !IRIX
C         end if !IRIX
C         write (chstr, 5100) knum(3) !IRIX
C 5100    format (i4) !IRIX
C         chdate(8:11) = chstr !IRIX
C         CALL TIME (CHSTR) !IRIX
C         CHDATE(15:19)=CHSTR !IRIX
C      else !IRIX
C         CHDATE(1:9)=CHSTR !IRIX
C         CALL TIME (CHSTR) !IRIX
C         CHDATE(13:17)=CHSTR !IRIX
C      end if !IRIX
C      end !IRIX
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MYCONT ()
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C  Inputs changes to Control Variables.
C  If your compiler doesn't accept NAMELIST, then you will have to modify this
C    subprogram to use only standard I/O statements.
C  Below is illustrated how you can make your own special-purpose
c    modifications to input data.  In the example below, if QUICK has been
C    input .TRUE., then other Control Variables are changed to force the
C    fastest possible (but crude) analysis.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
      INCLUDE 'lcmodel.inc'
      external ilen
      parameter (mch_line_long=mchsimul+12)
      CHARACTER LINE*(mch_line_long)
      logical is_sptype, standard_refs
      include 'nml_lcmodl.inc'
      include 'nml_lcmodel.inc'
      CHSUBP='MYCONT'
C     -------------------------------------------------------------------------
c     Preliminary default settings based on preliminary reading of Namelist.
c     SPTYPE = 'version5 will be used after 2nd Namelist read (to overwrite
c                          earlier settings)
c                        These settings are dependent on the Block data
c                          settings of CHRATO, NRATIO, etc.
c     Copy LCONTR to LCONTR_SCRATCH, since Cygwin g77 (& NAG90) cannot rewind
c       STDIN.
C     -------------------------------------------------------------------------
      open(lcontr_scratch, status='scratch')
 5104 format(a)
      do 104 jline = 1, 9999
         read(lcontr, 5104, end=105) line
         llen = ilen(line)
         write(lcontr_scratch, 5104) line(1:llen)
 104  continue
 105  rewind lcontr_scratch
      READ (LCONTR_SCRATCH, NML=LCMODL, end=801, err=802)
c     -------------------------------------------------------------------------
c     At long TE, delete ratio prior for GSH/(Gln+GSH), which can get larger
c        than normal.
c     -------------------------------------------------------------------------
      if (echot .gt. 49.) nratio = 12
c     -------------------------------------------------------------------------
c     For high fields, reset Lip* & MM* to sharper values.
c     Also remove ratio constraint on Tau (assuming rodent).
c     This assumes that HZPPPM or HIFMM will not be set in any of the *.inc
c       files below.
c     -------------------------------------------------------------------------
      if (hzpppm .gt. hifmm) then
         chsimu(1)= 'Lip13a @ 1.28 +- .01 FWHM= .05 < .07 +- .02
     1               AMP= 2.'
         chsimu(2)= 'Lip13b @ 1.28 +- .01 FWHM= .029 < .03 +- .02
     1               AMP= 2.'
         chsimu(3)= 'Lip13c @ 1.30 +- .01 FWHM= .029 < .03 +- .02
     1               AMP= 2.'
         chsimu(4)= 'Lip13d @ 1.26 +- .01 FWHM= .029 < .03 +- .02
     1               AMP= 2.'
         chsimu(5)= 'Lip09 @ .89 +- .02 FWHM= .05 < .07 +- .02
     1               AMP= 3.'
         chsimu(6)= 'MM09 @ .91 +- .02 FWHM= .05 < .06 +- .01
     1               AMP= 3.'
         chsimu(7)= 'Lip20 @ 2.04 +- .005 FWHM=.05 < .07 +- .02
     1               AMP=1.33  @ 2.25 FWHM=.05 AMP=.67
     2               @ 2.8 FWHM=.07 AMP=.87'
         chsimu(8)= 'MM20 @ 2.08 +- .005 FWHM=.05 < .06 +- .01
     1               AMP=1.33  @ 2.25 FWHM=.07 AMP=.33
     2               @1.95 FWHM=.05 AMP=.33  @ 3. FWHM=.07 AMP=.4'
         chsimu(9)= 'MM12 @ 1.21 +- .01 FWHM= .05 < .07 +- .02 AMP= 2.'
         chsimu(10)= 'MM14 @ 1.43 +- .02 FWHM= .06 < .07 +- .02 AMP= 2.'
         chsimu(11)= 'MM17 @ 1.67 +- .03 FWHM= .05 < .06 +- .02 AMP= 2.'
         nnorat= 1
         norato(1)= 'Tau'
      end if
      if (nunfil .lt. 64) call errmes (12, 4, chsubp)
      call remove_blank_start (sptype)
      if (nchgam .gt. 0   .and.   chgam .ne. ' '   .and.
     1    dkngam .gt. 0.   .and.   nchgam .le. mchgam   .and.
     2    sptype .eq. ' ') then
         if (index(filbas, chgam(1:nchgam)) .gt. 0   .and.
     1       hzpppm .ge. hzpgam(1)   .and.  hzpppm .le. hzpgam(2))
     2      dkntmn(1) = dkngam
      end if
      call toupper_lower (.false., sptype)
      is_sptype= .false.
      if (sptype(:3) .eq. 'csf') then
         is_sptype= .true.
         nnorat = 5
         norato(1) = 'Asp'
         norato(2) = 'GABA'
         norato(3) = 'Glc'
         norato(4) = 'Scyllo'
         norato(5) = 'Tau'
         nratio = 13
         reflac = .true.
         useglc = .true.
      end if
      if (sptype(:6) .eq. 'nulled') then
         is_sptype= .true.
         badref = .true.
         incsmx = 1
         namrel = 'Lip13a+Lip13b'
         nobasi = .true.
         nratio = 6
         nrefpk(2) = 2
         nsidmn = 1
         nsidmx = 1
         nsimul = 11
         ppmref(1,2) = 1.28
         ppmref(2,2) = 0.90
         vitro = .true.
      end if
      if (sptype(:5) .eq. 'tumor') then
c        ----------------------------------------------------------------------
c        BADREF=T will use every metabloite in Prel, but only those with
c           NUSE1 & CHUSE1 will be used to determine Prel FWHMST, which
c           determines NSIDES.  Must have NUSE1>0 to determine NSIDES.  Can
c           only have NUSE1=0 when lineshape is not used, as with liver, etc.,
c           where NOBASI=T.
c        CHRATO & NRATIO must be adjusted below if there is any change of
c           these in BLOCK DATA.
c        ----------------------------------------------------------------------
         is_sptype= .true.
         badref = .true.
         chrato(13) = 'NAAG/NAA = .15 +- .15'
         chuse1(1) = 'GPC'
         chuse1(2) = 'Cr'
         dkntmn(1) = .35
         namrel = 'GPC'
         nrefpk(2) = 4
         nratio = 13
         nuse1 = 2
         ppmend = .2
         ppmref(1,2) = 3.03
         ppmref(2,2) = 3.22
         ppmref(3,2) = 1.28
         ppmref(4,2) = 0.90
         rbackg(1) = 12.
         shifmn(2) =  -.07
         shifmx(2) = .07
      end if
      include 'muscle-1.inc'
      include 'liver-1.inc'
      include 'lipid-1.inc'
      if (sptype(:7) .eq. 'muscle-'   .or.
     1    sptype(:6) .eq. 'liver-'   .or.
     2    sptype(:7) .eq. 'breast-'   .or.
     3    sptype(:6) .eq. 'lipid-') then
         biglip = .true.
         ppm_truncate_max = 3.5
         ppm_truncate_min = 2.2
      end if
      if (.not.(is_sptype   .or.
     1          sptype(:8) .eq. 'version5'   .or.
     2          sptype(:9) .eq. 'version-5'   .or.
     3          sptype(:1) .eq. ' ')) call errmes (8, 4, chsubp)
c     -------------------------------------------------------------------------
c     Overwrite above defaults by rereading Namelist
c     -------------------------------------------------------------------------
      REWIND LCONTR_SCRATCH
      fwhmba_sav = fwhmba
      fwhmba = -.1
      READ (LCONTR_SCRATCH, NML=LCMODL, end=801, err=802)
c     -------------------------------------------------------------------------
c     Red herring; nothing will be done here.
      if (ldwfft .gt. 0) nlin = nlin + 1024
c     -------------------------------------------------------------------------
c     DOECC must be set before call to AVERAGE
c     -------------------------------------------------------------------------
      doecc_active = doecc
      if (biglip) doecc = forecc
      call remove_blank_start (sptype)
      call toupper_lower (.false., sptype)
      if (sptype(:8) .eq. 'version5'   .or.
     1    sptype(:9) .eq. 'version-5') then
         nratio = 0
         nsimul = 0
         ppmend = 1.0
         ppmst = 3.85
      end if
      if (iauto .eq. 1   .and.   ppmend .gt. 9998.) then
         if (echot .ge. 100.) then
            ppmend = 1.
            nuse1 = 4
         end if
c        ---------------------------------------------------------------------
c        No longer truncate CSI analyses with PPMEND.
c         if (max0(ndslic, ndrows, ndcols) .gt. 1) ppmend = 1.8
c        ---------------------------------------------------------------------
      end if
      if (ppmst .lt. -9998.) ppmst= 4.
      if (ppmend .gt. 9998.) ppmend= .2
      fwhmba_in_control = fwhmba .gt. 0.
      if (.not.fwhmba_in_control) then
         fwhmba = fwhmba_sav
      end if
      if (sptype(:7) .eq. 'muscle-'   .or.
     1    sptype(:6) .eq. 'liver-'   .or.
     2    sptype(:7) .eq. 'breast-'   .or.
     3    sptype(:6) .eq. 'lipid-') then
         if (ppmend .ge. -0.9) call errmes(17, 2, chsubp)
         if (ppmst .ge. 5.) then
            if (ppmst .le. 7.9) call errmes(18, 2, chsubp)
         else
            if ((sptype(:6) .eq. 'liver-')   .and.
     3          (ppmst .ge. 4.01   .or.   ppmst .le. 3.59))
     4         call errmes(19, 2, chsubp)
            if ((sptype(:7) .eq. 'breast-')   .and.
     3          (ppmst .le. 3.79   .or.   ppmst .ge. 4.01))
     4         call errmes(20, 2, chsubp)
            if ((sptype(:6) .eq. 'lipid-')   .and.
     2          (ppmst .le. 3.39   .or.   ppmst .ge. 4.01))
     3         call errmes(21, 2, chsubp)
         end if
      end if
      if ((sptype(:10) .eq. 'only-cho-1'   .or.
     1     sptype(:10) .eq. 'only-cho-2')   .and.
     2    (ppmst .le. 3.79   .or.   ppmst .ge. 4.01   .or.
     3     ppmend .ge. 2.81   .or.   ppmend .le. 2.59))
     4   call errmes(22, 2, chsubp)
c     -------------------------------------------------------------------------
c     Setup for no baseline.
c     -------------------------------------------------------------------------
      if (nobase) then
         alpbmx = alpbmn
         alpbst = alpbmn
         alpbpn = 1.0001 * alpbmn
         idgppm = -1
         nbackg = 0
         usemxb = .false.
      end if
c     -------------------------------------------------------------------------
c     Reset RBASMX & RSDGP3.
c     IDGPPM = -1 (default in Block Data) to skip fixed-DEGPPM series and only
c                  do free analysis.
c            = 0 (default in lipid-2 & PPMST<5) for analyses stressing flat
c                baselines (by allowing medium SSQ/SSQMIN), where there is no
c                small peak that can be lost by incorrect DEGPPM.  This can
c                also be used to input any values for RBASMX & RSDGP3; only
c                differences between IDGPPM=0 & >0 are settings below.
c            = 1 (default in breast-2, only-cho-1 & liver-2) for analyses
c                stressing good phasing and fit (at expense of flat
c                baselines), to find small peaks (like Cho).  Forces
c                fixed-DEGPPM series, even if baseline in unconstrained
c                solution is flat.
c            = 2 (default in muscle-2) (actually between 0 and 1) stresses
c                good fit, but does not force fixed-DEGPPM series if
c                unconstrained solution has flat baseline.
c     -------------------------------------------------------------------------
      if (idgppm .gt. 0) rsdgp3 = amin1(rsdgp3, 1.1)
      if (idgppm .eq. 1) rbasmx(1) = 0.
c     -------------------------------------------------------------------------
c     With dongle, omit "Data of:" allowing blank line if nothing is input (by
c       Hitachi).
c     -------------------------------------------------------------------------
      OWNOUT='Data of: ' // OWNER
       ownout = owner!Cyg
      fwhmst = amin1(fwhmst, fwhmmx)
c     -------------------------------------------------------------------------
c     USEGLC = T to add Glc to Preliminary Analysis
c     -------------------------------------------------------------------------
      if (useglc   .and.   max0(nuse1, nkeep) .lt. mmetab_extra) then
         nkeep = max0(1, nkeep + 1)
         chkeep(nkeep) = 'Glc'
         do 110 juse1 = 1, nuse1
            if (chuse1(juse1) .eq. 'Glc') go to 115
 110     continue
         nuse1 = nuse1 + 1
         chuse1(nuse1) = 'Glc'
 115     continue
      end if
c     -------------------------------------------------------------------------
c     REFLAC = T to use Lac for referencing and for the  Preliminary
c                Analysis, provided PPMEND < 1.33.
c     -------------------------------------------------------------------------
      if (reflac   .and.   ppmend .lt. 1.33   .and.
     1    nuse1 .lt. mmetab_extra) then
         do 120 juse1 = 1, nuse1
            if (chuse1(juse1) .eq. 'Lac') go to 125
 120     continue
         nuse1 = nuse1 + 1
         chuse1(nuse1) = 'Lac'
 125     dorefs(2) = .true.
         do 130 j = 1, nrefpk(2)
            if (abs(ppmref(j,2) - 1.33) .le. .01   .and.
     1          abs(hzref(j,2) - 3.6) .le. .1) go to 140
 130     continue
         nrefpk(2) = min0(mrefpk, max0(2, nrefpk(2) + 2))
         do 135 j = nrefpk(2), 3, -1
            ppmref(j,2) = ppmref(j-2,2)
            hzref(j,2) = hzref(j-2,2)
 135     continue
         ppmref(1,2) = 1.33
         ppmref(2,2) = 1.33
         hzref(1,2) = -3.6
         hzref(2,2) = 3.6
 140     continue
      end if
 
      IF (QUICK) THEN
         NDEGZ(2)=3
         NDGPPM(2)=1
         NSHIFT=1
         DOREFS(1)=.FALSE.
         NREFPK(2)=3
         PPMREF(1,2)=2.01
         PPMREF(2,2)=3.03
         PPMREF(3,2)=3.22
         HZREF(1,2)=0.
         HZREF(2,2)=0.
         HZREF(3,2)=0.
         DOFULL=.FALSE.
      END IF
c     -------------------------------------------------------------------------
c     Change CV's for calibration.
c     -------------------------------------------------------------------------
      if (ncalib .gt. 0) then
         nratio = 0
         dorefs(2) = .true.
         ncalib = min0(ncalib, mmetab)
         nuse1 = ncalib
         vitro = .true.
         dkntmn(2) = 99.
         standard_refs = .true.
         do 210 j = 1, ncalib
            chuse1(j) = chcali(j)
            standard_refs = standard_refs  .and.
     1           (chcali(j) .eq. 'Lac'  .or.
     2            chcali(j) .eq. 'NAA'  .or.
     3            chcali(j) .eq. 'Cr'  .or.
     4            chcali(j) .eq. 'Cre'  .or.
     5            chcali(j) .eq. 'GPC'  .or.
     6            chcali(j) .eq. 'PCh'  .or.
     7            chcali(j) .eq. 'Cho')
 210     continue
         if (standard_refs) then
            nrefpk(2) = 0
            do 220 j = 1, ncalib
               if (chcali(j) .eq. 'Lac') then
                  nrefpk(2) = nrefpk(2) + 2
                  ppmref(nrefpk(2)-1, 2) = 1.33
                  ppmref(nrefpk(2), 2) = 1.33
                  hzref(nrefpk(2)-1, 2) = -3.6
                  hzref(nrefpk(2), 2) = 3.6
               else if (chcali(j) .eq. 'NAA') then
                  nrefpk(2) = nrefpk(2) + 1
                  ppmref(nrefpk(2), 2) = 2.01
                  hzref(nrefpk(2), 2) = 0.
               else if (chcali(j) .eq. 'Cr'  .or.
     1                  chcali(j) .eq. 'Cre' ) then
                  nrefpk(2) = nrefpk(2) + 1
                  ppmref(nrefpk(2), 2) = 3.03
                  hzref(nrefpk(2), 2) = 0.
               else if (chcali(j) .eq. 'Cho'  .or.
     1                  chcali(j) .eq. 'GPC'  .or.
     1                  chcali(j) .eq. 'PCh' ) then
                  nrefpk(2) = nrefpk(2) + 1
                  ppmref(nrefpk(2), 2) = 3.22
                  hzref(nrefpk(2), 2) = 0.
               end if
 220        continue
         end if
      end if
C     -------------------------------------------------------------------------
c     Change CV's for Basis calibration when BASCAL=T
C     -------------------------------------------------------------------------
      if (bascal) then
         sddegz=amin1(sddegz, 3.)
         sddegp=amin1(sddegp, 1.)
         doecc=.false.
         absval=.false.
         nsimul = 0
         nratio = 0
         ndslic = 1
         ndrows = 1
         ndcols = 1
      end if
C     -------------------------------------------------------------------------
c     Set IKNTMN DKNTMN & ALPB* for special cases.
C     -------------------------------------------------------------------------
      IKNTMN=1
      IF (VITRO) IKNTMN=2
      IF (DKNTMN(IKNTMN) .LE. 0.) call errmes (25, 4, chsubp)
c     -------------------------------------------------------------------------
c     Bring DKNTMN into a reasonable range, so that ALPBMX is reasonable.
c     DKNTMN cannot be greater than the limits below, because
c        NBACKG >= NBCKMN is imposed later.
c     -------------------------------------------------------------------------
      if (nbckmn .lt. 4) call errmes (24, 4, chsubp)
      dkntmn(ikntmn) = amin1((ppmst-ppmend) / float(nbckmn - 3),
     1                       dkntmn(ikntmn))
C     -------------------------------------------------------------------------
C     This rescaling is necessary, because the scaling of alphaB with
C       knot spacing is still not invariant.
c     Assume that the ALPB* are originally set for DKNTMN_STANDARD.
C     -------------------------------------------------------------------------
      scale = (DKNTMN(IKNTMN) / DKNTMN_standard)**3
      ALPBMX = ALPBMX * scale
      ALPBMn = ALPBMn * scale
      ALPBst = ALPBst * scale
      ALPBpn = ALPBpn * scale
C     -------------------------------------------------------------------------
c     Fix phases for absolute-value spectra.
C     -------------------------------------------------------------------------
      if (absval) then
         fwhmst = sqrt(3.) * fwhmst
         fwhmmx = sqrt(3.) * fwhmmx
         degzer = 0.
         degppm = 0.
         sddegz = 0.
         sddegp = 0.
      end if
      sddegp_input = sddegp
c     -------------------------------------------------------------------------
c     Red herring; nothing will be done here.
      if (lwfft .gt. 0) npar = npar + 32
c     -------------------------------------------------------------------------
C     -------------------------------------------------------------------------
c     Set DOWS according to IAVERG
C     -------------------------------------------------------------------------
      if (iaverg .eq. 1   .or.   iaverg .eq. 4) then
         dows = .true.
      else if (iaverg .eq. 2) then
         dows = .false.
      else if (iaverg .ne. 0   .and.   iaverg .ne. 3   .and.
     1         iaverg .ne. 31   .and.   iaverg .ne. 32) then
         call errmes (13, 4, chsubp)
      end if
C     -------------------------------------------------------------------------
c     Check and order dimensions ND* I* *SK
C     -------------------------------------------------------------------------
      if (min0(ndrows, ndcols, ndslic, icolst, icolen, irowst,
     1         irowen, islice) .le. 0) call errmes (9, 4, chsubp)
      i1 = irowst
      irowst = min0(irowst, irowen)
      irowen = max0(i1, irowen)
      i1 = icolst
      icolst = min0(icolst, icolen)
      icolen = max0(i1, icolen)
      if (icolen .gt. ndcols   .or.   irowen .gt.  ndrows   .or.
     1    islice .gt. ndslic) call errmes (10, 4, chsubp)
      nvoxsk = min0(nvoxsk, mvoxsk)
      do 310 j = 1, nvoxsk
         if (min0(irowsk(j), icolsk(j)) .le. 0   .or.
     1       irowsk(j) .gt. ndrows   .or.
     2       icolsk(j) .gt. ndcols) call errmes (14, 4, chsubp)
 310  continue
C     -------------------------------------------------------------------------
c     Check L* & FIL* (except FILBAS, because test conversion only will not
c       have it).
C     -------------------------------------------------------------------------
      if (lps .gt. 0   .and.   filps .eq. ' ')
     1      call errmes (3, -4, chsubp)
      if (lcoord .gt. 0   .and.   filcoo .eq. ' ')
     1      call errmes (4, 4, chsubp)
      if (lcoraw .gt. 0   .and.   filcor .eq. ' ')
     1      call errmes (11, 4, chsubp)
      if (lcsv .gt. 0   .and.   filcsv .eq. ' ')
     1      call errmes (15, 4, chsubp)
      if (ltable .gt. 0   .and.   filtab .eq. ' ')
     1      call errmes (5, 4, chsubp)
      if (lraw .le. 0   .or.   filraw .eq. ' ')
     1      call errmes (6, 4, chsubp)
      if ((dows .or. doecc   .or.   unsupr)   .and.
     1    (lh2o .le. 0   .or.   filh2o .eq. ' '))
     2      call errmes (7, 4, chsubp)
      if (imethd .eq. 2   .and.
     1    (ipowrg .ne. 1   .and.   ipowrg .ne. 2))
     2   call errmes (23, 4, chsubp)
      linerr_mycont = linerr
      return
 801  call errmes (1, -4, chsubp)
 802  call errmes (2, -4, chsubp)
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine check_zero_voxels ()
c
c Go through all voxels RAW file and set ZERO_VOXEL(IVOXEL)=T for zero voxels.
c
      INCLUDE 'lcmodel.inc'
      CHARACTER FMTDAT*(MCHFmt), ID*(MCHID)
      LOGICAL bruker, seqacq
c
      NAMELIST /NMID/ ID, FMTDAT, TRAMP, VOLUME, SEQACQ, BRUKER
      CHSUBP='ZEROVX'
c
C     -------------------------------------------------------------------------
C     Open LRAW
c     (2cygwin does not left-shift below; others do)
C     -------------------------------------------------------------------------
C     IF (FILRAW .NE. ' ') OPEN (LRAW, FILE=FILRAW, STATUS='OLD',!sun
C    1                           err=803)!sun
      IF (FILRAW .NE. ' ') OPEN (LRAW, FILE=FILRAW, STATUS='OLD',!Cyg
     1                           err=803)!Cyg
C      IF (FILRAW .NE. ' ') OPEN (LRAW, FILE=FILRAW, STATUS='OLD',!OSF
C     1                           READONLY, err=803)!OSF
C      IF (FILRAW .NE. ' ') OPEN (LRAW, FILE=FILRAW, STATUS='OLD',!IRIX
C     1                           READONLY, err=803)!IRIX
C     -------------------------------------------------------------------------
C     Read time-domain data into DATAT.
C     -------------------------------------------------------------------------
      READ (LRAW, NML=NMID, err=804, end=804)
      IF (FMTDAT .EQ. ' ') CALL ERRMES (1, 4, CHSUBP)
      ivoxel = 0
      do 110 idslic = 1, ndslic
         do 120 idrow = 1, ndrows
            do 130 idcol = 1, ndcols
               ivoxel = ivoxel + 1
               if (ivoxel .gt. mvoxel) call errmes (2, 4, chsubp)
               zero_voxel(ivoxel) = .false.
               READ (LRAW, FMTDAT, err=805, end=805)
     1                 (DATAT(J),J=1,NUNFIL)
               do 150 j = 1, nunfil
                  if (real(datat(j))**2 + aimag(datat(j))**2 .gt. 0.)
     1                  go to 130
 150           continue
               zero_voxel(ivoxel) = .true.
 130        continue
 120     continue
 110  continue
      rewind lraw
      return
 803  call errmes (3, 4, chsubp)
 804  call errmes (4, 4, chsubp)
 805  call errmes (5, 4, chsubp)
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine average ()
c
c Called when IAVERG>=1.
c
c IAVERG  normalize  weight
c    1       Y         Y  (need DOWS=T) (as in phased-array)
c    2       N         Y
c    3       N         N
c    4       Y         N  (hard to imagine a use)
c IAVERG = 2 or 3 could be used for a series of Philips or Elscint scans.
c          2 forces DOWS=F.
c          3 sums DATAT (and, if present) H2OT and then does ECC or WS with
c            final sums and not individual "voxels".
c          31 or 32 is the same as IAVERG=3, except that only the odd or even
c            voxels are summed.  These can be used with Philips data to
c            sum only the "OFF" spectra.
c
c Averages (phased-array) spectra stored in CSI format in DATAT, then
c    puts the average in DATAT and sets NDCOLS, etc to single-voxel values.
c Uses AREAWA for weighting according to signal strength.  All signals are 1st
c    normalized, so that the weights (for normalized DATAT) are 1/VAR.
C Cannot do weighted average of H2OT: phases vary greatly among channels.
c    However, normalizing DATAT means that AREA_WATER=1 (see Notes of 110311).
c
      INCLUDE 'lcmodel.inc'
      external areawa
      parameter (mchannel=2048)
      logical avgh2o
      real channel_rms(mchannel), channel_signal(mchannel)
      chsubp = 'AVERAG'
c     ------------------------------------------------------------------------
c     PPMINC, NDATA, FNDATA & RADIAN are needed by AREAWA
c     ------------------------------------------------------------------------
      PPMINC = DELTAT * float(2 * nunfil) * HZPPPM
      IF (PPMINC .LE. 0.) CALL ERRMES (1, 4, CHSUBP)
      PPMINC=1./PPMINC
      PI=3.141592654D0
      RADIAN=PI/180.
      NDATA=2*NUNFIL
      fndata = float(ndata)
      voxel1 = .true.
      lraw_at_top = .true.
      lprint_sav = lprint
      lprint = 0
      nterm = nback(1) - nback(2) + 1
      if (nterm .lt. 20   .or.   min0(nback(1), nback(2)) .lt. 0)
     1       call errmes (2, 4, chsubp)
      jvoxel = 0
      nchannel_used = 0
      avgh2o = (iaverg .eq. 3   .or.   iaverg .eq. 31   .or.
     1          iaverg .eq. 32)   .and.   (doecc .or. dows)
c     -------------------------------------------------------------------------
c     DATAT_WORK: weighted DATAT will be accumulated here.
c     H2OF_WORK: Unweighted sum of H2OTs accumulated here when
c                 IAVERG=3, 31 or 32.
c     -------------------------------------------------------------------------
      do 150 junfil = 1, nunfil
         datat_work(junfil) = (0.,0.)
         h2of_work(junfil) = (0.,0.)
 150  continue
      sumwt = 0.
      do 210 idslic = 1, ndslic
         do 220 idrow = 1, ndrows
            do 230 idcol = 1, ndcols
               jvoxel = jvoxel + 1
               call mydata ()
               avgh2o = avgh2o   .and.   havh2o
c              ----------------------------------------------------------------
c              DOWS=T is set in MYCONT if IAVERG=1 or 4.
c              ----------------------------------------------------------------
               if (.not.havh2o   .and.
     1             (iaverg .eq. 1   .or.   iaverg .eq. 4))
     2            call errmes (3, 4, chsubp)
               voxel1 = .false.
               lraw_at_top = .false.
               if (idrow .lt. irowst   .or.
     2             idrow .gt. irowen   .or.
     3             idcol .lt. icolst   .or.
     4             idcol .gt. icolen   .or.
     5             idslic .ne. islice   .or.
     6             zero_voxel(jvoxel)) go to 230
               do 235 j = 1, nvoxsk
                  if (idrow .eq. irowsk(j)   .and.
     2                idcol .eq. icolsk(j)) go to 230
 235           continue
               if (iaverg .eq. 31   .and.   mod(jvoxel, 2) .eq. 0)
     1            go to 230
               if (iaverg .eq. 32   .and.   mod(jvoxel, 2) .eq. 1)
     1            go to 230
               nchannel_used = nchannel_used + 1
               if (nchannel_used .gt. mchannel)
     1            call errmes (4, 4, chsubp)
               if (iaverg .eq. 1   .or.   iaverg .eq. 4) then
                  channel_signal(nchannel_used) = areawa(1)
               else
                  channel_signal(nchannel_used) = 1.
               end if
               if (channel_signal(nchannel_used) .le. 0.)
     1            call errmes (5, 4, chsubp)
c              ---------------------------------------------------------------
c              Normalize DATAT
c              ---------------------------------------------------------------
               do 240 j = 1, nunfil
                  datat(j) = datat(j) / channel_signal(nchannel_used)
 240           continue
               if (iaverg .gt. 2) then
                  channel_rms(nchannel_used) = 1.
               else
                  term = getvar()
                  if (term .le. 0.) call errmes (6, 4, chsubp)
                  channel_rms(nchannel_used) = sqrt(term /
     1                                             float(2 * nterm - 4))
                  if (channel_rms(nchannel_used) .le. 0.)
     1                 call errmes (7, 4, chsubp)
               end if
               term = channel_rms(nchannel_used)**2
               if (term .le. 0.) call errmes (9, 4, chsubp)
               channel_wt = 1. / term
c              --------------------------------------------------------------
c              The following is the popular weighting proportional to signal
c                 discussed by Wright & Wald
c                 channel_wt = channel_signal(nchannel_used)**2
c              --------------------------------------------------------------
               sumwt = sumwt + channel_wt
               do 250 junfil = 1, nunfil
                  datat_work(junfil) = datat_work(junfil) +
     1                                 channel_wt * datat(junfil)
                  if (avgh2o) h2of_work(junfil) = h2of_work(junfil) +
     1                                            h2ot(junfil)
 250           continue
 230        continue
 220     continue
 210  continue
      if (iaverg .eq. 31   .or.   iaverg .eq. 32)   then
         if (nchannel_used .lt. ndrows * ndcols / 2)
     1      call errmes (8, 2, chsubp)
      else if (nchannel_used .lt. ndrows * ndcols) then
         call errmes (8, 2, chsubp)
      end if
      if (ldump(1)) then
         write (6, 5310) (j, channel_signal(j), channel_rms(j),
     1                     1. / channel_rms(j),  j = 1, nchannel_used)
 5310    format('ch#', 8x, 'areawa', 11x, 'rms', 4x, 'sgnl/noise' /
     2          (i3, 1p3e14.4))
      end if
c     ------------------------------------------------------------------------
c     SUMWT = 0 would most likely be caused by excluding all (non-zero)
c               voxels.
c     ------------------------------------------------------------------------
      if (sumwt .le. 0) call errmes (10, 4, chsubp)
      do 360 j = 1, nunfil
         datat(j) = datat_work(j) / sumwt
         if (avgh2o) h2ot(j) = h2of_work(j) / sumwt
 360  continue
      if (doecc   .and.   avgh2o) call ecc_truncate ()
      lprint = lprint_sav
      ndslic = 1
      ndrows = 1
      ndcols = 1
      islice = 1
      irowst = 1
      irowen = 1
      icolst = 1
      icolen = 1
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function getvar ()
c
c Computes variance of noise near end of DATAT, determined by NBACK(2), using
c    SSQ of linear regression (Draper & Smith, p 16).
c Must have already verified that NBACK(1)-NBACK(2) > 20
c Called (twice) from AVERAGE.
c
      INCLUDE 'lcmodel.inc'
      double precision dnterm, drn, sx, sxx, sy(2), sxy(2), syy(2),
     1                 yv(2)
c
      do 110 k = 1, 2
         sy(k) = 0.
         sxy(k) = 0.
         syy(k) = 0.
 110  continue
      sx = 0.
      sxx = 0.
      dnterm = 0.d0
      do 210 junfil = nunfil - nback(1), nunfil - nback(2)
         dnterm = dnterm + 1.d0
         sx = sx + dnterm
         sxx = sxx + dnterm**2
         yv(1) = dble(real(datat(junfil)))
         yv(2) = dble(aimag(datat(junfil)))
         do 220 k = 1, 2
            sy(k) = sy(k) + yv(k)
            sxy(k) = sxy(k) + dnterm * yv(k)
            syy(k) = syy(k) + yv(k)**2
 220     continue
 210  continue
      getvar = 0.
      drn = 1.d0 / dnterm
      do 310 k = 1, 2
c        ---------------------------------------------------------------------
c        From Draper & Smith, p 16.
c        ---------------------------------------------------------------------
         getvar = getvar + syy(k) - drn * sy(k)**2 -
     1                     (sxy(k) - drn * sx * sy(k))**2 /
     2                     (sxx - drn * sx**2)
 310  continue
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine restore_settings ()
c
c When VOXEL1 = T, saves settings that will be changed.
c When VOXEL1 = F, restores settings that have been changed.
c
      INCLUDE 'lcmodel.inc'
      save chrato_sav, degppm_sav, degzer_sav,
     1     fcalib_sav, fwhmst_sav,
     1     hzref_sav, lcoord_sav, lcoraw_sav,
     1     lprint_sav,
     1     ltable_sav, miter1_sav, ndegz_sav, ndgppm_sav,
     1     nlin_sav, nnolsh_sav, nnot1_sav, npar_sav,
     1     nrefpk_sav, ppmref_sav, sdmshf_sav,
     1     title_sav
      character chrato_sav(mmetab)*(264), title_sav*244
      integer ndegz_sav(2), ndgppm_sav(3), nrefpk_sav(2)
      real hzref_sav(mrefpk,2), ppmref_sav(mrefpk,2)
c
      if (voxel1) then
         do 102 j = 1, nratio
            chrato_sav(j) = chrato(j)
 102     continue
         degppm_sav = degppm
         degzer_sav = degzer
c        ----------------------------------------------------------------------
C        DGPPM*_ORIG will be used globally and not reset (or restored)
c        ----------------------------------------------------------------------
         dgppmn_orig = dgppmn
         dgppmx_orig = dgppmx
         fcalib_sav = fcalib
         fwhmst_sav = fwhmst
         lcoord_sav = lcoord
         lcoraw_sav = lcoraw
c         ldwfft_sav = ldwfft
         lprint_sav = lprint
         ltable_sav = ltable
c         lwfft_sav = lwfft
         miter1_sav = miter(1)
         nlin_sav = nlin
         nnolsh_sav = nnolsh
         nnot1_sav = nnot1
         npar_sav = npar
         sdmshf_sav = sdmshf
         do 110 jset = 1, 2
            ndegz_sav(jset) = ndegz(jset)
            ndgppm_sav(jset) = ndgppm(jset)
            nrefpk_sav(jset) = nrefpk(jset)
            do 120 jpeak = 1, min0(mrefpk, nrefpk(jset))
               hzref_sav(jpeak, jset) = hzref(jpeak, jset)
               ppmref_sav(jpeak, jset) = ppmref(jpeak, jset)
 120        continue
c           -------------------------------------------------------------------
C           SHIFM*_ORIG will be used globally and not reset (or restored)
c           -------------------------------------------------------------------
            shifmn_orig(jset) = shifmn(jset)
            shifmx_orig(jset) = shifmx(jset)
 110     continue
         title_sav = title
      else
         do 202 j = 1, nratio
            chrato(j) = chrato_sav(j)
 202     continue
c        ---------------------------------------------------------------------
c             DEGPPM & DEGZER = 0 at the end of an analysis, because of calls
c        to REPHAS.  So, if they are input, they must be restrored from 0,
c        unless they are non-zero here, in which case they have been set in
c        UPDATE_PRIORS for this analysis and should be left unchanged.
c        ---------------------------------------------------------------------
         if (amax1(abs(degppm), abs(degzer)) .lt. 1.e-5) then
            degppm = degppm_sav
            degzer = degzer_sav
         end if
         fcalib = fcalib_sav
         fwhmst = fwhmst_sav
         lcoord = lcoord_sav
         lcoraw = lcoraw_sav
c         ldwfft = ldwfft_sav
c         linerr = 0
         linerr = linerr_mycont
         linetc = 0
         lintbl = 0
         lprint = lprint_sav
         ltable = ltable_sav
c         lwfft = lwfft_sav
         miter(1) = miter1_sav
         nermes = linerr_mycont
         nlin = nlin_sav
         nnolsh = nnolsh_sav
         nnot1 = nnot1_sav
         npar = npar_sav
         sdmshf = sdmshf_sav
         do 210 jset = 1, 2
            ndegz(jset) = ndegz_sav(jset)
            ndgppm(jset) = ndgppm_sav(jset)
            nrefpk(jset) = nrefpk_sav(jset)
            do 220 jpeak = 1, min0(mrefpk, nrefpk(jset))
               hzref(jpeak, jset) = hzref_sav(jpeak, jset)
               ppmref(jpeak, jset) = ppmref_sav(jpeak, jset)
 220        continue
 210     continue
 
         do 230 j = linerr_mycont + 1, mmerm9
            nerror(j) = 0
 230     continue
         title = title_sav
      end if
      area_met_norm = 0.
      istago = 0
      nratio_used = 0
      wsdone = .false.
      do 320 j = 1, 2
         phitot(j) = 0.
 320  continue
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine update_priors ()
c
c      Use sample variances to update priors for DEGZER, DEGPPM & SHIFM*.
c      Only have to accumulate sums and sums of squares for SDDEGP & SHIFM*.
c      For SDDEGZ, must accumulate sums of sin & cos, then compute mean angle,
c and then compute variance from mean, since angle must be chosen to be
c closest to mean.
c SUM(J) & SUM2(J): J = 1 for shift
c                       2     DEGPPM
c                       3     cos(DEGZER) (no SUM2)
c                       4     sin(DEGZER)
c
      INCLUDE 'lcmodel.inc'
      parameter (msamples=4096)
      save degzer_sample, initial, nsamples, rppminc_shift,
     1     sddegp_min, sddegz_min, sum, sum2
      logical initial
      real degzer_sample(msamples), sd(2), smean(2),
     1     sum(4), sum2(3)
      data nsamples/0/, sum/4*0./, sum2/3*0./
      data initial/.true./
c     -------------------------------------------------------------------------
c     The following are somewhat arbitrary and may have to be changed.
c     -------------------------------------------------------------------------
      data rppminc_shift/5./, sddegp_min/1./, sddegz_min/3./
      chsubp = 'UPDPRI'
 
      if (nsamples .ge. msamples) return
 
      if (lcsi_sav_2 .eq. 13   .and.   initial) then
         initial = .false.
 5150    format (i5)
         read (13, 5150, err=150, end=160) nsamples
         read (13, 5160, err=150, end=150) sum, sum2,
     1      (degzer_sample(j), j = 1, nsamples)
 5160    format (1p5e16.6)
         go to 160
 150     call errmes (1, 4, chsubp)
      end if
 160  nsamples = nsamples + 1
      term = ppminc * float(ishifd)
      sum(1) = sum(1) + term
      sum2(1) = sum2(1) + term**2
 
      sum(2) = sum(2) + phitot(2)
      sum2(2) = sum2(2) + phitot(2)**2
 
      degzer_sample(nsamples) = phitot(1)
      term = cos(radian * phitot(1))
      sum(3) = sum(3) + term
      term = sin(radian * phitot(1))
      sum(4) = sum(4) + term
 
      if (lcsi_sav_2 .eq. 13) then
         rewind 13
         write (13, 5150) nsamples
         write (13, 5160) sum, sum2, (degzer_sample(j), j = 1, nsamples)
      end if
 
      if (nsamples .lt. max0(2, mnsamp)) return
 
c     ------------------------------------------------------------------------
c     j = 1 for shift; j = 2 for DEGPPM
c     ------------------------------------------------------------------------
      do 210 j = 1, 2
         smean(j) = sum(j) / float(nsamples)
         sd(j) = rsdsam(j) *
     1           sqrt(amax1(0.,
     1                      (sum2(j) - smean(j)**2 * float(nsamples)) /
     2                      float(nsamples - 1)))
 210  continue
      sd(1) = amax1(rppminc_shift * ppminc, sd(1))
c     ------------------------------------------------------------------------
c     Keep SHIFM* within bounds of SHIFM*, to avoid excessive shift ranges.
c     ------------------------------------------------------------------------
      do 220 j = 1, 2
         shifmn(j) = amin1(shifmx_orig(j),
     1                     amax1(smean(1) - sd(1), shifmn_orig(j)))
         shifmx(j) = amax1(shifmn_orig(j),
     1                     amin1(smean(1) + sd(1), shifmx_orig(j)))
 220  continue
      degppm = smean(2)
      dgppmx = amax1(degppm, dgppmx_orig)
      dgppmn = amin1(degppm, dgppmn_orig)
c     ------------------------------------------------------------------------
c     Do not allow SDDEGP to exceed SDDEGP_INPUT (presumably <=20) to avoid
c     destabilizing the analysis (Notes of 090201)
c     ------------------------------------------------------------------------
      sddegp = amin1(sddegp_input, amax1(sd(2), sddegp_min))
c     ------------------------------------------------------------------------
c     DEGZER
c     ------------------------------------------------------------------------
      degzer = atan2(sum(4), sum(3)) / radian
c     ------------------------------------------------------------------------
c     Make DEGZER between 0 and 360 (as DEGZER_SAMPLE already is).
c     ------------------------------------------------------------------------
      if (degzer .lt. 0.) degzer = degzer + 360.
      sum2(3) = 0.
      do 250 j = 1, nsamples
         term = degzer_sample(j) - degzer
         sum2(3) = sum2(3) + amin1(abs(term), abs(term + 360.),
     1                                        abs(term - 360.))**2
 250  continue
      sd(1) = rsdsam(3) * sqrt(sum2(3) / float(nsamples - 1))
      sddegz = amax1(sd(1), sddegz_min)
 
      if (lcsi_sav_1 .eq. 12) then
         rewind 12
         read (12, 5150) j
         write (12, 5160) degppm, degzer, dgppmn, dgppmx,
     1                    sddegp, sddegz, shifmn, shifmx
      end if
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine open_output ()
      INCLUDE 'lcmodel.inc'
      EXTERNAL ILEN
      save lchcol, lchrow, lchslic,
     1     lchcol_st, lchrow_st, lchslic_st,
     2     splcoo, splcor, splpri, splps, spltab
      CHARACTER CHSTR*40, chidcol*6, chidrow*6, chidslic*6,
     1          chinsert*36, ch256*256,
     2          splcoo(2)*(mchfil+1), splcor(2)*(mchfil+1),
     3          splpri(2)*(mchfil+1),
     4          splps(2)*(mchfil+1), spltab(2)*(mchfil+1)
      integer jtitle(2)
C
      include 'nml_lcmodl.inc'
      include 'nml_lcmodel.inc'
      CHSTR = ' '
      chsubp = 'OPENOU'
C     -------------------------------------------------------------------------
c     Get date & time info.
C     -------------------------------------------------------------------------
      CHDATE=' '
C      call idate (knum) !sun
C      call getdate(chdate) !sun
       call fdate(chdate)!Cyg
C      call y2k !OSF
C      call y2k !IRIX
C     -------------------------------------------------------------------------
c     Modify TITLE and filenames to include information on voxel coordinates.
c     This will be done unless the data set is a single-voxel one, i.e.,
c       unless SINGLE_VOXEL = T.
C     -------------------------------------------------------------------------
      if (.not.single_voxel) then
         if (voxel1) then
            if (lcoord .gt. 0) call split_filename(filcoo, 'coord',
     1                           'COORD', 'Coord', 5, splcoo)
            if (lcoraw .gt. 0) call split_filename(filcor, 'coraw',
     1                           'CORAW', 'Coraw', 5, splcor)
            if (lprint .gt. 0) call split_filename(filpri, 'print',
     1                           'PRINT', 'Print', 5, splpri)
            if (lps .gt. 0) call split_filename(filps, 'ps',
     1                                          'PS', 'Ps', 2, splps)
            if (ltable .gt. 0) call split_filename(filtab, 'table',
     1                           'TABLE', 'Table', 5, spltab)
            lchslic_st = icharst(chslic, len(chslic))
            lchrow_st = icharst(chrow, len(chrow))
            lchcol_st = icharst(chcol, len(chcol))
            lchcol = ilen(chcol)
            lchrow = ilen(chrow)
            lchslic = ilen(chslic)
         end if
         call chstrip_int6(idcol, chidcol, lidcol)
         call chstrip_int6(idrow, chidrow, lidrow)
         call chstrip_int6(idslic, chidslic, lidslic)
         if (ndslic .gt. 1) then
            ch256 = 'Slice#' // chidslic(:lidslic) //
     1              ' Row#' // chidrow(:lidrow) //
     2              ' Col#' // chidcol(:lidcol) // '  ' // title
         else
            ch256 = 'Row#' // chidrow(:lidrow) //
     1              ' Col#' // chidcol(:lidcol) // '  ' // title
         end if
         title = ch256
 
         if (lchslic_st .gt. 0) then
            chinsert = chslic(lchslic_st:lchslic) // chidslic(:lidslic)
         else
            chinsert = chidslic(:lidslic)
         end if
         if (lchrow_st .gt. 0) then
            lstr = ilen(chinsert)
            chinsert = chinsert(:lstr) // chrow(lchrow_st:lchrow)
         end if
         lstr = ilen(chinsert)
         chinsert = chinsert(:lstr) // chidrow(:lidrow)
         if (lchcol_st .gt. 0) then
            lstr = ilen(chinsert)
            chinsert = chinsert(:lstr) // chcol(lchcol_st:lchcol)
         end if
         lstr = ilen(chinsert)
         chinsert = chinsert(:lstr) // chidcol(:lidcol)
         linsert = ilen(chinsert)
 
         if (lcoord .gt. 0) then
            lsplit1 = ilen(splcoo(1))
            lsplit2 = ilen(splcoo(2))
            if (lsplit1 + linsert + lsplit2 .gt. mchfil)
     1         call errmes (1, 4, chsubp)
            filcoo = splcoo(1)(:lsplit1) // chinsert(:linsert) //
     1               splcoo(2)(:lsplit2)
         end if
         if (lcoraw .gt. 0) then
            lsplit1 = ilen(splcor(1))
            lsplit2 = ilen(splcor(2))
            if (lsplit1 + linsert + lsplit2 .gt. mchfil)
     1         call errmes (1, 4, chsubp)
            filcor = splcor(1)(:lsplit1) // chinsert(:linsert) //
     1               splcor(2)(:lsplit2)
         end if
         if (lprint .gt. 0) then
            lsplit1 = ilen(splpri(1))
            lsplit2 = ilen(splpri(2))
            if (lsplit1 + linsert + lsplit2 .gt. mchfil)
     1         call errmes (1, 4, chsubp)
            filpri = splpri(1)(:lsplit1) // chinsert(:linsert) //
     1               splpri(2)(:lsplit2)
         end if
         if (lps .gt. 0) then
            lsplit1 = ilen(splps(1))
            lsplit2 = ilen(splps(2))
            if (lsplit1 + linsert + lsplit2 .gt. mchfil)
     1         call errmes (1, 4, chsubp)
            filps = splps(1)(:lsplit1) // chinsert(:linsert) //
     1               splps(2)(:lsplit2)
         end if
         if (ltable .gt. 0) then
            lsplit1 = ilen(spltab(1))
            lsplit2 = ilen(spltab(2))
            if (lsplit1 + linsert + lsplit2 .gt. mchfil)
     1         call errmes (1, 4, chsubp)
            filtab = spltab(1)(:lsplit1) // chinsert(:linsert) //
     1               spltab(2)(:lsplit2)
         end if
      end if
      if (skip_voxel) then
c        ---------------------------------------------------------------------
c        Temporarily prevent opening output files and return.
c        LPS is left unchanged, because it would only be produced on
c          an error exit.
c        ---------------------------------------------------------------------
         lcoord = 0
         lcoraw = 0
         lprint = 0
         ltable = 0
         return
      end if
c     -------------------------------------------------------------------------
c     NLIN = 16384
c     -------------------------------------------------------------------------
      if (nratio_used .gt. 0) then
         nlin = 4 * nlin + 8192
      else
         nlin = 4 * nlin + 16384
      end if
c     -------------------------------------------------------------------------
c     Split TITLE into TITLE_LINE(1) & TITLE_LINE(2).
c     -------------------------------------------------------------------------
      call split_title ()
      JOWNER=ILEN(OWNOUT)
      JDATE=ILEN(CHDATE)
      JVERSI=INDEX(VERSIO,'Copyright')-2
      JTITLE(1)=ILEN(TITLE_line(1))
      JTITLE(2)=ILEN(TITLE_line(2))
      IF (LPRINT .GT. 0) THEN
         IF (FILPRI .NE. ' ') OPEN (LPRINT, FILE=FILPRI,
     1                              STATUS='UNKNOWN', err=802)
         if (nlines_title .eq. 1) then
            WRITE (LPRINT,5000) VERSIO,
     1                          TITLE_line(1)(1:JTITLE(1)),
     2                          OWNOUT(1:JOWNER), CHDATE(1:JDATE)
 5000       FORMAT (///10X,A//1X,A//1X,A//1X,A////)
         else
            WRITE (LPRINT,5002) VERSIO,
     1            (TITLE_line(j)(1:JTITLE(j)), j=1,2),
     2            OWNOUT(1:JOWNER), CHDATE(1:JDATE)
 5002       FORMAT (///10X,A//1X,A/1x,a//1X,A//1X,A////)
         end if
         if (mermes .eq. 8000) then
            mermes = mmerms
            WRITE (LPRINT, NML=LCMODL)
         else
            WRITE (LPRINT, NML=LCMODeL)
         END IF
      END IF
      IF (LCOORD .GT. 0) THEN
         IF (FILCOO .NE. ' ') OPEN (LCOORD, FILE=FILCOO,
     1                              STATUS='UNKNOWN', err=803)
         WRITE (LCOORD,5580) VERSIO,
     1            (TITLE_line(j)(1:JTITLE(j)), j=1,nlines_title)
C         write (lcoord, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!sun
C         write (lcoord, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!IRIX
C         write (lcoord, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!OSF
 5580    FORMAT (1X, A)
      END IF
      IF (LCORAW .GT. 0) THEN
         IF (FILCOR .NE. ' ') OPEN (LCORAW, FILE=FILCOR,
     1                              STATUS='UNKNOWN', err=805)
      END IF
      IF (LTABLE .GT. 0) THEN
         IF (FILTAB .NE. ' ') OPEN (LTABLE, FILE=FILTAB,
     1                              STATUS='UNKNOWN', err=804)
         WRITE (LTABLE,5580) VERSIO(1:JVERSI),
     1            (TITLE_line(j)(1:JTITLE(j)), j=1,nlines_title)
C         write (ltable, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!sun
C         write (ltable, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!IRIX
C         write (ltable, 5580) OWNOUT(1:JOWNER), CHDATE(1:JDATE)!OSF
      END IF
c     -------------------------------------------------------------------------
c     NPAR = 20
c     -------------------------------------------------------------------------
      if (ldwfft .gt. 0) then
         npar = 2 * npar + 40
      else
         npar = npar + 20
      end if
      return
 802  call errmes (2, 4, chsubp)
 803  call errmes (3, 4, chsubp)
 804  call errmes (4, 4, chsubp)
 805  call errmes (5, 4, chsubp)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE split_filename (filename, chtype1, chtype2, chtype3,
     1                           lchtype, split)
C
c Splits FILENAME into 2 parts into SPLIT for insertion of voxel
c   identifier between the 2 parts.
c CHTYPE1-3 are the 3 typical variations in the file type, e.g., ps PS Ps.
c LCHTYPE = length of strings chtype*; they must all have the same length.
c There are 3 cases for the multi-voxel identifier (CHINSERT) to be inserted.
c
      EXTERNAL ILEN
      character char*1, chtype1*5, chtype2*5, chtype3*5, filename*(*),
     1          split(2)*(*)
c
      lfile = ilen(filename)
      ichtype = lfile - lchtype + 1
      if (max0(index(filename(ichtype:), chtype1(:lchtype)),
     1         index(filename(ichtype:), chtype2(:lchtype)),
     2         index(filename(ichtype:), chtype3(:lchtype))) .eq. 1)
     3      then
c        ---------------------------------------------------------------------
c        CHTYPE is at the end of filename, as it normally should be.  Check
c          for two cases:
c        ---------------------------------------------------------------------
         islash_dot = ichtype - 1
         char = filename(islash_dot: islash_dot)
         if (char .eq. '/') then
c           ------------------------------------------------------------------
c           (1) Default LCMgui FILENAME, e.g., (...)/ps.  Split as
c               (...)/ | ps and insert (chinsert). to get
c               (...)/(chinsert).ps
c           ------------------------------------------------------------------
            split(1) = filename(:islash_dot)
            split(2) = '.' // filename(ichtype:)
            return
         end if
         if (char .eq. '.') then
c           ------------------------------------------------------------------
c           (2) Alternate LCMgui FILENAME, e.g., (...).ps.  Split as
c               (...) | .ps and insert _(chinsert) to get
c               (...)_(chinsert).ps
c           ------------------------------------------------------------------
            split(1) = filename(:islash_dot - 1) // '_'
            split(2) = filename(islash_dot:)
            return
         end if
      end if
c     ------------------------------------------------------------------------
c     (3) Exceptional FILENAME, e.g., (...).  Split as
c         (...) | ' ' and insert _(chinsert) to get
c         (...)_(chinsert)
c     ------------------------------------------------------------------------
      split(1) = filename(:lfile) // '_'
      split(2) = ' '
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function icharst(ch, lch)
c
c Returns index of first non-space charcter in CH
c         -1 if ch(1:lch) are all spaces
c Normally this should be called by ICHARST(CH, LEN(CH)), but the 2nd argument
c   could also be < LEN(CH).
c
      character ch(*)
      do 110 j = 1, lch
         if (ch(j) .ne. ' ') then
            icharst = j
            return
         end if
 110  continue
      icharst = -1
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine chstrip_int6(iarg, chi, leni)
c
c CHI = string with IARG value and no spaces of length LENI.
c     will be restricted to 6 characters, i.e., -99999 <= CHI value <= 999999
c     must be at least 6 characters long in calling program
c
      character chi*(*)
      i = max0(-99999, min0(999999, iarg))
      if (i .ge. 0   .and.   i .le. 9) then
         write(chi, 5001) i
         leni = 1
 5001    format(i1)
      else if (i .ge. -9   .and.   i .le. 99) then
         write(chi, 5002) i
         leni = 2
 5002    format(i2)
      else if (i .ge. -99   .and.   i .le. 999) then
         write(chi, 5003) i
         leni = 3
 5003    format(i3)
      else if (i .ge. -999   .and.   i .le. 9999) then
         write(chi, 5004) i
         leni = 4
 5004    format(i4)
      else if (i .ge. -9999   .and.   i .le. 99999) then
         write(chi, 5005) i
         leni = 5
 5005    format(i5)
      else
         write(chi, 5006) i
         leni = 6
 5006    format(i6)
      end if
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE split_title ()
C
C  Split TITLE into TITLE_LINE(1) & TITLE_LINE(2).
c  TITLE_LINE(1), after call to STRCHK, will not have more than 122 characters.
C
      INCLUDE 'lcmodel.inc'
      EXTERNAL ILEN
      ltitle=ilen(title)
c     -------------------------------------------------------------------------
c     STRCHK can at most double the number of characters (if all have to be
c       escaped).
c     -------------------------------------------------------------------------
      if (ntitle .eq. 1   .or.   ltitle .le. 61) then
         nlines_title = 1
         title_line(1) = title
         return
      end if
      iend_line1 = min0(ltitle,122)
c     -------------------------------------------------------------------------
c     NESCAPE1 = extra escape characters that will be added to BUFOUT [and
c                therefore will reduce the allowed length of TITLE_LINE(1)].
c     -------------------------------------------------------------------------
      nescape = 0
      nescape1 = 0
      do 110 i = 1, ltitle
         IF (title(I:I).EQ.'('
     1     .OR.title(I:I).EQ.')'
     2     .OR.title(I:I).EQ.'%'
     3     .OR.title(I:I).EQ.'\') then
c
c          A ' to correct error in emacs display
c
            nescape = nescape + 1
            if (i .le. iend_line1) nescape1 = nescape1 + 1
         end if
 110  continue
      if (ltitle + nescape .le. 122) then
         nlines_title = 1
         title_line(1) = title
         return
      end if
c     -------------------------------------------------------------------------
c     Go backwards on TITLE_LINE(1) to find a space for a break.
c     MAX_BACK < 122, because above test ensures (LTITLE + NESCAPE) > 122
c     -------------------------------------------------------------------------
      max_back = 244 - ltitle - nescape
      iend_line1 = iend_line1 - nescape1
      ibreak = iend_line1
      istart2 = iend_line1 + 1
      do 120 i = iend_line1, max0(4, iend_line1 - max_back), -1
         if (title(i:i) .eq. ' ') then
            ibreak = i - 1
            istart2 = i + 1
            go to 150
         end if
 120  continue
 150  nlines_title = 2
      title_line(1) = title(1:ibreak)
      title_line(2) = title(istart2:ltitle)
      return
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE LOADCH ()
C
C  Loads changes input in standard input into CHANGE for later output.
C
      INCLUDE 'lcmodel.inc'
      parameter (mch_line_long=mchsimul+12)
      EXTERNAL ILEN
      CHARACTER CH*(mCHLIN), CHLONG*(mch_line_long), DOLLAR*1,
     1          LINE*(mch_line_long), line_compact*(mch_line_long)
      LOGICAL ATEND
      CHSUBP='LOADCH'
      DOLLAR='$'
      IF (MERMES .GT. MMERMS) THEN
         CALL ERRMES (1, 3, CHSUBP)
         MERMES=MMERMS
      END IF
c     -------------------------------------------------------------------------
c     NLIN = 4 * 16384 - 1 = 65535
c     -------------------------------------------------------------------------
      if (lwfft .gt. 0) then
         nlin = 2 * nlin - 3
      else
         nlin = 4 * nlin - 1
      end if
      do 100 jpage = 1, 2
c        ----------------------------------------------------------------------
c        Put NCHLIN within the limits fixed below.
c        ----------------------------------------------------------------------
         nchlin(jpage) = min0(mchlin, max0(30, nchlin(jpage)))
         REWIND LCONTR_SCRATCH
         ATEND=.FALSE.
         LINCHG(JPAGE)=0
         DO 110 JLINE=1,9999
            READ (LCONTR_SCRATCH,5110, err=810, end=810) LINE
 5110       FORMAT (A)
            L=INDEX(LINE,DOLLAR)
            IF (L.GE.1 .AND. L.LE.mch_line_long-6) THEN
               IF((LINE(L+1:L+1).EQ.'L' .OR. LINE(L+1:L+1).EQ.'l') .AND.
     1            (LINE(L+2:L+2).EQ.'C' .OR. LINE(L+2:L+2).EQ.'c') .AND.
     1            (LINE(L+3:L+3).EQ.'M' .OR. LINE(L+3:L+3).EQ.'m') .AND.
     1            (LINE(L+4:L+4).EQ.'O' .OR. LINE(L+4:L+4).EQ.'o') .AND.
     1            (LINE(L+5:L+5).EQ.'D' .OR. LINE(L+5:L+5).EQ.'d') .AND.
     1            (LINE(L+6:L+6).EQ.'L' .OR. LINE(L+6:L+6).EQ.'l'))
     2                                                       GO TO 120
            END IF
 110     CONTINUE
         RETURN
C        ---------------------------------------------------------------------
C        $LCMODL has been found; $ is at position L.
C        ---------------------------------------------------------------------
 120     LEND=ILEN(LINE)
         DO 130 JLINE=1,MLINES
            IF (JLINE.GT.1 .OR. L+6.EQ.LEND) READ (LCONTR_SCRATCH, 5110,
     1                                             END=200) LINE
            IF (ilen(line).le.1 .AND. LINE(1:1).EQ.' ') GO TO 130
            call compact_string(line, line_compact, len_compact)
            if (ldump(3)) then
               j = ilen(line)
               write (6, 5130) line(1:j), len_compact,
     1                         line_compact(1:len_compact)
 5130          format (5x, a/ i3, 2x, a/)
            end if
            IF (MAX0(INDEX(line_compact,'TITLE='),
     1               INDEX(line_compact,'title='),
     2               INDEX(line_compact,'Title='),
     3               INDEX(line_compact,'OWNER='),
     4               INDEX(line_compact,'owner='),
     5               INDEX(line_compact,'Owner='),
     6               INDEX(line_compact,'KEY='),
     7               INDEX(line_compact,'key='),
     8               INDEX(line_compact,'Key=')) .GT. 0) GO TO 130
            if (index(line, '/.lcmodel/temp/') .gt. 0) then
               if (max0( index(line, 'filcoo='), index(line, 'filh2o='),
     1                   index(line, 'filpri='), index(line, 'filps='),
     2                   index(line, 'filraw='), index(line, 'filtab='),
     3                   index(line, 'filcor='), index(line, 'filcsv='),
     4                   index(line, 'lcsi_sav'),
     5                   index(line, 'filcsi_sav'))
     6             .gt. 0 )  go to 130
            end if
            L=INDEX(LINE,DOLLAR)
            IF (L.GE.1 .AND. L.LE.mch_line_long-3) THEN
               IF((LINE(L+1:L+1).EQ.'E' .OR. LINE(L+1:L+1).EQ.'e') .AND.
     1            (LINE(L+2:L+2).EQ.'N' .OR. LINE(L+2:L+2).EQ.'n') .AND.
     1            (LINE(L+3:L+3).EQ.'D' .OR. LINE(L+3:L+3).EQ.'d')) THEN
                  IF (L .le. 2) GO TO 200
                  CHLONG=LINE(1:L-1)
                  IF (ILEN(CHLONG) .LE. 1) GO TO 200
                  LINE=CHLONG
                  ATEND=.TRUE.
               END IF
            END IF
            LEND=ILEN(LINE)
            IF (LEND.le.1 .AND. LINE(1:1).EQ.' ') GO TO 130
            IF (LEND .LE. NCHLIN(JPAGE)) THEN
               NPARTS=1
            ELSE
               NPARTS=(LEND-NCHLIN(JPAGE)-1)/(NCHLIN(JPAGE)-3)+2
            END IF
            LINCHG(JPAGE)=LINCHG(JPAGE)+NPARTS
            IF (LINCHG(JPAGE) .GT. MLINES) GO TO 200
            LSTOP=MIN0(NCHLIN(JPAGE),LEND)
            CHANGE(LINCHG(JPAGE), jpage)=LINE(1:LSTOP)
            DO 150 JCHG=LINCHG(JPAGE)-1,LINCHG(JPAGE)-NPARTS+1,-1
               LSTART=LSTOP+1
               LSTOP=min0(lend, LSTART+NCHLIN(JPAGE)-4)
               CHANGE(JCHG, jpage)='   ' // LINE(LSTART:LSTOP)
 150        CONTINUE
            IF (ATEND) GO TO 200
 130     CONTINUE
C        ----------------------------------------------------------------------
C        Reverse order in CHANGE.
C        ----------------------------------------------------------------------
 200     K=LINCHG(JPAGE)
         DO 210 J=1,LINCHG(JPAGE)/2
            CH=CHANGE(K, jpage)
            CHANGE(K, jpage)=CHANGE(J, jpage)
            CHANGE(J, jpage)=CH
            K=K-1
 210     CONTINUE
 100  continue
c     -------------------------------------------------------------------------
c     NPAR = 2 * 20 + 10 = 50
c     -------------------------------------------------------------------------
      if (lwfft .gt. 0) then
         npar = 3 * npar + 40
      else
         npar = 2 * npar + 10
      end if
      IF (LPRINT .GT. 0) WRITE (LPRINT,5210)
     1                         (CHANGE(J,2),J=1,LINCHG(2))
 5210 FORMAT (/// 10X, 'Input changes to Control Variables'/ (1X,A))
      return
 810  call errmes (2, 4, chsubp)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine compact_string (str_in, str_out, len_out)
c
c
c Remove all white space before, after and within a string.
c
      EXTERNAL ILEN
      character str_in*(*), str_out*(*)
      len_in = ilen(str_in)
      len_out = len_in
      str_out(1:len_in) = str_in(1:len_in)
      do 110 jtest = len_in, 1, -1
         if (str_in(jtest:jtest) .eq. ' ') then
            do 120 k = jtest, len_out - 1
               str_out(k:k) = str_out(k+1:k+1)
 120        continue
            len_out = len_out - 1
         end if
 110  continue
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DCFFT_R (DATAT, FT, N, LDWFFT,
     1                  DWFFTC)
C
C  Input:
C    DATAT(J), J=1,N
C    LDWFFT = argument in last call to CALL DFFTCI (LDWFFT, DWFFTC)
C  Output:
C    FT = FFT of DATAT WITH rearrangement.  The FFT is normalized; i.e.,
C             the sum is divided by sqrt(N).  This is compatible with CFFT
C             and CFFTIN and Siemens FFTs.
C
      COMPLEX*16 cterm, DATAT(N), FT(N)
      DOUBLE PRECISION DWFFTC(4*N+15)
      IF (N .NE. LDWFFT) THEN
         CALL DFFTCI (N, DWFFTC)
         LDWFFT=N
      END IF
      CALL DF2TCF (N, DATAT, FT, DWFFTC)
C     ------------------------------------------------------------------------
C     Scale & rearrange
C     ------------------------------------------------------------------------
      FACT=1./SQRT(FLOAT(N))
      nunfil = n / 2
      DO 210 J=1,nunfil
         CTERM=FT(NUNFIL+J)
         FT(NUNFIL+J)=FT(J) * fact
         FT(J)=CTERM * fact
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ERRMES (NUMBER, ILEVEL, CHSUBP)
C  Prints error message and, if IABS(ILEVEL)=4 or 5, aborts (first calling
C    EXITPS if ILEVEL>0).
C  |ILEVEL| = 1 for information.
C           = 2 for a warning.
C           = 3 for a non-fatal error.
C           = 4 for a fatal error.
C           = 5 for an illogical stop.
C           < 0 means that the error message was called while EXITPS was being
C               executed.  Thus, it cannot be called again (to avoid an
C               infinite recursive loop) and the run must be stopped here
C               without producing a PS plot file.  The error messages can only
C               be obtained on unit LPRINT.
      INCLUDE 'lcmodel.inc'
      save ertype
      CHARACTER ERTYPE(5)*15
      DATA ERTYPE/'    INFORMATION', '        WARNING',
     1            'NON-FATAL ERROR', '    FATAL ERROR',
     2            ' ILLOGICAL STOP'/
      DO 110 JERR=1,LINERR
         IF (CHERR(JERR).EQ.CHSUBP .AND. IERRNO(JERR).EQ.NUMBER
     1       .AND. LEVERR(JERR).EQ.ILEVEL) THEN
            NERROR(JERR)=NERROR(JERR)+1
            GO TO 120
         END IF
  110 CONTINUE
c     ------------------------------------------------------------------------
c     Normal case.
c     ------------------------------------------------------------------------
      LINERR=LINERR+1
      NERROR(LINERR)=1
      CHERR(LINERR)=CHSUBP
      IERRNO(LINERR)=NUMBER
      LEVERR(LINERR)=ILEVEL
c     ------------------------------------------------------------------------
c     Illogical stop.
c     ------------------------------------------------------------------------
  120 IF (IABS(ILEVEL).LT.1 .OR. IABS(ILEVEL).GT.5) THEN
 5000    FORMAT (/' Illogical Stop In ERRMES.  ',A6,I3)
         IF (LPRINT .GT. 0) WRITE (LPRINT,5000) CHSUBP, NUMBER
         WRITE (*,5000) CHSUBP, NUMBER
         LINERR=LINERR+1
         NERROR(LINERR)=1
         CHERR(LINERR)='ERRMES'
         IERRNO(LINERR)=1
         LEVERR(LINERR)=5
         IF (ILEVEL .LT. 0) THEN
 5050       FORMAT ('***  ',A15,2X,A6,I3,/
     1              '     This error occurred before the plot could ',
     2              'be produced.'/
     3              '     See the Diagnostics list in the LCModel ',
     4              'Manual.  ***')
            IF (LPRINT .GT. 0) WRITE (LPRINT,5050) ERTYPE(IABS(ILEVEL)),
     1                                             CHSUBP, NUMBER
            WRITE (*,5050) ERTYPE(IABS(ILEVEL)), CHSUBP, NUMBER
            STOP
         END IF
         CALL EXITPS (.true.)
      END IF
c     ------------------------------------------------------------------------
c     Normal case.
c     ------------------------------------------------------------------------
 5100 FORMAT (1X,A15,2X,A6,I3,'.   (Check LCModel Manual).  ',36(2H**)/)
      IF (LPRINT .GT. 0) WRITE (LPRINT,5100) ERTYPE(IABS(ILEVEL)),
     1                                       CHSUBP, NUMBER
      NERMES=NERMES+1
c     ------------------------------------------------------------------------
c     Maximum error messages.
c     ------------------------------------------------------------------------
      IF (NERMES .GE. MERMES) THEN
 5800    FORMAT (//' Stopping after',I5,' diagnostic messages.')
         IF (LPRINT .GT. 0) WRITE (LPRINT,5800) MERMES
         LINERR=LINERR+1
         NERROR(LINERR)=1
         CHERR(LINERR)='ERRMES'
         IERRNO(LINERR)=2
         LEVERR(LINERR)=4
         IF (ILEVEL .LT. 0) THEN
            IF (LPRINT .GT. 0) WRITE (LPRINT,5050) ERTYPE(IABS(ILEVEL)),
     1                                             CHSUBP, NUMBER
            WRITE (*,5050) ERTYPE(IABS(ILEVEL)), CHSUBP, NUMBER
            STOP
         END IF
         CALL EXITPS (.true.)
      END IF
c     ------------------------------------------------------------------------
c     Abort on FATAL.
c     ------------------------------------------------------------------------
      IF (IABS(ILEVEL) .GE. 4) THEN
         IF (ILEVEL .LT. 0) THEN
            IF (LPRINT .GT. 0) WRITE (LPRINT,5050) ERTYPE(IABS(ILEVEL)),
     1                                             CHSUBP, NUMBER
            WRITE (*,5050) ERTYPE(IABS(ILEVEL)), CHSUBP, NUMBER
            STOP
         END IF
         CALL EXITPS (.true.)
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INITIA ()
C
C  Initialization at start of run.
C
      INCLUDE 'lcmodel.inc'
      EXTERNAL DIFF
      DOUBLE PRECISION DIFF
      CHSUBP='INITIA'
      PI=3.141592654D0
      RADIAN=PI/180.
      NDATA=2*NUNFIL
      FNDATA=FLOAT(NDATA)
c     -------------------------------------------------------------------------
c     NLIN = 65535 + 14 = 65549
c     -------------------------------------------------------------------------
      if (ldwfft .gt. 0) then
         nlin = 2 * nlin + 15
      else
         nlin = nlin + 14
      end if
      PPMINC=DELTAT*FNDATA*HZPPPM
      IF (PPMINC .LE. 0.) CALL ERRMES (1, 4, CHSUBP)
      PPMINC=1./PPMINC
      IF (PPMINC.LE.0. .OR. NUNFIL.GT.MUNFIL) THEN
         IF (LPRINT .GT. 0) THEN
            J=MUNFIL
            WRITE (LPRINT, 5110) NUNFIL, J
 5110       FORMAT (//' NUNFIL =', I6, ';   it cannot exceed', I6/)
         END IF
         CALL ERRMES (2, 4, CHSUBP)
      END IF
      LDATST=NINT((PPMCEN-PPMST)/PPMINC)+1+nunfil
      LDATEN=NINT((PPMCEN-PPMEND)/PPMINC)+1+nunfil
      NY=LDATEN-LDATST+1
      IF (LDATST.LE.0 .OR. LDATEN.GT.ndata) CALL ERRMES (3, 4, CHSUBP)
      IF (LDATST.GE.LDATEN .OR. NY.GT.MY) CALL ERRMES (4, 4, CHSUBP)
      DELPPM(1)=-PPMINC*FLOAT(LDATST-1-nunfil)
      ppm(1) = ppmst
      DO 110 JY=2,NY
         DELPPM(JY)=DELPPM(JY-1)-PPMINC
         ppm(jy) = ppm(jy - 1) - ppminc
  110 CONTINUE
C     -------------------------------------------------------------------------
c     Count & order PPMGAP (in the usual descending order in ppm)
c     LCY_SKIP(JY) = T if CY(JY) is to be skipped because it is in a gap.
C     -------------------------------------------------------------------------
      ngap = 0
      do 112 jgap = 1, mgap
         gmax = amax1(ppmgap(1, jgap), ppmgap(2, jgap))
         if (gmax .gt. 9.e36) go to 115
         ngap = ngap + 1
         gmin = amin1(ppmgap(1, jgap), ppmgap(2, jgap))
         ppmgap(1, jgap) = gmax
         ppmgap(2, jgap) = gmin
         if (jgap .gt . 1) then
            if (gmax .ge. ppmgap(2, jgap - 1))
     1         call errmes (9, 4, chsubp)
         end if
         do 114 jy = 1, ny
            lcy_skip(jy) = lcy_skip(jy)   .or.
     1                     (ppm(jy) .le. gmax   .and.
     2                      ppm(jy) .ge. gmin)
 114     continue
 112  continue
C     -------------------------------------------------------------------------
c     NYUSE = # points in Analysis Window actually used (i.e., with gaps
c             removed).
c     Compact array PPM to NYUSE array with the gaps removed.
C     -------------------------------------------------------------------------
 115  nyuse = 0
      do 116 jy = 1, ny
         if (lcy_skip(jy)) go to 116
         nyuse = nyuse + 1
         ppm(nyuse) = ppm(jy)
 116  continue
      if (nyuse .lt. 64) call errmes(8, 3, chsubp)
C     -------------------------------------------------------------------------
C     Convert SDMSHF from ppm to radians/s.
C     -------------------------------------------------------------------------
      SDMSHF=SDMSHF*2.*PI*HZPPPM
      ALOG2=ALOG(2.)
      TOFWHM=SQRT(4.*ALOG2)/(PI*HZPPPM)
c     -------------------------------------------------------------------------
c     NPAR = 2 * 50 + 20 = 120
c     -------------------------------------------------------------------------
      if (ldwfft .gt. 0) then
         npar = 4 * npar + 60
      else
         npar = 2 * npar + 20
      end if
      PRECIS=1.D-6
      DO 120 J=1,100
         IF (DIFF(1.D0+PRECIS, 1.D0) .LE. 0.D0) GO TO 125
         PRECIS=.1D0*PRECIS
  120 CONTINUE
      CALL ERRMES (5, 5, CHSUBP)
  125 PRECIS=100.D0*PRECIS
      IF (LPRINT .GT. 0) WRITE (LPRINT,5125) PPMINC, NYuse, PRECIS
 5125 FORMAT (///10X, ' Other quantities'/
     1           ' PPMINC =', 1PE10.3/
     2           ' NY =', I5/
     3           ' PRECIS =', E9.2)
      IF (PRECIS .GE. .99D-10) CALL ERRMES (6, 2, CHSUBP)
      IF (PRECIS .GE. .99D-5) CALL ERRMES (7, 4, CHSUBP)
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DATAIN ()
C
      INCLUDE 'lcmodel.inc'
      parameter (mdevx=400)
      save havlin, wdone
      logical havlin, wdone
C      CHARACTER CHLI(MDEV)*6!OSF
C      INTEGER IGETAD!OSF
      INTEGER LDEV(MDEV), ldevx(mdevx)
       EXTERNAL ILEN!Cyg
C     EXTERNAL ILEN!sun
C      EXTERNAL ILEN!IRIX
C      EXTERNAL ILEN, IGETAD!OSF
      data havlin/.false./, wdone/.false./
      CHSUBP='DATAIN'
C     -------------------------------------------------------------------------
C     Get time-domain data, scale them, and get non-H2O-suppressed time-domain
C       data.
C     -------------------------------------------------------------------------
      if (iaverg .le. 0) CALL MYDATA ()
C     -------------------------------------------------------------------------
C     ISTAGO = 1 when DATAT is loaded with time-domain data.
C     -------------------------------------------------------------------------
      ISTAGO=1
      if (skip_voxel) return
c     ------------------------------------------------------------------------
c     Change spectral data to absolute value.  The zero-filling below is
c       necesary to avoid huge high-frequency oscillations in spectrum.
c       Presumably zero-filling brings in enough info, so that these symmetric
c       parts can be discarded.
c     The 2.0 * CABS below is to compensate for discarding the
c       symmetric/antisymmetric 2nd half in points NUNFIL+1--NDATA.
c     1st time-data point (the integral over the spectrum) must be halved,
c       to remove the offset produced by halving the number of points.  See
c       Notes of 010819.
c     ------------------------------------------------------------------------
      if (absval) then
         do 130 j = nunfil+1, ndata
            datat(j) = 0.
 130     continue
         call cfft (datat, dataf, ndata, lwfft, wfftc)
         do 135 j = 1, ndata
            dataf(j) = 2. * cabs(dataf(j))
 135     continue
         call cfftin (dataf, datat, ndata, lwfft, wfftc)
         datat(1) = .5 * datat(1)
      end if
C
C     ----------------------------------------------------------------------
C     Test license.
C     LI = host identification number initially; later this is modifed by
C          OWNER information.
C     LLINE = T if the run is to be aborted because of no license, and this is
c               not the Test Data.
c     COUNTC used instead of GETID for routine name for getting ID.
c     NOLINE = T if there is no license; used instead of NOLIC for security.
C     ----------------------------------------------------------------------
      CALL FTDATA(0)
c     -----------------------------------------------------------------------
c     Skip license tests with Intel Windows version & dongle.
c     -----------------------------------------------------------------------
       if (nlin .gt. 0) go to 200!Cyg
c     -----------------------------------------------------------------------
c     Check for Master KEY.
c     -----------------------------------------------------------------------
      if (key(1) .eq. 210387309) go to 200
      LLINE=.FALSE.
      NOLINE=.FALSE.
c     -------------------------------------------------------------------------
c     The numerical year & month below should be updated to foil backdating.
c     -------------------------------------------------------------------------
      IF (KNUM(3) .LT. 100) KNUM(3)=KNUM(3)+2000
      if (knum(3) .lt. 2018   .or.
     1    (knum(3) .eq. 2018   .and.   knum(2) .lt. 7)) go to 195
c     -------------------------------------------------------------------------
c     GETFFT produces a segmentation fault if called more than once (although
c        NPAR & NLIN were verified to be correct.
c     Therefore, HAVLIN=T (set below) if Linux valid license was found, and
c        this allows skipping a 2nd test.
c     -------------------------------------------------------------------------
      if (havlin) go to 200
c     -------------------------------------------------------------------------
c     Initialize LDEVX, since some may be commented out in DEVX*.INC
c     -------------------------------------------------------------------------
      do 250 j = 1, mdevx
         ldevx(j) = -2
 250  continue
      ndevx = 999999
c     ------------------------------------------------------------------------
c     The following statement will never be reached with the Cyg flag, but
c        the start of the IF must be there to avoid diagnostics with Sun,
c        etc.  Should use INCLUDEs instead.
c     -------------------------------------------------------------------------
C     if (nlin .gt. 0) then!Linux
       if (nlin .gt. 0) then!Cyg
C      if (nlin .lt. 0) then!sun
C      if (nlin .lt. 0) then!IRIX
C      if (nlin .lt. 0) then!OSF
c        --------------------------------------------------------------------
c        LETT = STATUS for ERRMES from GETFFT (GETLIC) = 0 for valid license.
c        NPAR = 120 = NSTART in GETLIC
c        NLIN = 65549 = ARG1 in GETLIC
c        --------------------------------------------------------------------
         len_owner = 0
C        call getfft(npar, nlin, lett, owner, len_owner, knum3, !Linux
C    1               knum2, idevx, version_lcm)!Linux
         if (len_owner .gt. 0) then
c           -----------------------------------------------------------------
c           The following 2 lines are a clumsy way to properly fill OWNER with
c              blanks, avoiding the output of a binary character in the
c              Namelist output of OWNER.
c           -----------------------------------------------------------------
            ownout = OWNER(1:len_owner) // ' '
            owner = ownout
            OWNOUT='Data of: ' // OWNER(1:len_owner)
            JOWNER=ILEN(OWNOUT)
            JDATE=ILEN(CHDATE)
            if (ltable .gt. 0)
     1         write (ltable, 5410) OWNOUT(1:JOWNER), CHDATE(1:JDATE)
 5410       format(1x, a)
            IF (LCOORD .GT. 0)
     1         write (lcoord, 5410) OWNOUT(1:JOWNER), CHDATE(1:JDATE)
         end if
C        include 'devx_linux.inc'!Linux
         do 410 j = 1, ndevx
            if (idevx .eq. ldevx(j)) lett = 40
 410     continue
         if (lett .eq. 0) then
            if (knum(3) .gt. knum3   .or.
     1          (knum(3) .eq. knum3   .and.   knum(2) .gt. knum2))
     2           lett = 50
         end if
         havlin = lett .eq. 0
         if (havlin) go to 200
      else
         NDEV=1
         L=0
C         CALL countc (L)!sun
C         CALL countc (L)!IRIX
         LDEV(1)=L
C        ------------------------------------------------------------------
C        Initialize LDEV1 to only output NDEV LDEV1 values in One-Page Output.
C        ------------------------------------------------------------------
         DO 185 JDEV=2,MDEV
            LDEV1(JDEV)=-1
 185     CONTINUE
C        ------------------------------------------------------------------
C        OSF IGETAD can return NDEV>1 LDEV values.
C        ------------------------------------------------------------------
C         MDEVAR=MDEV!OSF
C         NDEV=IABS(IGETAD (CHLI, MDEVAR))!OSF
C         DO 186 JDEV=1,NDEV!OSF
C            LDEV(JDEV)=77*ICHAR(CHLI(JDEV)(1:1))+!OSF
C     1                 7*ICHAR(CHLI(JDEV)(2:2))+!OSF
C     2                 26*ICHAR(CHLI(JDEV)(3:3))+!OSF
C     3                 65536*ICHAR(CHLI(JDEV)(4:4))+!OSF
C     4                 256*ICHAR(CHLI(JDEV)(5:5))+!OSF
C     5                 ICHAR(CHLI(JDEV)(6:6))!OSF
C  186    CONTINUE!OSF
C         NDEV=MAX0(1, NDEV)!OSF
C
C         include 'devx.inc'!IRIX
C         include 'devx.inc'!OSF
C         include 'devx.inc'!sun
         DO 187 JDEV=1,NDEV
            LDEV1(JDEV) = IGETP (34481 + IABS(LDEV(JDEV)), 35)
c           -------------------------------------------------------------
c           Reject license if LDEV1 is on blacklist in LDEVX
c           -------------------------------------------------------------
            do 1875 jdevx = 1, ndevx
               if (LDEV1(JDEV) .eq. ldevx(jdevx)) go to 195
 1875       continue
            LDEV(JDEV)=LDEV1(JDEV)
C        ----------------------------------------------------------------
C        Modify LDEV(JDEV) with OWNER info.
C        ----------------------------------------------------------------
            DO 188 J=1,ILEN(OWNER)
               LDEV(JDEV)=LDEV(JDEV)-(J+9)*ICHAR(OWNER(J:J))
 188        CONTINUE
            LDEV(JDEV)=IABS(LDEV(JDEV))
C        ----------------------------------------------------------------
C        Test LDEV(JDEV) for lifetime license.
C        ----------------------------------------------------------------
            K = IGETP (LDEV(JDEV) + 8829, 59)
            DO 189 J = 1, MKEY
               IF (KEY(J) .EQ. K) GO TO 200
 189        CONTINUE
C           -------------------------------------------------------------
C           Test host identification number for license for next 25 months.
C           KNUM(2) = month
C           KNUM(3) = year
C           -------------------------------------------------------------
            KJ = 100*(KNUM(3) - 1997)
            KM = KNUM(2) - 1
            DO 191 JPAR = 1, 25
               KM = KM + 1
               IF (KM .GT. 12) THEN
                  KM = KM - 12
                  KJ = KJ + 100
               END IF
               K = IGETP (LDEV(JDEV) + (KJ + KM)*3678, 41)
               DO 193 J = 1, MKEY
                  IF (KEY(J) .EQ. 0) GO TO 191
                  IF (KEY(J) .EQ. K) GO TO 200
 193           CONTINUE
 191        CONTINUE
 187     CONTINUE
      end if
C     ----------------------------------------------------------------------
C     There is no license.  Check if this is demo data.
C     ----------------------------------------------------------------------
 195  NOLINE=.TRUE.
      PPM1=3.4
      PPM2=2.8
      SUML=0.
      DO 9210 J=NINT((4.65-PPM1)/PPMINC)+1+nunfil,
     1     NINT((4.65-PPM2)/PPMINC)+1+nunfil
         SUML=AMAX1(SUML,CABS(DATAF(J)))
 9210 CONTINUE
C
      PPM1=2.8
      PPM2=2.2
      SUMM=0.
      DO 9220 J=NINT((4.65-PPM1)/PPMINC)+1+nunfil,
     1     NINT((4.65-PPM2)/PPMINC)+1+nunfil
         SUMM=AMAX1(SUMM,CABS(DATAF(J)))
 9220 CONTINUE
C
      PPM1=2.2
      PPM2=1.8
      SUMR=0.
      DO 9230 J=NINT((4.65-PPM1)/PPMINC)+1+nunfil,
     1     NINT((4.65-PPM2)/PPMINC)+1+nunfil
         SUMR=AMAX1(SUMR,CABS(DATAF(J)))
 9230 CONTINUE
C
      TEST=SUML-SUMM
      LLINE=ABS(TEST).LE.0.
      IF (.NOT.LLINE) THEN
         TEST=(SUMR-SUMM)/TEST
C         write (*, 5230) test!delete
C5230     format (1pe15.6)!delete
C         if (nunfil .gt. 0) stop!delete
C         LLINE=MOD(NINT(3.*TEST),2) .EQ. 0
          LLINE=TEST.LT.2.16 .OR. TEST.GT.2.17
      END IF
C     -------------------------------------------------------------------------
c     NOLINE = T if there is no license.
C     LLINE = T if the run is to be aborted because of no license, and this is
c               not the Test Data.
c     fndata = 0. causes the run to be aborted in MYBASI.  This will sabotage
c                 the run, even if the abort is stopped.
C     -------------------------------------------------------------------------
      if (lline) fndata = 0.
 200  RETURN
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER FUNCTION ICYCLE_r (J, NDATA)
C
C  Computes the subscript from J for a rearranged array of length NDATA.
c  Simply sets extensions beyond array to the endpoint.
C
      ICYCLE_R=max0(1,min0(ndata,j))
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER FUNCTION ICYCLE (J, NDATA)
C
C  Computes the subscript from J for a cyclic array of length NDATA.
C  J > -NDATA must be satisfied.
C
      CHARACTER CHSUBP*6
      CHSUBP='ICYCLE'
      K=J-1+NDATA
      IF (K .LT. 0) CALL ERRMES (1, 4, CHSUBP)
      ICYCLE=MOD(K,NDATA)+1
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MYDATA ()
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C  Read the NUNFIL complex time-domain data into DATAT.
C  Scale DATAT.
C  Read the NUNFIL complex time-domain non-H2O-suppressed data into H2OT.
c  If DOECC, then do ECC.
C
C  If your compiler doesn't accept NAMELIST, then you will have to modify
C    this subprogram to use only standard I/O statements.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
      INCLUDE 'lcmodel.inc'
      save BRUKER_h2o, BRUKER_raw, fmtdat_h2o, fmtdat_raw, id,
     1     scale_h2o, scale_raw, seqacq_h2o, seqacq_raw
      CHARACTER FMTDAT*(MCHFmt), FMTDAT_h2o*(MCHFmt),
     1          FMTDAT_raw*(MCHFmt), ID*(MCHID), seq*5
      LOGICAL bruker, BRUKER_h2o, BRUKER_raw,
     1        seqacq, seqacq_h2o, seqacq_raw
      namelist /seqpar/ hzpppm, echot, seq
      NAMELIST /NMID/ ID, FMTDAT, TRAMP, VOLUME, SEQACQ, BRUKER
      CHSUBP='MYDATA'
      if (voxel1) then
c        ---------------------------------------------------------------------
c        If not present in seqpar, then echot_raw = -1., and seq_raw = ' '.
c        ---------------------------------------------------------------------
         hzpppm_sav=hzpppm
         hzpppm=-1.
         seq=' '
         echot=-1.
         read (lraw, nml=seqpar, err=102, end=102)
 102     rewind lraw
         if (hzpppm .gt. 0.) then
            if (abs(hzpppm_sav / hzpppm - 1.) .gt. .05) call errmes
     1                                                  (4, 2, chsubp)
         end if
         hzpppm=hzpppm_sav
         echot_raw=echot
         seq_raw=seq
      end if
      if (voxel1   .or.   lraw_at_top) then
c        ---------------------------------------------------------------------
c        Read NMID.
c        ---------------------------------------------------------------------
         tramp = 1.
         volume = 1.
         id = ' '
         FMTDAT=' '
         BRUKER=.FALSE.
         SEQACQ=.FALSE.
         READ (LRAW, NML=NMID, err=807, end=807)
         IF (FMTDAT .EQ. ' ') CALL ERRMES (1, 4, CHSUBP)
         FMTDAT_RAW=fmtdat
         BRUKER_RAW=bruker
         SEQACQ_RAW=seqacq
      end if
      READ (LRAW, FMTDAT_RAW, err=808, end=808) (DATAT(J),J=1,NUNFIL)
      IF (LPRINT .GT. 0) then
         if (nlines_title .eq. 1) then
            WRITE (LPRINT,5110) TITLE_line(1), ID
 5110       FORMAT (1X,A//' Data set ID = ',A////)
         else
            WRITE (LPRINT,5112) TITLE_line, ID
 5112       FORMAT (1X,A/1x,a//' Data set ID = ',A////)
         end if
      end if
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     TRAMP and VOLUME are the transmitter amplitude and VOI used for scaling
C       the spectra to get absolute concentrations, as discussed in MRM.  If
C       you do not have this data, then you can leave TRAMP and VOLUME out of
C       NAMELIST NMID, and leave them at their default values of
C       TRAMP=VOLUME=1; the spectrum data will remain unscaled (with SCALE=1.)
C       below.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
      if (voxel1) then
         IF (AMIN1(TRAMP,VOLUME,fcalib) .LE. 0.)
     1      CALL ERRMES (2, 4, CHSUBP)
         SCALE_RAW=fcalib*TRAMP/VOLUME
      end if
      if (skip_voxel) go to 190
      DO 110 JDATA=1,NUNFIL
         DATAT(JDATA)=DATAT(JDATA)*SCALE_RAW
  110 CONTINUE
      if (fwhmsm .gt. 0.) then
         RSD = PI* fwhmsm * deltat * hzpppm / sqrt(2.*alog(2.))
         DO 112 JUNFIL=1,NUNFIL
            DATAt(JUNFIL)=DATAT(JUNFIL)*
     1                    EXP(-.5*(RSD*FLOAT(JUNFIL-1))**2)
 112     CONTINUE
      endif
 
      IF (SEQACQ_RAW) THEN
         CALL SEQTOT (DATAT, DATAF, NUNFIL, LWFFT, WFFTC)
      ELSE IF (BRUKER_RAW) THEN
         DO 115 JDATA=1,NUNFIL
            DATAT(JDATA)=CONJG(DATAT(JDATA))
  115    CONTINUE
      END IF
      if (smtail) call smooth_tail(datat)
 190  havh2o = .false.
      IF (filh2o .ne. ' ') THEN
         if (voxel1) then
C           -------------------------------------------------------------------
C           Read time-domain non-H2O-suppressed data into H2OT.
C           -------------------------------------------------------------------
             OPEN (LH2O, FILE=FILH2O, STATUS='OLD', err=800)!Cyg
C           OPEN (LH2O, FILE=FILH2O, STATUS='OLD', err=800) !sun
C            OPEN (LH2O, FILE=FILH2O, STATUS='OLD', READONLY, err=800)!OSF
C            OPEN (LH2O, FILE=FILH2O, STATUS='OLD', READONLY, err=800)!IRIX
         end if
         if (voxel1   .or.   lraw_at_top) then
            FMTDAT=' '
            BRUKER=.FALSE.
            SEQACQ=.FALSE.
            READ (LH2O, NML=NMID, err=809, end=809)
            IF (FMTDAT .EQ. ' ') CALL ERRMES (3, 4, CHSUBP)
            FMTDAT_H2O=fmtdat
            BRUKER_H2O=bruker
            SEQACQ_H2O=seqacq
         end if
         READ (LH2O, FMTDAT_H2O, err=810, end=810) (H2OT(J),J=1,NUNFIL)
         havh2o = .true.
         if (voxel1) then
            if (amin1(tramp, volume) .le. 0.) then
c              ----------------------------------------------------------------
c              This could cause an error in water-scaling, if water-reference
c                & suppresed TRAMP & VOLUME do not match.
c              ----------------------------------------------------------------
               if (dows) call errmes (11, 2, chsubp)
               scale_h2o = scale_raw
            else
               SCALE_H2O=TRAMP/VOLUME
            end if
         end if
         if (skip_voxel) go to 800
         DO 210 JDATA=1,NUNFIL
            H2OT(JDATA)=H2OT(JDATA)*SCALE_H2O
 210     CONTINUE
         IF (SEQACQ_H2O) THEN
            CALL SEQTOT (H2OT, DATAF, NUNFIL, LWFFT, WFFTC)
         ELSE IF (BRUKER_H2O) THEN
            DO 215 JDATA=1,NUNFIL
               H2OT(JDATA)=CONJG(H2OT(JDATA))
  215       CONTINUE
         END IF
         if (smtail) call smooth_tail(h2ot)
         if (fwhmsm .gt. 0.) then
            RSD = PI* fwhmsm * deltat * hzpppm / sqrt(2.*alog(2.))
            DO 220 JUNFIL=1,NUNFIL
               h2ot(JUNFIL)=h2ot(JUNFIL)*
     1                      EXP(-.5*(RSD*FLOAT(JUNFIL-1))**2)
 220        CONTINUE
         endif
c        ----------------------------------------------------------------------
c        UNSUPR = T to use H2OT as DATAT.  Cannot be used with phased arrays,
c                   because DOECC=F would sum incoherently.
c        ----------------------------------------------------------------------
         if (unsupr) then
            doecc = .false.
            do 250 JDATA=1,NUNFIL
               datat(jdata) = H2OT(JDATA)
 250        CONTINUE
         end if
      END IF
      if (havh2o   .and.   biglip   .and.   iaverg .gt. 0   .and.
     1    .not.doecc) then
         call phase_with_max_real ()
      else if (doecc_active   .and.   .not.doecc) then
         doecc_active = .false.
         sddegz = amin1(999., 333. * sddegz)
      end if
      IF (DOECC   .and.   havh2o   .and.   (iaverg .ne. 3   .and.
     1     iaverg .ne. 31   .and.   iaverg .ne. 32))
     2   call ecc_truncate ()
 800  if (.not.havh2o   .and.   (doecc .or. dows)) then
         call errmes (5, 3, chsubp)
         doecc = .false.
         dows = .false.
      end if
      if (.not.havh2o   .and.   unsupr) call errmes (12, 4, chsubp)
      RETURN
 807  call errmes (7, 4, chsubp)
 808  call errmes (8, 4, chsubp)
 809  call errmes (9, 4, chsubp)
 810  call errmes (10, 4, chsubp)
      END
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine  phase_with_max_real ()
c
c Phase (zero-order)  DATAT & H2OT so that integral from ppmst_phalip &
c    ppmend_phalip of real part of H2OF is max.
c Assumes that H2OT is positive (as with lipid spectra).
c
      include 'lcmodel.inc'
      complex cfact, cfact_best, delta_cfact
      logical sddegz_done
      save sddegz_done
      data sddegz_done/.false./
      chsubp = 'PHALIP'
c     ------------------------------------------------------------------------
c     Account for increased uncertainty in DEGZER compared to using ECC.
c     ------------------------------------------------------------------------
      if (.not.sddegz_done) then
         sddegz_done = .true.
         sddegz = 2. * sddegz
      end if
c     ------------------------------------------------------------------------
c     PPMINC, NDATA, FNDATA & RADIAN initialized in AVERAGE.
c     ------------------------------------------------------------------------
      do 150 jdata = nunfil + 1, ndata
         h2ot(jdata) = cmplx(0., 0.)
 150  continue
      call csft_r (h2ot, h2of, ndata)
c     ------------------------------------------------------------------------
c     Increment is 1 degree
c     ------------------------------------------------------------------------
      delta_cfact = cexp(cmplx(0., radian))
      cfact = 1.
      cfact_best = cfact
      ldatst = max0(1, NINT((PPMCEN - PPMST_phalip)/PPMINC)+1+nunfil)
      LDATEN = min0(ndata,
     1              nint((PPMCEN - PPMEND_phalip)/PPMINC)+1+nunfil)
      if (ldaten .le. ldatst) call errmes (1, 4, chsubp)
      sum_best = -1.e30
      ldeg = -1
      do 210 jdeg = 0, 359
         sum = 0.
         do 250 jdata = ldatst, ldaten
            sum = sum + real(h2of(jdata) * cfact)
 250     continue
         if (sum .ge. sum_best) then
            sum_best = sum
            cfact_best = cfact
            ldeg = jdeg
         end if
         cfact = cfact * delta_cfact
  210  continue
      do 310 junfil = 1, nunfil
         datat(junfil) = datat(junfil) * cfact_best
         h2ot(junfil) = h2ot(junfil) * cfact_best
 310  continue
c      write (6, 5310) ldeg
c 5310 format ('zero-order phase =', i4)
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smooth_tail(cdatat)
c
c Smooth tail to get rid of oscillations at end with Siemens (and sometimes
c   Bruker) time-domain data.
c CDATAT = DATAT or H2OT
c
      include 'lcmodel.inc'
      complex cdatat(munfil)
      real out_imag(munfil), out_real(munfil), work_in(munfil)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common out_imag, out_real, work_in
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 110 j = 1, nunfil
         work_in(j) = real(cdatat(j))
 110  continue
      call smooth_tail_2(work_in, out_real,
     1                   munfil, nunfil, lprint, voxel1)
 
      do 210 j = 1, nunfil
         work_in(j) =  aimag(cdatat(j))
 210  continue
      call smooth_tail_2(work_in, out_imag,
     1                   munfil, nunfil, lprint, voxel1)
 
      do 310 j = 1, nunfil
         cdatat(j) = cmplx(out_real(j), out_imag(j))
 310  continue
C      write(6, 9310) (cdatat(j), j = 1, nunfil)!DELETE
C 9310 format(1p2e15.6)!DELETE
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine smooth_tail_2(work_in, out,
     1                         munfil, nunfil, lprint, voxel1)
c
c Average over next neighbors as long as zig-zag pattern holds.
c j = point at which pattern is broken; out(j) is not computed.
c
      logical voxel1
      real out(munfil), work_in(munfil)
c
      do 110 j = nunfil - 1, 2, -1
         if ((work_in(j + 1) - work_in(j)) *
     1       (work_in(j) - work_in(j - 1)) .gt. 0.) go to 210
         out(j) = .5 * work_in(j) +
     1            .25 * (work_in(j + 1) + work_in(j - 1))
 110  continue
      j = 1
 210  if (j .lt. nunfil - 1) then
         out(nunfil) = out(nunfil - 1)
      else
         out(nunfil) = work_in(nunfil)
      end if
      do 220 junfil = 1, j
         out(junfil) = work_in(junfil)
 220  continue
      if (lprint. gt. 0   .and.   voxel1)
     1      write(lprint, 5210) nunfil - j
 5210 format(/i4, ' points in tail smoothed')
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ecc_truncate ()
      include 'lcmodel.inc'
      real awork(munfil)
      complex h2ot_ecc(munfil)
      common awork, h2ot_ecc
c
      chsubp = 'ECC'
      PPMINC2 = DELTAT * float(nunfil) * HZPPPM
      IF (PPMINC2 .LE. 0.) CALL ERRMES (1, 4, CHSUBP)
      PPMINC2=1./PPMINC2
      nunfil_half = nunfil / 2
      if (ppm_truncate_max .lt. 0.   .or.
     1    ppm_truncate_min .gt. ppmh2o) then
         do 130 junfil = 1, nunfil
            h2ot_ecc(junfil) = h2ot(junfil)
 130     continue
      else
         lmin = nunfil_half + 1 + (ppmcen - ppm_truncate_max) / ppminc2
         lmax = nunfil_half + 1 + (ppmcen - ppm_truncate_min) / ppminc2
         if (lmin .lt. 3   .or.   lmax .gt. ndata -2   .or.
     1       lmin .ge. lmax) call errmes (2, 4, chsubp)
         call csft_r(h2ot, h2of_work, nunfil)
         do 240 j = 1, nunfil
            awork(j) = cabs(h2of_work(j))
 240     continue
         absmin = 1.e30
         ltruncate = 0
         do 250 j = lmin, lmax
c           ------------------------------------------------------------------
c           TEST = least-squares line at j (Hildebrand, p 295)
c           ------------------------------------------------------------------
            test = awork(j-2) + awork(j-1) + awork(j) +
     1             awork(j+1) + awork(j+2)
            if (test .lt. absmin) then
               absmin = test
               ltruncate = j
            end if
 250     continue
         ppm_truncate = ppmcen - (ltruncate - nunfil_half -1) * ppminc2
         if (lprint .gt. 0) write (lprint, 5250) ppm_truncate
 5250    format (/'Truncated at ', f6.3, ' ppm for ECC')
         do 270 j = ltruncate + 1, nunfil
            h2of_work(j) = (0.,0.)
 270     continue
         call csftin_r (h2of_work, h2ot_work, h2ot_ecc, nunfil)
      end if
C     -------------------------------------------------------------------------
C     Klose's eddy-current correction.
C     -------------------------------------------------------------------------
      DO 310 JUNFIL=1,NUNFIL
         XX=REAL(H2OT_ecc(JUNFIL))
         YY=AIMAG(H2OT_ecc(JUNFIL))
         ABSSQ=XX*XX+YY*YY
         IF (ABSSQ .gt. 0.) THEN
C           -------------------------------------------------------------------
C           Zero points can occur when SEQACQ=T.  This will only be at the
C             end of the t-range, and is therefore not serious.
C           -------------------------------------------------------------------
            DATAT(JUNFIL)=DATAT(JUNFIL)*CMPLX(XX,-YY)/SQRT(ABSSQ)
         END IF
 310  CONTINUE
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION IGETP (ISTART, NITER)
C
C  Uses RANDOM to encode ISTART and produce IGETP (max. 9 digits) for hostid
C    tests.
C
C  RANDOM PRODUCES A PSEUDORANDOM REAL ON THE OPEN INTERVAL (0.,1.).
C  DIX (IN DOUBLE PRECISION) will BE INITIALIZED TO A WHOLE NUMBER
C      BETWEEN 1.D0 AND 2147483646.D0 BEFORE THE FIRST CALL TO RANDOM
C      AND NOT CHANGED BETWEEN SUCCESSIVE CALLS TO RANDOM.
C  BASED ON L. SCHRAGE, ACM TRANS. ON MATH. SOFTWARE 5, 132 (1979).
C-----------------------------------------------------------------------
C
C  PORTABLE RANDOM NUMBER GENERATOR
C   USING THE RECURSION
C    DIX = DIX*A MOD P
C
      DOUBLE PRECISION A,P,DIX,B15,B16,XHI,XALO,LEFTLO,FHI,K
C
C  7**5, 2**15, 2**16, 2**31-1
 
      DATA A/16807.D0/,B15/32768.D0/,B16/65536.D0/,P/2147483647.D0/
C
      IF (ISTART .EQ. 8829) THEN
C        ----------------------------------------------------------------------
C        Special return for hostid=0
C        ISTART = hostid + 8829 must be in call from LCModel.
C        ----------------------------------------------------------------------
         IGETP = 0
         RETURN
      ENDIF
      DIX = DBLE(MAX0(1, MOD(IABS(ISTART), 2147483647)))
      DO 110 JITER = 1, MAX0(1, NITER)
C
C        GET 15 HI ORDER BITS OF DIX
C
         XHI = DIX / B16
         XHI = XHI - DMOD(XHI,1.D0)
C
C        GET 16 LO BITS IF DIX AND FORM LO PRODUCT
C
         XALO=(DIX-XHI*B16)*A
C
C        GET 15 HI ORDER BITS OF LO PRODUCT
C
         LEFTLO = XALO/B16
         LEFTLO = LEFTLO - DMOD(LEFTLO,1.D0)
C
C        FORM THE 31 HIGHEST BITS OF FULL PRODUCT
C
         FHI = XHI*A + LEFTLO
C
C        GET OVERFLO PAST 31ST BIT OF FULL PRODUCT
C
         K = FHI/B15
         K = K - DMOD(K,1.D0)
C
C        ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C        THE PARENTHESES ARE ESSENTIAL
C
         DIX = (((XALO-LEFTLO*B16) - P) + (FHI-K*B15)*B16) + K
C
C        ADD P BACK IN IF NECESSARY
C
         IF (DIX .LT. 0.D0) DIX = DIX + P
  110 CONTINUE
C
C     MULTIPLY BY 1/(2**31-1)
C
      DIX=DIX*4.656612875D-10
      IGETP=1.D9*DIX
      RETURN
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MYBASI (LSTAGE)
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C  Most of this subprogram you should not have to change (or understand).
C    Only those parts set off by the usual lines of "v" and "^" characters are
C    intended for you to possibly change.
C
C  Read the basis set on unit LBASIS (file FILBAS).  There are 3 read
C    statements:
C      (1) NAMELIST BASIS1 (the main header).
C      (2) NAMELIST BASIS (a header for the basis spectrum to be read next)
C      (3) BASISF the basis spectrum.
C      Item (1) is read only once. Items (2) and (3) are read consecutively in
C        a loop until there are no more spectra left.
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
C  Input and scale (also to unit concentrations) frequency-domain basis
C    vectors obtained with zero-filling.
C  Put their time-domain vectors into BASIST(JDATA,JMETAB).
C  LSTAGE = 1 for initial analyses with restricted model.
C         = 2 for analyses with full model.
C
      INCLUDE 'lcmodel.inc'
      external areaba, ilen
      save hzpppm_basis
      PARAMETER (MCHIDB=80)
      CHARACTER file_basout*262, FMTBAS*(MCHFMT), ID*(MCHID),
     1          IDBASI*(MCHIDB), METABO*(MCHMET), seq*5
      COMPLEX BASISF, offset_end, offset_start
      double precision dix
      LOGICAL corr_areaba, dows_now, encryp, ppmmet_in_gap(2),
     1        US1FUL(MMETAB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON BASISF(8*MDATA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      namelist /seqpar/ fwhmba, hzpppm, echot, seq
      NAMELIST /BASIS1/ IDBASI, FMTBAS, BADELT, NDATAB
      NAMELIST /BASIS/ ID, METABO, CONC, TRAMP, VOLUME, ISHIFT
      CHSUBP='MYBASI'
      IF (LSTAGE .EQ. 1) THEN
C        ----------------------------------------------------------------------
C        Abort if no license and no test data; FNDATA = 0. (from DATAIN)
c           causes abort.
C        ----------------------------------------------------------------------
         IF (AMIN1(PPMINC, FNDATA, TOFWHM) .LE.0. .OR.
     1       LDATST.GT.LDATEN .OR. LWFFT.LT.0)
     2       CALL ERRMES (max0(1, lett), 4, 'DEMO  ')
C        ----------------------------------------------------------------------
C        LETT > 0 at this stage is only possible with linux with Test Data.
C        ----------------------------------------------------------------------
         if (lett .gt. 0) call errmes (lett, 1, 'DEMO  ')
      END IF
 
      if ((scafwh   .or.  imethd .eq. 3)   .and.
     1    (.not.nobasi .or. badref)) call errmes (40, 5, chsubp)
C     -------------------------------------------------------------------------
c     RT2_SCALE: SDRT2 & EXRT2 will be scaled by
c                RT2_SCALE = sqrt(HZPPPM / HZPPPM_RT2_REF), assuming
c                that 1/T2 is roughly proportional to sqrt(HZPPPM).
c     Input values of SDRT2 will then pertain to HZPPPM_RT2_REF, which is
c       set to 2T below.
C     -------------------------------------------------------------------------
      RT2_SCALE = sqrt(hzpppm / 85.16)
      IF (LSTAGE .EQ. 1) THEN
         DO 105 JUSE1=1,NUSE1
            US1FUL(JUSE1)=.FALSE.
  105    CONTINUE
      END IF
c     -------------------------------------------------------------------------
c     NONNEG below must be later reset for all JPAR>NMETAB and when BASIST
c            correspond to Taylor expansion terms when TAYLOR=T
c     LSHAPE below must be later reset when convolution with Lineshape is
c            being suppressed.
c     -------------------------------------------------------------------------
      do 107 jmetab = 1, mmetab
         nonneg(jmetab) = .true.
         lshape(jmetab) = .true.
 107  continue
      NMETAB=0
      IF (LPRINT .GT. 0) THEN
         IF (LSTAGE .EQ. 1) THEN
            WRITE (LPRINT, 5105)
 5105       FORMAT (/////20X, 'Basis Spectra used for the ',
     1              'Preliminary Analysis')
         ELSE
            WRITE (LPRINT, 5106)
 5106       FORMAT (/////20X, 'Basis Spectra used for the ',
     1              'Final Analysis')
         END IF
      END IF
c     -------------------------------------------------------------------------
c     NOBASI=T to skip BASIS file and only use simulated spectra.
c     -------------------------------------------------------------------------
      if (nobasi) then
         if (bascal   .or.   ncalib .gt. 0   .or.   nsimul .le. 0)
     1      call errmes (27, 4, chsubp)
         scasim = .false.
c        ----------------------------------------------------------------------
c             Set DESDSH if it has not been input positive or set positive in
c        Block Data.  Of course, PPMINC (not PPMINC_BASIS) must be used here.
c        ----------------------------------------------------------------------
         if (desdsh .le. 0.) then
            if (rincsh .le. 0.   .or. sdshmn .le. 0.   .or.
     1          sdshmx .lt. sdshmn) call errmes (26, 4, chsubp)
            desdsh = amin1(sdshmx, amax1(sdshmn, rincsh * ppminc))
         end if
         IF (LPRINT .GT. 0) WRITE (LPRINT, 5110)
 5110    FORMAT (/' No.   Metabolite',3X,
     1            'Expec[delta(1/T2)](1/s)', 3X, 'SDdelta(1/T2)]',
     2            3X, 'Shift(Pts)', 3X, 'SD[Shift](ppm)   ID')
         go to 300
      end if
C     -------------------------------------------------------------------------
C     BASISF(JDATA) = unshifted, unscaled frequency-domain basis vector with
C                      zero-filling and marker peaks removed.
c     (2cygwin does not left-shift below; others do)
C     -------------------------------------------------------------------------
      if (lbasis .le. 0   .or.   filbas .eq. ' ')
     1      call errmes (33, 4, chsubp)
C     IF (FILBAS .NE. ' ') OPEN (LBASIS, FILE=FILBAS, STATUS='OLD',!sun
C    1                           err=810)!sun
      IF (FILBAS .NE. ' ') OPEN (LBASIS, FILE=FILBAS, STATUS='OLD',!Cyg
     1                           err=810)!Cyg
C      IF (FILBAS .NE. ' ') OPEN (LBASIS, FILE=FILBAS, STATUS='OLD',!OSF
C     1                           READONLY, err=810)!OSF
C      IF (FILBAS .NE. ' ') OPEN (LBASIS, FILE=FILBAS, STATUS='OLD',!IRIX
C     1                           READONLY, err=810)!IRIX
C
c     -------------------------------------------------------------------------
c     BASCAL=T for analyzing the CHBCAL basis spectrum as the data.
c       Since the basis spectra are also divided by CONC, the LCModel
c         concentration should be 1.0, e.g., if CHBCAL=PCh or NAAG and
c         CHCALI(1)=GPC or NAA.
c       When BASCAL=T, anything can be submitted as FILRAW, since it is not
c         read.  With LCMgui, one can use the "Other" fidType and enter
c         HZPPPM, etc., by hand.
c     -------------------------------------------------------------------------
      bascal = bascal .and. lstage.eq.1
      if (bascal) then
         if (ncalib.le.0 .or. chbcal.eq.' ') call errmes (12, 4, chsubp)
         do 108 jcalib = 1, ncalib
            if (chbcal .eq. chcali(jcalib)) call errmes (13, 4, chsubp)
 108     continue
      end if
      if (lstage .eq. 1  .and.  .not.bascal) then
         hzpppm_sav=hzpppm
         hzpppm=-1.
         seq=' '
         echot=-1.
         fwhmba_sav = fwhmba
         fwhmba = -1.
         read (lbasis, nml=seqpar, err=110, end=110)
         go to 120
 110     rewind lbasis
 120     hzpppm_basis = hzpppm
         hzpppm=hzpppm_sav
         if (amin1(echot,echot_raw) .gt. 0.) then
            test=abs(echot - echot_raw)
            if (test .ge. 3.001) call errmes (10, 3, chsubp)
         end if
 
         if (echot .le. temm   .and.   ppmend .gt. .6   .and.
     1       (sptype(:1) .eq. ' '   .or.
     2        sptype(:5) .eq. 'tumor'   .or.
     3        sptype(:6) .eq. 'nulled')) then
c           -----------------------------------------------------------------
c           Modify PPMMET so that MMs are also included, even though MM09
c             (and its ratio priors) is excluded, since especially MM20 is
c             important for NAA, NAAG, Glx, etc.
c           Settings below must be coordinated with BLOCK DATA. <<<<<<<<<<<<
c           -----------------------------------------------------------------
            call errmes (41, 2, chsubp)
C           ppmmet(1, 45) = 1.41 !MM12 @ 1.21; FWHM=.2
            ppmmet(2, 45) = 1.21
            ppmmet(3, 45) = 1.21
            ppmmet(4, 45) = 1.01
C           ppmmet(1, 46) = 1.63 !MM14 @ 1.43; FWHM=.2
            ppmmet(2, 46) = 1.43
            ppmmet(3, 46) = 1.43
            ppmmet(4, 46) = 1.23
C           ppmmet(1, 47) = 1.84 !MM17 @ 1.67; FWHM=.17
            ppmmet(2, 47) = 1.67
            ppmmet(3, 47) = 1.67
            ppmmet(4, 47) = 1.50
C           ppmmet(1, 48) = 2.43 !MM20 @ 2.08 & 2.25 & 1.95; FWHM=.18
C           ppmmet(2, 48) = 1.81 !PPMEND often 1.8
            ppmmet(3, 48) = 2.43
            ppmmet(4, 48) = 1.81
         end if
 
         call toupper_lower (.true., seq)
         call toupper_lower (.true., seq_raw)
         if ( (seq .eq. 'PRESS'  .and. seq_raw .eq. 'STEAM')  .or.
     1        (seq .eq. 'STEAM'  .and. seq_raw .eq. 'PRESS') )
     2        call errmes (11, 3, chsubp)
         if (fwhmba_in_control .or. fwhmba.le.0.) fwhmba = fwhmba_sav
c        ----------------------------------------------------------------------
c        For Lorentzians, FWHM_absv/FWHM_real = sqrt(3) = 1.73.  Use 2.0 below
c        to give more freedom to lineshape to fit strange (apparently)
c        Gaussian windowed Friedman spectra when ABSVAL=T.
c        ----------------------------------------------------------------------
         if (.not.fwhmba_in_control  .and.  absval) fwhmba = 2. * fwhmba
      end if
c
      READ (LBASIS, NML=BASIS1, err=820, end=820)
C
      IF (LPRINT .GT. 0)  WRITE (LPRINT, 5210) IDBASI
 5210 FORMAT (/' Basis set ID = ',A//' No.   Metabolite',3X,
     1        'Expec[delta(1/T2)](1/s)', 3X, 'SD[delta(1/T2)]', 3X,
     2        'Shift(Pts)', 3X, 'SD[Shift](ppm)   ID')
c     -----------------------------------------------------------------------
c     ENCRYP = T to decrypt basis spectra
c     -----------------------------------------------------------------------
      encryp = badelt .lt. 0.
      badelt = abs(badelt)
      IF (AMIN1(DELTAT,BADELT).LE.0. .OR. MIN0(NDATA,NDATAB).LE.0  .or.
     1    ndatab .gt. mdata) CALL ERRMES (1, 4, CHSUBP)
c     ------------------------------------------------------------------------
c     Compare bandwidths.
c     IERROR_BW = 0 if the BWs are equal (and BASIS_BW does not need to be
c                   changed);
C                 1 if BASIS_BW != DATA_BW, which means that the BASIS_BW
c                                           must be converted to the DATA_BW;
c                 2 if BASIS_BW < DATA_BW and 90% of the BASIS_BW is needed
c                   for the Analysis Window, which means that there will very
c                   likely be distortions near the edge.
c                 4 if BASIS_BW < range needed by PPMEND.
c     ------------------------------------------------------------------------
      if (hzpppm_basis .gt. 0.) then
         ppminc_basis = 1. / (badelt * float(ndatab) * hzpppm_basis)
      else
         ppminc_basis = 1. / (badelt * float(ndatab) * hzpppm)
      end if
c     -------------------------------------------------------------------------
c     Set DESDSH if it has not been input positive or set positive in Block
c       Data.
c     -------------------------------------------------------------------------
      if (desdsh .le. 0.) then
         if (rincsh .le. 0.   .or. sdshmn .le. 0.   .or.
     1       sdshmx .lt. sdshmn) call errmes (26, 4, chsubp)
         desdsh = amin1(sdshmx, amax1(sdshmn, rincsh * ppminc_basis))
      end if
      ierror_bw = 0
      IF (ABS(1.-DELTAT/BADELT).GT.1.E-4 .and. .not.bascal) then
         if (BADELT .GT. DELTAT) then
            DATEN = (PPMCEN-PPMEND)/ppmINC_basis + 1.
            if (daten .gt. .45*float(ndatab)) then
               ierror_bw = 2
               if (nint(daten) .ge. ndatab / 2) ierror_bw = 4
               if (lstage .eq. 1   .and.  ierror_bw .ge. 2)
     1            call errmes (2, ierror_bw, chsubp)
            end if
         end if
      END IF
      corr_areaba = .false.
      DO 210 JMETAB=1,999
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C     TRAMP and VOLUME are the transmitter amplitude and VOI used for scaling
C       the spectra to get absolute concentrations, as discussed in MRM.  If
C       you do not have this data, then you can input TRAMP=VOLUME=1. in
C       NAMELIST BASIS, and the basis spectra will only be normalized to unit
C       concentration (with SCALE=1./CONC) below.  See the User's Manual for
C       the other input.
C
C        Don't forget the END=300 in your READ statement.
C
         READ (LBASIS, NML=BASIS, END=300, err=830)
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C
         if (metabo .eq. 'GSH'   .and.   gshgua   .and.
     1       nomit .lt. mpmet) then
            nomit = nomit + 1
            chomit(nomit) = 'Gua'
         end if
C        ----------------------------------------------------------------------
C        Omit this spectrum from the basis set, according to the
C          specifications in CHCALI, CHKEEP, PPMMET, CHUSE1, SYNUS1, or CHOMIT.
C        ----------------------------------------------------------------------
         if (ncalib .gt. 0) then
            do 211 j = 1, ncalib
               if (chcali(j) .eq. metabo) go to 216
 211        continue
            if (bascal  .and.  chbcal .eq. metabo) then
               nmetab = nmetab + 1
               go to 245
            end if
            go to 210
         end if
         DO 213 JOMIT=1,NOMIT
            IF (CHOMIT(JOMIT) .EQ. METABO) GO TO 210
  213    CONTINUE
         DO 214 J=1,NKEEP
            IF (CHKEEP(J) .EQ. METABO) GO TO 216
  214    CONTINUE
         DO 215 J=1,MPMET
            IF (CHPMET(J) .EQ. METABO) THEN
               IF ((PPMMET(1,J).GT.PPMST .OR. PPMMET(2,J).LT.PPMEND)
     1             .AND.
     2             (PPMMET(3,J).GT.PPMST .OR. PPMMET(4,J).LT.PPMEND))
     3               GO TO 210
               ppmmet_in_gap(1) = .false.
               ppmmet_in_gap(2) = .false.
               do 2152 jgap = 1, ngap
                  ppmmet_in_gap(1) = ppmmet_in_gap(1)   .or.
     1              (ppmmet(1, j) .le. ppmgap(1, jgap)   .and.
     1               ppmmet(1, j) .ge. ppmgap(2, jgap))   .or.
     1              (ppmmet(2, j) .le. ppmgap(1, jgap)   .and.
     1               ppmmet(2, j) .ge. ppmgap(2, jgap))
                  ppmmet_in_gap(2) = ppmmet_in_gap(2)   .or.
     1              (ppmmet(3, j) .le. ppmgap(1, jgap)   .and.
     1               ppmmet(3, j) .ge. ppmgap(2, jgap))   .or.
     1              (ppmmet(4, j) .le. ppmgap(1, jgap)   .and.
     1               ppmmet(4, j) .ge. ppmgap(2, jgap))
 2152          continue
               if (ppmmet_in_gap(1)   .and.   ppmmet_in_gap(2))
     1            go to 210
            END IF
  215    CONTINUE
 216     IF (LSTAGE .EQ. 1) THEN
            DO 217 JUSE1=1,NUSE1
               IF (US1FUL(JUSE1)) GO TO 217
               IF (CHUSE1(JUSE1) .EQ. METABO) THEN
                  US1FUL(JUSE1)=.TRUE.
                  GO TO 219
               END IF
               DO 218 JSYN=1,MPMET
                  IF (SYNUS1(1,JSYN).EQ.' ' .OR. SYNUS1(2,JSYN).EQ.' ')
     1                 GO TO 217
                  IF ((SYNUS1(1,JSYN).EQ.CHUSE1(JUSE1) .AND.
     1                 SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2                (SYNUS1(2,JSYN).EQ.CHUSE1(JUSE1) .AND.
     3                 SYNUS1(1,JSYN).EQ.METABO)) THEN
                     US1FUL(JUSE1)=.TRUE.
                     GO TO 219
                  END IF
  218          CONTINUE
  217       CONTINUE
            do 2182 j = 1, nnot1
               if (chnot1(j) .eq. metabo) go to 210
 2182       continue
c           -------------------------------------------------------------------
c           BADREF = T will include every metabolite into the Preliminary
c                      Analysis, unless it has been eliminated above.  The
c                      included metabolites will have CHNOLS; i.e., no
c                      lineshape convolution.
c           -------------------------------------------------------------------
            if (badref   .and.   nnolsh .lt. MMETAB_extra) then
c              ----------------------------------------------------------------
c              Must not let a SYNUS1 of another METABO in, since one could fit
c                the peak and the other create havoc.  This is mainly if
c                it is line-shape-enabled, but even if not, there is no sense
c                in ever having two SYNUS1 together in a Preliminary Analysis.
c              ----------------------------------------------------------------
               do 2184 KMETAB=1,NMETAB
                  DO 2185 JSYN=1,MPMET
                     IF (SYNUS1(1,JSYN).EQ.' ' .OR.
     1                   SYNUS1(2,JSYN).EQ.' ') GO TO 2184
                     IF ((SYNUS1(1,JSYN).EQ.NACOMB(kmetab) .AND.
     1                    SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2                   (SYNUS1(2,JSYN).EQ.NACOMB(kmetab) .AND.
     3                    SYNUS1(1,JSYN).EQ.METABO)) go to 210
 2185             CONTINUE
 2184          CONTINUE
               nnolsh = nnolsh + 1
               chnols(nnolsh) = metabo
               go to 219
            end if
            GO TO 210
         else
            do 2188 j = 1, nnot2
               if (chnot2(j) .eq. metabo) go to 210
 2188       continue
         END IF
  219    IF (NMETAB .GE. MMETAB) CALL ERRMES (5, 4, CHSUBP)
         NMETAB=NMETAB+1
         NACOMB(NMETAB)=METABO
         table_top(nmetab) = .true.
         DO 230 J=1,NSDSH
            IF (CHSDSH(J) .EQ. METABO) THEN
               SDSHIF(NMETAB)=ALSDSH(J)
               GO TO 232
            END IF
  230    CONTINUE
         SDSHIF(NMETAB)=DESDSH
  232    DO 234 J=1,NSDT2
            IF (CHSDT2(J) .EQ. METABO) THEN
               SDRT2(NMETAB)=ALSDT2(J)
               GO TO 236
            END IF
  234    CONTINUE
         SDRT2(NMETAB)=DESDT2
  236    DO 238 J=1,NEXT2
            IF (CHEXT2(J) .EQ. METABO) THEN
               EXRT2(NMETAB)=ALEXT2(J)
               GO TO 240
            END IF
  238    CONTINUE
         EXRT2(NMETAB)=DEEXT2
  240    EXRT2(NMETAB) = EXRT2(NMETAB) * RT2_SCALE
         SDRT2(NMETAB) = SDRT2(NMETAB) * RT2_SCALE
         IF (LPRINT .GT. 0) THEN
            IF (LSTAGE .EQ. 1) THEN
               WRITE (LPRINT,5234) NMETAB, METABO, ISHIFT, ID
 5234          FORMAT (1X, I3, 7X, A6, I57, 20X, A20)
            ELSE
               WRITE (LPRINT,5235) NMETAB, METABO, EXRT2(NMETAB),
     1                             SDRT2(NMETAB), ISHIFT,
     2                             SDSHIF(NMETAB), ID
 5235          FORMAT (1X, I3, 7X, A6, 1PE26.4, E18.4, I13, E17.2, 3X,
     1                 A20)
            END IF
         END IF
C        ----------------------------------------------------------------------
C        Convert SDSHIF from ppm to radians/s.
C        ----------------------------------------------------------------------
         SDSHIF(NMETAB)=SDSHIF(NMETAB)*2.*PI*HZPPPM
C        ----------------------------------------------------------------------
C        NCOMPO(JCONC) = no. of metabolites for concentration JCONC
C                      = 1, except for combinations.
C        LCOMPO(J,JCONC) = the subscript for the Jth metabolite in
C                            concentration JCONC, J=1,NCOMPO(JCONC),
C                            JCONC=1,NCONC.
C        ----------------------------------------------------------------------
         NCOMPO(NMETAB)=1
         LCOMPO(1,NMETAB)=NMETAB
C        ----------------------------------------------------------------------
C        BASISF(JDATA) = frequency-domain basis spectrum (COMPLEX).  BASISF
C                          is multiplied by SCALE=TRAMP/(VOLUME*CONC), but not
C                          rearranged and not shifted.
C                        By also dividing by CONC, the solution is already in
C                          concentration units (e.g. mM).
C        ISHIFT = referencing shift for frequency-domain spectrum NMETAB; a
C                 positive ISHIFT shifts the spectrum ISHIFT grid points to
C                 the left, i.e., toward larger ppm.
C        BASIST(JDATA,NMETAB) = inverse FFT of BASISF (COMPLEX).
C        ----------------------------------------------------------------------
C
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv For user vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C
C        Read full basis array into BASISF.
C
 245     READ (LBASIS, FMTBAS, err=831, end=831) (BASISF(J),J=1,NDATAb)
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For user ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
         if (encryp) then
            dix = 1499.d0
            do 247 j = 1, ndatab
               basisf(j) = -basisf(j) * exp(-20. * random(dix) + 10.)
 247        continue
         end if
C
         IF (AMIN1(TRAMP,VOLUME,CONC) .LE. 0.) CALL ERRMES (6, 4,
     1                                                      CHSUBP)
         SCALE=TRAMP/(VOLUME*CONC)
C        ----------------------------------------------------------------------
C        Temporarily put the scaled, but unshifted, frequency-domain data into
C          BASIST(*,NMETAB).
C        ----------------------------------------------------------------------
         DO 250 J=1,NDATAb
            BASIST(J,NMETAB)=BASISF(J)*SCALE
  250    CONTINUE
C        ----------------------------------------------------------------------
c        Take absolute value of basis spectra.  This must be done here, and
c          not below, because full zero-filling of the time data is needed to
c          preserve full information in the absolute value spectrum.
c        The factor of 2. * CABS is to compensate for discarding the
c          symmetric/antisymmetric 2nd half in points NUNFIL+1--NDATA.
c        1st time-data point (the integral over the spectrum) must be halved,
c          to remove the offset produced by cutting off the 2nd half of
c          time data.  See Notes of 010819.
C        ----------------------------------------------------------------------
         if (absval) then
            do 256 j = 1, ndatab
               basisf(j) = 2. * cabs(BASIST(J,NMETAB))
 256        continue
            CALL CFFTIN (BASISF, BASIST(1,NMETAB), NDATAb,
     1           LWFFT, WFFTC)
            BASIST(1,NMETAB) = .5 * BASIST(1,NMETAB)
            do 258 j = (1 + ndatab/2), ndatab
               BASIST(j,NMETAB) = 0.
 258        continue
            call cfft (basist(1,nmetab), basist(1,nmetab), ndatab,
     1           lwfft, wfftc)
         end if
C        ----------------------------------------------------------------------
C        Put shifted frequency-domain data into BASISF.
C        ----------------------------------------------------------------------
         jshift = ishift + nint((4.65 - ppmcen) / ppminc_basis)
         DO 260 J=1,NDATAb
            jshift = jshift + 1
            BASISF(J)=BASIST(ICYCLE(JSHIFT,NDATAb), NMETAB)
  260    CONTINUE
 
c        ----------------------------------------------------------------------
c        AREABA_ORIG_BASISF = T in the normal case; else, AREABA with
c           transformed BASISF can be inaccurate, e.g., with basis spectra
c           that are truncated because of large AQ and have a peak that is
c           difficult to integrate.  (See Notes of 120929.)
c        Much code is repeated from below, to allow (inaccurate) case to be
c           tested below when AREABA_ORIG_BASISF=F.
c        ----------------------------------------------------------------------
         if (.not.areaba_orig_basisf) go to 269
         dows_now = dows   .and.   lstage .eq. 2   .and.   havh2o  .and.
     1              .not.bascal   .and.
     2              fcalib .gt. .999999   .and.   fcalib .lt. 1.000001
         if (.not.(dows_now .or. (scasim .and. nsimul.gt.0))) go to 269
         if (wsmet .ne. metabo) then
            DO 262 JSYN=1,MPMET
               IF (SYNUS1(1,JSYN).EQ.' ' .OR. SYNUS1(2,JSYN).EQ.' ')
     1              GO TO 269
               IF ((SYNUS1(1,JSYN).EQ.wsmet .AND.
     1              SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2              (SYNUS1(2,JSYN).EQ.wsmet .AND.
     3              SYNUS1(1,JSYN).EQ.METABO)) go to 264
 262        CONTINUE
         end if
 264     area_met_norm = areaba(basisf, ppminc_basis, ndatab/2)
         corr_areaba = .true.
 
c        ----------------------------------------------------------------------
c        Check (and correct for) inconsistent BWs or field strengths.
c        ----------------------------------------------------------------------
 269     if (bascal) then
            ndata_freq = ndatab
         else
c           ------------------------------------------------------------------
c           Choose NDATA_FREQ so that bw_basis_new = bw_data.
c           If HZPPPM_BASIS is available, then the bandwidth will be adjusted
c             to compensate for the difference between HZPPPM & HZPPPM_BASIS.
c           The following reasoning is used:
c           The basis spectra must be adjusted so that 2 conditions are
c             matched:
c             (1) Basis bandwidth must be changed to data bandwidth.
c                 This is done by truncating or zero-filling the basis
c                   spectra.  The spacing of the basis spectra remains fixed
c                   in this step at 1/(ndatab*badelt*hzpppm_basis).
c                 RNDATA_FREQ = the number of points in DATAT such that the
c                               total time range of DATAT equals that of
c                               BASIST; i.e., such that DATAF has the same
c                               spacing as BASISF
c                             = NDATAB*BADELT*HZPPPM_BASIS/(DELTAT*HZPPPM)
c                 Zero-filling or truncating BASISF to NDATA_FREQ points and
c                   doing an FFTINV produces a BASIST with the same spacing as
c                   DATAT.  This part is independent of NDATA.
c             (2) Basis frequency spacing must be changed to data spacing.
c                   BASIST is zero-filled or truncated to NDATA points, so
c                   that BASISF matches the spacing of DATAF (as calculated in
c                   Step (1)).
c                 This could cause significant errors due to inconsistencies
c                   in the scaling of the basis spectra caused by (differing
c                   effects of truncation).  There should be a warning to be
c                   consistent with NUNFIL*DELTAT with all in vivo data.  In
c                   practice, could partially reduce this effect with FWHMSM
c                   in MakeBasis, especially if AUTOSC=T, but it is nice to
c                   have sharp basis spectra, with plenty of possibility of
c                   lineshape convolution to correct for lineshape.
c           ------------------------------------------------------------------
            if (hzpppm_basis .gt. 0.) then
               rndata_freq = float(ndatab) * badelt * hzpppm_basis /
     1                       (deltat * hzpppm)
               if (lstage .eq. 1   .and.   nmetab .eq. 1) then
c                 -------------------------------------------------------------
c                 Test consistency of field strengths
c                 Variations of GE field strengths among sites seems to be
c                   only  about 1 part in 3000, corresponding to TEST = .00033
c                 However, Danielsen's 1.5 T Siemens was about 0.995 * 1.5T.
c                 Furthermore, the BW correction means that the only
c                   differences will be in the coupling.  So, the threshold
c                   for warning was raised to .05 (from .005).
c                 -------------------------------------------------------------
                  test=abs(hzpppm_basis/hzpppm - 1.)
                  if (test .ge. .2) then
                     call errmes (9, 4, chsubp)
                  else if (test .ge. .05) then
                     call errmes (9, 2, chsubp)
                   end if
               end if
            else
               rndata_freq = float(ndatab) * badelt / deltat
            end if
c           -------------------------------------------------------------------
c           NDATA_FREQ must be even
c           -------------------------------------------------------------------
            ndata_freq = 2 * nint(.5 * rndata_freq)
c           -------------------------------------------------------------------
c           Leave NDATA_FREQ=NDATA if their ratio is within BWTOLR.
c           -------------------------------------------------------------------
            if (abs(1. - rndata_freq / float(ndata)) .le. bwtolr)
     1         ndata_freq = ndata
            if (ndata_freq .gt. mdata) call errmes (35, 4, chsubp)
c           -------------------------------------------------------------------
c           Correct for different lengths of FFTs.  The FFTs must be
c              multiplied by sqrt(N), where N is the number of points.  The
c              factor is sqrt(ndata_freq/ndata).  (See also Notes of 120929.)
c           -------------------------------------------------------------------
            if (corr_areaba) then
               area_met_norm = area_met_norm *
     1            sqrt(float(ndata_freq) / float(ndata))
               corr_areaba = .false.
            end if
            if (ndata_freq .ne. ndatab) then
c              ----------------------------------------------------------------
c              Change bw_basis to bw_data, and put result temporarily into
c                BASIST.  This part is independent of NDATA.  It is Step (1)
c                described after 260 above.
c              ----------------------------------------------------------------
               if (ndata_freq .gt. ndatab) then
                  do 270 j = 1, ndata_freq
                     basist(j, nmetab) = 0.
 270              continue
               end if
c              ----------------------------------------------------------------
c              1st into 1st half
c              ----------------------------------------------------------------
               nhalf = min0(ndata_freq, ndatab) / 2
               do 272 j = 1, nhalf
                  basist(j, nmetab) = basisf(j)
 272           continue
c              ----------------------------------------------------------------
c              Now into 2nd half
c              ----------------------------------------------------------------
               jin = ndatab
               do 274 jout = ndata_freq, ndata_freq - nhalf + 1, -1
                  basist(jout, nmetab) = basisf(jin)
                  jin = jin - 1
 274           continue
c              ----------------------------------------------------------------
c              Now everything back into BASISF
c              ----------------------------------------------------------------
               do 275 j = 1, ndata_freq
                  basisf(j) = basist(j, nmetab)
 275           continue
            end if
         end if
C        ----------------------------------------------------------------------
C        Do an inverse FFT of BASISF to put the time-domain data
C          into BASIST(*,NMETAB).
C        ----------------------------------------------------------------------
         if (bascal  .and.  chbcal .eq. metabo) then
            CALL CFFTIN (BASISF, datat, NDATA_freq,
     1                LWFFT, WFFTC)
            nmetab = nmetab - 1
            go to 210
c            go to 9999
         else
            CALL CFFTIN (BASISF, BASIST(1,NMETAB), NDATA_freq,
     1           LWFFT, WFFTC)
         end if
C        ----------------------------------------------------------------------
c        The following was used from 5.2-3B to 6.0-2C based on the incorrect
c          reasoning.  See Notes 040118.
C        ----------------------------------------------------------------------
c         basist(1, nmetab) = basist(1, nmetab) * float(ndata) /
c     1                                           float(ndata_freq)
C        ----------------------------------------------------------------------
c        Zero-fill BASIST.  This is (the final) Step (2) described after 260
c          above.
C        ----------------------------------------------------------------------
         if (ndata_freq .lt. ndata) then
            do 280 j = ndata_freq + 1, ndata
               basist(j, nmetab) = 0.
 280        continue
         end if
C        ----------------------------------------------------------------------
c        BASOUT = T to write out modified BASIS file.
c        ABSVAL = T is not permitted when BASOUT=T
C        ----------------------------------------------------------------------
         if (basout   .and.   absval) call errmes (14, 3, chsubp)
         if (basout   .and.   ndata_freq .eq. ndata   .and.
     1       ndata_freq .eq. ndatab) call errmes (15, 1, chsubp)
         if (basout   .and.   (ndata_freq .ne. ndata   .or.
     1                         ndata_freq .ne. ndatab)   .and.
     2       lstage .eq. 2   .and.   .not.absval) then
            if (nmetab .eq. 1) then
               lfilbas = ilen(filbas)
               lastdot = 0
               lastslash = 0
               do 282 j = 1, lfilbas
                  if (filbas(j:j) .eq. '.') lastdot = j
                  if (filbas(j:j) .eq. '/') lastslash = j
 282           continue
               if (lastdot .gt. lastslash) then
                  file_basout = filbas(1:lastdot-1) // '-new-bw' //
     1                          filbas(lastdot:lfilbas)
               else
                  file_basout = filbas(1:lfilbas) // '-new-bw'
               end if
               OPEN (21, FILE=file_basout, STATUS='UNKNOWN', err=832)
               write (21, nml=seqpar)
               j = min0(ilen(idbasi),mchidb-3)
               idbasi(j+2:j+3) = 'BW'
               badelt_save = badelt
               badelt = deltat
               ndatab_save = ndatab
               ndatab = ndata
               write (21, nml=basis1)
               badelt = badelt_save
               ndatab = ndatab_save
            end if
            conc=1.
            tramp=1.
            volume=1.
            ishift = nint((ppmcen - 4.65) / ppminc)
            write (21, nml=basis)
            call cfft (basist(1,nmetab), basisf, ndata, lwfft, wfftc)
            write (21, FMTBAS) (BASISF(J),J=1,NDATA)
         end if
c        --------------------------------------------------------------------
c        Suppress convolution with Lineshape (and broaden with a Gaussian
c          with FWHM=FWHMST) according to CHNOLS, when LSTAGE=1.
c        This is mainly for Lac & Ala when simulated lipids are being used.
c        Broadening from CHANGE-RAW.f is used below
c        --------------------------------------------------------------------
         if (lstage .eq. 1   .and.   nnolsh .gt. 0) then
            do 285 j = 1, nnolsh
               if (metabo .eq. chnols(j)) then
                  call set_lshape_false ()
                  go to 289
               end if
 285        continue
         end if
 289     continue
c        --------------------------------------------------------------------
c        Try scaling, if needed.
c        (The 2 statements containing 9999 can be uncommented to test the
c          SCASIM=T scaling using BASCAL with CHCALI(1)='SimCr' & CHBCAL='Cr'.)
c        --------------------------------------------------------------------
c9999     continue
         dows_now = dows   .and.   lstage .eq. 2   .and.   havh2o  .and.
     1              .not.bascal   .and.
     2              fcalib .gt. .999999   .and.   fcalib .lt. 1.000001
         if (.not.(dows_now .or. (scasim .and. nsimul.gt.0))) go to 210
         if (wsmet .ne. metabo) then
            DO 290 JSYN=1,MPMET
               IF (SYNUS1(1,JSYN).EQ.' ' .OR. SYNUS1(2,JSYN).EQ.' ')
     1              GO TO 210
               IF ((SYNUS1(1,JSYN).EQ.wsmet .AND.
     1              SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2              (SYNUS1(2,JSYN).EQ.wsmet .AND.
     3              SYNUS1(1,JSYN).EQ.METABO)) go to 295
 290        CONTINUE
         end if
c        --------------------------------------------------------------------
c        Get BASISF for the new BW.
c        Scaling by sqrt(ndata_freq/ndata) leads to absurd results in extreme
c           cases.
c        --------------------------------------------------------------------
 295     if (ndata_freq .ne. ndata)
     1      call cfft (basist(1, nmetab), basisf, ndata, lwfft, wfftc)
c        --------------------------------------------------------------------
c        This call to AREABA is superceded.  It will normally be skipped,
c           because AREA_MET_NORM will have been calculated above in the
c           normal case that the original basis spectrum instead of the
c           transformed spectrum is integrated above, because
c           AREABA_ORIG_BASISF = T.
c        --------------------------------------------------------------------
         if (area_met_norm .le. 0.)
     1          area_met_norm = areaba(basisf, ppminc, nunfil)
C        ----------------------------------------------------------------------
c        Try water scaling
C        ----------------------------------------------------------------------
         if (dows_now   .and.   area_met_norm .gt. 0.)
     1      call water_scale ()
  210 CONTINUE
c     -------------------------------------------------------------------------
c     NSIMUL <= MMETAB = # extra model spectra to be synthesized
c     MGAU = 20 = max # of Gaussian components in a synthesized model spectrum
c     CHSIM(JSIMUL)*6 = name
c     SIFWMN(JGAU,JSIMUL) = FWHM of Gaussian JGAU
c                         < 0 for SIFWMN=FWHMBA and EXRT2 SDRT2 SDSH = DE*
c     SIFWEX(JSIMUL) = Expected FWHM of Gaussian JGAU=1  This determines EXRT2
c                      (of course, for the whole simulated spectrum, not just
C                       JGAU=1)
c     SIFWSD(JSIMUL) = SD of FWHM of Gaussian JGAU=1.  This determines SDRT2
c                      (of course, for the whole simulated spectrum, not just
C                       JGAU=1)
c     NGAU(JSIMUL) = # of Gaussian components in a synthesized model spectrum
c                    JSIMUL
c     SISDSH(JSIMUL) = SD of shift (ppm) (as ALSDSH & DESDSH)
c     SIAMP(JGAU, JSIMUL) = amplitude
c     SIPPM(JGAU, JSIMUL) = chemical shift, except for IMETHD=2:
 
c     IMETHD = 2 to replace lineshape regularizor with a linewidth
c                regularizor  minimizing
c                c_j/c^0_j * (gamma^0_j + RLRNTZ*gamma_j)**IPOWRG, where
c                c_j = CONC parameter (as usual)
c                gamma_j = 1/RT2 (as usual)
c                gamma^0_j = pi * HZPPPM * SIFWMN(1, JSIMUL)
c                Single-peaked simulations are assumed.
c                c^0_jsimul = SIAMP(2, JSIMUL).  This is signaled by:
c                SIPPM(2, JSIMUL) >= 999, which causes NGAU=1 to be used.
c                So, 2nd component is: @999. FWHM=1. AMP=c^0_jsimul
 
c     IMETHD = 3 to broaden by multiplying DATAT by factors of form
c                exp(-PARNLN(J) * T ** POWER(K)), where model NMETAB has
c                POWER(1), ... POWER(NPOWER(NMETAB)) powers.
c
 
c     SIDUMP(JSIMUL) = T to abort and produce a plot of the real part of
c                        simulated model spectrum JSIMUL.
c
c     -------------------------------------------------------------------------
  300 if (nsimul .gt. 0) then
         if (nsimul .gt. mmetab - nmetab) call errmes (16, 4, chsubp)
         if (lstage .eq. 1) call parse_chsimu ()
         if (gauss_rt2   .and.   lprint .gt. 0) write(lprint, 5305)
 5305    format(' No.   Metabolite', 9x, 'EXRT2', 6x,
     1          'SIFWEX  input', 3x, 'SIFWMN', 14x, 'SDRT2', 6x,
     1          'SIFWSD   input', 3x, 'SD[Shift](ppm)')
         do 310 jsimul = 1, nsimul
            metabo = chsim(jsimul)
            if (ncalib .gt. 0) then
               do 312 j = 1,ncalib
                  if (chcali(j) .eq. metabo) go to 316
 312           continue
               go to 310
            end if
c           -------------------------------------------------------------------
c           Omit spectrum according to CHOMIT.
c           -------------------------------------------------------------------
            DO 303 JOMIT=1,NOMIT
               IF (CHOMIT(JOMIT) .EQ. METABO) GO TO 310
 303        CONTINUE
            DO 314 J=1,NKEEP
               IF (CHKEEP(J) .EQ. METABO) GO TO 316
 314        CONTINUE
            DO 315 J=1,MPMET
               IF (CHPMET(J) .EQ. METABO) THEN
                  IF ((PPMMET(1,J).GT.PPMST .OR. PPMMET(2,J).LT.PPMEND)
     1                 .AND.
     2                (PPMMET(3,J).GT.PPMST .OR. PPMMET(4,J).LT.PPMEND))
     3                 GO TO 310
               END IF
 315        CONTINUE
c           -------------------------------------------------------------------
c           Omit spectrum according to SIPPM and [PPMST,PPMEND]
c           Could set
c             DISTMX = .5 * SIFWEX(JSIMUL)
c             since it could intrude.  However, user should be sensible enough
c             to include the whole peak.
c           -------------------------------------------------------------------
            if (chksim) then
               distmx = 0.
               do 305 jgau = 1, ngau(jsimul)
                  do 307 jgap = 1, ngap
                     if (sippm(jgau,jsimul) .le. ppmgap(1, jgap)   .and.
     1                   sippm(jgau,jsimul) .ge. ppmgap(2, jgap))
     2                  go to 310
 307               continue
                  if (sippm(jgau,jsimul) - distmx .lt. ppmst   .and.
     1                sippm(jgau,jsimul) + distmx .gt. ppmend) go to 316
 305           continue
               go to 310
            end if
 316        if (lstage .eq. 1) then
c              ----------------------------------------------------------------
c              Omit spectrum according to CHNOT1 & CHUSE1
c              ----------------------------------------------------------------
               do 3182 j = 1, nnot1
                  if (chnot1(j) .eq. metabo) go to 310
 3182          continue
               DO 317 JUSE1=1,NUSE1
                  IF (US1FUL(JUSE1)) GO TO 317
                  IF (CHUSE1(JUSE1) .EQ. METABO) THEN
                     US1FUL(JUSE1)=.TRUE.
                     GO TO 329
                  END IF
                  DO 318 JSYN=1,MPMET
                     IF (SYNUS1(1,JSYN).EQ.' ' .OR.
     1                   SYNUS1(2,JSYN).EQ.' ') GO TO 317
                     IF ((SYNUS1(1,JSYN).EQ.CHUSE1(JUSE1) .AND.
     1                    SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2                    (SYNUS1(2,JSYN).EQ.CHUSE1(JUSE1) .AND.
     3                    SYNUS1(1,JSYN).EQ.METABO)) THEN
                        US1FUL(JUSE1)=.TRUE.
                        GO TO 329
                     END IF
 318              CONTINUE
 317           CONTINUE
c              ----------------------------------------------------------------
c              BADREF = T will include every metabolite into the Preliminary
c                         Analysis, unless it has been eliminated above.  The
c                         included metabolites will have CHNOLS; i.e., no
c                         lineshape convolution.
c              Arbitrarily set SIFWSD; it will not be used (but will
c                still be checked below when LSTAGE=2).
c              ----------------------------------------------------------------
               if (badref) then
c                 -------------------------------------------------------------
c                      Must not let a SYNUS1 of another METABO in, since one
c                 could fit the peak and the other create havoc.  This is
c                 mainly if it is line-shape-enabled, but even if not,
c                 there is no sense in ever having two SYNUS1 together in a
c                 Preliminary Analysis.
c                 -------------------------------------------------------------
                  do 3184 KMETAB=1,NMETAB
                     DO 3185 JSYN=1,MPMET
                        IF (SYNUS1(1,JSYN).EQ.' ' .OR.
     1                       SYNUS1(2,JSYN).EQ.' ') GO TO 3184
                        IF ((SYNUS1(1,JSYN).EQ.NACOMB(kmetab) .AND.
     1                       SYNUS1(2,JSYN).EQ.METABO)  .OR.
     2                      (SYNUS1(2,JSYN).EQ.NACOMB(kmetab) .AND.
     3                       SYNUS1(1,JSYN).EQ.METABO)) go to 310
 3185                CONTINUE
 3184             CONTINUE
                  do 322 j = 1, nnot2
                     if (chnot2(j) .eq. metabo)
     1                  sifwsd(jsimul) = .5 * sifwex(jsimul)
 322              continue
                  go to 329
               end if
               GO TO 310
            else
c              ----------------------------------------------------------------
c              LSTAGE=2.  Omit spectrum according to CHNOT2
c              ----------------------------------------------------------------
               do 325 j = 1, nnot2
                  if (chnot2(j) .eq. metabo) go to 310
 325           continue
               if (scafwh) then
c                 ------------------------------------------------------------
c                 Scale SIFW* by FWHMST from LSTAGE=1
c                 ------------------------------------------------------------
                  term = amax1(fwhmst, fwhmmn)
                  sifwmn(1, jsimul) = sifwmn(1, jsimul) * term
                  sifwex(jsimul) = sifwex(jsimul) * term
                  sifwsd(jsimul) = sifwsd(jsimul) * term
               end if
            END IF
 329        nmetab = nmetab + 1
            if (nmetab .gt. mmetab) call errmes (25, 4, chsubp)
            NACOMB(NMETAB)=METABO
c           ----------------------------------------------------------------
c           Use Lineshape only if simulating basis spectrum or if specified
c             by CHLSHA & NLSHAP (or if NOBASI=T ???)
c           ----------------------------------------------------------------
            lshape(nmetab) = sifwmn(1, jsimul) .le. 0.
            table_top(nmetab) = lshape(nmetab)
            if (.not.lshape(nmetab)) then
               do 330 j = 1, nlshap
                  if (chlsha(j) .eq. metabo) lshape(nmetab) = .true.
 330           continue
            end if
c           -------------------------------------------------------------------
c           Compute SDRT2 & EXRT2 from the 3 SIFW*, assuming additivity of
c             variances.
c           *_EXTRA are the ranges needed; SDRT2 should be about HZ_EXTRA/3.
c           EXRT2_MIN, SDRT2_MIN & SDSHIF_MIN are defined here rather than
c             input.  They are to avoid ill-conditioning due to too strong
c             priors.
c           -------------------------------------------------------------------
            if (sifwmn(1, jsimul) .le. 0.) then
c              ----------------------------------------------------------------
c              Special case for simulating basis spectra, using FWHMBA and
c                DE* for broadening & shift priors.
c              In this case, only CHSIM, SIPPM & SIAMP are used from CHSIMU;
c                SISDSH, SIFWEX & SIFWSD are not used (and must have SIFWMN<0).
c              ----------------------------------------------------------------
               exrt2(nmetab) = deext2 * RT2_SCALE
               sdrt2(nmetab) = desdt2 * RT2_SCALE
               sdshif(nmetab) = desdsh
            else
               if (imethd .eq. 3) then
c                 -------------------------------------------------------------
c                 Special use of SIFWSD = NPOWER
c                 -------------------------------------------------------------
                  npower(nmetab) = nint(sifwsd(jsimul))
                  if (.not.nobasi   .or.   npower(nmetab) .lt. 1
     1                .or.   npower(nmetab) .gt. mmpowr)
     2                 call errmes (39, 4, chsubp)
               else
c                 -------------------------------------------------------------
c                 SIFWMN (and all other SI*) input (not simulating with
c                    FWHMBA, etc)
c                 -------------------------------------------------------------
                  if (sifwmn(1,jsimul) .ge. sifwex(jsimul)  .or.
     1               sifwsd(jsimul) .le. 0.) call errmes (18, 4, chsubp)
c                 -------------------------------------------------------------
c                 GAUSS_RT2 = T only for old computation of EXRT2 & SDRT2
c                              based (incorrectly) on Gaussian broadening,
c                              instead of Lorentzian broadening.  This will be
c                              the default, for backward compatibility.
c                 -------------------------------------------------------------
                  if (gauss_rt2) then
                     fwhm_ex = sqrt(sifwex(jsimul)**2 -
     1                              sifwmn(1,jsimul)**2)
                  else
                     fwhm_ex = sifwex(jsimul) - sifwmn(1,jsimul)
                  end if
                  exrt2(nmetab) = pi * hzpppm * fwhm_ex
                  EXRT2_MIN = .5
                  if (exrt2(nmetab) .lt. EXRT2_MIN) then
                     call errmes (19, 2, chsubp)
                     exrt2(nmetab) = EXRT2_MIN
                  end if
                  if (gauss_rt2) then
c                 -------------------------------------------------------------
c                 GAUSS_RT2 led to excessive SIFWSD.
c                 The following formula is the solution of the quadratic
c                   equation requiring that (EXRT2 + SDRT2) produce a
c                   broadened peak with FWHM = SIFWEX + SIFWSD; i.e.,
c                 sqrt{SIFWMN**2 + [(EXRT2+SDRT2)/pi*HZPPPM]**2} = SIFWEX +
c                                                                  SIFWSD
c                 It was also checked numerically
c                 -------------------------------------------------------------
                     sdrt2(nmetab) = pi * hzpppm *
     1                               sqrt((sifwex(jsimul) +
     2                                     sifwsd(jsimul))**2 -
     3                              sifwmn(1,jsimul)**2) - exrt2(nmetab)
                  else
                     sdrt2(nmetab) = pi * hzpppm * sifwsd(jsimul)
                  end if
                  SDRT2_MIN = .25
                  if (sdrt2(nmetab) .lt. SDRT2_MIN) then
                     call errmes (20, 2, chsubp)
                     sdrt2(nmetab) = SDRT2_MIN
                  end if
               end if
c              ----------------------------------------------------------------
c              SISDSH is only necessary when LSTAGE=2 is used.  With BADREF=T,
c                it may be that it is only used when LSTAGE=1.
c              ----------------------------------------------------------------
               if (sisdsh(jsimul) .le. 0.  .and.  lstage .eq. 2)
     1              call errmes (21, 4, chsubp)
               sdshif(nmetab) = sisdsh(jsimul)
c              ----------------------------------------------------------------
c              Might have to change SDSHIF_MIN for other nuclei.
c              ----------------------------------------------------------------
               SDSHIF_MIN = ppminc / 10.
               if (sdshif(nmetab) .lt. SDSHIF_MIN   .and.
     1             lstage .eq. 2) then
                  if (sptype(:7) .ne. 'muscle-')
     1               call errmes (21, 2, chsubp)
                  sdshif(nmetab) = SDSHIF_MIN
               end if
            end if
 
            IF (LPRINT .GT. 0) THEN
               IF (LSTAGE .EQ. 1) THEN
                  WRITE (LPRINT,5334) NMETAB, METABO
 5334             FORMAT (1X, I3, 7X, A6)
               ELSE
                  if (imethd .eq. 3) then
 5337                format(1X, I3, 7X, A6, i4, f9.3)
                     write (lprint, 5337) NMETAB, METABO,
     1                                    npower(nmetab), sdshif(nmetab)
                  else
                     if (gauss_rt2) then
                        write(lprint, 5336) NMETAB, METABO,
     1                     EXRT2(NMETAB),
     2                     exrt2(nmetab)/(pi*hzpppm) + sifwmn(1,jsimul),
     3                     sifwex(jsimul), sifwmn(1, jsimul),
     4                     sdrt2(nmetab), sdrt2(nmetab)/(pi*hzpppm),
     5                     sifwsd(jsimul), sdshif(nmetab)
 5336                   format(1X, I3, 7X, A6, f14.3, f12.3, f7.3, f9.3,
     1                         f19.3, f12.4, f8.4, f17.3)
                     else
                        WRITE (LPRINT,5335) NMETAB, METABO,
     1                  EXRT2(NMETAB), SDRT2(NMETAB), SDSHIF(NMETAB)
 5335                   FORMAT (1X, I3, 7X, A6, 1PE26.4, E18.4, 13x,
     1                          E17.2)
                     end if
                  end if
               END IF
            END IF
C           -------------------------------------------------------------------
C           Convert SDSHIF from ppm to radians/s.
C           -------------------------------------------------------------------
            SDSHIF(NMETAB)=SDSHIF(NMETAB)*2.*PI*HZPPPM
C           -------------------------------------------------------------------
C           NCOMPO(JCONC) = no. of metabolites for concentration JCONC
C                         = 1, except for combinations.
C           LCOMPO(J,JCONC) = the subscript for the Jth metabolite in
C                             concentration JCONC, J=1,NCOMPO(JCONC),
C                             JCONC=1,NCONC.
C           -------------------------------------------------------------------
            NCOMPO(NMETAB)=1
            LCOMPO(1,NMETAB)=NMETAB
c           -------------------------------------------------------------------
c           Simulate (real) spectrum, using strategy from
c             real-spectra-only/make-raw.f
c           SUM = sum over spectrum, for direct comparison with PPMINC*sum in
c                 AREA_MET_NORM
c           -------------------------------------------------------------------
            if (ngau(jsimul) .le. 0   .or.
     1          ngau(jsimul) .gt. mgau) call errmes (17, 4, chsubp)
 
            do 340 jdata = 1, ndata
               basisf(jdata) = (0.,0.)
 340        continue
            expmax = 2. * alog(rrange)
c            sum = 0.
            kpower = 0
            if (imethd .eq. 2   .and.    ngau(jsimul) .ne. 2)
     1         call errmes (40, 4, chsubp)
            do 350 jgau = 1, ngau(jsimul)
               if (sifwmn(1, jsimul) .le. 0.) then
c                 -------------------------------------------------------------
c                 Special case for simulating basis spectra.
c                 -------------------------------------------------------------
                  if (nobasi   .or.   imethd .eq. 2   .or.
     1               imethd .eq. 3) call errmes (37, 4, chsubp)
                  rsd = 2. * sqrt(2. * alog(2.)) / fwhmba
               else
                  if (sifwmn(jgau, jsimul) .le. 0.)
     1                  call errmes (24, 4, chsubp)
                  rsd = 2. * sqrt(2. * alog(2.)) / sifwmn(jgau, jsimul)
               end if
               if (imethd .eq. 2   .and.  jgau .eq. 2) then
                  if (sippm(2, jsimul) .le. 998.   .or.
     1                siamp(2, jsimul) .le. 0.   .or.   .not.nobasi)
     2               call errmes (36, 4, chsubp)
                  if (lstage .eq. 1) then
                     conc_expect(nmetab) = siamp(2, jsimul)
                  else
                     if (fconc_expect .le. 0.)
     1                   call errmes (38, 4, chsubp)
                     conc_expect(nmetab) = siamp(2, jsimul) *
     1                                     fconc_expect
                  end if
                  rt2min(nmetab) = pi * hzpppm * sifwmn(1, jsimul)
                  go to 350
               else if (imethd .eq. 3) then
                  if (sippm(jgau, jsimul) .gt. 998.) then
                     kpower = kpower + 1
                     fract_power_sd(kpower, nmetab) =
     1                  sifwmn(jgau, jsimul)
                  end if
               else
                  if (abs(siamp(jgau,jsimul)) .le. 0.   .or.
     1                abs(sippm(jgau,jsimul)) .gt. 998.)
     2               call errmes (22, 4, chsubp)
               end if
               if (sippm(jgau, jsimul) .le. 998.) then
                  anorm = siamp(jgau,jsimul) * rsd / sqrt(2. * pi)
                  xppm = ppmcen + nunfil * ppminc
                  do 360 jdata = 1, ndata
                     expon = (rsd * (xppm - sippm(jgau,jsimul)))**2
                     if (expon .lt. expmax) then
                        term = anorm * exp(-.5 * expon)
                        basisf(jdata) = basisf(jdata) + term
c                        sum = sum + term
                     end if
                     xppm = xppm - ppminc
 360              continue
               end if
 350        continue
c           -------------------------------------------------------------------
c           Arrange (as BASISF normally is, i.e., "not rearranged").
c           -------------------------------------------------------------------
            do 370 j = 1, nunfil
               cterm(1) = basisf(nunfil + j)
               basisf(nunfil + j) = basisf(j)
               basisf(j) = cterm(1)
 370        continue
            if (absval) then
c              ----------------------------------------------------------------
c              Take absolute value of basis spectra and compute BASIST.
c              ----------------------------------------------------------------
               do 372 j = 1, ndata
                  basisf(j) = cabs(BASISf(J))
 372           continue
               CALL CFFTIN (BASISF, BASIST(1,NMETAB), NDATA,
     1                      LWFFT, WFFTC)
               BASIST(1,NMETAB) = .5 * BASIST(1,NMETAB)
            else
               call  CFFTin (basisf, basist(1, nmetab), NDATA, LWFFT,
     1                       WFFTC)
            end if
c           -------------------------------------------------------------------
c           Replace 2nd half with zeroes.
c           See comments after 250 explaining that the factor of 2 below
c           is needed, apparently because the 2nd half of the
c           symmetric/antisymmetric part is zeroed.
c                This factor of 2 also gave correct calibration between a
c           simulated Cr and a basis Cr, both with ABSVAL=T & F.
c           -------------------------------------------------------------------
            do 373 junfil = 1, nunfil
               basist(junfil, nmetab) = 2. * basist(junfil, nmetab)
               basist(nunfil + junfil, nmetab) = (0.,0.)
 373        continue
c           -------------------------------------------------------------------
c                Correct possible remaining offset, apparently due to the
c           finiteness of the frequency range.  It is bigger with broad
c           simulated peaks.
c           -------------------------------------------------------------------
            call cfft_r (basist(1, nmetab), basisf, ndata, lwfft, wfftc)
            offset_start = (0., 0.)
            offset_end = (0., 0.)
            do 375 j = 1, 5
               offset_start = offset_start + basisf(j)
               offset_end = offset_end + basisf(ndata - j + 1)
 375        continue
            if (abs(real(offset_start)) .lt. abs(real(offset_end))) then
               basist(1, nmetab) = basist(1, nmetab) - .2 *
     1                             sqrt(float(ndata)) * offset_start
            else
               basist(1, nmetab) = basist(1, nmetab) - .2 *
     1                             sqrt(float(ndata)) * offset_end
            end if
c           -------------------------------------------------------------------
c           For more realistic LSTAGE=1 (and initial FWHM and lineshape),
c             exponentially broaden model spectrum according to EXRT2
c           SCAFWH assumes that there is Gaussian broadening in Prel (to
c             determine FWHMST and broadness of models.
c           -------------------------------------------------------------------
            if (lstage .eq. 1   .and.   .not.lshape(nmetab)   .and.
     1          .not.scafwh   .and.   imethd .ne. 3) then
               factor = exp(-deltat * exrt2(nmetab))
               term = 1.
               do 380 junfil = 1, nunfil
                  basist(junfil, nmetab) = term * basist(junfil, nmetab)
                  term = term * factor
 380           continue
            end if
c           ------------------------------------------------------------------
c           Scale spectrum to be consistent with Basis Set.  If total SIAMP is
c             the number of visible protons, then CONC should be in mM.
c           Correct normalization (to 5 places) was verified by checking that
c           PPMINC*SUM = SIAMP(1, JSIMUL) when NGAU=1, by evaluating SUM above.
c           ------------------------------------------------------------------
            if (scasim) then
               if (area_met_norm .le. 0.) then
c                 ------------------------------------------------------------
c                      This error message will be called once for each
c                 simulated metabolite that cannot be scaled in the
c                 Final Analysis.
c                 ------------------------------------------------------------
                  if (.not.dofull) call errmes (34, 3, chsubp)
               else
c                  factor = siamp(1,jsimul) / (ppminc * sum)
c                  if (lprint .gt. 0) then
c                     write (lprint, 5380) factor
c 5380                format (/'The following normalization check ',
c     1                        'should equal 1:', 1pe12.4)
c                  end if
                  do 385 junfil = 1, nunfil
                     basist(junfil, nmetab) = area_met_norm *
     1                                        basist(junfil, nmetab)
 385              continue
               end if
            end if
c           -------------------------------------------------------------------
c           Dump model spectrum for viewing absolute value.
c           -------------------------------------------------------------------
            if (sidump(jsimul)) then
               do 390 jdata = 1, ndata
                  datat(jdata) = basist(jdata, nmetab)
 390           continue
c               open (22, file='/tmp/tmp/dump.raw', status='unknown')
c               write (22, 9370)
c 9370          format (' $NMID fmtdat=''(1p2e15.6)'' $END')
c               write (22, 9375) (datat(j), j=1,ndata)
c 9375          format (1p2e15.6)
               istago = 1
               call errmes (23, 4, chsubp)
            end if
c           -------------------------------------------------------------------
c           Add 2 additional BASIST for 1st-order Taylor terms for shift &
c             broadening for Preliminary Analysis.
c           -------------------------------------------------------------------
            if (lstage .eq. 1   .and.   sitayl(jsimul)) then
               if (nmetab + 2 .gt. mmetab) call errmes (36, 4, chsubp)
               rtime = 0.
               do 410 junfil = 1, nunfil
                  cterm(1) = rtime * basist(junfil, nmetab)
                  basist(junfil, nmetab + 1) = cterm(1)
                  basist(junfil, nmetab + 2) = cmplx(-aimag(cterm(1)),
     1                                               real(cterm(1)))
                  basist(junfil + nunfil, nmetab + 1) = (0.,0.)
                  basist(junfil + nunfil, nmetab + 2) = (0.,0.)
                  rtime = rtime + deltat
 410           continue
               if (jsimul .lt. 10) then
                  write (metabo, 5410) jsimul
 5410             format ('Brod-', i1)
               else
                  write (metabo, 5415) jsimul
 5415             format ('Brod', i2)
               end if
               NMETAB=NMETAB+1
               table_top(nmetab) = .false.
               NACOMB(NMETAB)=METABO
               NCOMPO(NMETAB)=1
               LCOMPO(1,NMETAB)=NMETAB
               nonneg(nmetab) = .false.
 
               if (jsimul .lt. 10) then
                  write (metabo, 5420) jsimul
 5420             format ('Shif-', i1)
               else
                  write (metabo, 5425) jsimul
 5425             format ('Shif', i2)
               end if
               NMETAB=NMETAB+1
               table_top(nmetab) = .false.
               NACOMB(NMETAB)=METABO
               NCOMPO(NMETAB)=1
               LCOMPO(1,NMETAB)=NMETAB
               nonneg(nmetab) = .false.
            end if
 310     continue
      end if
      IF (NMETAB .LE. 0) CALL ERRMES (7, 4, CHSUBP)
      IF (NMETAB.LT.NUSE1 .AND. LSTAGE.EQ.1   .and.   .not.omit_chless)
     1     CALL ERRMES (8, 2, CHSUBP)
c     -------------------------------------------------------------------------
c     Make CPRIOR_SHIFT for constraining shifts of a group from its mean.
c     -------------------------------------------------------------------------
      if (lstage .eq. 2) call make_cgroup_shift()
c     -------------------------------------------------------------------------
c     Water-scaling when NOBASI=T
c     -------------------------------------------------------------------------
      if (nobasi   .and.   dows   .and.   lstage .eq. 2   .and.
     1    havh2o  .and.
     2    fcalib .gt. .999999   .and.   fcalib .lt. 1.000001) then
         area_met_norm = 1.
         call water_scale ()
      end if
      if (basout   .and.   ndata_freq .ne. ndata   .and.
     1    lstage .eq. 2   .and.   .not.absval) then
c        ---------------------------------------------------------------------
c        Intel Windows version outputs Namelists that cannot be read; so do
c           not bother to fix Cyg Namelist.
c        ---------------------------------------------------------------------
C        call fix_g77_namelist(21) !sun
         close (21)
      end if
      CLOSE (LBASIS)
      return
 810  call errmes (28, 4, chsubp)
 820  call errmes (29, 4, chsubp)
 830  call errmes (30, 4, chsubp)
 831  call errmes (31, 4, chsubp)
 832  call errmes (32, 4, chsubp)
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine make_cgroup_shift()
c
      INCLUDE 'lcmodel.inc'
      external ilen
c
      chsubp = 'GRPSHF'
      nrow_group_shift = 0
      if (ngrsh .gt. mgroup_shift) then
         call errmes (1, 3, chsubp)
         ngrsh = mgroup_shift
      end if
      do 110 jgroup = 1, ngrsh
         nin_group = 0
         do 130 jmetab = 1, nmetab
            lchgr = ilen(chgrsh(jgroup))
            if (lchgr .le. 0) go to 130
            if (index(nacomb(jmetab), chgrsh(jgroup)(:lchgr)) .eq. 1)
     1           nin_group = nin_group + 1
 130     continue
         if (nin_group .le. 1) go to 110
         if (sdgrsh(jgroup) .le. 0.) call errmes(2, 4, chsubp)
c        ----------------------------------------------------------------------
c        Might have to change SDSHIF_MIN  for other nuclei.
c        ----------------------------------------------------------------------
         SDSHIF_MIN = ppminc / 50.
         if (sdgrsh(jgroup) .lt. SDSHIF_MIN) then
            call errmes (3, 2, chsubp)
            sdgrsh(jgroup) = SDSHIF_MIN
         end if
         term = 1. / float(nin_group)
         do 140 kmetab = 1, nmetab
            lchgr = ilen(chgrsh(jgroup))
            if (lchgr .le. 0) go to 140
            if (index(nacomb(kmetab), chgrsh(jgroup)(:lchgr)) .ne. 1)
     1           go to 140
            nrow_group_shift = nrow_group_shift + 1
            lmetab_shift_prior(nrow_group_shift) = kmetab
C           -------------------------------------------------------------------
C           Convert SD from ppm to radians/s.
C           -------------------------------------------------------------------
            sdgroup_shift_row(nrow_group_shift) = sdgrsh(jgroup) * 2. *
     1                                            PI * HZPPPM
            do 150 j = 1, mmetab
               cgroup_shift(nrow_group_shift, j) = 0.
 150        continue
            do 160 jmetab = 1, nmetab
               lchgr = ilen(chgrsh(jgroup))
               if (lchgr .le. 0) go to 160
               if (index(nacomb(jmetab), chgrsh(jgroup)(:lchgr)) .ne. 1)
     1              go to 160
               cgroup_shift(nrow_group_shift, jmetab) = -term
 160        continue
            cgroup_shift(nrow_group_shift, kmetab) =
     1           cgroup_shift(nrow_group_shift, kmetab) + 1.
 140     continue
 110  continue
c     -------------------------------------------------------------------------
c     Dump matrix of group-shift priors (if IPDUMP >= 3)
c     -------------------------------------------------------------------------
      if (min0(nrow_group_shift, lprint) .gt. 0) then
         if (ipdump .ge. 3) then
            write (lprint, 5210) (nacomb(j), j = 1, nmetab)
 5210       format (//20x, 'Prior matrix for group shifts'//
     1             (8x, (10(6x, a6))))
            do 210 j = 1, nrow_group_shift
               write (lprint, 5220) nacomb(lmetab_shift_prior(j)),
     1                  (cgroup_shift(j, jmetab), jmetab = 1, nmetab)
 5220          format (/2x, a6, 1p10e12.3 / (8x, 1p10e12.3))
 210        continue
         else
            write (lprint, 5230) (nacomb(lmetab_shift_prior(jrow)),
     1           jrow = 1, nrow_group_shift)
 5230       format (//' Group-shift priors used for:'/(a6))
         end if
      end if
      end
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine parse_chsimu ()
c
c     Parse CHSIMU input strings into component parts shown below:
c     CHSIMU must be of the form (all in one record):
c CHSIM @ SIPPM +- SISDSH FWHM= SIFWMN < SIFWEX +- SIFWSD AMP= SIAMP
C       @ SIPPM           FWHM= SIFWMN                    AMP= SIAMP
c
      INCLUDE 'lcmodel.inc'
      external ilen
      character chreturn*(mchmet)
c
      chsubp = 'PARSIM'
      do 210 jsimul = 1, nsimul
         ngau(jsimul) = 1
         len_chsimu = ilen(chsimu(jsimul))
         istart = 1
         call get_field ('@', 1, 1, 0, chsim(jsimul), freturn,
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 0
            go to 800
         end if
         call get_field ('+-', 2, 2, 0, chreturn, sippm(1, jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 1
            go to 800
         end if
         call get_field ('FWHM=', 5, 2, 0, chreturn, sisdsh(jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 2
            go to 800
         end if
         call get_field ('<', 1, 2, 0, chreturn, sifwmn(1, jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 3
            go to 800
         end if
         call get_field ('+-', 2, 2, 0, chreturn, sifwex(jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 4
            go to 800
         end if
         call get_field ('AMP=', 4, 2, 0, chreturn, sifwsd(jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 5
            go to 800
         end if
         call get_field ('@', 1, 2, 1, chreturn, siamp(1, jsimul),
     1                   istart, len_chsimu, chsimu(jsimul))
         if (istart .le. 0) then
            ierr = 6
            go to 800
         end if
         if (istart .gt. len_chsimu) go to 210
         do 310 jgau = 2, mgau
            ngau(jsimul) = jgau
            call get_field ('FWHM=', 5, 2, 0, chreturn,
     1                      sippm(jgau, jsimul),
     1                      istart, len_chsimu, chsimu(jsimul))
            if (istart .le. 0) then
               ierr = 7
               go to 800
            end if
            call get_field ('AMP=', 4, 2, 0, chreturn,
     1                      sifwmn(jgau, jsimul),
     1                      istart, len_chsimu, chsimu(jsimul))
            if (istart .le. 0) then
               ierr = 8
               go to 800
            end if
            call get_field ('@', 1, 2, 1, chreturn,
     1                      siamp(jgau, jsimul),
     1                      istart, len_chsimu, chsimu(jsimul))
            if (istart .le. 0) then
               ierr = 9
               go to 800
            end if
            if (istart .gt. len_chsimu) go to 210
 310     continue
c        ----------------------------------------------------------------------
c        Error -- NGAU > MGAU
c        ----------------------------------------------------------------------
         if (lprint .gt. 0) write (lprint, 5210)
     1        chsimu(jsimul)(1:132), chsimu(jsimul)(133:264),
     2        chsimu(jsimul)(265:396), chsimu(jsimul)(397:528), istart
 5210    format ('Incorrect CHSIMU follows:', / a132 / a132 / a132 /
     1           a132 / 'ISTART =', i3)
         call errmes (0, 4, chsubp)
 210  continue
      return
 800  if (lprint .gt. 0) write (lprint, 5210)
     1     chsimu(jsimul)(1:132), chsimu(jsimul)(133:264),
     2     chsimu(jsimul)(265:396), chsimu(jsimul)(397:528), istart
      call errmes (100*ierr + jsimul, 4, chsubp)
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine set_lshape_false ()
      INCLUDE 'lcmodel.inc'
      lshape(nmetab) = .false.
      if (fwhmba .ge. fwhmst) return
      fwhm_extra = sqrt(fwhmst**2 - fwhmba**2)
      rsd = pi * fwhm_extra * deltat * hzpppm /
     1      sqrt(2. * alog(2.))
      expmax = 2. * alog(rrange)
      do 110 jdata = 1, ndata
         expon = (rsd * float(jdata - 1))**2
         if (expon .lt. expmax) then
            basist(jdata, nmetab) = exp(-.5 * expon) *
     1                              basist(jdata, nmetab)
         else
            basist(jdata, nmetab) = (0.,0.)
         end if
 110  continue
      return
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine water_scale ()
      INCLUDE 'lcmodel.inc'
      chsubp = 'WSCALE'
      if (iaverg .eq. 1   .or.   iaverg .eq. 4) then
         area_water = 1.
      else
         area_water = areawa(2)
      end if
      if (lprint .gt. 0) write (lprint, 5110) area_water
 5110 format (//'Area of unsuppressed water peak =', 1pe13.5)
      if (area_water .le. 0.) then
         call errmes (1, 3, chsubp)
         go to 800
      end if
      if (amin1(atth2o, wconc) .le. 0.) then
         call errmes (2, 3, chsubp)
         go to 800
      end if
      water_norm = area_water / (2. * atth2o * wconc)
      fcalib = area_met_norm / water_norm
      wsdone = .true.
      if (lprint .gt. 0) write (lprint, 5250) fcalib
 5250 format ('FCALIB =', 1pe15.5/)
      do 260 j = 1, nunfil
         datat(j) = fcalib * datat(j)
 260  continue
      do 270 jy = 1, ny
         cy(jy) = cy(jy) * fcalib
 270  continue
      rmsamp = rmsamp * fcalib
      if (imethd .eq. 2) then
         do 310 j = 1, nmetab
            conc_expect(j) = conc_expect(j) * fcalib
 310     continue
      end if
 800  return
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function areawa(istage)
C
C  Compute integral of water peak.
C  ISTAGE != 2 to avoid using AREAW2 when calling from AVERAGE, which can
c              cause big errors (and often an abort) with very weak water
c              signal (as sometimes in bone).  In this case, whole spectrum
c              (including lipids) is used for weighting in AVERAG.
c
      INCLUDE 'lcmodel.inc'
      external areaw2
      double precision sx, sy, sxy, sxx, term1, xterm, yterm
      chsubp = 'AREAWA'
      if (iareaw .eq. 2   .and.   istage .eq. 2) then
         areawa = areaw2()
         return
      end if
c     ------------------------------------------------------------------------
c     Log-linear regression for water integral.
c     ------------------------------------------------------------------------
      npts = nwsend - nwsst + 1
      if (nwsst .lt. 1   .or.   nwsend .gt. nunfil   .or.
     1    npts .lt. 10) then
         call errmes (1, 3, chsubp)
         go to 800
      end if
      sx = 0.d0
      sy = 0.d0
      sxy = 0.d0
      sxx = 0.d0
      do 110 j = nwsst, nwsend
         xterm = dble(float(j))
         yterm = cabs(h2ot(j))
         if (yterm .le. 0.d0) then
            call errmes(4, 3, chsubp)
            go to 800
         end if
         yterm = alog(sngl(yterm))
c         write (LPRINT, 9110) j, alog10(cabs(h2ot(j)))
c 9110    format (i4, 1pe15.7)
         sx = sx + xterm
         sy = sy + yterm
         sxy = sxy + xterm * yterm
         sxx = sxx + xterm**2
 110  continue
      term1 = dble(float(nwsend - nwsst + 1)) * sxx
      denom = term1 - sx**2
      if (abs(denom) .lt. 1.e-10 * sngl(term1)) then
c        ---------------------------------------------------------------------
c        The difference of the two nonnegative terms has lost too much
c          precision.
c        ---------------------------------------------------------------------
         call errmes (2, 3, chsubp)
         go to 800
      end if
      rnum = sxx * sy - sx * sxy
c     ------------------------------------------------------------------------
c     sqrt(ndata) below is needed, because the amplitude is an inverse FFT of
c       the spectrum, and this inverse transform is divided by sqrt(ndata).
c       The value at zero time is thus the sum of all spectral points divided
c       by sqrt(ndata). The metabolite spectral peak, on the other hand is
c       simply a sum of the points (without any division by sqrt(ndata)).
c     0.5 is needed to get good agreement with water-scaling-2 (analyzing
c       water as Scyllo or Cr).  This may be related to the 0.5 that must
c       multiply DATAT(1) when ABSVAL=T (discussed above and in Notes of
c       990830 & 001203).
c     PPMINC is nesessary to get an area, rather than just the FFT sum.
c     ------------------------------------------------------------------------
      expmax = alog(rrange)
      if (abs(rnum) .ge. expmax * abs(denom)) then
         call errmes (3, 3, chsubp)
         go to 800
      end if
      areawa = .5 * ppminc * exp(rnum / denom) * sqrt(float(2 * nunfil))
      return
 800  areawa = -1.
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function areaw2()
c
c AREAW2 = area under (unsuppressed) water peak in H2OF.
c Uses peak phasing routine of MakeBasis.
c Uses integration routine of AREABA, which is based on MakeBasis routine.
      INCLUDE 'lcmodel.inc'
 
      chsubp = 'AREAW2'
c     ------------------------------------------------------------------------
c     H2OF_WORK = smoothed water spectrum.
c     No zero-filling means only half the number of points (N) in the FT,
c        which is normalized by dividing by the sqrt(N).  So must divide the
c        FTs by sqrt(2), as they would be with the original NDATA points.
c     ------------------------------------------------------------------------
      nunfil_half = nunfil / 2
      ppminc2 = ppminc * 2.
      RSD=2.*PI*SDSMOO(4)/(PPMINC2*float(nunfil))
      DO 120 JUNFIL=1,NUNFIL
         h2ot_work(JUNFIL)=H2Ot(JUNFIL)*
     1                     EXP(-.5*(RSD*FLOAT(JUNFIL-1))**2)
 120  CONTINUE
      call csft_r(h2ot_work, h2of_work, Nunfil)
      call csft_r(h2ot, h2of, Nunfil)
c     ----------------------------------------------------------------------
c     No zero-filling means only half the number of points (N) in the FT,
c        which is normalized by dividing by the sqrt(N).  So must divide the
c        FTs by sqrt(2), as they would be with the original NDATA points.
c     ----------------------------------------------------------------------
      rsqrt2 = sqrt(.5)
      do 130 junfil = 1, nunfil
         h2of_work(junfil) = h2of_work(junfil) * rsqrt2
         h2of(junfil) = h2of(junfil) * rsqrt2
 130  continue
c     ------------------------------------------------------------------------
c     ppmh2o_act = true position of water peak.
c     ------------------------------------------------------------------------
      if (ppm_water_range  .le. 0.) then
         ppmh2o_corr = ppmh2o
      else
         kystrt = max0(1,
     1            nint((PPMcen-PPMH2O-ppm_water_range)/PPMINC2)+
     2                NUNFIL_half+1)
         kyend = min0(nunfil,
     1           nint((PPMcen-PPMH2O+ppm_water_range)/PPMINC2)+
     2                 NUNFIL_half+1)
         lmax = 0
         rmax = -rrange
         do 140 j = kystrt, kyend
            if (cabs(h2of(j)) .gt. rmax) then
               rmax = cabs(h2of(j))
               lmax = j
            end if
 140     continue
         ppmh2o_corr = ppmcen - float(lmax - 1 - nunfil_half) * ppminc2
         if (lprint .gt. 0) write (lprint, 5140) ppmh2o_corr
 5140    format (/'Corrected PPMH2O =', f8.4)
      end if
      ly = (ppmcen - ppmh2o_corr) / ppminc2 + nunfil_half + 1
      kystrt = NINT((PPMcen-PPMH2O_corr-HWDWAT(1))/PPMINC2)+
     1         NUNFIL_half + 1
      KYEND=NINT((PPMcen-PPMH2O_corr+HWDWAT(1))/PPMINC2)+
     1      NUNFIL_half + 1
      NYPEAK=KYEND-KYSTRT+1
      if (kystrt .lt. 3  .or.  kyend .gt. nunfil-2   .or.
     1    nypeak .le. 0) call errmes (1, 4, chsubp)
      if (.not.havh2o) call errmes (2, 4, chsubp)
c 9140 format ('ly, lmax =' 2i6)
c      write (6, 9140) ly, lmax
c     -------------------------------------------------------------------------
c     H2OF_WORK is smoothed water spectrum on input to GETPHA and phased
c               (unsmoothed) spectrum (with PPMINC2 spacing) on return.
c     -------------------------------------------------------------------------
      call getpha (KYSTRT, KYEND, h2oF, h2of_work, nunfil, RADIAN,
     1            NYPEAK, rwork, rwork2,
     2            degzer_calc)
      if (lprint .gt. 0) write (lprint, 5210) degzer_calc
 5210 format (/'Zero-order phase correction for unsuppressed water =',
     1        f6.1/)
      kystrt = NINT((PPMcen-PPMH2O_corr-HWDWAT(2))/PPMINC2)+
     1         NUNFIL_half + 1
      KYEND=NINT((PPMcen-PPMH2O_corr+HWDWAT(2))/PPMINC2)+
     1         NUNFIL_half + 1
      NYPEAK=KYEND-KYSTRT+1
      if (kystrt .lt. 3  .or.  kyend .gt. nunfil-2   .or.
     1    nypeak .le. 0) call errmes (1, 4, chsubp)
      nwndo = nint(ppmbas(2) / ppminc)
      call integrate (h2of_work, ppminc2,
     1                area_water,
     2                kyend, kystrt, ly, nunfil, nwndo)
      areaw2 = area_water
      return
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine getpha (KYSTRT, KYEND, dataf, dataw, nunfil, radian,
     1                   nypeak, yorig, yinterp,
     2                   degzer_calc)
c
C  Take optimal 0-order phase correction as the one producing the minimum sum
c    of absolute differences of the spectral values equally distant from the
c    maximum in the smoothed real spectrum between KYSTRT and KYEND.  The sum
c    is normalized by dividing by abs(max - min).
c   DATAW : must contain smoothed spectrum (with PPMINC2 spacing) in necessary
c           region.
c
      parameter (ninterp = 10)
      CHARACTER CHSUBP*6
      COMPLEX CFACT, CFACT_BEST, CFINC, cterm, DATAF(nunfil),
     1        DATAW(nunfil)
      REAL yinterp(nunfil), yorig(nunfil)
      CHSUBP='GETPHA'
 
c     -------------------------------------------------------------------------
c     Start main loop through 0-order phase angles.
c     -------------------------------------------------------------------------
      CFINC=CEXP(CMPLX(0.,RADIAN))
      CFACT=(1.,0.)
      cfact_best = cfact
      ABSMIN=1.E30
      kystrt_best = 999999
      kyend_best = -999999
      DO 110 JDEG=1,360
         CFACT=CFACT*CFINC
c        ----------------------------------------------------------------------
c        YORIG = real part of spectrum with trial 0-order phase.
c        ----------------------------------------------------------------------
         do 120 jy = max0(1, kystrt-2-nypeak),
     1               min0(nunfil, kyend+2+nypeak)
            yorig(jy) = real( dataw(jy) * cfact )
 120     continue
c        ----------------------------------------------------------------------
c        LY is at max in a smoothed YORIG
c        ----------------------------------------------------------------------
         rmin = 1.e37
         RMAX=-1.E37
         ly = 999999
         DO 210 JY=KYSTRT,KYEND
            TERM=yorig(JY)
            IF (TERM .GT. RMAX) THEN
               RMAX=TERM
               LY=JY
            END IF
            IF (TERM .lt. Rmin) THEN
               RMin=TERM
            END IF
 210     CONTINUE
         if (rmax + rmin   .le.   0.) go to 110
c        ----------------------------------------------------------------------
c        Center test region at LY
c        ----------------------------------------------------------------------
         KYSTRT_new=LY-NYPEAK/2
         KYEND_new=LY+NYPEAK/2
         IF (KYSTRT_new.Le.2 .OR. KYEND_new.Ge.NUNFIL-1)
     1        CALL ERRmes (1, 4, CHSUBP)
         NYPEAK=KYEND_new-KYSTRT_new+1
c        ----------------------------------------------------------------------
c        YINTERP = The linearly interpolated YORIG at NINTERP interpolations
c                  on each side to allow for the fact that the peak on the
c                  discrete grid may not be at the center of symmetry.
c        ----------------------------------------------------------------------
         fract = 1.
         rinterp = ninterp
         delta_fract = .5 / rinterp
         do 310 jinterp = 1, 2 * ninterp + 1
            fract = fract - delta_fract
            fract2 = 1 - fract
            if (jinterp .eq. ninterp) fract = 1.
            if (jinterp .lt. ninterp) then
               do 312 jy = kystrt_new - 1, kyend_new + 1
                  yinterp(jy) = fract * yorig(jy) + fract2 * yorig(JY-1)
 312           continue
            else
               do 314 jy = kystrt_new - 1, kyend_new + 1
                  yinterp(jy) = fract * yorig(jy) + fract2 * yorig(JY+1)
 314           continue
            end if
            rmax = yorig(LY)
            rmin = rmax
            sum = 0.
            DO 320 JYREL=1,NYPEAK/2
               rmax = amax1(rmax, yinterp(LY+JYREL), yinterp(LY-JYREL))
               rmin = amin1(rmin, yinterp(LY+JYREL), yinterp(LY-JYREL))
               sum = sum + sqrt(float(jyrel)) *
     1                     abs(yinterp(LY+JYREL) - yinterp(LY-JYREL))
 320        continue
            span = rmax - rmin
            if (span .le. 0.) go to 110
            acrit = sum / span
            IF (ACRIT .LT. ABSMIN) THEN
               kystrt_best = kystrt_new
               kyend_best = kyend_new
               ABSMIN=ACRIT
               DEGZER_CALC=FLOAT(JDEG)
               CFACT_BEST=CFACT
            END IF
 310     continue
  110 CONTINUE
      kystrt = kystrt_best
      kyend = kyend_best
      nypeak = kyend - kystrt + 1
      DO 410 J=1,NUNFIL
         DATAW(J)=DATAF(J)*CFACT_BEST
  410 CONTINUE
      END
 
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE INTEGRATE (DATAF, ppminc2,
     2                      RINTEG,
     3                      KYEND, KYSTRT, ly, Nunfil, NWNDO)
C
C  Crude integration of peak (with PPMINC2 spacing) between KYSTRT and KYEND
C    with the baseline taken as the line running from the mean over the 2
c    bordering regions NWNDO wide.
C
      COMPLEX DATAF(Nunfil)
C     -------------------------------------------------------------------------
c     Restricting the integration range by moving KY* inward to 1 beyond the
c       minimum on each side led to erratic results for hennig/9911/press270,
c       where there were eddy-current distortions with one side sharply
c       dropping to min near the peak.
c     So, the integration range is shrunk symmetrically to 1 point beyond the
c       max distance (over both sides) of the min from the peak.
C     -------------------------------------------------------------------------
      rmin = 1.e+30
      lmin = kystrt
      do 102 jy = kystrt, ly - 1
         term = real(dataf(jy))
         if (term .lt. rmin) then
            rmin = term
            lmin = jy
         end if
 102  continue
      ldist_left = ly - lmin
c
      rmin = 1.e+30
      lmin = kyend
      do 104 jy = ly + 1, kyend
         term = real(dataf(jy))
         if (term .lt. rmin) then
            rmin = term
            lmin = jy
         end if
 104  continue
      ldist_right = lmin - ly
      ldist_max = max0( ldist_left, ldist_right)
      kystrt = ly - ldist_max
      kyend = ly + ldist_max
c
      COUNT=0.
      AVG_LEFT=0.
      DO 110 JY=MAX0(1,KYSTRT-NWNDO),KYSTRT-1
         COUNT=COUNT+1.
         AVG_LEFT=AVG_LEFT+REAL(DATAF(JY))
  110 CONTINUE
      AVG_LEFT=AVG_LEFT/COUNT
      COUNT=0.
      AVG_RIGHT=0.
      DO 120 JY=KYEND+1,MIN0(Nunfil,KYEND+NWNDO)
         COUNT=COUNT+1.
         AVG_RIGHT=AVG_RIGHT+REAL(DATAF(JY))
  120 CONTINUE
      avg_right = avg_right / count
C     -------------------------------------------------------------------------
C     Subtracting AVG from every data point in the peak range is
C       equivalent to using an unbiased estimate of a line through the
C       points in the bordering regions.  It is more robust to outliers
C       due to neighboring peaks than least squares.
C     -------------------------------------------------------------------------
      avg = .5 * (avg_left + avg_right)
      RINTEG=0.
      DO 130 JY=KYSTRT,KYEND
         RINTEG=RINTEG+REAL(DATAF(JY))
  130 CONTINUE
      rinteg = rinteg - float(kyend - kystrt + 1) * avg
      rinteg = rinteg * ppminc2
c 9130 format('rinteg =', 1pe13.5/ 'ly =', i5/ 'kystrt =', i5/
c     1       'kyend =', i5/ 'avg_right =', e13.5/ 'avg_left=', e13.5/
c     2       'nwndo =', i5/ 'count =', e13.5/)
c      write (6, 9130) rinteg, ly, kystrt, kyend, avg_right,
c     1                     avg_left, nwndo, count
c 9140 format (i5, 1pe13.5)
c      write (6, 9140) (j, real(dataf(j)), j = kystrt, kyend)
c      if (rinteg .le. 0.) stop
      END
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function areaba(basisf, ppminc_arg, nunfil_arg)
C
C  Compute normalized integral of metabolite peak.  This area is
c    normalized by dividing by N1HMET & ATTMET.  So, it should correspond to
c    1mM of protons.
c
C
      INCLUDE 'lcmodel.inc'
      external icycle
      complex basisf(*)
      chsubp = 'AREABA'
      areaba = 0.
      ndata_arg = 2 * nunfil_arg
      if (sptype(:10) .eq. 'mega-press') then
         do 120 jy = 1, ndata_arg
            basisf(jy) = -basisf(jy)
 120     continue
      end if
      ly = nint((ppmcen - wsppm) / ppminc_arg) + 1
c     ----------------------------------------------------------------------
c     LY = corrected index of metabolite peak
c     ----------------------------------------------------------------------
      if (r_areaba .gt. 0.) then
         nyhalf = nint(r_areaba * desdsh / ppminc_arg)
         kystrt = ly - nyhalf
         kyend = ly + nyhalf
         if (kystrt .le. -nunfil_arg   .or.
     1       kyend .ge. nunfil_arg) then
            call errmes (1, 3, chsubp)
            go to 800
         end if
         rmax = -rrange
         do 140 jy = kystrt, kyend
            term = real(basisf(icycle(jy, ndata_arg)))
            if (term .gt. rmax) then
               rmax = term
               ly = jy
            end if
 140     continue
         term = ppmcen - (ly - 1) * ppminc_arg
         if (lprint .gt. 0) write (lprint, 5140) term
 5140    format (/'Corrected WSPPM =', f8.4)
      end if
 
      nwndo = nint(ppmbas(1) / ppminc_arg)
      if (ldump(5)) write (lprint, 9210) nwndo, ppmbas(1), ppminc_arg
 9210 format ('NWNDO, PPMBAS(1), PPMINC =', i4, 1p2e15.7)
      hwdsca = .5 * rfwbas * fwhmba
c     -------------------------------------------------------------------------
c          When ABSVAL=T, increase HWDSCA by another factor of 1.5, because of
c     extremely broad absolute-value tails.  Even with Cr, cannot go any
c     higher, since FWHMBA has already been doubled.  For example, with
c     original FWHMBA=.04 and RFWBAS=10, we would get
c     HWDSCA = .5 * 10 * 1.5 * 2 * .04 = .60, which would be more than halfway
c     to CrCH2.
c     -------------------------------------------------------------------------
      if (absval) hwdsca = hwdsca * 1.5
      nyhalf = nint(hwdsca / ppminc_arg)
      kystrt = ly - nyhalf
      kyend = ly + nyhalf
      if (kystrt - nwndo .le. -nunfil_arg   .or.
     1    kyend + nwndo .ge. nunfil_arg
     2    .or.   min0(nyhalf, nwndo) .le. 0) then
         call errmes (1, 3, chsubp)
         go to 800
      end if
      rmin = 1.e+30
      lmin = kystrt
      do 210 jy = kystrt, ly - 1
         term = real(basisf(icycle(jy, ndata_arg)))
         if (ldump(5)) write (lprint, 9250) jy,
     1                       REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
         if (term .lt. rmin) then
            rmin = term
            lmin = jy
         end if
 210  continue
      ldist_left = ly - lmin
      if (ldump(5)) write (lprint, 9211) lmin, ldist_left, ly
 9211 format ('LMIN, LDIST, LY =', 3i5)
 
      rmin = 1.e+30
      lmin = kyend
      do 220 jy = ly + 1, kyend
         term = real(basisf(icycle(jy, ndata_arg)))
         if (ldump(5)) write (lprint, 9250) jy,
     1                        REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
         if (term .lt. rmin) then
            rmin = term
            lmin = jy
         end if
 220  continue
      ldist_right = lmin - ly
      if (ldump(5)) write (lprint, 9211) lmin, ldist_right, ly
      ldist_max = max0( ldist_left, ldist_right)
      kystrt = ly - ldist_max
      kyend = ly + ldist_max
c
      COUNT=0.
      AVG_LEFT=0.
      DO 230 JY = KYSTRT-NWNDO, KYSTRT-1
         COUNT=COUNT+1.
 
         if (ldump(5)) write (lprint, 9250) jy,
     1                        REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
         if (count .gt. 0.) AVG_LEFT=AVG_LEFT+
     1                               REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
  230 CONTINUE
      AVG_LEFT=AVG_LEFT/COUNT
      if (ldump(5)) write (lprint, 9240) nint(count), avg_left
      COUNT=0.
      AVG_RIGHT=0.
      DO 240 JY = KYEND+1, KYEND+NWNDO
         COUNT=COUNT+1.
 
         if (ldump(5)) write (lprint, 9250) jy,
     1                       REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
         AVG_RIGHT=AVG_RIGHT+REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
  240 CONTINUE
      if (count .gt. 0.) avg_right = avg_right / count
      if (ldump(5)) write (lprint, 9240) nint(count), avg_right
 9240 format (i4, 1pe15.7/)
C     -------------------------------------------------------------------------
C     Subtracting AVG from every data point in the peak range is
C       equivalent to using an unbiased estimate of a line through the
C       points in the bordering regions.  It is more robust to outliers
C       due to neighboring peaks than least squares.
C     -------------------------------------------------------------------------
      avg = .5 * (avg_left + avg_right)
      AREA_MET=0.
      DO 250 JY=KYSTRT,KYEND
         AREA_MET=AREA_MET+REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
 
         if (ldump(5)) write (lprint, 9250) jy,
     1                       REAL(BASISF(ICYCLE(JY, NDATA_ARG)))
 9250    format (i5, 1pe15.7)
 
  250 CONTINUE
      if (lprint .gt. 0) write (lprint, 9251) area_met
 9251 format ('Uncorrected AREA_MET =', 1pe13.5)
      area_met = area_met - float(kyend - kystrt + 1) * avg
      if (n1hmet .le. 0   .or.   attmet .le. 0.) then
         call errmes (2, 3, chsubp)
         go to 800
      end if
c     -------------------------------------------------------------------------
c     PPMINC_ARG is necessary to get area, rather than just the FFT sum.
c     -------------------------------------------------------------------------
      areaba = ppminc_arg * area_met / (float(n1hmet) * attmet)
      if (lprint .gt. 0) write (lprint, 5250) areaba
 5250 format ('Normalized area of reference Basis singlet =', 1pe15.5/)
c     -------------------------------------------------------------------------
c     If AREABA<=0, then neither water-scaling nor scaling of simulated
c       spectra will be possible.
c     -------------------------------------------------------------------------
      if (areaba .le. 0.) call errmes (3, 3, chsubp)
 800  if (sptype(:10) .eq. 'mega-press') then
         do 820 jy = 1, ndata_arg
            basisf(jy) = -basisf(jy)
 820     continue
      end if
      return
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE COMBIS ()
C
C  Compute NCOMPO and LCOMPO for combinations of metabolites.
C
      INCLUDE 'lcmodel.inc'
      LOGICAL ATEND
      CHSUBP='COMBIS'
      ncombi = min0(ncombi, mpmet)
c     -------------------------------------------------------------------------
c     Special (tedious) case to eliminate combinations of 2's if all 3 of Cho
c       GPC & PCh and CHCOMB=Cho+GPC+PCh are present.
c     -------------------------------------------------------------------------
      ncho = 0
      do 50 jmetab = 1, nmetab
         nacom2(jmetab) = ' '
         if (nacomb(jmetab) .eq. 'Cho'   .or.
     1       nacomb(jmetab) .eq. 'GPC'   .or.
     2       nacomb(jmetab) .eq. 'PCh') ncho = ncho + 1
 50   continue
      if (ncho .ge. 3) then
         do 52 jcombi = 1, ncombi
            if (CHCOMB(JCOMBI) .eq. 'Cho+GPC+PCh') go to 54
 52      continue
         go to 100
 54      do 55 jtry = 1, 3
            do 56 jcombi = 1, ncombi
               if (CHCOMB(JCOMBI) .eq. 'GPC+PCh'   .or.
     1             CHCOMB(JCOMBI) .eq. 'GPC+Cho'   .or.
     2             CHCOMB(JCOMBI) .eq. 'PCh+Cho') go to 57
 56         continue
            go to 100
 57         ncombi = ncombi - 1
            do 58 j = jcombi, ncombi
               chcomb(j) = chcomb(j + 1)
 58         continue
 55      continue
      end if
 100  JCONC=NMETAB
      DO 110 JCOMBI=1,NCOMBI
         JCONC=JCONC+1
         IF (JCONC .GT. MCONC) CALL ERRMES (1, 4, CHSUBP)
         NACOMB(JCONC)=CHCOMB(JCOMBI)
         NACOM2(JCONC)=CHCOM2(JCOMBI)
         table_top(jconc) = .true.
         NCOMPO(JCONC)=0
         ISTART=1
         DO 120 KMETAB=1,NMETAB
            LENGTH=INDEX(CHCOMB(JCOMBI)(ISTART:),'+')-1
            ATEND=LENGTH .EQ. -1
            IF (ATEND) THEN
C              ---------------------------------------------------------------
C              No more + characters.
C              ---------------------------------------------------------------
               LENGTH=INDEX(CHCOMB(JCOMBI)(ISTART:),' ')-1
               IF (LENGTH .EQ. -1) LENGTH=LEN(CHCOMB(JCOMBI))+1-ISTART
            END IF
            IF (LENGTH .GT. MCHMET) THEN
               IF (LPRINT .GT. 0) WRITE (LPRINT,5120) CHCOMB(JCOMBI)
 5120          FORMAT (' Incorrect CHCOMB =',A)
               CALL ERRMES (2, 4, CHSUBP)
            END IF
            DO 130 JMETAB=1,NMETAB
            IF (INDEX(NACOMB(JMETAB),' ')-1 .EQ. LENGTH) THEN
               IF (CHCOMB(JCOMBI)(ISTART:ISTART+LENGTH-1) .EQ.
     1             NACOMB(JMETAB)(1:LENGTH)) GO TO 135
               END IF
  130       CONTINUE
            JCONC=JCONC-1
            GO TO 110
  135       NCOMPO(JCONC)=NCOMPO(JCONC)+1
            ISTART=ISTART+LENGTH+1
            IF (NCOMPO(JCONC) .GT. MCOMPO) THEN
               IF (LPRINT .GT. 0) WRITE (LPRINT,5120) CHCOMB(JCOMBI)
               CALL ERRMES (3, 4, CHSUBP)
            END IF
            LCOMPO(NCOMPO(JCONC),JCONC)=JMETAB
            table_top(jconc) = table_top(jconc) .and. table_top(jmetab)
            IF (ATEND) GO TO 110
  120    CONTINUE
         JCONC=JCONC-1
  110 CONTINUE
      NCONC=JCONC
      jline = 0
      do 210 jconc = 1, nconc
         if (table_top(jconc)) then
            jline = jline + 1
            iconc_line_table(jline) = jconc
         end if
 210  continue
      ntable_top = jline
      do 220 jconc = 1, nconc
         if (.not.table_top(jconc)) then
            jline = jline + 1
            iconc_line_table(jline) = jconc
         end if
 220  continue
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE STARTV (ipass)
C
C  Crude fast analysis to get starting values for fwhm (FWHMST), DEGZER,
C    and DEGPPM, and to determine ISHIFD and NSIDES.
C
      INCLUDE 'lcmodel.inc'
      save degppm_sav_startv, degzer_sav_startv, ishfst, lwidth
      DOUBLE PRECISION SDSHBS
      INTEGER INTSHF(MSHIFT,2), ISHFST(2), LREFPK(MREFPK), LSHFMN(3),
     1        lshfmn_orig(2), LSHFMX(3), lshfmx_orig(2), LWIDTH(2)
      LOGICAL CNVRG1, dosub_bas_ccf
      real cpy2(mdata)
C      COMPLEX CYFIT(MY)!CT1
C      LOGICAL LTEST!CT1
C      COMMON /BLT1/ CYFIT, LTEST!CT1
      CHSUBP='STARTV'
      do 101 j = 1, 2
         ishfst(j) = 0
         lwidth(j) = 999999
 101  continue
      lshfbs_fixshf = 99999
c     -------------------------------------------------------------------------
c     PPMSHF cannot shift Analysis Window out of spectrum.
c     Max negative shift would shift first point of entire rearranged
c       spectrum to LDATST.
c     Max positive shift would shift last point of entire rearranged
c       spectrum to LDATEN.
c     Otherwise, the shift would require spectral data that do not exist.
c     -------------------------------------------------------------------------
      if (abs(ppmshf) .lt. 1.e5) then
         if ((ppmshf .le. 0.   .and.
     1       -ppmshf .gt. ppminc * float(ldatst))   .or.
     2       (ppmshf .gt. 0.   .and.
     3        ppmshf .gt. ppminc * float(ndata - ldaten))) then
            call errmes (26, 4, chsubp)
         else
            LGRID=1
            ISHFST(1)=PPMSHF/PPMINC
            LSHFMN(3)=ISHFST(1)
            LSHFMX(3)=ISHFST(1)
            LGRDMN=9999
            GO TO 300
         end if
      END IF
c     -----------------------------------------------------------------------
c     FIXSHF = T (used for muscle-4) for initial CCF with PPMREF to fix
c                the referencing shift.  Prel is only used for initial
c                phasing; its optimal shift is ignored.
c     PPMSHF (with FIXSHF=T) will fix the referencing shift at PPMSHF.
c     -----------------------------------------------------------------------
      if (fixshf) then
         rfwhcc = 0.
         fwhh2o = 0.
      end if
      IF (.NOT.(DOREFS(1) .OR. DOREFS(2))) CALL ERRMES (1, 4, CHSUBP)
      LSHFMN(3)=99999
      LSHFMX(3)=-99999
      DO 105 JRFSET=1,2
         IF (.NOT.DOREFS(JRFSET)) GO TO 105
C        ----------------------------------------------------------------------
C        Reduce NREFPK(JRFSET) if a PPMREF is out of analysis window.
C        ----------------------------------------------------------------------
         KREFPK=NREFPK(JRFSET)
         DO 106 J=1,NREFPK(JRFSET)
            IF (J.GT.KREFPK .OR. JRFSET.EQ.1) GO TO 108
            TEST=PPMREF(J,JRFSET)+HZREF(J,JRFSET)/HZPPPM
            IF (PPMST.LT.TEST .OR. PPMEND.GT.TEST) THEN
               KREFPK=KREFPK-1
               DO 107 K=J,KREFPK
                  PPMREF(K,JRFSET)=PPMREF(K+1,JRFSET)
                  HZREF(K,JRFSET)=HZREF(K+1,JRFSET)
  107          CONTINUE
            END IF
  106    CONTINUE
  108    NREFPK(JRFSET)=KREFPK
         LSHFMN(JRFSET)=NINT(SHIFMN(JRFSET)/PPMINC)
         LSHFMX(JRFSET)=NINT(SHIFMX(JRFSET)/PPMINC)
         LSHFMN_orig(JRFSET)=NINT(SHIFMN_orig(JRFSET)/PPMINC)
         LSHFMX_orig(JRFSET)=NINT(SHIFMX_orig(JRFSET)/PPMINC)
         IF (LSHFMN(JRFSET) .GT. LSHFMX(JRFSET)) CALL ERRMES (2, 4,
     1                                                        CHSUBP)
C        ----------------------------------------------------------------------
C        DATAF = (temporarily) Gaussian-smoothed FFT of zero-filled
C                (time-domain) DATAT (COMPLEX).
C        DATAF(1+nunfil): corresponds to PPMCEN; i.e., the spectrum is
c                         rearranged.
c                  /main/projects/lcm/lcmodel/src/misc/test-fft.f verifies
c                    that the 1st point is the zero-frequency transform (in
c                    the old non-rearranged version).
C        ----------------------------------------------------------------------
         RSD=2.*PI*SDSMOO(JRFSET)/(PPMINC*FNDATA)
         DO 110 JUNFIL=1,NUNFIL
            DATAF(JUNFIL)=DATAT(JUNFIL)*
     1                    EXP(-.5*(RSD*FLOAT(JUNFIL-1))**2)
            DATAF(NUNFIL+JUNFIL)=(0.,0.)
  110    CONTINUE
         CALL CFFT_r (DATAF, DATAF, NDATA, LWFFT,
     1              WFFTC)
C        ----------------------------------------------------------------------
C        CPY2 = smoothed power spectrum.
C        ----------------------------------------------------------------------
         DO 120 JDATA=1,NDATA
            cpy2(JDATA)=REAL(DATAF(JDATA))**2+AIMAG(DATAF(JDATA))**2
  120    CONTINUE
C        ----------------------------------------------------------------------
C        ISHFST = starting value (in grid points) for referencing shift in the
C                 spectrum.
C        Compute cross-correlation function (CCF) with a sum of unit reference
C          delta-functions at PPMREF(J,JRFSET).
C        Ordinarily, NREFPK(1)=1 and PPMREF(1,1)=4.65; i.e., the water peak is
C          used, despite the possible distortion in the peak position due to
C          water suppression.
C        Typically, NREFPK(2)=3 and PPMREF(*,2)=2.01, 3.03, 3.22; i.e., the
C          usual landmarks are used to extend the range of the grid search in
C          case the water peak has been strongly displaced by the water
C          suppression.
C        A positive shift shifts the data spectrum to the left (larger ppm).
C          Thus SHIFMN and SHIFMX specify the allowed range of the correction
C          to the apparent ppm postion of the data spectrum; i.e.,
C          (apparent ppm)+SHIFMN < true ppm < (apparent ppm)+SHIFMX.
C          SHIFMN <= SHIFMX always.
C          Usually SHIFMN<0, SHIFMX>0
C        The same sign convention applies to all other shift parameters, such
C          as LSHFM*, ISH*, etc.
C        ----------------------------------------------------------------------
C        DOSUB_BAS_CCF = T to subtract rough baseline from smoothed power
c                          spectrum.
c        Set DOSUB_BAS_CCF = F if one of the reference peaks is within
c          PPM_WATER_TOL (typically 0.08) of 4.68 (to avoid baseline
c          subtraction at the huge water peak).
C          PPM_WATER_TOL < 0 will negate this test.
c        NBAS_CCF <= 0 will omit baseline subtraction.
c        CPY = smoothed power spectrum with baseline subtracted.
c        Limits hardcoded for (input values) RFWHMST_CCF & NBAS_CCF below.
C        ----------------------------------------------------------------------
         IF (NREFPK(JRFSET).LE.0 .OR. NREFPK(JRFSET).GT.MREFPK) CALL
     1         ERRMES (3, 4, CHSUBP)
         dosub_bas_ccf = .true.
         DO 130 J=1,NREFPK(JRFSET)
            LREFPK(J)=NINT((PPMCEN-PPMREF(J,JRFSET)-
     1                      HZREF(J,JRFSET)/HZPPPM)/PPMINC)+1+nunfil
            dosub_bas_ccf = dosub_bas_ccf   .and.
     1                  abs(PPMREF(J,JRFSET) - 4.68) .gt. ppm_water_tol
  130    CONTINUE
         dosub_bas_ccf = dosub_bas_ccf   .and.   rfwhmst_ccf .gt. 0.
     1                   .and.   nbas_ccf .gt. 0
         if (dosub_bas_ccf) then
            rfwhmst_ccf = amin1(rfwhmst_ccf, 4.)
            nbas_ccf = min0(nbas_ccf, 40)
            lsub_start = max0(1, nint(rfwhmst_ccf * fwhmst / ppminc))
         else
            nbas_ccf = 0
            lsub_start = 0
         end if
         do 132 jdata = 1, ndata
            bas_sum = 0.
            do 133 jsub = lsub_start, lsub_start + nbas_ccf - 1
               bas_sum = bas_sum +
     1                   cpy2(icycle_r(jdata + jsub, ndata)) +
     2                   cpy2(icycle_r(jdata - jsub, ndata))
 133        continue
            if (dosub_bas_ccf) then
               cpy(jdata) = amax1(0., cpy2(jdata) -
     1                                bas_sum / float(2 * nbas_ccf))
            else
               cpy(jdata) = cpy2(jdata)
            end if
 132     continue
C        ----------------------------------------------------------------------
C        ISHFST = shift with max. correlation with the NREFPK(JRFSET)
C                 reference delta-functions.
C        ----------------------------------------------------------------------
  140    CORRMX=-RRANGE
         CORRMN=RRANGE
         DO 150 JSHIFT=LSHFMN(JRFSET),LSHFMX(JRFSET)
            CCF=0.
            DO 160 JREFPK=1,NREFPK(JRFSET)
               CCF=CCF+CPY(ICYCLE_r(LREFPK(JREFPK)+JSHIFT, NDATA))
  160       CONTINUE
            CORRMN=AMIN1(CORRMN, CCF)
            IF (CCF .GE. CORRMX) THEN
               CORRMX=CCF
               ISHFST(JRFSET)=JSHIFT
            END IF
  150    CONTINUE
c        ---------------------------------------------------------------------
c        Only make this test if original SHIFM* are being used, not if
c          restricted SHIFM* from multi-voxel priors for SHIFM*.
c        ---------------------------------------------------------------------
         IF ((ISHFST(JRFSET).EQ.LSHFMN(JRFSET) .OR.
     1        ISHFST(JRFSET).EQ.LSHFMX(JRFSET))   .and.
     2       lshfmn(jrfset) .eq. lshfmn_orig(jrfset)   .and.
     2       lshfmx(jrfset) .eq. lshfmx_orig(jrfset)) THEN
            IF (JRFSET .EQ. 1) THEN
C              ----------------------------------------------------------------
C              An extreme shift SHIFM*(1) is being called for.  The water peak
C                is probably very distorted or weak.
C              ----------------------------------------------------------------
               CALL ERRMES (4, 2, CHSUBP)
            ELSE
               IF (ISHFST(JRFSET) .EQ. LSHFMX(JRFSET)) THEN
C                 -------------------------------------------------------------
C                 The extreme shift SHIFMX(2) is being called for.
C                 -------------------------------------------------------------
                  CALL ERRMES (5, 2, CHSUBP)
               ELSE
                  IF (NREFPK(JRFSET) .GT. MAX0(1,NRF2MN)) THEN
C                    ----------------------------------------------------------
C                    The CCF has probably been distorted by incomplete water
C                      suppression.  Decrement NREFPK(2).  This assumes that
C                      the elements of PPMREF(*,2) are increasing in ppm.
C                    ----------------------------------------------------------
                     CALL ERRMES (6, 1, CHSUBP)
                     NREFPK(2)=NREFPK(2)-1
                     GO TO 140
                  ELSE
C                    ----------------------------------------------------------
C                    All but NRF2MN peaks have been removed, and the extreme
C                      SHIFMN(2) is still being called for.
C                    ----------------------------------------------------------
                     CALL ERRMES (7, 2, CHSUBP)
                  END IF
               END IF
            END IF
         END IF
C        ----------------------------------------------------------------------
C        FWHMCC = FWHM of CCF peak.  It will be used to set limits on the grid
C                 search for the shift.
C        CORRMN is subtracted from CCF as a background.
C        ----------------------------------------------------------------------
         HALFMX=.5*(CORRMX-CORRMN)
         DO 170 JSHIFT=ISHFST(JRFSET)+1,ISHFST(JRFSET)+NUNFIL/2
            CCF=-CORRMN
            DO 180 JREFPK=1,NREFPK(JRFSET)
               CCF=CCF+CPY(ICYCLE_R(LREFPK(JREFPK)+JSHIFT, NDATA))
  180       CONTINUE
            IF (CCF .LE. HALFMX) GO TO 190
  170    CONTINUE
         CALL ERRMES (8, 2, CHSUBP)
  190    DO 200 KSHIFT=ISHFST(JRFSET)-1, ISHFST(JRFSET)-NUNFIL/2, -1
            CCF=-CORRMN
            DO 210 JREFPK=1,NREFPK(JRFSET)
               CCF=CCF+CPY(ICYCLE_R(LREFPK(JREFPK)+KSHIFT, NDATA))
  210       CONTINUE
            IF (CCF .LE. HALFMX) GO TO 220
  200    CONTINUE
         CALL ERRMES (9, 2, CHSUBP)
C        ----------------------------------------------------------------------
C        The 8(ln2)*SDSMOO**2 roughly corrects FWHMCC for the convolution with
C          the Gaussian.  This tends to be an undercorrection; so, too large
C          an SDSMOO would cause too wide a grid and too long computation
c          times.
C          SDSMOO about 0.02 seems reasonable.
C        ----------------------------------------------------------------------
  220    FWHMCC=SQRT(AMAX1(0.,
     1                     (PPMINC*FLOAT(JSHIFT-KSHIFT))**2-
     2                     8.*ALOG2*SDSMOO(JRFSET)**2))
C     -------------------------------------------------------------------------
C        The starting grid will normally not extend beyond RFWHCC*FWHMCC on
C          either side of grid center.  Only a multimodal water peak can
C          extend this.
C     -------------------------------------------------------------------------
         LCC=NINT(RFWHCC*FWHMCC/PPMINC)
         LH2OMN=LCC
         LH2OMX=LCC
         IF (NREFPK(JRFSET).EQ.1 .AND.
     1       ABS(PPMREF(1,JRFSET) - 4.7) .LT. .1) THEN
C           -------------------------------------------------------------------
C           The residual water peak may be multimodal.  Search over the ppm
C             range FWHH2O (centered at LH2OPK) to see if the power spectrum
C             attains values of at least THRESH=HALFMX*2*FH2OMX; i.e.,
C             FH2OMX=.25 will set the threshold at the quarter-maximum (a
C             reasonable value).
C           LH2OMN = no. of points toward the left where the power spectrum
C                    exceeds the threshold.
C           LH2OMX = no. of points toward the right ...
C           -------------------------------------------------------------------
            LH2OPK=LREFPK(1)+ISHFST(JRFSET)
            THRESH=HALFMX*2.*FH2OMX+CORRMN
            DO 230 JSHIFT=1,NINT(.5*FWHH2O/PPMINC)
               IF (CPY(ICYCLE_R(LH2OPK-JSHIFT,NDATA)) .GT. THRESH)
     1               LH2OMN=JSHIFT
               IF (CPY(ICYCLE_R(LH2OPK+JSHIFT,NDATA)) .GT. THRESH)
     1               LH2OMX=JSHIFT
  230       CONTINUE
         END IF
         LSHFMN(3)=MIN0(LSHFMN(3),
     1                  MAX0(LSHFMN(JRFSET), ISHFST(JRFSET)-LH2OMN))
         LSHFMX(3)=MAX0(LSHFMX(3),
     1                  MIN0(LSHFMX(JRFSET), ISHFST(JRFSET)+LH2OMX))
 5230    FORMAT (' Peak in CCF at', I4, 5X, 'Range =', 2I5)
         IF (LPRINT .GT. 0) WRITE (LPRINT,5230) ISHFST(JRFSET),
     1                                          -LH2OMN, LH2OMX
         LWIDTH(JRFSET)=LH2OMX+LH2OMN
  105 CONTINUE
      IF (LWIDTH(1) .LT. LWIDTH(2)) THEN
         LGRID=1
      ELSE
         LGRID=2
      END IF
      LGRDMN=max0(1, NINT(RFWHST*FWHMST/PPMINC))
 300  if (fixshf) then
         lshfbs_fixshf = lshfmn(3)
         lgrid = 1
         ishfst(1) = lshfmn(3)
         lgrdmn = 9999
      end if
C     -------------------------------------------------------------------------
C     Get DATAF = frequency-domain data shifted by ISHFST(LGRID).
C     -------------------------------------------------------------------------
      CALL FTDATA (ISHFST(LGRID))
C     -------------------------------------------------------------------------
C     Get starting estimates for the 0- and 1st-order phase corrections.
c     First restore current values of DEGPPM & DEGZER if in 2nd pass, since
c        they are zeroed at the end of STARTV and are needed in PHASTA to set
c        EXDEGP & EXDEGZ.
C     -------------------------------------------------------------------------
      if (ipass .eq. 1) then
         degzer_sav_startv = degzer
         degppm_sav_startv = degppm
      else
         degzer = degzer_sav_startv
         degppm = degppm_sav_startv
      end if
      CALL PHASTA ()
C
C     -------------------------------------------------------------------------
C     Set up for calls to PLINLS and search for starting values.
C     -------------------------------------------------------------------------
      NSIDES=0
      CALL SETUP (1)
      IF (NSHIFT .GE. MSHIFT) THEN
          NSHIFT=MSHIFT
          CALL ERRMES (10, 3, CHSUBP)
      END IF
      NDEGZ(2)=MAX0(NDEGZ(2),1)
      DDEGZ=360./FLOAT(NDEGZ(2))
      IF (NDGPPM(2) .LE. 1) THEN
         NDGPPM(2)=1
         DDGPPM=DGPPMX-DGPPMN
      ELSE
         DDGPPM=(DGPPMX-DGPPMN)/FLOAT(NDGPPM(2)-1)
         IF (DDGPPM .LE. 0.) CALL ERRMES (11, 4, CHSUBP)
      END IF
      DGPPST=DEGPPM
C     -------------------------------------------------------------------------
C     PARBES(*,1) will be updated to the best so far in the grid search below.
C     PARBES(*,2) must be initialized here because it overwrites PARNLN when
C                 REPHAS is called in the grid search below, and PARNLN must
C                 be properly initialized for PLINLS.
C     -------------------------------------------------------------------------
      PARBES(LSHIST,2)=PARNLN(LSHIST)
      PARBES(LRT2ST,2)=PARNLN(LRT2ST)
      SDBEST(1)=DRANGE
 5310 FORMAT (//' Preliminary search for starting phases ',
     1              'and referencing shift with starting shifts ',
     2              'between', I5, ' and', I5, ' and mesh of', I5,
     3              ' points')
      IF (LPRINT .GT. 0) WRITE (LPRINT,5310) LSHFMN(3), LSHFMX(3),
     1                                       LGRDMN
      CNVRG1=.FALSE.
C     -------------------------------------------------------------------------
C     Start of main loop for search for starting values.
C     Original CY (with no rephasing) is used.
C     -------------------------------------------------------------------------
      DO 310 JSHIFT=1,NSHIFT
         SDSHBS=DRANGE
         KSHFBS=ISHIFD
         DO 320 JDEGZ=1,NDEGZ(2)
            TEST=AMOD(ABS(DEGZER-EXDEGZ),360.)
            IF (AMIN1(TEST,360.-TEST).GT.4.*SDDEGZ .AND. .NOT.FXDEGZ)
     1            GO TO 328
            DEGPPM=DGPPST
            DO 330 JDGPPM=1,NDGPPM(2)
               IF (ABS(DEGPPM-EXDEGP).GT.4.*SDDEGP .AND. .NOT.FXDEGP)
     1               GO TO 338
               IF (LPRINT .GT. 0) WRITE (LPRINT,5330) ISHIFD
 5330          FORMAT (////' Starting shift =', I5, ' points')
               PARNLN(LPHAST)=DEGZER*RADIAN
               PARNLN(LPHAST+1)=DEGPPM*RADIAN
               PARNLN(LSHIST)=0.D0
C              ----------------------------------------------------------------
C              Here the peak broadening is done by multiplying the time-domain
C                data by exp{-[PARNLN(LRT2ST)*t]**2}, t in s, and
C                therefore PARNLN(LRT2ST) in radians/s.
C              FWHMST is in ppm.
C              ----------------------------------------------------------------
               IF (FWHMST .LE. 0.) CALL ERRMES (12, 4, CHSUBP)
               PARNLN(LRT2ST)=FWHMST/TOFWHM
C               LTEST=.TRUE.!CT1
               CALL PLINLS (1, IERROR)
               CNVRG1=CNVRG1 .OR. IERROR.EQ.1
C!CT1
C               PARNLN(LPHAST)=DEGZER*RADIAN!CT1
C               PARNLN(LPHAST+1)=DEGPPM*RADIAN!CT1
C               PARNLN(LSHIST)=0.D0!CT1
C               PARNLN(LRT2ST)=FWHMST/TOFWHM!CT1
C               DO 910 JY=1,NY!CT1
C                  CY(JY)=CYFIT(JY)!CT1
C  910          CONTINUE!CT1
C               CALL PLINLS (1, IERROR)!CT1
C               IF (LTEST) STOP!CT1
C!CT1
               IF (STDDEV .LT. SDSHBS) THEN
                  SDSHBS=STDDEV
C                 -------------------------------------------------------------
C                 Subtraction in the next line is necessary, because ISHIFD
C                   shifts the data and PARNLN(LSHIST) shifts the model.
C                 -------------------------------------------------------------
                  KSHFBS=ISHIFD-NINT(SNGL(PARNLN(LSHIST))*DELTAT*FNDATA/
     1                                 (2.*PI))
                  IF (LPRINT .GT. 0) WRITE (LPRINT,5331) KSHFBS,
     1                                 FLOAT(KSHFBS)*PPMINC
 5331             FORMAT (/' ******** Best data shift for this ',
     1                     'starting shift so far = ', I5,
     2                     ' points =', 1PE10.2, ' ppm')
                  IF (STDDEV .LT. SDBEST(1)) THEN
                     CALL SAVBES (1)
                     LSHFBS=KSHFBS
                     IF (LPRINT .GT. 0) WRITE (LPRINT,5332) LSHFBS,
     1                                    FLOAT(LSHFBS)*PPMINC
 5332                FORMAT (/' ************************ Best data ',
     1                        'shift for all starting shifts so far = ',
     2                        I5, ' points =', 1PE10.2, ' ppm')
                  END IF
               END IF
  338          TEST=DEGPPM+DDGPPM
               IF (TEST .LE. DGPPMX) THEN
                  DEGPPM=TEST
               ELSE
C                 -------------------------------------------------------------
C                 New DEGPPM would exceed DGPPMX.  Set DEGPPM to an endpoint;
C                   which endpoint is determined by how far TEST overshot
C                   DGPPMX.
C                 -------------------------------------------------------------
                  IF (TEST-DGPPMX .LE. .5*DDGPPM) THEN
                     DEGPPM=DGPPMX
                  ELSE
                     DEGPPM=DGPPMN
                  END IF
               END IF
  330       CONTINUE
  328       DEGZER=AMOD(DEGZER+DDEGZ,360.)
  320    CONTINUE
         INTSHF(JSHIFT,1)=MIN0(KSHFBS,ISHIFD)
         INTSHF(JSHIFT,2)=MAX0(KSHFBS,ISHIFD)
C        ----------------------------------------------------------------------
C        Use grid point in grid range that is farthest from the intervals
C          already covered.  If this max. distance is less than LGRDMN, then
C          exit.
C        ----------------------------------------------------------------------
         LMAX=LGRDMN-1
         DO 340 L=LSHFMN(3),LSHFMX(3)
            LMIN=NDATA
            DO 350 KSHIFT=1,JSHIFT
               LMIN=MIN0(LMIN, MAX0(INTSHF(KSHIFT,1)-L,
     1                              L-INTSHF(KSHIFT,2)))
  350       CONTINUE
            IF (LMIN .GT. LMAX) THEN
               LMAX=LMIN
               ISHFUS=L
            END IF
  340    CONTINUE
         IF (LMAX .LT. LGRDMN) GO TO 400
         CALL SHIFTD (ISHFUS)
  310 CONTINUE
      CALL ERRMES (13, 1, CHSUBP)
  400 IF (SDBEST(1) .GE. DRANGE) CALL ERRMES (14, 4, CHSUBP)
C     -------------------------------------------------------------------------
C     Load best values into PARNLN.
C     -------------------------------------------------------------------------
      CALL SAVBES (2)
      IF (LSHFBS .NE. ISHIFD) THEN
         CALL SHIFTD (LSHFBS)
         CALL REPHAS ()
C        ----------------------------------------------------------------------
C        Analysis with optimal shift to get best starting values with
C          optimally shifted and rephased data.
C        ----------------------------------------------------------------------
         PARNLN(LSHIST)=0.D0
 5400    FORMAT (////' Analysis with optimal data shift to get optimal',
     1               ' starting values')
         IF (LPRINT .GT. 0) WRITE (LPRINT,5400)
         MITER(1)=MITER(1)*2
         CALL PLINLS (1, IERROR)
         CNVRG1=CNVRG1 .OR. IERROR.EQ.1
         IF (OBJECT .GE. DRANGE) CALL ERRMES (15, 4, CHSUBP)
         CALL SAVBES (1)
         CALL SAVBES (2)
         LSHFBS=LSHFBS-NINT(SNGL(PARNLN(LSHIST))*DELTAT*FNDATA/(2.*PI))
      END IF
      IF (.NOT.CNVRG1) CALL ERRMES (16, 2, CHSUBP)
      FWHMST=AMAX1(0., TOFWHM*SNGL(PARBES(LRT2ST,1)))
      fwhmst_full = fwhmst
      IF (.not.DOFULL) THEN
c        ----------------------------------------------------------------------
c        Set ISHIFD=LSHFBS, just for final output as "Data shift" (without
c          shifting the data).
c        ----------------------------------------------------------------------
         ishifd = lshfbs
      else
         if (fwhmst .gt. fwhmmx) then
            fwhmst = fwhmmx
c           ------------------------------------------------------------------
c           With IMETHD=3, FWHMMX is set fairly low.  It determines the max
c              starting FWHM of the (e.g., Voigt) power functions.
c           Similarly, FWHMST determines the SIFWMN when SCAFWH=T
c           ------------------------------------------------------------------
            if (imethd .ne. 3   .and.   .not.scafwh)
     1         call errmes (25, 2, chsubp)
         end if
         IF (LSHFBS .NE. ISHIFD) THEN
            PARBES(LPHAST,2)=PARBES(LPHAST,2)+PHITOT(1)
            PARBES(LPHAST+1,2)=PARBES(LPHAST+1,2)+PHITOT(2)
            CALL SHIFTD (LSHFBS)
         END IF
         if (fixshf   .and.   lshfbs .ne. lshfbs_fixshf) then
            lshfbs = lshfbs_fixshf
            PARBES(LPHAST,2)=PARBES(LPHAST,2)+PHITOT(1)
            PARBES(LPHAST+1,2)=PARBES(LPHAST+1,2)+PHITOT(2)
            call shiftd(lshfbs)
         end if
         CALL REPHAS ()
         DEGZER=0.
         DEGPPM=0.
         if (scafwh   .or.   imethd .eq. 3) then
C           -------------------------------------------------------------------
C           SCAFWH assumes that there is no lineshape in Final
C           -------------------------------------------------------------------
            nsides = 0
            nside2 = 0
            incsid = 1
C           -------------------------------------------------------------------
C           FCONC_EXPECT = scale factor to multiply CONC_EXPECT to get
c                          ALPHAS in constant range (only used with IMETHD=2)
C           -------------------------------------------------------------------
            fconc_expect = 0.
            do 450 jmetab = 1, nmetab
               fconc_expect = fconc_expect + solbes(jmetab, 2)
 450        continue
            if (imethd .eq. 3) fwhmst = amax1(fwhmst, fwhmmn)
         else
C           -------------------------------------------------------------------
C           Compute NSIDES and INCSID, after checking their limits.
C           -------------------------------------------------------------------
            IF (NSIDMX .GT. MSIDES) THEN
               NSIDMX=MSIDES
               CALL ERRMES (17, 3, CHSUBP)
            END IF
            IF (NSIDMN.LE.0 .OR. NSIDMN.GT.NSIDMX) CALL ERRMES (18, 4,
     1                                                          CHSUBP)
            IF (INCSMX.LE.0 .OR. INCSMX.GT.MINCSD) THEN
               CALL ERRMES (19, 3, CHSUBP)
               INCSMX=MIN0(MINCSD,MAX0(1,INCSMX))
            END IF
            NSIDES=NINT(.5*FWHMST*RFWHM/PPMINC)
            IF (NSIDES .LT. NSIDMN) THEN
               CALL ERRMES (20, 1, CHSUBP)
               NSIDES=NSIDMN
            END IF
            INCSID=(NSIDES-1)/NSIDMX +1
            NSIDES=NINT((FLOAT(NSIDES)+.001)/FLOAT(INCSID))
            IF (INCSID .GT. INCSMX) THEN
               INCSID=INCSMX
               NSIDES=NSIDMX
               if (.not.nobasi) CALL ERRMES (21, 3, CHSUBP)
            END IF
C           -------------------------------------------------------------------
C           The data with ppm < PPMCEN run from NUNFIL+1 to NDATA
C           -------------------------------------------------------------------
            INCDIM=MIN0(LDATST-1,(ndata-LDATEN-1))/NSIDES
            IF (INCDIM .LE. 0) CALL ERRMES (22, 4, CHSUBP)
            IF (INCSID .GT. INCDIM) THEN
               INCSID=INCDIM
               CALL ERRMES (23, 3, CHSUBP)
            END IF
            IF (INCSID .GT. 4) then
               CALL ERRMES (24, 2, CHSUBP)
            else if (INCSID .GT. 1) then
               CALL ERRMES (24, 1, CHSUBP)
            end if
         end if
         IF (LPRINT .GT. 0) WRITE (LPRINT,5410) LSHFBS,
     1         FLOAT(LSHFBS)*PPMINC, PHITOT(1)/RADIAN,
     2         PHITOT(2)/RADIAN, fwhmst_full, NSIDES, INCSID
 5410     FORMAT (////' Starting values for final analysis'/
     1                ' Shift =', I5, ' points =', F7.4, ' ppm'/
     2                ' Phi(0) =', F8.1, ' deg'/
     3                ' Phi(1) =', F8.2, ' deg/ppm'/
     4                ' FWHM =', F6.3, ' ppm'/
     5                ' NSIDES =', I3,/
     6                ' INCSID =', I2///)
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FTDATA (ISHIFT)
C
C  DATAF(JDATA) = unshifted frequency-domain data (finally) (COMPLEX).
C  DATAF(1+nunfil): corresponds to PPMCEN; i.e., the spectrum is rearranged.
C  CY(JY) = shifted (by ISHIFT) frequency-domain data in analysis window.
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='FTDATA'
C     -------------------------------------------------------------------------
C     DATAF = zero-filled time-domain data (temporarily).
C           = unshifted frequency-domain data (finally).
C     -------------------------------------------------------------------------
      DO 120 JUNFIL=1,NUNFIL
         DATAF(JUNFIL)=DATAT(JUNFIL)
         DATAF(JUNFIL+NUNFIL)=(0.,0.)
  120 CONTINUE
      CALL CFFT_r (DATAF, DATAF, NDATA, LWFFT,
     1           WFFTC)
C     -------------------------------------------------------------------------
C     CY(JY) = frequency-domain data in window to be analyzed (COMPLEX).
C     -------------------------------------------------------------------------
      CALL SHIFTD (ISHIFT)
C     -------------------------------------------------------------------------
C     RMSAMP will be used in GBACKG to scale background, so that the columns
C       of the Hessian will be of similar magnitude for better performance in
C       PNNLS.
C     -------------------------------------------------------------------------
      RMSAMP=0.
      DO 150 JY=1,NY
         if (.not.lcy_skip(jy))
     1      RMSAMP=RMSAMP+REAL(CY(JY))**2+AIMAG(CY(JY))**2
  150 CONTINUE
      RMSAMP=SQRT(RMSAMP/FLOAT(NYuse))
      IF (RMSAMP .LE. 0.) then
         istago = 0
         CALL ERRMES (1, 4, CHSUBP)
      end if
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SHIFTD (ISHIFT)
C
C  Shift CY ISHIFT points from original unshifted position.
C  ISHIFD = ISHIFT = current shift.
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='SHIFTD'
      JDATA=LDATST+ISHIFT
      DO 150 JY=1,NY
         CY(JY)=DATAF(ICYCLE_r(JDATA,NDATA))
         JDATA=JDATA+1
  150 CONTINUE
C     -------------------------------------------------------------------------
C     ISTAGO = 2 when CY is loaded with frequency-domain data, and can
C                therefore be plotted in EXITPS.
C     -------------------------------------------------------------------------
      ISTAGO=2
      IF (PHITOT(1).NE.0. .OR. PHITOT(2).NE.0.) THEN
         PHITOT(1)=0.
         PHITOT(2)=0.
      END IF
      ISHIFD=ISHIFT
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SETUP (LSTAGE)
C
C  Set up variables for calls to PLINLS.
C  LSTAGE = 1 for initial analyses with restricted model.
C         = 2 for analyses with full model.
C
      INCLUDE 'lcmodel.inc'
      LOGICAL LDELTA, LERROR
      CHSUBP='SETUP'
c     -------------------------------------------------------------------------
c     If SCAFWH=T, use lineshape in Prel so that FWHMST can be used to
c        determine # of knots in Final.  All NMETAB in Prel will be with
c        lineshape; so must set CHUSE1 carefully and BADREF=F.  For FWHMST to
c        be useful, only sharpest models should be used in Prel.
c     If SCAFWH=T, FWHMST also determines the FWHM specs of the simulated
c        basis spectra.
c     -------------------------------------------------------------------------
      if ((scafwh   .or.   imethd .eq. 3)   .and.   lstage .eq. 1)
     1   then
         do 120 j = 1, nmetab
            lshape(j) = .true.
 120     continue
      end if
      if (nobase) go to 200
C     -------------------------------------------------------------------------
C     NBACKG is set so that there will be a knot every
C       FWHMST_full*RBACKG(LSTAGE) ppm.
C     -------------------------------------------------------------------------
      if (lstage .eq. 1) fwhmst_full = fwhmst
      TEST=AMAX1(DKNTMN(IKNTMN),fwhmst_full*RBACKG(LSTAGE))
      IF (TEST .LE. 0.) CALL ERRMES (1, 4, CHSUBP)
      NBACKG=NINT((delppm(1) - delppm(ny))/TEST)+3
c     -------------------------------------------------------------------------
c     Each gap can introduce 2 extra (external) knots (although large gaps
c        will reduce the # knots).
c     -------------------------------------------------------------------------
      IF (NBACKG .GT. MBACKG - 2 * ngap) THEN
         CALL ERRMES (2, 3, CHSUBP)
         NBACKG=MBACKG - 2 * ngap
      END IF
C     -------------------------------------------------------------------------
C     Get BACKGR(JY,JBACKG) = background terms (REAL).
c     Compute REGB = (temporarily) Squared baseline regularizor
C     -------------------------------------------------------------------------
      CALL GBACKG ()
C
C      DO 140 IROW=1,NBACKG!
C         WRITE (LPRINT,5150) (REGB(IROW,J),J=1,NBACKG) !DELETE
C 5150    FORMAT (/(1X,1P10E13.5)) !DELETE
C  140 CONTINUE!
C     -------------------------------------------------------------------------
C     Eigenvalue decomposition of squared regularizor (not the fastest way,
C       but simple and stable).  Because of the lack of boundary conditions,
C       REGB is only of rank NBACKG-2.
C     The IMSL routine EVCSF normalizes the eigenvectors in the infinity-norm.
C       This destroys the orthognality of the matrix of eigenvectors.  The
C       IMSL eigenvectors would have to be renormalized in the Euclidean norm.
C       Use EISPACK routines instead.
C     CPY = vector of eigenvalues (temporarily).
C     WFFTC = matrix of eigenvectors (temporarily).
C     -------------------------------------------------------------------------
C      CALL EVCSF (NBACKG, REGB, MBACKG, CPY, WFFTC, MBACKG)!
C      DO 9145 J=1,NBACKG!
C         RNORM=0.!
C         IPOINT=(J-1)*MBACKG!
C         DO 9146 K=1,NBACKG!
C            RNORM=RNORM+WFFTC(IPOINT+K)**2!
C 9146    CONTINUE!
C         RNORM=1./SQRT(RNORM)!
C         DO 9147 K=1,NBACKG!
C            WFFTC(IPOINT+K)=WFFTC(IPOINT+K)*RNORM!
C 9147    CONTINUE!
C 9145 CONTINUE!
      LWFFT=0
      CALL EIGVRS (MBACKG, NBACKG, REGB, CPY, WFFTC, DAMAT, DAMAT(1,2),
     1             IERROR)
      IF (IERROR .NE. 0) CALL ERRMES (4, 4, CHSUBP)
C      WRITE (LPRINT,5150) (CPY(J),J=1,NBACKG)!
C     -------------------------------------------------------------------------
C     REGB = regularizor
C          = square root of above squared regularizor in REGB
C     -------------------------------------------------------------------------
      DO 150 ICOL=1,NBACKG
         DO 160 IROW=1,NBACKG
            REGB(IROW,ICOL)=0.
  160    CONTINUE
  150 CONTINUE
      DO 170 ICOL=1,NBACKG
         IPOINT=(ICOL-1)*MBACKG
         IF (CPY(ICOL) .GT. 0.) THEN
            FACT=SQRT(CPY(ICOL))
            DO 180 JCOL=1,NBACKG
               TERM=WFFTC(IPOINT+JCOL)*FACT
               DO 181 IROW=1,NBACKG
                  REGB(IROW,JCOL)=REGB(IROW,JCOL)+WFFTC(IPOINT+IROW)*
     1                                            TERM
  181          CONTINUE
  180       CONTINUE
         END IF
  170 CONTINUE
C     -------------------------------------------------------------------------
C     Scale REGB for normalized B-splines.
C     The check for a zero denominator below was just done in GBACKG.
C     -------------------------------------------------------------------------
      FNORM=(FLOAT(NBACKG-3)/(DELPPM(1)-DELPPM(NY)))**4
      DO 183 ICOL=1,NBACKG
         DO 184 IROW=1,NBACKG
            REGB(IROW,ICOL)=REGB(IROW,ICOL)*FNORM
  184    CONTINUE
  183 CONTINUE
C      DO 9160 IROW=1,NBACKG !DELETE
C         IF (LPRINT .GT. 0) WRITE (LPRINT,5150) (REGB(IROW,J),
C     1                                           J=1,NBACKG) !
C 9160 CONTINUE !
C      DO 9195 IROW=1,NBACKG!
C         DO 9196 ICOL=1,NBACKG!
C            TERM=0.!
C            DO 9197 J=1,NBACKG!
C               TERM=TERM+REGB(J,IROW)*REGB(J,ICOL)!
C 9197       CONTINUE!
C            CPY(ICOL)=TERM!
C 9196    CONTINUE!
C         WRITE (LPRINT,5150) (CPY(ICOL),ICOL=1,NBACKG)!
C 9195 CONTINUE!
C      IF (NBACKG .GE. 0) STOP!
C     -------------------------------------------------------------------------
C     Initialize arrays for PLINLS.
C     First initialize linear parameters.
C     -------------------------------------------------------------------------
 200  NLIN=NMETAB+NBACKG
      IF (NLIN .GT. MPAR) CALL ERRMES (5, 4, CHSUBP)
      DO 210 JPAR=1,NLIN
         SOLUTN(JPAR)=0.D0
  210 CONTINUE
C     -------------------------------------------------------------------------
C     Initialize lineshape parameters.
C     PARNLN = normalized Gaussian with fwhm=FWHMST.
C     -------------------------------------------------------------------------
      NSIDE2=2*NSIDES
      NNONL=NSIDE2
      IF (NNONL .GT. MNONL) CALL ERRMES (6, 4, CHSUBP)
      RNORM=1.
      LDELTA=FWHMST .LT. .5*PPMINC
      rexp = rrange
      IF (.NOT.LDELTA) REXP=-4.*ALOG2*(FLOAT(INCSID)*PPMINC/FWHMST)**2
      DO 230 JSIDE=1,NSIDES
         IF (LDELTA) THEN
            TERM=0.
         ELSE
            TERM=EXP(REXP*FLOAT(JSIDE)**2)
         END IF
         PARNLN(NSIDES+1-JSIDE)=TERM
         PARNLN(NSIDES+JSIDE)=TERM
         RNORM=RNORM+2.*TERM
  230 CONTINUE
      RNORM=1./RNORM
      DO 250 JNONL=1,NNONL
         PARNLN(JNONL)=PARNLN(JNONL)*RNORM
         DPARMQ(JNONL)=DFLDMQ
  250 CONTINUE
      IF (LSTAGE .EQ. 2   .and.   imethd .ne. 2)
     1    THEN
C        ----------------------------------------------------------------------
C        DGAUSS = a 1-sided normalized Gaussian profile with a fwhm of FWHMBA
C                 ppm.
C        ----------------------------------------------------------------------
         IF (FWHMBA .LE. 0.) CALL ERRMES (7, 4, CHSUBP)
         REXP=-4.*ALOG2*(FLOAT(INCSID)*PPMINC/FWHMBA)**2
         RNORM=1.
         DGAUSS(0)=1.D0
         DO 255 J=1,NSIDE2+4
            TERM=EXP(REXP*FLOAT(J)**2)
            DGAUSS(J)=TERM
            RNORM=RNORM+2.*TERM
  255    CONTINUE
         RNORM=1./RNORM
         DO 256 J=0,NSIDE2+4
            DGAUSS(J)=DGAUSS(J)*RNORM
  256    CONTINUE
      END IF
C     IF (LSTAGE .EQ. 2) THEN!DELETE
C         WRITE (LPRINT,9250) (PARNLN(J),J=1,NNONL)!DELETE
C 9250    FORMAT (//' Starting Gaussian profile'/(1X, 1P10E13.3))!
C         IF (LPRINT .GT. 0) STOP!DELETE
C     END IF!DELETE
C     -------------------------------------------------------------------------
C     Initialize phase parameters.
C     -------------------------------------------------------------------------
      LPHAST=NNONL+1
      NNONL=NNONL+2
      IF (NNONL .GT. MNONL) CALL ERRMES (8, 4, CHSUBP)
      PARNLN(NNONL-1)=DEGZER*RADIAN
      PARNLN(NNONL)=DEGPPM*RADIAN
      DPARMQ(NNONL-1)=DDGZMQ(LSTAGE)*RADIAN
      DPARMQ(NNONL)=DDGPMQ(LSTAGE)*RADIAN
C     -------------------------------------------------------------------------
C     Initialize shift parameters.
C     -------------------------------------------------------------------------
      IF (LSTAGE .EQ. 1) THEN
         ALPHAB=0.D0
         ALPHAS=0.D0
         NEXPON=1
C        ----------------------------------------------------------------------
C        Set (max. step size)/FSTPMQ to RSHFMQ*(initial estimate of FWHM).
C        ----------------------------------------------------------------------
         DMARQ=RSHFMQ*FWHMST
C        ----------------------------------------------------------------------
C        Above DMARQ is in ppm; convert this to radians/s.
C        ----------------------------------------------------------------------
         DMARQ=DMARQ*2.*PI*HZPPPM
         LSHIST=NNONL+1
         NNONL=NNONL+NEXPON
         IF (NNONL .GT. MNONL) CALL ERRMES (11, 4, CHSUBP)
         PARNLN(lshist)=0.D0
         DPARMQ(lshist)=DMARQ
      ELSE
         ALPHAB=ALPBST
         IF (NSIDES.LE.0   .AND.   ALPSST.gt.0.   .and.   imethd .ne. 2)
     1       ALPSST=0.
         ALPHAS=ALPSST
         IF (ALPHAB.LE.0.D0 .OR. ALPHAS.LT.0.D0) CALL ERRMES (10, 4,
     1                                                        CHSUBP)
         NEXPON=NMETAB
         LSHIST=NNONL+1
         NNONL=NNONL+NEXPON
         IF (NNONL .GT. MNONL) CALL ERRMES (11, 4, CHSUBP)
c        ----------------------------------------------------------------------
c        Set (max. step size)/FSTPMQ to RSDSMQ * SDSHIF
c        ----------------------------------------------------------------------
         jmetab = 0
         DO 260 JNONL=LSHIST,NNONL
            PARNLN(JNONL)=0.D0
            jmetab = jmetab + 1
            DPARMQ(JNONL)= rsdsmq * sdshif(jmetab)
 260     CONTINUE
      END IF
C     -------------------------------------------------------------------------
C     Initialize broadening parameters.
C     -------------------------------------------------------------------------
      LRT2ST=NNONL+1
      NNONL=NNONL+NEXPON
      IF (NNONL .GT. MNONL) CALL ERRMES (12, 4, CHSUBP)
      IF (NSIDES .GT. 0   .or.   imethd .eq. 2) THEN
         JMETAB=0
         DO 270 JNONL=LRT2ST,NNONL
            JMETAB=JMETAB+1
            PARNLN(JNONL)=EXRT2(JMETAB)
 
C            write(6, 9270) jnonl, parnln(jnonl) * tofwhm!<<< = 4, .04995
c 9270       format('FWHM(', i2, ') =', 1pe10.3)
 
            DPARMQ(JNONL)=rrt2mq * EXRT2(JMETAB)
  270    CONTINUE
      ELSE
C        ----------------------------------------------------------------------
C        NSIDES = 0. Use FWHMST for initial estimate.  In the usual case that
C          LSTAGE=1, FWHMST is input data. When LSTAGE=2, it is from the
C          initial analysis with LSTAGE=1.
C        ----------------------------------------------------------------------
         IF (LSTAGE .EQ. 1) THEN
C           -------------------------------------------------------------------
C           Here the peak broadening is done by multiplying the time-domain
C             data by exp{-[PARNLN(LRT2ST)*t]**2}, t in s, and therefore
C             PARNLN(LRT2ST) in radians/s
C           FWHMST is in ppm.
C           -------------------------------------------------------------------
            IF (FWHMST .LE. 0.) CALL ERRMES (13, 4, CHSUBP)
            PARNLN(LRT2ST)=FWHMST/TOFWHM
         ELSE
            if (imethd .eq. 3) then
               call setup3 ()
            else
C              ----------------------------------------------------------------
C              Here the peak broadening will be done by multiplying the
C                time-domain data by exp(-PARNLN*t), t in s and therefore
C                PARNLN in 1/s.
C              A Lorenzian has 1/T2 = PI*fwhm.
C              FWHMST is in ppm.
C              ----------------------------------------------------------------
               PARNLN(LRT2ST)=PI*FWHMST*HZPPPM
            end if
         END IF
         DPARMQ(LRT2ST)=RRT2MQ*PARNLN(LRT2ST)
         if (imethd .ne. 3) then
C           -------------------------------------------------------------------
C           Set (max. step size)/FSTPMQ to RRT2MQ*(initial estimate).
C           -------------------------------------------------------------------
            DO 275 JNONL=LRT2ST+1,NNONL
               PARNLN(JNONL)=PARNLN(LRT2ST)
               DPARMQ(JNONL)=DPARMQ(LRT2ST)
 275        CONTINUE
         end if
      END IF
      DO 277 JNONL=LRT2ST,NNONL
         IF (PARNLN(JNONL) .LE. 0.D0) CALL ERRMES (14, 4, CHSUBP)
  277 CONTINUE
      NPAR=NLIN+NNONL
      IF (NPAR .GT. MPAR) CALL ERRMES (15, 4, CHSUBP)
C     -------------------------------------------------------------------------
C     Set NONNEG, except for
c       JPAR=1,NMETAB, which were set in MYBASI
c       JPAR=NMETAB+1,NMETAB+NBACKG, which must be set in GBACKG.
C     NONNEG for delta(1/T2) will be reset in PASTEP to constrain delta(1/T2)
C       to be nonnegative.
C     -------------------------------------------------------------------------
      DO 290 JPAR=NLIN+1,NPAR
         NONNEG(JPAR)=.FALSE.
  290 CONTINUE
      IF (NSIDES .GT. 0) THEN
C        ----------------------------------------------------------------------
C        Set up 2nd-order regularizor for lineshape coefficients with zero
C          boundary conditions.
C          The rows start from the center and work outward.  The 1st 3 rows
C            correspond to the equality constraint that has eliminated the
C            center point as a parameter; these 1st 3 also have extra
C            constants (2, -1, -1) on the rhs of the regularizor.
C        ----------------------------------------------------------------------
         NCOLRF=NSIDE2
         NREGF=NSIDE2+3
         DO 310 IROW=4,NREGF
            DO 320 ICOL=1,NCOLRF
               REGF(IROW,ICOL)=0.D0
  320       CONTINUE
  310    CONTINUE
         DO 330 ICOL=1,NCOLRF
            REGF(1,ICOL)=2.D0
            REGF(2,ICOL)=-1.D0
            REGF(3,ICOL)=-1.D0
  330    CONTINUE
         REGF(1,NSIDES)=3.D0
         REGF(1,NSIDES+1)=3.D0
         IF (NSIDES .GE. 2) REGF(2,NSIDES-1)=0.D0
         REGF(2,NSIDES)=-3.D0
         REGF(3,NSIDES+1)=-3.D0
         IF (NSIDES .GE. 2) REGF(3,NSIDES+2)=0.D0
         IROW=3
         DO 340 JCENTR=2,NSIDES+1
            IROW=IROW+1
            ICOL=NSIDES+JCENTR
            REGF(IROW,ICOL-1)=1.D0
            IF (JCENTR .LE. NSIDES) REGF(IROW,ICOL)=-2.D0
            IF (JCENTR .LT. NSIDES) REGF(IROW,ICOL+1)=1.D0
            IROW=IROW+1
            ICOL=NSIDES+1-JCENTR
            REGF(IROW,ICOL+1)=1.D0
            IF (ICOL .GE. 1) REGF(IROW,ICOL)=-2.D0
            IF (ICOL .GT. 1) REGF(IROW,ICOL-1)=1.D0
  340    CONTINUE
      END IF
C     -------------------------------------------------------------------------
C     First initialize SDREF with PMQACT=0 to avoid possible stalling later in
C       PLINLS.
C     -------------------------------------------------------------------------
      INISOL=.TRUE.
      SDREF=DRANGE
      CALL SOLVE (LSTAGE, .FALSE., 0.D0, .FALSE., LERROR)
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE setup3 ()
c
c Called when IMETHD=3 and LSTAGE=2.
c Set starting values for PARNLN & DPARMQ
c
      INCLUDE 'lcmodel.inc'
      double precision ddtime, dtime
      real pstart(mmpowr)
      chsubp = 'SETUP3'
      start1 = alog(pi * hzpppm * fwhmst)
      start2 = 2. * alog(fwhmst / tofwhm)
      ddtime = dble(deltat)
      if (mpower .gt. mmpowr   .or.   mpower .lt. 1   .or.
     1    amin1(fmain_power, fother_power) .le. -.001   .or.
     2    amax1(fmain_power, fother_power) .ge. 1.001)
     3   call errmes (1, 4, chsubp)
      do 110 jpower = 1, mpower
         if (power(jpower) .lt. .999d0   .or.
     1       power(jpower) .gt. 2.001d0) call errmes (2, 4, chsubp)
         dtime = 0.d0
         do 120 jdata = 1, ndata
            tpower(jpower, jdata) = dtime ** power(jpower)
            dtime = dtime + ddtime
 120     continue
c        ---------------------------------------------------------------------
c        Geometric mean
c        ---------------------------------------------------------------------
         term = power(jpower)
         pstart(jpower) = exp((2. - term) * start1 +
     1                        (term - 1.) * start2)
 110  continue
      ncoeff_power = 0
c     ------------------------------------------------------------------------
c     NNONL was already incremented by NMETAB in SETUP.  Set it back and then
c        increment it here.
c     ------------------------------------------------------------------------
      nnonl = lrt2st - 1
      lpowen(0) = nnonl
      do 210 jmetab = 1, nmetab
         ncoeff_power = ncoeff_power + 1
         nnonl = nnonl + 1
c        ---------------------------------------------------------------------
c        FMAIN_POWER = attenuation factor for starting values of PARNLN for
c                      1st broadening term (typically 0.6)
c        FOTHER_POWER = attenuation factor for starting values of PARNLN for
c                       other broadening terms (typically 0.1)
c        ---------------------------------------------------------------------
         parnln(nnonl) = fmain_power * pstart(1)
         dparmq(nnonl) = rpowmq * pstart(1)
c        ----------------------------------------------------------------------
c        COEFF_POWER_SD = prior SD of power coefficient (analogous to SDRT2).
c        All prior means are zero, to suppress excessive broadening.  Typical
c        values for FRACT_POWER_SD might be 0.5 for component 1 and somewhat
c        less for others.
c        ----------------------------------------------------------------------
         coeff_power_sd(1, jmetab) = fract_power_sd(1, jmetab) *
     1                               pstart(1)
         do 230 jpower = 2, npower(jmetab)
            ncoeff_power = ncoeff_power + 1
            nnonl = nnonl + 1
            if (nnonl .gt. mnonl   .or.
     1        ncoeff_power .gt. mcoeff_power) call errmes (3, 4, chsubp)
            parnln(nnonl) = fother_power * pstart(jpower)
            dparmq(nnonl) = rpowmq * pstart(jpower)
            coeff_power_sd(jpower, jmetab) =
     1         fract_power_sd(jpower, jmetab) * pstart(jpower)
 230     continue
         lpowen(jmetab) = nnonl
 210  continue
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PHASTA ()
C
C  Get starting estimates for the 0- and 1st-order phase corrections for the
C    data.
C  These can be specified in 5 ways:
C    (1) NDGPPM(1) > 1 will cause a grid search over NDGPPM(1) points in
C                      DEGPPM=[DGPPMN,DGPPMX] and NDEGZ(1) points in
c                      DEGZER=[0,180] for the minimum of the
c                      |vertical distance|**IPOWPH
c                      traveled by the real part of the smoothed spectrum.
c                      With IPOWPH about 6, the big jumps in the center of the
c                      dispersion part will count very strongly and help to
c                      eliminate these.
c                      This criterion will be
c                      weakened by noise and also if the residual water peak
c                      is incoherent (as in a strongly water-suppressed
c                      spectrum).  However, final results do not seem to be
c                      significantly affect by PHASTA.
c                      Until 6.2-0G, IPOWPH=1.
c                      Until 6.1-4K, the
c                      sum of the squares of the integrals over the real and
c                      imaginary parts was maximized. This now seems like
c                      nonsense (especially if the water tails intrudes).
c                      Also, this criterion would be identical for 90, 180 &
c                      270-deg phase changes.
C    (2) NDGPPM(1) = 1 will set DEGPPM=DEGPPM with DEGZER set as in (1).
C    (3) NDGPPM(1) < 1 will cause input values of DEGZER or DEGPPM to be used.
C    (4) SDDEGZ<45 will cause the input value of DEGZER to be used;
C        SDDEGZ < min[DDGZMQ(1),DDGZMQ(2)]/3 will cause DEGZER to be fixed to
C                 its input value throughout the analysis.
C    (5) SDDEGP < 10 will cause the input value of DEGPPM to be used;
C        SDDEGP < min[DDGPMQ(1),DDGPMQ(2)]*.4 will cause DEGPPM to be fixed to
C                 its input value throughout the analysis.
C  The criteria in (1) and (2) can be poor if incomplete water suppression
C    causes a significant part of the phased real part of the spectrum to be
C    negative.
c
      INCLUDE 'lcmodel.inc'
      complex data_ph1(my), data_zero(mdata)
      double precision dist, distmn
      LOGICAL STARTP, STARTZ
      CHSUBP='PHASTA'
      SDDEGP=AMAX1(SDDEGP,0.)
      SDDEGZ=AMAX1(SDDEGZ,0.)
C     -------------------------------------------------------------------------
C     START* = T means that DEGZER or DEGPPM will be used as starting values.
C     FXDEG* = T means that DEGZER or DEGPPM will be fixed at their input
C                values throughout the analysis.
C     -------------------------------------------------------------------------
      FXDEGZ=SDDEGZ .LT. AMIN1(DDGZMQ(1),DDGZMQ(2))/3.
      STARTZ=SDDEGZ.LT.45. .OR. FXDEGZ .OR. NDGPPM(1).LE.0
      IF (FXDEGZ) NDEGZ(2)=1
      FXDEGP=SDDEGP .LT. AMIN1(DDGPMQ(1),DDGPMQ(2))*.4
      STARTP=SDDEGP.LT.10. .OR. NDGPPM(1).LE.0 .OR. FXDEGP
      IF (FXDEGP) NDGPPM(2)=1
C     -------------------------------------------------------------------------
C     DEGZER and DEGPPM are the input expectation values; save these in EXDEGZ
C       and EXDEGP.
C     -------------------------------------------------------------------------
      DEGZER=AMOD(DEGZER,360.)
      EXDEGZ=DEGZER
      EXDEGP=DEGPPM
C     -------------------------------------------------------------------------
C     Put EXDEGP in range of DGPPM*.
C     -------------------------------------------------------------------------
      if (exdegp .lt. dgppmn) then
         exdegp = dgppmn
         call errmes (1, 3, chsubp)
      else if (exdegp .gt. dgppmx) then
         exdegp = dgppmx
         call errmes (1, 3, chsubp)
      end if
      IF (.NOT.STARTZ .OR. .NOT.STARTP) THEN
C        ----------------------------------------------------------------------
C        Grid search for starting values for DEGZER and/or DEGPPM.
C        ----------------------------------------------------------------------
         IF (STARTP) NDGPPM(1)=1
         IF (NDGPPM(1) .GT. 1) THEN
            DELTPP=(DGPPMX-DGPPMN)/FLOAT(NDGPPM(1)-1)
            DEGPPM=DGPPMN-DELTPP
         ELSE
            DELTPP=0.
         END IF
         if (startz) ndegz(1) = 1
         if (ndegz(1) .gt. 1) then
            deltz = 180. / float(ndegz(1))
         else
            deltz = 0.
         end if
C        ----------------------------------------------------------------------
C        DATAF_ZERO = (temporarily) Gaussian-smoothed FFT of zero-filled
C                     (time-domain) DATAT (COMPLEX).
C        DATA_ZERO(1+nunfil): corresponds to PPMCEN; i.e., the spectrum is
c                         rearranged.
c        Use SDSMOO(3) for SD of smoothing.
C        ----------------------------------------------------------------------
         RSD=2.*PI*SDSMOO(3)/(PPMINC*FNDATA)
         DO 110 JUNFIL=1,NUNFIL
            DATA_ZERO(JUNFIL)=DATAT(JUNFIL)*
     1                    EXP(-.5*(RSD*FLOAT(JUNFIL-1))**2)
            DATA_ZERO(NUNFIL+JUNFIL)=(0.,0.)
  110    CONTINUE
         CALL CFFT_r (DATA_ZERO, DATA_ZERO, NDATA, LWFFT,
     1              WFFTC)
C        ----------------------------------------------------------------------
c        Start of main loop for phase optimization.
C        ----------------------------------------------------------------------
         distmn = 1.0d300
         sumbes = 9.
         degppm_best = 999.
         degzer_best = 999.
         do 210 jdgppm = 1, ndgppm(1)
            DEGPPM=DEGPPM+DELTPP
            CTERM(1)=CEXP(CMPLX(0.,RADIAN*DELPPM(1)*DEGPPM))
            CTERM(2)=CEXP(CMPLX(0.,-RADIAN*PPMINC*DEGPPM))
            jy = 0
            DO 220 JDATA=LDATST,LDATEN
               jy = jy + 1
               data_ph1(jy) = data_zero(jdata) * cterm(1)
               CTERM(1)=CTERM(1)*CTERM(2)
 220        continue
            degzer_try = degzer - deltz
            do 230 jdegz = 1, ndegz(1)
               degzer_try = degzer_try + deltz
               cterm(3) = cexp(cmplx(0., radian * degzer_try))
               dist = 0.
               sum = 0.
               dterm(1) = real(data_ph1(1) * cterm(3))
               do 240 jy = 2, ny
                  dterm(2) = real(data_ph1(jy) * cterm(3))
                  if (.not.lcy_skip(jy)) then
                     dist = dist + dabs(dterm(2) - dterm(1))**ipowph
                     sum = sum + sngl(dterm(2))
                  end if
                  dterm(1) = dterm(2)
 240           continue
               if (dist .le. distmn) then
                  distmn = dist
                  sumbes = sum
                  degppm_best = degppm
                  degzer_best = degzer_try
               end if
 230        continue
 210     continue
         if (.not.startp) degppm = degppm_best
         if (.not.startz) degzer = degzer_best
c        ---------------------------------------------------------------------
c        SUMBES < 0 implies that spectrum is negative of properly phased
c        one.  This may be incorrect if there is a big offset or if
c        DEGPPM is so big that some peaks are negative and some positive
c        (but that would require a bigger DEGPPM than usually allowed).
c        ---------------------------------------------------------------------
         if (ndegz(1) .gt. 1   .and.   sumbes .lt. 0.)
     1         degzer = degzer - 180.
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GBACKG ()
C
C  Load BACKGR(JY,JBACKG) = background terms (REAL).
C  Background is a cubic B-spline with NBACKG equally spaced knots on the
C    (-DELPPM)-axis; i.e., the grid points are positive and increasing (for
C    computational convenience).
C  Background is scaled by RMSAMP (rms amplitude of frequency-domain data) for
C    similar magnitudes of Hessian columns and for near data-independence of
C    ALPHAB and ALPHAS values.
C  Knots 2 and (NBACKG-1) correspond to PPMST and PPMEND; i.e., there is 1
C    knot beyond each endpoint (to allow full freedom at the endpoints).
C  Neither speed nor precision are important in this subprogram.
C
      INCLUDE 'lcmodel.inc'
      real ppmmax(mgap + 1), ppmmin(mgap + 1)
      CHSUBP='GBACKG'
      IF (NBACKG .LT. nbckmn) THEN
         IF (MBACKG .LT. nbckmn) CALL ERRMES (2, 4, CHSUBP)
         NBACKG=nbckmn
      END IF
      DO 110 JBACKG=1,MBACKG
         DO 120 JY=1,NY
            BACKGR(JY,JBACKG)=0.
  120    CONTINUE
         do 130 icol = 1, mbackg
            regb(jbackg, icol) = 0.
 130     continue
  110 CONTINUE
      DELTA=(DELPPM(1)-DELPPM(NY))/FLOAT(NBACKG-3)
      IF (DELTA .LE. 0.) CALL ERRMES (3, 4, CHSUBP)
      region_min = 4. * delta
      FNORM=RMSAMP/(6.*DELTA**4)
c     -------------------------------------------------------------------------
c     NBACKG will be recomputed (and may change if NGAP>0) below.
c     PPMMAX(J) = the left (high) ppm value for the Jth background region.
c     PPMMIN(J) = the right (low) ppm value for the Jth background region.
c     NREGION = # of ppm-regions with separate baselines (which may be
c               reduced).
c     -------------------------------------------------------------------------
      do 210 jgap = 1, ngap
         ppmmax(jgap + 1) = ppmgap(2, jgap)
         ppmmin(jgap) = ppmgap(1, jgap)
 210  continue
      ppmmax(1) = delppm(1) + ppmcen
      ppmmin(ngap + 1) = delppm(ny) + ppmcen
      nregion = ngap + 1
      nsep = 0
      if (ngap .gt. 0) then
c        ----------------------------------------------------------------------
c        PPMSEP: the baseline should be separated into separate baselines with
c           separate free boundaries on each side of PPMSEP.
c        ----------------------------------------------------------------------
         do 230 jsep = 1, mgap
            if (ppmsep(jsep) .ge. ppmst - region_min   .or.
     1          ppmsep(jsep) .le. ppmend + region_min) go to 240
            nsep = nsep + 1
            ppmsep(nsep) = ppmsep(jsep)
            if (lprint .gt. 0) write (lprint, 5232) ppmsep(nsep)
 5232       format (//'PPMSEP =', f7.3/)
 230     continue
 240     continue
         if (ldump(2)) then
            write (6,5230) nbackg
 5230       format (/'NBACKG =', i3)
            write (6, 5240) ((ppmgap(j, k), j=1,2), k=1,ngap)
 5240       format ('PPMGAP =', 5(2f8.2, 3x))
            write (6,5241) (ppmmax(j), ppmmin(j), j=1,nregion)
 5241       format ('PPMMNX =',  5(2f8.2, 3x))
            write (6, 5242) (ppmsep(j), j=1,nsep)
 5242       format('PPMSEP =', 10f8.2)
         end if
      end if
c     -------------------------------------------------------------------------
c     Merge all baseline regions without PPMSEP between them.
c     -------------------------------------------------------------------------
      nregion_old = nregion
      do 250 kregion = 1, nregion_old - 1
         do 260 jregion = 1, nregion - 1
            do 270 jsep = 1, nsep
               if (ppmmin(jregion) .ge. ppmsep(jsep)   .and.
     1             ppmmax(jregion + 1) .le. ppmsep(jsep)) go to 260
 270        continue
            call merge_right(jregion, ppmmin, ppmmax, nregion)
c            write (6, 9270) jregion
c 9270       format (/'Right-merging region', i2)
c            write (6,5241) (ppmmax(j), ppmmin(j), j=1,nregion)
            go to 250
 260     continue
 250  continue
c     -------------------------------------------------------------------------
c     If a region is too small to support its own baseline (of 4 knots), then
c        merge it with a neighboring region (for background only; the gaps
c        remain).  Merging will be done over lowest-priority PPMSEP(J) (i.e.,
c        with the highest J).  Typical usage is to input PPMSEP(1)=4.65, when
c        the feet of the water peak are strong (especially when they are
c        dispersive, pointing in opposite directions).
c     -------------------------------------------------------------------------
      nregion_old = nregion
      do 280 kregion = 1, nregion_old
         do 290 jregion = 1, nregion
            test = ppmmax(jregion) - ppmmin(jregion)
            if (test .ge. region_min) go to 290
c           -------------------------------------------------------------------
c           Region JREGION is too small; merge it.
c           From the above merges, there must already be a PPMSEP between all
c              regions.  Merge over the lowest priority PPMSEP.
c           -------------------------------------------------------------------
            if (jregion .ge. nregion) then
               jsep_right = 0
            else
               do 310 jsep_right = 1, nsep
                  if (ppmmin(jregion) .ge. ppmsep(jsep_right)   .and.
     1                ppmmax(jregion + 1) .le. ppmsep(jsep_right))
     2               go to 320
 310           continue
 320           continue
            end if
            if (jregion .le. 1) then
               jsep_left = 0
            else
               do 330 jsep_left = 1, nsep
                  if (ppmmin(jregion - 1) .gt. ppmsep(jsep_left)   .and.
     1                ppmmax(jregion) .lt. ppmsep(jsep_left)) go to 340
 330           continue
 340           continue
            end if
c            write (6, 9340) jregion, jsep_right, jsep_left,
c     1                      ppmmax(jregion), ppmmin(jregion)
c 9340       format ('JREGION =', i2, 5x, 'JSEP_RIGHT =', i2, 5x,
c     1              'JSEP_LEFT =', i2, 5x, 'PPMMAX =', 1pe9.2, 5x,
c     2              'PPMMIN =', e9.2)
            if (jsep_left .gt. jsep_right) then
               call merge_left(jregion, ppmmin, ppmmax, nregion)
            else if (jsep_right .gt. 0) then
               call merge_right(jregion, ppmmin, ppmmax, nregion)
            end if
            go to 280
 290     continue
 280  continue
c     -------------------------------------------------------------------------
c     Compute BACKGR & REGB for each region.  REGB will be block-diagonal,
c        with a block for each Region.
c     -DELPPM axis is used
c     X = -DELPPM = PPMCEN - PPM increases with JBACKG
c     -------------------------------------------------------------------------
      nbackg = 0
      do 410 jregion = 1, nregion
         nbackg_use = max0(nbckmn, 3 +
     1                nint((ppmmax(jregion) - ppmmin(jregion)) / delta))
         nbackg_start = nbackg + 1
         nbackg_end = nbackg + nbackg_use
         delta_use = (ppmmax(jregion) - ppmmin(jregion)) /
     1               float(nbackg_use - 3)
         xknot = ppmcen - ppmmax(jregion) - 2. * delta_use
         DO 420 JBACKG=1,NBACKG_use
            nbackg = nbackg + 1
c           -------------------------------------------------------------------
c           This should have been avoided with SETUP 2
c           -------------------------------------------------------------------
            IF (NBACKG .gt. mbackg) CALL ERRMES (4, 5, CHSUBP)
            XKNOT=XKNOT+DELTA_USE
c            XMIN=amax1(XKNOT, ppmcen - ppmmax(jregion)) - 2. * delta_use
            XMIN=XKNOT - 2. * delta_use
c            XMAX=amin1(XKNOT+2.*DELTA_USE,
c     1                 ppmcen - ppmmin(jregion) + delta)
            XMAX=XKNOT+2.*DELTA_USE
            XPPM_KNOT=PPMCEN-XKNOT
            NONNEG(NMETAB+NBACKG)=XPPM_KNOT.LE.PPMPOS(1) .AND.
     1                            XPPM_KNOT.GE.PPMPOS(2)
            DO 430 JY=1,NY
               X=-DELPPM(JY)
               xppm = delppm(jy) + ppmcen
c              ----------------------------------------------------------------
c              This baseline is exclusively for PPM(JY) in this region.
c              PPM & knots only coincide exactly for PPMST & PPMEND.
c                 Tolerances below are mainly for these 2 cases.
c              ----------------------------------------------------------------
               if (xppm .gt. ppmmax(jregion) + .001*delta_use   .or.
     1             xppm .lt. ppmmin(jregion) - .001*delta_use) go to 430
               IF (X .LE. XMIN) GO TO 430
               IF (X .GE. XMAX) GO TO 420
               BACKGR(JY,NBACKG)=(DIM(X,XMIN)**3 -
     1                            4.*DIM(X,XMIN+DELTA_USE)**3 +
     2                            6.*DIM(X,XKNOT)**3 -
     3                            4.*DIM(X,XMAX-DELTA_USE)**3 +
     4                            DIM(X,XMAX)**3) * FNORM
 430        CONTINUE
 420     CONTINUE
         if (ldump(2)) then
            write (6, 5420) jregion
 5420       format (/'Region =', i2)
            write (6,5241) (ppmmax(j), ppmmin(j), j=1,nregion)
            write (6,5422) nbackg_use, delta_use, delta_use/delta
 5422       format ('NBACKG_USE =', i3/ 'DELTA_USE =', 1pe11.4/
     1              'delta ratio =', e11.4)
         end if
C        ----------------------------------------------------------------------
C        Compute regularizor for cubic B-splines with equally spaced knots and
C           no boundary conditions.
C        For simplicity, no use is made of the banded structure; the full
C           matrices are used with many 0 elements.
c        Fill upper triangle of REGB was zeroed above (DO 130).
C        NBACKG_use >= 6 was imposed above.
C        ----------------------------------------------------------------------
         DO 440 IROW=nbackg_start, NBACKG_end
            REGB(IROW,IROW)=16.
            IF (IROW+1 .LE. NBACKG_end) REGB(IROW,IROW+1)=-9.
            IF (IROW+3 .LE. NBACKG_end) REGB(IROW,IROW+3)=1.
 440     CONTINUE
         REGB(nbackg_start,nbackg_start)=2.
         REGB(NBACKG_END, NBACKG_END)=2.
         REGB(nbackg_start, nbackg_start + 1)=-3.
         REGB(NBACKG_END-1, NBACKG_END)=-3.
         REGB(nbackg_start + 1, nbackg_start + 1)=8.
         REGB(NBACKG_END-1,NBACKG_END-1)=8.
         REGB(nbackg_start + 1, nbackg_start + 2)=-6.
         REGB(NBACKG_END-2, NBACKG_END-1)=-6.
         REGB(nbackg_start + 2, nbackg_start + 2)=14.
         REGB(NBACKG_END-2, NBACKG_END-2)=14.
         DO 450 IROW=nbackg_start,NBACKG_END
            DO 460 ICOL=IROW+1,NBACKG_END
               REGB(ICOL,IROW)=REGB(IROW,ICOL)
 460        CONTINUE
 450     CONTINUE
 410  continue
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE merge_right (lregion, ppmmin, ppmmax, nregion)
C
C  Assumes that lregion < nregion
c
      character chsubp*6
      real ppmmax(*), ppmmin(*)
      chsubp = 'MERGRT'
      if (lregion .ge. nregion) call errmes (1, 5, chsubp)
      nregion = nregion - 1
      ppmmin(lregion) = ppmmin(lregion + 1)
      do 110 jregion = lregion + 1, nregion
         ppmmax(jregion) = ppmmax(jregion + 1)
         ppmmin(jregion) = ppmmin(jregion + 1)
 110  continue
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE merge_left (lregion, ppmmin, ppmmax, nregion)
C
C  Assumes that 1 < lregion
c
      character chsubp*6
      real ppmmax(*), ppmmin(*)
      chsubp = 'MERGLF'
      if (lregion .le. 1) call errmes (1, 5, chsubp)
      nregion = nregion - 1
      ppmmin(lregion - 1) = ppmmin(lregion)
      do 110 jregion = lregion, nregion
         ppmmax(jregion) = ppmmax(jregion + 1)
         ppmmin(jregion) = ppmmin(jregion + 1)
 110  continue
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine check_chless ()
c
c .TRUE. if STARTV is to be repeated with all CHLESS omitted from Prel.
c The concentrations [SOLBES(JMETAB,2)] with all NACOMB containing the unique
c    prefixes CHLESS(JLESS) and CHMORE will be summed.
c    OMIT_CHLESS = [CONCEN corresponding to CHLESS(JLESS) exceeds that for
c                   RLESMO*CHMORE].
c    Then Prel will be repeated with all CHLESS components omitted.
c Typically, CHLESS='L09' 'L20' and CHMORE='L13' to avoid Lip09 replacing
c    Lip13 in Prel.
c RLESMO: OMIT_CHLESS=T if CONC_LESS(JLESS) / CONC_MORE > RLESMO (default 0
c         will set OMIT_CHLESS=F).
c NCHLES > 0 to allow OMIT_CHLESS
c
      INCLUDE 'lcmodel.inc'
      external ilen
      integer lless(mmetab)
      real conc_less(mmetab)
      chsubp = 'LESSMO'
      omit_chless = .false.
      if (chmore .eq. ' '   .or.  rlesmo .le. 0.   .or.
     1    nchles .le. 0) return
      if (nchles .gt. mmetab) call errmes (1, 4, chsubp)
      do 110 jless = 1, nchles
         if (chless(jless) .eq. ' ')  return
         if (chless(jless) .eq. chmore) then
            call errmes (2, 3, chsubp)
            return
         end if
         conc_less(jless) =-1.
         lless(jless) = ilen(chless(jless))
 110  continue
      conc_more = -1.
      lmore = ilen(chmore)
      lnot1 = nnot1
      do 210 jmetab = 1, nmetab
         do 220 jless = 1, nchles
            if (index(nacomb(jmetab), chless(jless)(:lless(jless)))
     1          .eq. 1) then
               conc_less(jless) = amax1(0., conc_less(jless)) +
     1                            solbes(jmetab, 2)
               lnot1 = lnot1 + 1
               if (lnot1 .gt. mmetab_extra) then
                  call errmes (3, 3, chsubp)
                  return
               end if
               chnot1(lnot1) = nacomb(jmetab)
            end if
 220     continue
         if (index(nacomb(jmetab), chmore(:lmore)) .eq. 1)
     1      conc_more = amax1(0., conc_more) + solbes(jmetab, 2)
 210  continue
      do 310 jless = 1, nchles
         if (conc_less(jless) .gt. 0.) go to 320
 310  continue
      return
 320  do 350 jless = 1, nchles
         if (conc_less(jless) .ge. conc_more * rlesmo) go to 400
 350  continue
      return
 400  omit_chless = .true.
      call errmes (4, 1, chsubp)
      nnot1 = lnot1
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE TWOReg ()
C
C  Find best unconstrained solution.
c  If unusual conditions occur, do a series of NDEGPPM3 analyses with fixed
c     DEGPPM.
c  IDGPPM = -1 (default) to not do any DEGPPM-fixed analyses, only standard
c              free analysis.
c            0 (default in lipid-2 & PPMST<5) for analyses stressing flat
c              baselines (by allowing medium SSQ/SSQMIN), where there is no
c              small peak that can be lost by incorrect DEGPPM.  This can
c              also be used to input any values for RBASMX & RSDGP3; only
c              differences between IDGPPM=0 & >0 are settings in MYCONT.
c            1 (default in breast-2, only-cho-1 & liver-2) for analyses
c              stressing good phasing and fit (at expense flat baselines),
c              to find small peaks (like Cho).  Forces fixed-DEGPPM series,
c              even if baseline in unconstrained solution is flat.
c            2 (default in muscle-2) (actually between 0 and 1) stresses
c              good fit, but does not force fixed-DEGPPM series if
c              unconstrained solution is flat.
c  To switch on DEGPPM-analyses, only need to input IDGPPM>=0 (and check
c     RBASMX)
c
      INCLUDE 'lcmodel.inc'
      external r_base_sol_big
      logical big_base_free(2), lterm, r_base_sol_big
      CHSUBP='TWOREG'
c     -----------------------------------------------------------------------
c     NDEGPPM3_USED = 0 for free analysis.
c     -----------------------------------------------------------------------
      ndegppm3_used = 0
      call tworg2 (0, .false.)
      degppm_free = (phitot(2) + parbes(lphast + 1, 2)) / radian
      big_base_free(1) = r_base_sol_big(1)
      big_base_free(2) = r_base_sol_big(2)
      lterm = degppm_free .lt. dgppmn   .or.   degppm_free .gt. dgppmx
     1        .or.   big_base_free(1)
c     ------------------------------------------------------------------------
c     If there are no unusual conditions or IDGPPM=-1, return and skip
c        series with fixed DEGPPM.
c     ------------------------------------------------------------------------
      if (fxdegp   .or.   .not.lterm   .or.   idgppm .lt. 0) return
      call errmes (12, 1, chsubp)
      call tworeg_sav ()
      degp_degp(0) = degppm_free
      if (ddegp3 .le. 0.   .or.   mdegp3 .le. 2)
     1   call errmes (13, 4, chsubp)
c     ------------------------------------------------------------------------
c     It is necessary to start from near DEGPPM of the free DEGPPM analysis,
c        even if that is outside range of DGPPM*.
c     ------------------------------------------------------------------------
      dgppmn_use = amin1(degppm_free, dgppmn)
      dgppmx_use = amax1(degppm_free, dgppmx)
      if (degppm_free .lt. dgppmn   .or.   degppm_free .gt. dgppmx)
     1   call errmes (16, 1, chsubp)
      ndegppm3 = nint((dgppmx_use - dgppmn_use) / ddegp3)
      if (ndegppm3 .gt. mdegp3) then
         call errmes (14, 2, chsubp)
         ndegppm3 = mdegp3
      end if
      ddegp3_use = (dgppmx_use - dgppmn_use) / float(ndegppm3)
      degzer_free = (phitot(1) + parbes(lphast, 2)) / radian
      exdegp_orig = exdegp
      sddegp_orig = sddegp
c     -----------------------------------------------------------------------
c     Start with DEGPPM of previous unconstrained analysis (not EXDEGP)
c     -----------------------------------------------------------------------
      exdegp = degppm_free
      sddegp = 0.
      fxdegp = .true.
      ssqbes_degp = rrange
      do 110 jdgppm = 1, ndegppm3
         exdegp = exdegp + ddegp3_use
c        ---------------------------------------------------------------------
c        6*SD used (instead of original 4) in the two statements below
c        to allow for prostate with SDDEGP=30 & DGPPM*=30-150
c        ---------------------------------------------------------------------
         IF (ABS(exdegp-exdegp_orig).GT.6.*sddegp_orig) GO TO 110
         if (exdegp .gt. dgppmx_use + .5 * ddegp3_use) go to 200
         ndegppm3_used = ndegppm3_used + 1
         call tworg1 ()
 110  continue
c     ------------------------------------------------------------------------
c     Start decreasing DEGPPM from DEGPPM & DEGZER of solution with free
c        DEGPPM
c     ------------------------------------------------------------------------
 200  parbes(lphast, 2) = degzer_free * radian - phitot(1)
      call rephas()
      exdegp = degppm_free
      do 210 jdgppm = 1, ndegppm3
         exdegp = exdegp - ddegp3_use
          IF (ABS(exdegp-exdegp_orig).GT.6.*sddegp_orig) GO TO 210
         if (exdegp .lt. dgppmn_use - .5 * ddegp3_use) go to 300
         ndegppm3_used = ndegppm3_used + 1
         call tworg1 ()
 210  continue
c     ------------------------------------------------------------------------
c     Regenerate solution with minimum SSQ
c     SSQBES(7) is from initial analysis with free DEGPPM
c     SSQBES(4) is from min-SSQ DEGPPM-constrained
c     ------------------------------------------------------------------------
 300  if (ndegppm3_used .le. 0) go to 700
      if (ssqbes(7) .le. ssqbes(4)   .and.
     1    .not.big_base_free(2)) then
         call errmes (17, 1, chsubp)
         call savbes(-7)
      else
c        ---------------------------------------------------------------------
c        Select min-SSQ solution if it has R_BASE_SOL_BIG(2) = F.
C        No! Comment this out.  Always look for smaller spline distance.
c        ---------------------------------------------------------------------
         call savbes(-4)
c         if (r_base_sol_big(2)) then
c           ------------------------------------------------------------------
c           Try to find solution with the min spline distance that has an
c              SSQ < (min SSQ) * RSDGP3
c           RSDGP3 < 1. will force the above min-SSQ solution to be used,
c              even though its R_BASE_SOL_BIG(2) = T
c           RSDGP3 very large will force min spline distance solution, even
c              if its SSQ is very high.
c           ------------------------------------------------------------------
            bound = ssqbes(4) * rsdgp3
            distmn = rrange
            ldistmn = -1
            do 350 j = 0, ndegppm3_used
               if (ssq_degp(j) .lt. bound) then
                  if (dist_degp(j) .le. distmn) then
                     distmn = dist_degp(j)
                     ldistmn = j
                  end if
               end if
 350        continue
            if (ldistmn .ge. 0) call savbes(-ldistmn - 7)
c         end if
      end if
c     ------------------------------------------------------------------------
c     These have to be restored, for use in FINOUT.  This may not be optimal
c        for SOLVE in FINOUT, but removing these produced identical solution
c        & error estimates.
c     ------------------------------------------------------------------------
 700  exdegp = exdegp_orig
      fxdegp = .false.
      sddegp = sddegp_orig
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE tworg1 ()
c
c  Find best solution for a fixed DEGPPM
c
      INCLUDE 'lcmodel.inc'
      CHSUBP='TWOREG'
      degp_degp(ndegppm3_used) = exdegp
      IF (LPRINT .GT. 0) WRITE (LPRINT,5110) exdegp
 5110 FORMAT (////25x, 'Analysis with fixed DEGPPM =', f8.2)
      parbes(lphast + 1, 2) = exdegp * radian - phitot(2)
      call rephas()
      call tworg2(ndegppm3_used, .true.)
      call tworeg_sav ()
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE tworeg_sav ()
c
c  Save current solution, which is in *(ndegppm3_used) and compute baseline
c  criteria.
c
      INCLUDE 'lcmodel.inc'
      call savbes (ndegppm3_used + 7)
      ssq_degp(ndegppm3_used) = ssqbes(2)
      alpb_degp(ndegppm3_used) = alpbbs(2)
      dist_degp(ndegppm3_used) = 0.
      do 120 j = nmetab + 1, nmetab + nbackg - 1
         dist_degp(ndegppm3_used) = dist_degp(ndegppm3_used) +
     1                       dabs(solbes(j, 2) - solbes(j + 1, 2))
 120  continue
      if (ssq_degp(ndegppm3_used) .le. ssqbes_degp) then
         ssqbes_degp = ssq_degp(ndegppm3_used)
         call savbes(4)
      end if
      end
c
c
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE tworg2 (jpass, fixed_degppm_series)
C
C  Call TWORG3 to find ALPHAB and ALPHAS so that PRMNMX(1,1)<PROB1<PRMNMX(2,1).
c  ENDPHA = T to do fixed-phase analysis after analysis for ALPHAB & ALPHAS
c             (previously the final analysis)
c  USEMXB = T choose final solution with the largest ALPHAB
c  USEMXB=F was not tested (assuming that it will never be used).
C
      INCLUDE 'lcmodel.inc'
      external ldegmx
      double precision parbes_phase_sav(2)
      LOGICAL fixed_degppm_series, fixed_phase,
     1        fxdegp_orig, fxdegz_orig, LDEGMX,
     1        ldegmx_sav, PREJOK, usesol
      CHSUBP='TWOREG'
C     -------------------------------------------------------------------------
C     Repeated calls to SETUP produce slightly slower execution and more
c        PNNLS.
c     Repeated calls to "Preliminary full..." below speed up the analysis.
C     -------------------------------------------------------------------------
      if (jpass .le. 0) then
         CALL SETUP (2)
C        ----------------------------------------------------------------------
C        When doing a series of fixed DEGPPM, rephasing has already set the
c           following to zero, but SETUP sets them.  Set them back to zero.
C        ----------------------------------------------------------------------
         if (fixed_degppm_series) then
            parnln(lphast) = 0.d0
            parnln(lphast + 1) = 0.d0
            parbes(lphast, 2) = 0.d0
            parbes(lphast + 1, 2) = 0.d0
         end if
      end if
C     -------------------------------------------------------------------------
C     Get Starting Solution with ALPBST and ALPSST.
C     -------------------------------------------------------------------------
 5240 FORMAT (////20X, 'Preliminary full analysis with alphaB =',
     1                 1PE9.2E2, ' and alphaS =', E9.2E2)
      IF (LPRINT .GT. 0) WRITE (LPRINT,5240) ALPHAB, ALPHAS
      if (imethd .eq. 1   .and.   nratio .le. 0) miter(2) = 2 * miter(2)
      CALL PLINLS (2, IERROR)
      IF (DMAX1(OBJECT,SDREF) .GE. DRANGE) CALL ERRMES (1, 4, CHSUBP)
      CALL SAVBES (1)
      CALL SAVBES (2)
c     -------------------------------------------------------------------------
c     Special quick exit, e.g., for huge Roche-2 analyses, where above values
c       of ALPHA* can be used as final values, and where huge NY makes FSHSSQ
c       inaccurate and PRMNMX much too narrow (because of so much baseline and
c       few data in sharp peaks).
c     -------------------------------------------------------------------------
      if (imethd .eq. 1) then
         NINFL(1)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                    imethd)
         NEXTR(1)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                    imethd)
         if (nratio .gt. 0) then
            call conc_prior ()
            if (nratio_used .gt. 0) then
               if (ninfl(1) .gt. 2) alphas = alphas * rdalpb**2
               CALL PLINLS (2, IERROR)
               if (object .lt. drange) then
                  CALL SAVBES (1)
                  CALL SAVBES (2)
                  NINFL(2)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                             THRLIN, imethd)
                  NEXTR(2)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                             THRLIN, imethd)
               end if
            end if
         end if
         do 210 jrepha = 2, mrepha(1)
            IF (.NOT.LDEGMX(1)) go to 800
            parbes(lphast, 2) = parbes(lphast, 2) * frepha
            parbes(lphast + 1, 2) = parbes(lphast + 1, 2) * frepha
            call rephas()
            CALL ERRMES (8, 1, CHSUBP)
            CALL PLINLS (2, IERROR)
            if (object .lt. drange) then
               CALL SAVBES (1)
               CALL SAVBES (2)
               NINFL(2)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                          THRLIN, imethd)
               NEXTR(2)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                          THRLIN, imethd)
            end if
 210     continue
         IF (LDEGMX(1)) CALL ERRMES (18, 2, CHSUBP)
         go to 800
      end if
      IF (RDALPB.LE.1. .OR. ALPBST.LT.SNGL(.99999*ALPBMN) .OR.
     1    ALPBST.GT.SNGL(1.00001*ALPBMX)) CALL ERRMES (2, 4, CHSUBP)
C     -------------------------------------------------------------------------
C     There is no point in varying ALPHAS when NSIDES=1; the peak shape does
c        not change with ALPHAS.
C     -------------------------------------------------------------------------
      if (nsides .le. 1   .and.   imethd .ne. 2) mdalpb = 0
      if (nratio .gt. 0   .and.   jpass .le. 0) then
C        ----------------------------------------------------------------------
c        Set up CPRIOR, matrix for priors of concentration ratios, after doing
c          preliminary Reference Solution and then Regula Falsi search with
c          ALSBMN to get CONC estimates for weighting the priors.
C        ----------------------------------------------------------------------
         SSQAIM=0.
         PENBES(2)=DRANGE
         CALL RFALSI (1, 1, .TRUE., ALPBBS(2), ALPSMN, -1./RRANGE, 0.,
     1                RRANGE, RRANGE, 1., PREJOK, PREJ1)
         IF (PREJOK .OR. (ALPHAB.GE..99999*ALPBMX .AND. MDALPB.LE.0))
     1       CALL SAVBES (2)
         call conc_prior ()
      end if
      exdegz_orig = exdegz
      exdegp_orig = exdegp
      sddegz_orig = sddegz
      sddegp_orig = sddegp
      fxdegz_orig = fxdegz
      fxdegp_orig = fxdegp
      SSQAIM=0.
      alpbbs(3) = -1.
c     -----------------------------------------------------------------------
c     The following is just to avoid compile diagnostic about possible
c     undefined variable
      ldegmx_sav = .false.
c     -----------------------------------------------------------------------
      DO 310 JREPHA=1,MREPHA(1)
         call tworg3(iabs(jrepha))
         ldegmx_sav = ldegmx(1)
c        --------------------------------------------------------------------
c        Do "final" solution with fixed phases.
c        ENDPHA = F by default to avoid making "normal" analyses (including
c                   test run) from taking 50% longer.
c        --------------------------------------------------------------------
         parbes_phase_sav(1) = parbes(lphast, 2)
         parbes_phase_sav(2) = parbes(lphast + 1, 2)
         fixed_phase = (endpha   .and.
     1                  (.not.ldegmx_sav   .or.
     2                   jrepha . eq. mrepha(1)))   .or.
     3                 (ldegmx_sav   .and.   jrepha . eq. mrepha(1))
c     4                 .or.   idgppm .ge. 0
         if (fixed_phase) then
            if (lprint .gt. 0) write(lprint, 5380)
 5380       format(//20x,
     1             'Repeat analysis with fixed phases from above')
            call rephas ()
            exdegz = phitot(1) / radian
            exdegp = phitot(2) / radian
            sddegz = 0.
            sddegp = 0.
            fxdegz = .true.
            fxdegp = .true.
            call tworg3(iabs(jrepha))
            exdegz= exdegz_orig
            exdegp= exdegp_orig
            sddegz= sddegz_orig
            sddegp= sddegp_orig
            fxdegz= fxdegz_orig
            fxdegp= fxdegp_orig
         end if
c        --------------------------------------------------------------------
c        Save solution with the max ALPHAB with SAVBES(3)
c        --------------------------------------------------------------------
         if (usemxb   .and.   alpbbs(2) .ge. alpbbs(3)) call savbes(3)
         if (.not.ldegmx_sav) go to 700
         if (jrepha .lt. mrepha(1)) then
c           ------------------------------------------------------------------
c           This will reduce the rephasing corrections by a factor of FREPHA.
c           Version 6.0-2C and before use FREPHA=1.
c           A reasonable value would be 0.5, which will hopefully avoid the
c             oscillations between two extreme phase sets that occurred in
c             previous versions.
c           ------------------------------------------------------------------
            if (.not.fixed_phase) then
               parbes(lphast, 2) = parbes_phase_sav(1) * frepha
               parbes(lphast + 1, 2) = parbes_phase_sav(2) * frepha
            end if
         end if
         CALL ERRMES (8, 1, CHSUBP)
 310  continue
      CALL ERRMES (9, 1, CHSUBP)
c     ------------------------------------------------------------------------
c     Regenerate solution with the max ALPHAB, if no stable phases have been
c     found.
c     ------------------------------------------------------------------------
 700  if (usemxb   .and.  ldegmx_sav) call savbes(-3)
  800 RETURN
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tworg3(jrepha)
c
C  Find ALPHAB and ALPHAS so that PRMNMX(1,1)<PROB1<PRMNMX(2,1).
c  Shares CALL ERRMES numbers with TWOREG.
C
      INCLUDE 'lcmodel.inc'
      LOGICAL PREJOK, usesol
 
      chsubp='TWOREG'
      DO 312 J=1,2
         NEXTR(J)=0
         NINFL(J)=-1
 312  CONTINUE
      PENBES(2)=DRANGE
      CALL RFALSI (1, 1, .TRUE., ALPbbs(2), ALPSMN, -1./RRANGE, 0.,
     1             RRANGE, RRANGE, 1., PREJOK, PREJ1)
      if (mdalpb .le. 0   .and.  daimbs .lt. rrange) then
         usesol = prejok   .or.   ALPHAB.GE..99999*ALPBMX
         if (.not.usesol   .and.   useany) call errmes (10, 2, chsubp)
         IF (usesol .OR.   useany) THEN
            CALL SAVBES (2)
            NINFL(1)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                       imethd)
            NEXTR(1)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                       imethd)
         END IF
      end if
      DO 320 JDALPB=1,MDALPB
         IF (LPRINT .GT. 0) WRITE (LPRINT,5320) JREPHA, JDALPB
 5320    FORMAT (//////10X, 'Phase Pair', I2, 5X, 'Decrease', I2,
     1           ' of alphaB/alphaS')
         CALL RFALSI (3, 1, .FALSE., ALPBBS(1)/RDALPB,
     1                ALPSBS(1)*RDALPB, SSQREF, 0., RRANGE, RRANGE,
     2                1., PREJOK, PREJ1)
         IF (PREJOK) THEN
            IF (PENBES(1) .LT. PENBES(2)) THEN
               CALL SAVBES(2)
               NINFL(1)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                          THRLIN, imethd)
               NEXTR(1)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                          THRLIN, imethd)
            ELSE IF (PENBES(1) .GE. RPENMX*PENBES(2)) THEN
               GO TO 330
            END IF
            IF (NEXTRE(PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN, imethd)
     1          .LE. 1) GO TO 330
         END IF
         IF (ALPBBS(1).LE.ALPBMN .OR. ALPSBS(1).GE.ALPSMX) GO TO 330
 320  CONTINUE
      IF (MDALPB .GT. 0) CALL ERRMES (3, 2, CHSUBP)
 330  IF (PENBES(2) .GE. DRANGE) CALL ERRMES (4, 4, CHSUBP)
      IF (MDALPB.LE.0 .OR. USINFL) THEN
         ITEST=NINFL(1)-2
      ELSE
         ITEST=NEXTR(1)-1
      END IF
      IF (ITEST.GT.0 .AND. ALPBBS(2).GT.ALPBMN .AND.
     1    .NOT.VITRO  .AND.   ALPSBS(2).LT.ALPSMX) THEN
         IF (LPRINT .GT. 0) WRITE (LPRINT,5340)
 5340    FORMAT (//////20X, 'Step 2: Increase alphaS with fixed',
     1                 ' alphaB')
         IF (AMAX1(PRMNMX(2,1),PRMNMX(2,3)) .GT. PRMNMX(1,2)) CALL
     1         ERRMES (5, 4, CHSUBP)
         CALL SSRANG (2)
         DO 340 JNONL=1,NNONL
            PARNLN(JNONL)=PARBES(JNONL,2)
 340     CONTINUE
         CALL RFALSI (2, 2, .FALSE., ALPBBS(2), ALPSBS(2),
     1      SSQBES(2), 1., RRANGE, RRANGE, RALINC, PREJOK, PREJ1)
         IF (.NOT.PREJOK) THEN
            CALL ERRMES (6, 2, CHSUBP)
            GO TO 800
         END IF
         if (skip_step3) then
            call savbes(2)
            NINFL(2)=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                       THRLIN, imethd)
            NEXTR(2)=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                       THRLIN, imethd)
            go to 800
         end if
         if (imethd .eq. 2   .and.   accept_step2) call savbes(2)
         IF (MDALPB.LE.0 .OR. USINFL) THEN
            ITEST=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                    imethd)-2
         ELSE
            ITEST=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS, THRLIN,
     1                    imethd)-1
         END IF
         IF (ITEST .LE. 0) THEN
C           -------------------------------------------------------------------
C           Test for special exit from Step 2 with 2-inflection-point or
C              unimodal lineshape and PREJ1<PRMNMX(2,3).
C           -------------------------------------------------------------------
            IF (PREJ1 .GT. PRMNMX(2,3)) THEN
               IF (LPRINT .GT. 0) WRITE (LPRINT,5350)
 5350          FORMAT (//////20X, 'Step 3A: Decrease alphaB and ',
     1                       'alphaS')
               CALL SSRANG (3)
               CALL RFALSI (3, 3, .FALSE., ALPBBS(2), ALPSBS(1),
     1                      SSQREF, 0., SSQBES(1), 1., 1./RALINC,
     2                      PREJOK, PREJ1)
            END IF
         ELSE
            IF (LPRINT .GT. 0) WRITE (LPRINT,5360)
 5360       FORMAT (//////20X, 'Step 3B: Decrease alphaB with ',
     1                    'fixed alphaS')
            CALL SSRANG (3)
            CALL RFALSI (1, 3, .FALSE., ALPBBS(2), ALPSBS(1),
     1                   -1./RRANGE, 0., SSQBES(1), 1., 1./RALINC,
     2                   PREJOK, PREJ1)
c           ------------------------------------------------------------------
c           When IMETHD=2, use ALPBMN in Step 3, even if PREJ is not down to
c              PRMNMX(2,3) to avoid using results of Step 1.
c           ------------------------------------------------------------------
            if (imethd .eq. 2   .and.   accept_alpbmn   .and.
     1          is_alpbmn) then
               call savbes(1)
               prejok = .true.
            end if
         END IF
         IF (PREJOK) THEN
            NEXTR_test=NEXTRE (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                         THRLIN, imethd)
            NINFL_test=INFLEC (PARBES(1,1), NSIDE2, DPY, DGAUSS,
     1                         THRLIN, imethd)
            if (ninfl_test .ge. ninfl(1)  .and.
     1          nextr_test .ge. nextr(1)   .and.   imethd .ne. 2) then
               IF (LPRINT .GT. 0) WRITE (LPRINT,5370)
 5370          FORMAT (//////' Steps 2 and 3 did not decrease ',
     1              'the no. of inflections or extrema; ',
     2              'they will not be used.'//)
               GO TO 800
            END IF
            CALL SAVBES (2)
            NINFL(2)=ninfl_test
            NEXTR(2)=NEXTR_test
         ELSE
            CALL ERRMES (7, 1, CHSUBP)
         END IF
      END IF
 800  return
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      logical function r_base_sol_big(istage)
c
c ISTAGE = 1 for testing free solution.
c          2 for testing min-SSQ solution from fixed-DEGPPM series.
c .TRUE. if the ratio of the height difference in the baseline to the height
c    difference in the solution > RBASMX(ISTAGE)
c Uses the values currently in *BES(2).
c
      INCLUDE 'lcmodel.inc'
      logical lerror
      chsubp = 'RBASOL'
c     ------------------------------------------------------------------------
c     Statements down thru 165 (except SAVBES) are from FINOUT.  Skipping
c        the 1st call to SOLVE below caused a FATAL
c     ------------------------------------------------------------------------
      call savbes(6)
      DO 150 JPAR=1,NLIN
         SOLUTN(JPAR)=SOLBES(JPAR,2)
c         NONNEG(JPAR)=SOLUTN(JPAR) .EQ. 0.D0
  150 CONTINUE
      DO 160 JNONL=1,NNONL
         PARNLN(JNONL)=PARBES(JNONL,2)
  160 CONTINUE
      ALPHAB=ALPBBS(2)
      ALPHAS=ALPSBS(2)
      CALL SOLVE (2, .TRUE., PMQBES(2), .FALSE., LERROR)
      IF (LERROR) CALL ERRMES (1, 4, CHSUBP)
C     -------------------------------------------------------------------------
C     Now call SOLVE to get the following:
C     YFITRE(JY,0) = fit
C     BACKRE(JY) = background
C     -------------------------------------------------------------------------
      DO 165 JPAR=1,NLIN
         SOLUTN(JPAR)=SOLBES(JPAR,2)
  165 CONTINUE
      CALL SOLVE (2, .FALSE., 0.D0, .TRUE., LERROR)
      IF (LERROR) CALL ERRMES (2, 4, CHSUBP)
      call savbes(-6)
      solmax = -rrange
      basmin = rrange
      basmax = -rrange
      lmax = 0
      lmin = 0
      do 210 j = 1, nyuse
         solmax = amax1(solmax, abs(yfitre(j, 0) - backre(j)))
         if (backre(j) .le. basmin) then
            basmin = backre(j)
            lmin = j
         end if
         if (backre(j) .ge. basmax) then
            basmax = backre(j)
            lmax = j
         end if
 210  continue
c     ------------------------------------------------------------------------
c     Count distances from both (1) & (nyuse) to any internal extremum.  This
c        properly counts extra distance to a bowed baseline.
c     ------------------------------------------------------------------------
c      base_dist = 0.
c      if (lmin .ne. 1   .and.   lmin .ne. nyuse)
c     1   base_dist = abs(backre(1) - basmin) +
c     2               abs(backre(nyuse) - basmin)
c      if (lmax .ne. 1   .and.   lmax .ne. nyuse)
c     1   base_dist = base_dist + abs(backre(1) - basmax) +
c     2                           abs(backre(nyuse) - basmax)
c      if (base_dist .eq. 0.) base_dist = basmax - basmin
      base_dist = basmax - basmin
      if (solmax .le. 0.) then
         call errmes (3, 2, chsubp)
         r_base_sol_big = .true.
      else
         r_base_sol_big = base_dist / solmax .gt. rbasmx(istage)
 9210    format ('base_dist =', 1pe11.3, 3x, 'solmax =', e11.3, 3x,
     1           'ratio =',  e11.3, 3x, 'ratio limit =', e11.3/)
         if (lprint .gt. 0) write (lprint, 9210) base_dist, solmax,
     1                         base_dist / solmax, rbasmx(istage)
      end if
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine conc_prior()
c
c     Make CPRIOR, the matrix of priors for the CONC ratios in CHRATI.
c     CHRATI must contain exactly one Metabolite Name in the numerator,
c followed by a slash, follwed by a CONC sum (with no whitespace).
c     CHRATW can contain a further CONC sum.  These will not go into the
c constraint, but will be used to compute CSUM, for the weights.  For example,
c might have CHRATW='Lip09' when the denominator in CHRATI is MM09, which is
c very difficult to resolve from Lip09.
c     In many cases, a prior will simply be skipped if an error condition is
c found.
c
c     The prior is
c c_num / c_sum = exrati +- sdrati
c c_num - c_sum * exrati = 0 +- c_sum * sdrati
c where c_num is a single metabolite;
c       c_sum can be a sum of several metabolites.
c       For weighting, c_sum is estimated from the CONC values of the
c         Preliminary Full Analysis with ALP*ST.
c         CHRATW can contain further metabolites added to c_sum for weighting.
c
      INCLUDE 'lcmodel.inc'
      character chterm*(mchmet)
      logical denom_absent
      CHSUBP='CONCPR'
c     -------------------------------------------------------------------------
c     NRATIO = input number of priors
c     NRATIO_USED = actual number of priors that will be used.
c     -------------------------------------------------------------------------
      if (nratio .gt. mmetab) call ERRMES (1, 4, CHSUBP)
      call parse_prior ()
      nratio_used = 0
      do 110 jratio = 1, nratio
c        ----------------------------------------------------------------------
c        Increment NRATIO_USED as though this try will be accepted.
c           NRATIO_USED will have to be decremented at 300 if this try is not
c           accepted.
c        ----------------------------------------------------------------------
         nratio_used = nratio_used + 1
         do 120 j = 1, mmetab
            cprior(nratio_used, j) = 0.
 120     continue
         islash = index(chrati(jratio), '/')
         if (islash .le. 1) then
            IF (LPRINT .GT. 0) WRITE (LPRINT,5120) CHrati(jratio)
 5120       FORMAT (' Missing slash in the ratio in CHRATO: ',A)
            call errmes (2, 4, chsubp)
         end if
c        ----------------------------------------------------------------------
c        Parse numerator, which must have exactly one Metabolite Name.
c        ----------------------------------------------------------------------
         LENGTH=islash - 1
         IF (LENGTH .GT. MCHMET) THEN
            IF (LPRINT .GT. 0) WRITE (LPRINT,5130) CHrati(jratio)
 5130       FORMAT (' Incorrect numerator in CHRATI =',A)
            CALL ERRMES (3, 4, CHSUBP)
         END IF
         chterm = chrati(jratio)(1:length)
         do 125 jomit = 1, nnorat
            if (chterm .eq. norato(jomit)) go to 300
 125     continue
         do 130 jmetab = 1, nmetab
            if (chterm .eq. nacomb(jmetab)) then
               if (exrati(jratio) .lt. 0.   .or.
     1             sdrati(jratio) .le. 0.) call ERRMES (4, 4, CHSUBP)
               lmetab_prior(nratio_used) = jmetab
               cprior(nratio_used, jmetab) = 1.
               go to 150
            end if
 130     continue
         go to 300
c        ----------------------------------------------------------------------
c        Parse denominator.
c        ----------------------------------------------------------------------
 150     csum = 0.
         denom_absent = .true.
         call parse_sum(exrati(jratio), chrati(jratio)(islash+1:),
     1                  mchratio - islash, nratio_used, csum,
     2                  denom_absent)
c        ----------------------------------------------------------------------
c        Parse CHRATW for extra metabolites to be included for weighting.
c        ----------------------------------------------------------------------
         call parse_sum(0., chratw(jratio), mchratio, nratio_used, csum,
     1                  denom_absent)
c        ----------------------------------------------------------------------
c        Skip prior if metabolites in denominator & weight terms are not in
c          analysis.
c        ----------------------------------------------------------------------
         if (denom_absent) go to 300
         if (csum .le. 0.) then
c           -------------------------------------------------------------------
c           CSUM from denominator was 0.  Use sum over all CONCs times a small
c             factor to produce a weight, which is all that CSUM affects.
c             CSUM = 0 implies that the numerator CONC should also be zero.
c             Therefore, the weight should be strong, and FCSUM small, e.g.,
c             FCSUM = 1.e-2 or 1.e-3
c           -------------------------------------------------------------------
            dterm(1) = 0.d0
            do 160 jmetab = 1, nmetab
               dterm(1) = dterm(1) + solbes(jmetab, 1)
 160        continue
            if (fcsum .le. 0.) call ERRMES (5, 4, CHSUBP)
            csum = fcsum * sngl(dterm(1))
c           -------------------------------------------------------------------
c           If all CONCs are 0, then skip this constraint
c           -------------------------------------------------------------------
            if (csum .eq. 0.) then
               call ERRMES (6, 2, CHSUBP)
               go to 300
            end if
         end if
         sqrtwt = 1. / (csum * sdrati(jratio))
         do 170 jmetab = 1, nmetab
            cprior(nratio_used, jmetab) = cprior(nratio_used, jmetab) *
     1                                    sqrtwt
 170     continue
c        ----------------------------------------------------------------------
c        Move *(JRATIO) down to *(NRATIO_USED) for later output.
c        ----------------------------------------------------------------------
         chrato(nratio_used) = chrato(jratio)
         exrati(nratio_used) = exrati(jratio)
         sdrati(nratio_used) = sdrati(jratio)
         sqrtwt_ratio_used(nratio_used) = sqrtwt
         go to 110
 300     nratio_used = nratio_used - 1
 110  continue
c     -------------------------------------------------------------------------
c     Dump prior matrix (if IPDUMP >= 3) or active CHRATO
c     -------------------------------------------------------------------------
      if (min0(nratio_used, lprint) .gt. 0) then
         if (ipdump .ge. 3) then
            write (lprint, 5210) (nacomb(j), j = 1, nmetab)
 5210       format (//20x, 'Prior matrix for concentration ratios'//
     1             (8x, (10(6x, a6))))
            do 210 jratio = 1, nratio_used
               write (lprint, 5220) nacomb(lmetab_prior(jratio)),
     1              (cprior(jratio, jmetab), jmetab = 1, nmetab)
 5220          format (/2x, a6, 1p10e12.3 / (8x, 1p10e12.3))
 210        continue
         else
            write (lprint, 5230) (chrato(jratio)(1:132),
     1                            chrato(jratio)(133:264),
     2                            jratio = 1, nratio_used)
 5230       format (//10x, 'CHRATO ratio priors used' // (a132))
         end if
      end if
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine parse_prior ()
c
c     Parse CHRATO input strings into CHRATI, CHRATW, EXRATI & SDRATI.
c     CHRATO must be of form:
c CHRATI = EXRATI +- SDRATI +WT= CHRATW
c
      INCLUDE 'lcmodel.inc'
      external ilen
      character chreturn*(mchratio)
c
      chsubp = 'PPRIOR'
      do 210 jratio = 1, nratio
         len_chrato = ilen(chrato(jratio))
         istart = 1
         call get_field ('=', 1, 1, 0,
     1                   chrati(jratio), freturn,
     2                   istart, len_chrato, chrato(jratio))
         if (istart .le. 0) then
c           -------------------------------------------------------------------
c           This can be caused by an empty string; i.e., no input for CHRATO.
c           -------------------------------------------------------------------
            ierr = 0
            go to 800
         end if
         call get_field ('+-', 2, 2, 0,
     1                   chreturn, exrati(jratio),
     2                   istart, len_chrato, chrato(jratio))
         if (istart .le. 0) then
            ierr = 1
            go to 800
         end if
         call get_field ('+WT=', 4, 2, 1,
     1                   chreturn, sdrati(jratio),
     2                   istart, len_chrato, chrato(jratio))
         if (istart .gt. len_chrato) go to 210
         if (istart .lt. 0) then
            ierr = 2
            go to 800
         end if
         call get_field (' ', 0, 1, 2,
     1                   chratw(jratio), freturn,
     2                   istart, len_chrato, chrato(jratio))
         if (istart .le. 0) then
            ierr = 3
            go to 800
         end if
 210  continue
      return
 800  if (lprint .gt. 0) write (lprint, 5810)
     1     chrato(jratio)(1:132), chrato(jratio)(133:264), istart
 5810 format ('Incorrect CHRATO follows:', / a132 / a132 /
     1        'ISTART =', i3)
      call errmes (100*ierr + jratio, 4, chsubp)
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_field (chseparator, len_chseparator, ifield_type,
     2                      iatend, chreturn, freturn,
     3                      istart, len_string_in, string_in)
c
c     Extracts field from STRING_IN into
c CHRETURN (if IFIELD_TYPE = 1)
C FRETURN  (if IFIELD_TYPE = 2)
c
c     STRING_IN, the string being searched is always a full string in the
c calling routine (e.g., CHRATO).  Substrings of it are explicitly defined
c here in GET_FIELD with ISTART, IEND_FIELD, etc.  This avoids possible
c compiler problems with substrings as actual arguments and INDEX.
c
      implicit none
      character chseparator*(*), chreturn*(*), chsubp*6, fmt*7,
     1          string_in*(*)
      integer i, iatend, iend_field, ifield_type, istart, istart_field,
     3        len_chseparator, len_string_in
      real freturn
c
      chsubp = 'GFIELD'
      iend_field = 0
c     -------------------------------------------------------------------------
c     Locate field in STRING_IN
c     -------------------------------------------------------------------------
      if (len_string_in .lt. istart) call errmes(1, 5, chsubp)
      if (iatend .eq. 0) then
c        ----------------------------------------------------------------------
c        CHSEPARATOR must be present.
c        ----------------------------------------------------------------------
         iend_field = index(string_in(istart:len_string_in),
     1                      chseparator(:len_chseparator)) - 2 + istart
         if (iend_field .lt. istart) then
c           -------------------------------------------------------------------
c           Error -- separator absent.
c           -------------------------------------------------------------------
            istart = -1
            return
         end if
      else if (iatend .eq. 1) then
c        ----------------------------------------------------------------------
c        CHSEPARATOR may be present, but the field may also be at the end.
c        ----------------------------------------------------------------------
         iend_field = index(string_in(istart:len_string_in),
     1                      chseparator(:len_chseparator)) - 2 + istart
         if (iend_field .lt. istart) iend_field = len_string_in
      else if (iatend .eq. 2) then
c        ----------------------------------------------------------------------
c        The field is at the end; CHSEPARATOR is not present.
c        ----------------------------------------------------------------------
         iend_field = len_string_in
      else
         call errmes(2, 5, chsubp)
      end if
c     -------------------------------------------------------------------------
c     Remove leading white space.
c     -------------------------------------------------------------------------
      do 210 istart_field = istart, iend_field
         if (string_in(istart_field:istart_field) .ne. ' ') go to 220
 210  continue
c     -------------------------------------------------------------------------
c     Error -- entire field is blank.
c     -------------------------------------------------------------------------
      istart = -2
      return
c     -------------------------------------------------------------------------
c     Increment ISTART (since IEND_FIELD is going to be changed below).
c     -------------------------------------------------------------------------
 220  istart = iend_field + len_chseparator + 1
c     -------------------------------------------------------------------------
c     Remove trailing white space.
c     -------------------------------------------------------------------------
      istart = iend_field + len_chseparator + 1
      i = iend_field
      do 230 iend_field = i, istart_field, -1
         if (string_in(iend_field:iend_field) .ne. ' ') go to 240
 230  continue
 240  if (ifield_type .eq. 1) then
         chreturn = string_in(istart_field:iend_field)
      else if (ifield_type .eq. 2) then
         i = iend_field - istart_field + 1
         if (i .le. 9) then
            write (fmt, 5240) i
 5240       format ('(f', i1, '.0)')
         else if (i .le. 99) then
            write (fmt, 5250) i
 5250       format ('(f', i2, '.0)')
         else
            istart = -3
            return
         end if
         read (string_in(istart_field:iend_field), fmt, err=810,
     1         end=810) freturn
      else
         istart = -4
         return
      end if
      return
 810  istart = -5
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine parse_sum(exrati_arg, substring, len_substring,
     1                     lratio, csum, denom_absent)
c
c     Parse CHRATI or CHRATW, and compute (unweighted) CPRIOR elements and
c CSUM.  CPRIOR will have to be multiplied by the weight (1/SD*CSUM).
C     If no metabolites are present from the denominator, then CPRIOR will
c have all zeroes, except for the numerator, which will then have an
c expectation of zero.
c     If CSUM=0, then the original constraint with the denominator terms will
c hold, but the weight (presumably large) will have to be set in calling
c program.
c
      INCLUDE 'lcmodel.inc'
      character chterm*(mchmet), substring*(*)
      logical atend, denom_absent, wildcard
      CHSUBP='PARSUM'
      istart = 1
c      write (lprint, 9000) substring
c 9000 format (a)
      do 110 jterm = 1, mmet_ratio - 1
         LENGTH=INDEX(substring(ISTART:len_substring),'+')-1
c         write (lprint, 9110) length
c 9110    format (' length =', i3)
         ATEND=LENGTH .EQ. -1
         IF (ATEND) LENGTH=INDEX(substring(ISTART:len_substring),' ')-1
C        ----------------------------------------------------------------------
C        Blank in position 1 or right after last "+"
C        ----------------------------------------------------------------------
         IF (LENGTH .EQ. 0) return
c        ----------------------------------------------------------------------
c        CHRATI has (MCHMET+1)*mmet_ratio + 1.  Thus, there should
c          be room for at least one blank at the end, if <= MMET_RATIO
c          terms are used.  If not, abort (probably because there are
c          > MMET_RATIO terms.
c        ----------------------------------------------------------------------
         IF (LENGTH .EQ. -1) then
            IF (LPRINT .GT. 0) WRITE (LPRINT,5110)
     1        substring(:len_substring)
 5110       FORMAT (' Too long a sum in CHRATO, ending as follows:'/
     1              a)
            call errmes (1, 4, chsubp)
         end if
         IF (LENGTH .GT. MCHMET) THEN
            IF (LPRINT .GT. 0) WRITE (LPRINT,5120)
     1        substring(:len_substring)
 5120       FORMAT
     1           (' Too long a Metabolite Name in CHRATO ',
     2            'ending as follows:' / A)
            CALL ERRMES (2, 4, CHSUBP)
         END IF
         wildcard = substring(istart+length-1:istart+length-1) .eq. '*'
         do 130 jmetab = 1, nmetab
            factor = 0.
            if (wildcard) then
               if (length .eq. 1) then
                  factor = 1.
               else
                  if (substring(istart:istart+length-2) .eq.
     1                nacomb(jmetab)(1:length-1)) factor = 1.
               end if
            else
               chterm = substring(istart:istart+length-1)
               if (chterm .eq. nacomb(jmetab)) factor = 1.
               if (chterm .eq. 'totCho'   .and.
     1             (nacomb(jmetab) .eq. 'Cho'   .or.
     2              nacomb(jmetab) .eq. 'GPC'   .or.
     3              nacomb(jmetab) .eq. 'PCh')) factor = 1.
               if (chterm .eq. 'totCr'   .and.
     1             (nacomb(jmetab) .eq. 'Cr'   .or.
     2              nacomb(jmetab) .eq. 'Cre'   .or.
     3              nacomb(jmetab) .eq. 'PCr')) factor = 1.
               if (chterm .eq. 'totNAA'   .and.
     1             (nacomb(jmetab) .eq. 'NAA'   .or.
     2              nacomb(jmetab) .eq. 'NAAG')) factor = 1.
               if (chterm .eq. 'Big3') then
                  if (nacomb(jmetab) .eq. 'Cr'   .or.
     1                nacomb(jmetab) .eq. 'Cre'   .or.
     2                nacomb(jmetab) .eq. 'PCr'   .or.
     3                nacomb(jmetab) .eq. 'NAA'   .or.
     4                nacomb(jmetab) .eq. 'NAAG') then
                     factor = 1.
                  else if (nacomb(jmetab) .eq. 'Cho'   .or.
     1                     nacomb(jmetab) .eq. 'GPC'   .or.
     2                     nacomb(jmetab) .eq. 'PCh') then
                     factor = 3.
                  end if
               end if
            end if
            if (factor .gt. 0.) then
               denom_absent = .false.
               csum = csum + factor * sngl(solbes(jmetab, 1))
               if (exrati_arg .gt. 0.) cprior(lratio, jmetab) =
     1                                 cprior(lratio, jmetab) -
     2                                 factor * exrati_arg
            end if
 130     continue
         if (atend) return
         istart = istart + length + 1
 110  continue
      IF (LPRINT .GT. 0) WRITE (LPRINT,5130) substring(:len_substring)
 5130 FORMAT (' Too many metabolites in CHRATI or CHRATW ',
     1        'ending as follows:' / A)
      call errmes (3, 4, chsubp)
      end
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      LOGICAL FUNCTION LDEGMX (IDEGMX)
C
C  LDEGMX = T if the max. phase correction to the data would exceed
C             DEGMAX(IDEGMX)
C
      INCLUDE 'lcmodel.inc'
      PHACOR(ppmarg) = ABS(PHIZER + (ppmarg - ppmcen) * PHIONE)
      CHSUBP='LDEGMX'
      IF ((IDEGMX.NE.1 .AND. IDEGMX.NE.2) .OR. RADIAN.LE.0.) CALL
     1      ERRMES (1, 5, CHSUBP)
      IF (DEGMAX(IDEGMX) .LE. 0.) CALL ERRMES (2, 4, CHSUBP)
      PHIZER=PARBES(LPHAST,2)
      PHIONE=PARBES(LPHAST+1,2)
      ppmsig_max = amax1(ppmsig(1), ppmsig(2))
      ppmsig_min = amin1(ppmsig(1), ppmsig(2))
      ppmmax = amin1(ppmsig_max, ppm(1))
      ppmmin = amax1(ppmsig_min, ppm(nyuse))
      LDEGMX=AMAX1(PHACOR(ppmmax),PHACOR(ppmmin)) .GT.
     1       DEGMAX(IDEGMX)*RADIAN
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER FUNCTION NEXTRE (PARNL, NSIDE2, DPY, DGAUSS, THRLIN,
     1                         imethd)
C  NEXTRE = no. of extrema, assuming one zero point outside of each boundary
C           and that the point between PARNL(NSIDE2/2) and PARNL(NSIDE2/2+1)
C           is 1-(sum over the NSIDE2 PARNL).
C  Lineshape is first smoothed with DGAUSS.
C  An extremum is not counted if the absolute values of it and one of its
C     neighboring extrema are less than
C     THRLIN*(maximum of the absolute value of all points); i.e., if it is
C     just a notch near the noisy Baseline.
C
      DOUBLE PRECISION PARNL(NSIDE2), DPY(NSIDE2+3,2), DSUM,
     1                 DGAUSS(0:NSIDE2+2), DTHR
      if (imethd .eq. 2) then
         nextre = 99
         return
      end if
      DPY(1,1)=0.D0
      DPY(NSIDE2+3,1)=0.D0
      DSUM=1.D0
      DO 110 JPAR=1,NSIDE2
         DSUM=DSUM-PARNL(JPAR)
  110 CONTINUE
      IMIDDL=NSIDE2/2+2
      JPAR=0
      DO 120 JPY=2,NSIDE2+2
         IF (JPY .EQ. IMIDDL) THEN
            DPY(JPY,1)=DSUM
        ELSE
            JPAR=JPAR+1
            DPY(JPY,1)=PARNL(JPAR)
         END IF
  120 CONTINUE
      DTHR=0.D0
      DO 150 JPY=1,NSIDE2+3
         DSUM=0.D0
         DO 160 KPY=1,NSIDE2+3
            DSUM=DSUM+DPY(KPY,1)*DGAUSS(IABS(KPY-JPY))
  160    CONTINUE
         DPY(JPY,2)=DSUM
         DTHR=DMAX1(DTHR, DABS(DSUM))
  150 CONTINUE
      DTHR=DTHR*THRLIN
C      WRITE (6,5150) DTHR, (DPY(J,2),J=1,NSIDE2+3)!DELETE
C 5150 FORMAT (/' Threshold =',F7.4/(1X,10F12.4))!DELETE
      KEXTR=1
      DPY(1,1)=0.D0
      DO 210 JPAR=2,NSIDE2+2
         IF ((DPY(JPAR,2)-DPY(JPAR-1,2))*(DPY(JPAR+1,2)-DPY(JPAR,2))
     1       .LT. 0.D0) THEN
            KEXTR=KEXTR+1
            DPY(KEXTR,1)=DABS(DPY(JPAR,2))
         END IF
  210 CONTINUE
      DPY(KEXTR+1,1)=0.D0
      NEXTRE=0
      DO 220 JEXTR=2,KEXTR
         IF (DMAX1(DPY(JEXTR,1), DMIN1(DPY(JEXTR+1,1),
     1                                 DPY(JEXTR-1,1))) .GT. DTHR)
     2         NEXTRE=NEXTRE+1
  220 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER FUNCTION INFLEC (PARNL, NSIDE2, DPY, DGAUSS, THRLIN,
     1                         imethd)
C
C  INFLEC = no. of inflection points assuming 2 zero points outside of each
C           boundary and that the point between PARNL(NSIDE2/2) and
C           PARNL(NSIDE2/2+1) is 1-(sum over the NSIDE2 PARNL).
C  Lineshape is first smoothed with DGAUSS.
C  Inflections are only counted if smoothed point exceeds THRLIN*(maximum of
C    the abolute values of all points).
C
      DOUBLE PRECISION PARNL(NSIDE2), DPY(NSIDE2+5,2), DSUM, DEL,
     1                 DELOLD, DGAUSS(0:NSIDE2+4), DTHR
      if (imethd .eq. 2) then
         inflec = 99
         return
      end if
      DPY(1,1)=0.D0
      DPY(2,1)=0.D0
      DPY(NSIDE2+4,1)=0.D0
      DPY(NSIDE2+5,1)=0.D0
      DSUM=1.D0
      DO 110 JPAR=1,NSIDE2
         DSUM=DSUM-PARNL(JPAR)
  110 CONTINUE
      IMIDDL=NSIDE2/2+3
      JPAR=0
      DO 120 JPY=3,NSIDE2+3
         IF (JPY .EQ. IMIDDL) THEN
            DPY(JPY,1)=DSUM
         ELSE
            JPAR=JPAR+1
            DPY(JPY,1)=PARNL(JPAR)
         END IF
  120 CONTINUE
      DTHR=0.D0
      DO 150 JPY=1,NSIDE2+5
         DSUM=0.D0
         DO 160 KPY=1,NSIDE2+5
            DSUM=DSUM+DPY(KPY,1)*DGAUSS(IABS(KPY-JPY))
  160    CONTINUE
         DPY(JPY,2)=DSUM
         DTHR=DMAX1(DTHR, DABS(DSUM))
  150 CONTINUE
      DTHR=DTHR*THRLIN
      INFLEC=0
      DELOLD=DPY(3,2)
      DO 210 JPAR=3,NSIDE2+4
         DEL=DPY(JPAR-1,2)-2.D0*DPY(JPAR,2)+DPY(JPAR+1,2)
         IF (DEL*DELOLD.LT.0.D0 .AND. DABS(DPY(JPAR,2)).GT.DTHR)
     1     INFLEC=INFLEC+1
         DELOLD=DEL
  210 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SSRANG (IRANGE)
C
C  Set parameters for Regula Falsi search.
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='SSRANG'
      IF (IRANGE.LT.1 .OR. IRANGE.GT.3) CALL ERRMES (1, 5, CHSUBP)
      SSQMIN=FSHSSQ (PRMNMX(1,IRANGE),
     1               0, NYuse, FLOAT(NDFREF), SSQREF,
     2                      LPRINT, RRANGE)
      SSQMAX=FSHSSQ (PRMNMX(2,IRANGE),
     1               0, NYuse, FLOAT(NDFREF), SSQREF,
     2                      LPRINT, RRANGE)
 5210 FORMAT (//' SSQREF =',1PE14.6,'   SSQMIN =',E14.6,'   SSQMAX =',
     1        E14.6,'   PREJMN =',0PF7.4,'   PREJMX =',F7.4//)
      IF (LPRINT .GT. 0) WRITE (LPRINT,5210) SSQREF, SSQMIN, SSQMAX,
     1                            PRMNMX(1,IRANGE), PRMNMX(2,IRANGE)
      SSQAIM=.5*(SSQMIN+SSQMAX)
      IF (SSQMIN .GE. SSQMAX) CALL ERRMES (2, 4, CHSUBP)
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RFALSI (IALPHA, IRANGE, LREPHA, ALPHB, ALPHS, ASSQLO,
     1                   AALPLO, ASSQHI, AALPHI, AALPHA, PREJOK, PREJ1)
C
C  Regula falsi search for an ALPHA with SSQ between SSQMIN and SSQMAX, where
C    ALPHAB = ALPHB*ALPHA when IALPHA=1 or 3
C    ALPHAS = ALPHS*ALPHA when IALPHA=2 or 3
C
C  IALPHA = 1 keeping ALPHAS fixed.
C         = 2 keeping ALPHAB fixed.
C         = 3 for normal search varying ALPHAB and ALPHAS.
c
c  IRANGE: PRMNMX(*, IRANGE) will be used.
C
C  LREPHA = T to rephase and get Reference Solution.
C
C  Needs starting values from previous analysis and call to SAVBES(2).
C
      INCLUDE 'lcmodel.inc'
      DOUBLE PRECISION ALPBOL, ALPBUS, alpha_alphab, alpha_alphas,
     1                 ALPHB, ALPHS, ALPSOL, ALPSUS, PENLTY
      LOGICAL LREPHA, PREJOK
C      COMPLEX CYFIT(MY)!CT2
C      LOGICAL LTEST!CT2
C      COMMON /BLT1/ CYFIT, LTEST!CT2
 
      CHSUBP='RFALSI'
      IF (MIN0(IRANGE,IALPHA).LT.1 .OR. MAX0(IRANGE,IALPHA).GT.3) CALL
     1      ERRMES (1, 5, CHSUBP)
      IF (AMIN1(RALIMN,RALINC) .LE. 1.) CALL ERRMES (2, 4, CHSUBP)
      is_alpbmn = .false.
      ALPBUS=DMIN1(ALPBMX,DMAX1(ALPBMN,ALPHB))
      ALPSUS=DMIN1(ALPSMX,DMAX1(ALPSMN,ALPHS))
      SSQLO=ASSQLO
      ALPLO=AALPLO
      SSQHI=ASSQHI
      ALPHI=AALPHI
      ALPHA=AALPHA
      IF (LREPHA) THEN
C        ----------------------------------------------------------------------
C        Rephase data and get Reference Solution.
C        ----------------------------------------------------------------------
         IF (LPRINT .GT. 0) WRITE (LPRINT,5110)
 5110    FORMAT (//////20X, 'Reference Solution for rephased data')
         CALL REPHAS ()
         ALPHAB=ALPBMN
         ALPHAS=ALPSMN
C         LTEST=.TRUE.!CT2
         CALL PLINLS (2, IERROR)
         IF (OBJECT .GE. DRANGE) CALL ERRMES (3, 4, CHSUBP)
C!CT2
C         DO 2110 JNONL=1,NNONL!CT2
C            PARNLN(JNONL)=PARBES(JNONL,2)!CT2
C 2110    CONTINUE!CT2
C         DO 2120 JY=1,NY!CT2
C            CY(JY)=CYFIT(JY)!CT2
C 2120    CONTINUE!CT2
C         CALL PLINLS (2, IERROR)!CT2
C         IF (LTEST) STOP!CT2
C!CT2
         INISOL=.FALSE.
         SSQREF=DSSQ
         NDFREF=NDF
         CALL SSRANG (IRANGE)
      END IF
      DAIMBS=RRANGE
      NINFLE=99
      PREJOK=.FALSE.
      ALPHAB=ALPBUS
      ALPHAS=ALPSUS
      DO 310 JTRY=1,MFNDAL
         ALPBOL=ALPHAB
         ALPSOL=ALPHAS
         IF (IALPHA .EQ. 1) THEN
            ALPHAB=DMAX1(ALPBMN,DMIN1(ALPBMX,ALPBUS*ALPHA))
c           ------------------------------------------------------------------
c           ALPHA must be adjusted (for next ALPHA interplolation) to agree
c              with any limit imposed by ALPBMN or ALPBMX.
c           The factors .99999 & 1.00001 allow for rounding error in
c              calculation of ALPHA.
c           ------------------------------------------------------------------
            alpha = alphab / alpbus
            IF (JTRY .GT. 1) THEN
               IF (DMIN1(ALPBOL,ALPHAB) .GE. .99999*ALPBMX) THEN
                  if (.not.nobase) CALL ERRMES (4, 1, CHSUBP)
                  GO TO 700
               ELSE IF (DMAX1(ALPBOL,ALPHAB) .LE. 1.00001*ALPBMN) THEN
                  CALL ERRMES (5, 1, CHSUBP)
                  is_alpbmn = .true.
                  GO TO 700
               END IF
            END IF
         ELSE IF (IALPHA .EQ. 2) THEN
            ALPHAS=DMAX1(ALPSMN,DMIN1(ALPSMX,ALPSUS*ALPHA))
c           ------------------------------------------------------------------
c           ALPHA must be adjusted (for next ALPHA interplolation) to agree
c              with any limit imposed by ALPSMN or ALPSMX.
c           ------------------------------------------------------------------
            alpha = alphas / alpsus
            IF (JTRY .GT. 1) THEN
               IF (DMIN1(ALPSOL,ALPHAS) .GE. .99999*ALPSMX) THEN
                  CALL ERRMES (6, 1, CHSUBP)
                  GO TO 700
               ELSE IF (DMAX1(ALPSOL,ALPHAS) .LE. 1.00001*ALPSMN) THEN
                  CALL ERRMES (7, 2, CHSUBP)
                  GO TO 700
               END IF
            END IF
         ELSE
            ALPHAB=DMAX1(ALPBMN,DMIN1(ALPBMX,ALPBUS*ALPHA))
            ALPHAS=DMAX1(ALPSMN,DMIN1(ALPSMX,ALPSUS*ALPHA))
c           ------------------------------------------------------------------
c           The limits imposed by ALP*MN or ALP*MX produce two ALPHAs, at least
c              as close to 1 (on log scale) as the original ALPHA.  Readjust
c              ALPHA (for interpolation later) to the one furthest from 1 (on
c              log scale), i.e., closest to the original ALPHA.
c           ------------------------------------------------------------------
            alpha_alphab = alphab/alpbus
            alpha_alphas = alphas/alpsus
            if (dabs(dlog(alpha_alphab)) .gt.
     1          dabs(dlog(alpha_alphas))) then
               alpha = alpha_alphab
            else
               alpha = alpha_alphas
            end if
            IF (JTRY .GT. 1) THEN
               IF (DMIN1(ALPBOL,ALPHAB).GE. .99999*ALPBMX .AND.
     1             DMIN1(ALPSOL,ALPHAS).GE. .99999*ALPSMX) THEN
C                 -------------------------------------------------------------
C                 PREJOK = T for special case that both ALPHAs are max, and
C                            SSQ<SSQMAX nevertheless (or, below, when
C                            unimodality has already been achieved during
C                            increase of ALPHAS).
C                 -------------------------------------------------------------
                  PREJOK=SSQ .LE. SSQMAX
                  CALL ERRMES (8, 1, CHSUBP)
                  GO TO 700
               END IF
               IF (DMAX1(ALPBOL,ALPHAB).LE. 1.00001*ALPBMN .AND.
     1              DMAX1(ALPSOL,ALPHAS).LE. 1.00001*ALPSMN) THEN
                  CALL ERRMES (9, 2, CHSUBP)
                  GO TO 700
               END IF
            END IF
         END IF
         IF (LREPHA) THEN
            CALL PLINLS (2, IERROR)
         ELSE
            CALL PLINLS (3, IERROR)
         END IF
         SSQ=DSSQ
         IF (OBJECT .LT. DRANGE) THEN
            NINFLE=INFLEC (PARNLN, NSIDE2, DPY, DGAUSS, THRLIN, imethd)
            NEXT=NEXTRE (PARNLN, NSIDE2, DPY, DGAUSS, THRLIN, imethd)
            IF (MDALPB.LE.0 .OR. USINFL) THEN
               ITEST=NINFLE-2
            ELSE
               ITEST=NEXT-1
            END IF
            IF (ITEST.LE.0 .AND. SSQ.LE.SSQMAX) THEN
C              ----------------------------------------------------------------
C              The peak is smooth and the fit is good; so set
C              PREJOK = T for special exit to stop increasing ALPHAS
C                         if only ALPHAS is being increased (IALPHA=2) or
C                         if ALPHAB=ALPBMX during the normal search varying
C                           ALPHAB & ALPHAS (IALPHA=3), or when NSIDES=1
c                           (i.e. only ALPHAB will be varied).
C              ----------------------------------------------------------------
               PREJOK=IALPHA.EQ.2 .OR.
     1                (ALPHAB.GE. .99999*ALPBMX   .and.
     2                 (IALPHA.EQ.3   .or.   nsides .le. 1))
C               if (ALPHAB.GE.ALPBMX .AND. IALPHA.EQ.3) call errmes (99,
C     1                                                   1, CHSUBP)
            END IF
            if (imethd .eq. 2   .and.   ialpha .eq. 2   .and.
     1          alphas .ge. .99999 * alpsmx   .and.
     2          ssq .le. ssqmax) then
               prejok = .true.
               call errmes (18, 1, chsubp)
            end if
            RTERM=ABS(SSQ-SSQAIM)
            IF (RTERM .LT. DAIMBS) THEN
               DAIMBS=RTERM
               CALL SAVBES (1)
            END IF
         ELSE
            CALL ERRMES (10, 3, CHSUBP)
            ALPHA=ALPHA*RALINC
            GO TO 310
         END IF
         IF (LPRINT .GT. 0) THEN
            IF (SDREF**2 .GT. 0.D0) THEN
               PENB=PENLTY (ALPHAB, 0.D0, SOLUTN, PARNLN)/SDREF**2
               PENS=PENLTY (0.D0, ALPHAS, SOLUTN, PARNLN)/SDREF**2
               IF (PENS .GT. 0.) THEN
                  RPEN=PENB/PENS
               ELSE
                  RPEN=RRANGE
               END IF
               PEN=(OBJECT-DSSQ)/SDREF**2
               IF (SNGL(ALPHAB) .LT. ALPBPN) PEN=PEN+PNALPB
            ELSE
               PENB=RRANGE
               PENS=RRANGE
               RPEN=RRANGE
               PEN=RRANGE
            END IF
            WRITE (LPRINT,5310) JTRY, ALPHAB, ALPHAS,
     1                           NINFLE, SSQLO, SSQ, SSQHI,
     2                           PENB, PENS, RPEN, NEXT,
     3                           PEN, OBJECT-DSSQ
 5310       FORMAT (/' ITER =',I2,'     ALPHAB =', 1PE10.3, 3X,
     1      'ALPHAS =', E10.3, 6X, 'Ninfle =', I2, 6X, '(SSQ', E14.6,
     2      ' <', E14.6,' <',E14.6,')'/
     3      14X, 'penB/penS =', E9.2, ' /', E9.2, ' =', E9.2, 3X,
     4      'Nextre =', I2, /
     5      14X, 'Penalty =', E12.4, 24X,
     6      'Penalty*SDREF**2 (w/o prior) =', E10.2)
            if (iter_dump .eq. jtry   .and.
     1           dabs(alphab - alphab_dump) .le. 1.d-3 * alphab   .and.
     2           dabs(alphas - alphas_dump) .le. 1.d-3 * alphas) then
               call savbes(1)
               call savbes(2)
               call finout ()
               stop
            end if
         END IF
         IF (PREJOK) THEN
C           -------------------------------------------------------------------
C           Special exit when either of the above 2 criteria (before LPRINT)
c              are satisfied.
C           -------------------------------------------------------------------
            CALL SAVBES (1)
            GO TO 700
         END IF
         ALPOLD=ALPHA
         IF (SSQ .LE. AMAX1(SSQREF,SSQLO)) THEN
            IF (SSQ .GE. SSQREF) THEN
C              ----------------------------------------------------------------
C              This error is probably due to a local min. in a preceding
C                solution.
C              ----------------------------------------------------------------
               CALL ERRMES (11, 1, CHSUBP)
            ELSE
C              ----------------------------------------------------------------
C              This error is probably due to a local min. in the Reference
C                Solution.  Redefine the present solution as the Reference
C                Solution.
C              ----------------------------------------------------------------
               CALL ERRMES (12, 1, CHSUBP)
               SSQREF=SSQ
               NDFREF=NDF
               CALL SSRANG (IRANGE)
            END IF
            SSQLO=SSQ
            ALPLO=ALPHA
            IF (ALPHI .LT. RRANGE) THEN
               ALPHA=AMIN1(.5*(ALPHA+ALPHI), ALPHA*RALINC)
            ELSE
               ALPHA=ALPHA*RALINC
            END IF
            GO TO 320
         END IF
         IF (SSQ .GE. SSQHI) THEN
C           -------------------------------------------------------------------
C           This error is probably due to round-off or to a local min. in this
C             solution.
C           -------------------------------------------------------------------
            ALPHA=AMAX1(.5*(ALPHA+ALPLO), ALPHA/RALINC)
            CALL ERRMES (13, 1, CHSUBP)
            GO TO 320
         END IF
         IF (SSQ .GT. SSQMAX) THEN
            IF (SSQHI-SSQAIM .LE. SSQAIM-SSQLO) THEN
C              ----------------------------------------------------------------
C              SSQLO is too far away for interpolation.  Extrapolate using
C                SSQHI.
C              ----------------------------------------------------------------
               ALPHA=AMAX1(ALPHA/RALINC,
     1                 ALPHI-(ALPHI-ALPHA)*(SSQHI-SSQAIM)/(SSQHI-SSQ))
C              ----------------------------------------------------------------
C              If extrapolation goes back behind ALPLO, then take midpoint.
C              ----------------------------------------------------------------
               IF (ALPHA .LE. ALPLO) ALPHA=.5*(ALPOLD+ALPLO)
            ELSE IF (SSQLO .LT. 0.) THEN
               ALPHA=ALPHA/RALINC
            ELSE
C              ----------------------------------------------------------------
C              Normal interpolation.
C              ----------------------------------------------------------------
               ALPHA=ALPLO+(ALPHA-ALPLO)*(SSQAIM-SSQLO)/(SSQ-SSQLO)
            END IF
            ALPHI=ALPOLD
            SSQHI=SSQ
         ELSE IF (SSQ .LT. SSQMIN) THEN
            IF (SSQHI .GE. RRANGE) THEN
               IF (SSQLO .LT. 0.) THEN
                  ALPHA=ALPHA*RALINC
               ELSE
C                 -------------------------------------------------------------
C                 Limited extrapolation.
C                 -------------------------------------------------------------
                  test = ALPLO + (ALPHA - ALPLO) *
     1                           (SSQAIM - SSQLO) / (SSQ - SSQLO)
                  alpha = amin1(test, xtrpmx * alpha)
               END IF
            ELSE IF (SSQAIM-SSQLO.LE.SSQHI-SSQAIM .AND. SSQLO.GT.0.)
     1                 THEN
C              ----------------------------------------------------------------
C              SSQHI is too far away for interpolation.  Extrapolate using
C                SSQLO.
C              ----------------------------------------------------------------
               ALPHA=ALPLO+(ALPHA-ALPLO)*(SSQAIM-SSQLO)/(SSQ-SSQLO)
               IF (ALPHA .GE. ALPHI) THEN
C                 -------------------------------------------------------------
C                 Limit extrapolation beyond ALPHI.
C                 Remove ALPHI as upper bound.  It could have been a local
C                   min. that has since been improved.
C                 -------------------------------------------------------------
                  ALPHA=AMIN1(ALPHA,ALPOLD*RALINC)
                  ALPHI=RRANGE
                  SSQHI=RRANGE
               END IF
            ELSE
C              ----------------------------------------------------------------
C              Normal interpolation.
C              ----------------------------------------------------------------
               ALPHA=ALPHA+(ALPHI-ALPHA)*(SSQAIM-SSQ)/(SSQHI-SSQ)
            END IF
            ALPLO=ALPOLD
            SSQLO=SSQ
         ELSE
            GO TO 700
         END IF
  320    IF (ALPHA/ALPOLD.GT.1. .AND. ALPHA/ALPOLD.LT.RALIMN) THEN
C           -------------------------------------------------------------------
C           This too small an increase in ALPHA could be due to too low an
C             upper bound; remove it.
C           -------------------------------------------------------------------
            SSQHI=RRANGE
            ALPHI=RRANGE
         ELSE IF (ALPOLD/ALPHA.GT.1. .AND. ALPOLD/ALPHA.LT.RALIMN) THEN
C           -------------------------------------------------------------------
C           This too small a decrease in ALPHA could be due to too high a
C             lower bound; remove it.
C           -------------------------------------------------------------------
            SSQLO=-1./RRANGE
            ALPLO=0.
         END IF
  310 CONTINUE
      CALL ERRMES (14, 2, CHSUBP)
  700 IF (DAIMBS .GE. RRANGE) CALL ERRMES (15, 4, CHSUBP)
      IF (SSQREF*FLOAT(NDFREF) .LE. 0.) THEN
         CALL ERRMES (16, 3, CHSUBP)
         PREJ1=-1.
      ELSE
         PREJ1=FISHNI(AMAX1(0.,(SSQBES(1)-SSQREF))*FLOAT(NYuse-NDFREF)/
     1                (SSQREF*FLOAT(NDFREF)), FLOAT(NDFREF),
     2                FLOAT(NYuse-NDFREF), LPRINT)
         IF (LPRINT .GT. 0) WRITE (LPRINT,5710) PREJ1
 5710    FORMAT (/' Probability to reject =',F7.4)
         PREJOK=PREJOK .OR. (PREJ1.GE..99*PRMNMX(1,IRANGE) .AND.
     1                       PREJ1.LE.1.01*PRMNMX(2,IRANGE))
C        ----------------------------------------------------------------------
C        PREJOK = F together with an alpha = alpha_max is probably harmless,
C                   e.g., when ALPHAB=ALPBMX and ALPHAS=ALPSMN at start of
C                   grid scan.
C        ----------------------------------------------------------------------
         IF (.NOT.PREJOK .AND. ALPHAB.LT..99999*ALPBMX .AND.
     1       ALPHAB.GT.1.00001*ALPBMN) CALL ERRMES (17, 2, CHSUBP)
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DOUBLE PRECISION FUNCTION PENLTY (ALPB, ALPS, SOL, PARNL)
C
C  PENLTY = regularizor penalty for ALPHAB=ALPB and ALPHAS=ALPS multiplied by
C           weighting factor SDREF.
C         = 0 if SDREF has not yet been computed.
C  REGF is assumed to be for the lineshape coefficients with the center
C       coefficient eliminated by the equality constraint.  Therefore, the 1st
C       3 rows are assumed to have the constants 2, -1, -1 on the rhs, rather
C       than 0.
C
      INCLUDE 'lcmodel.inc'
      DOUBLE PRECISION ALPB, ALPS, PARNL(MNONL), SOL(MPAR)
      CHSUBP='PENLTY'
      PENLTY=0.D0
      IF (SDREF .GE. DRANGE) RETURN
      IF (ALPB .GT. 0.D0) THEN
         DO 110 IROW=1,NBACKG
            DTERM(1)=0.D0
            DO 120 ICOL=1,NBACKG
               DTERM(1)=DTERM(1)+REGB(IROW,ICOL)*SOL(NMETAB+ICOL)
  120       CONTINUE
            PENLTY=PENLTY+(ALPB*DTERM(1))**2
  110    CONTINUE
      END IF
      IF (ALPS .GT. 0.D0) THEN
         if (imethd .eq. 2) then
            do 205 jmetab = 1, nmetab
               jnonl = lrt2st - 1 + jmetab
c              ---------------------------------------------------------------
c              CONC_EXPECT(JMETAB) = expectation of relative CONC
c              PI * HZPPPM * RT2MIN(JMETAB) = (Lorentzian) 1/T2 of
c                 of initial Gaussian (sharpness limit), SIFWMN.
c              RLRNTZ = extra factor (>1) to penalize Lorentzian broadening
c                       compared to component with initial Gaussian.
c              ---------------------------------------------------------------
               dterm(1) = alps * sol(jmetab) / conc_expect(jmetab)
               if (ipowrg .eq. 1) then
                  penlty = penlty + (dterm(1) * (rt2min(jmetab) +
     1                                        rlrntz * parnl(jnonl)))**2
               else
                  dterm(2) = (rt2min(jmetab) + rlrntz * parnl(jnonl))**2
                  penlty = penlty + (dterm(1) * dterm(2))**2
               end if
c               write(lprint, 9205) sdref, alps, sol(jmetab),
c     1            conc_expect(jmetab), rt2min(jmetab), parnl(jnonl),
c     2            dterm(1), dterm(2), penlty, jnonl, jmetab, nmetab
c 9205          format(1p9e10.2, 3i4)
 
 205        continue
         else
            DO 210 IROW=1,NREGF
               IF (IROW .EQ. 1) THEN
                  DTERM(1)=-2.D0
               ELSE IF (IROW .LE. 3) THEN
                  DTERM(1)=1.D0
               ELSE
                  DTERM(1)=0.D0
               END IF
               DO 220 ICOL=1,NSIDE2
                  DTERM(1)=DTERM(1)+REGF(IROW,ICOL)*PARNL(ICOL)
 220           CONTINUE
               PENLTY=PENLTY+(ALPS*DTERM(1))**2
 210        CONTINUE
         END IF
      END IF
      PENLTY=PENLTY*SDREF**2
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE REPHAS ()
C
C  Rephases CY using phases in PARBES(*,2)
C  Loads PARBES(*,2) into PARNLN
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='REPHAS'
      PHITOT(1)=PHITOT(1)+PARBES(LPHAST,2)
      PHITOT(2)=PHITOT(2)+PARBES(LPHAST+1,2)
      if (lprint .gt. 0) write(lprint, 5110) phitot(1)/radian,
     1                                       phitot(2)/radian
 5110 format(/'Rephasing to', f7.1, ' deg;   ', f8.2, ' deg/ppm'/)
      CTERM(1)=CEXP(CMPLX(0.,SNGL(PARBES(LPHAST,2)+
     1                            PARBES(LPHAST+1,2)*DELPPM(1))))
      CTERM(2)=CEXP(CMPLX(0.,-SNGL(PARBES(LPHAST+1,2)*PPMINC)))
      JDATA=LDATST
      DO 110 JY=1,NY
         CY(JY)=CY(JY)*CTERM(1)
         CTERM(1)=CTERM(1)*CTERM(2)
  110 CONTINUE
      PARBES(LPHAST,2)=0.D0
      PARBES(LPHAST+1,2)=0.D0
      DO 120 JNONL=1,NNONL
         PARNLN(JNONL)=PARBES(JNONL,2)
  120 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL FUNCTION FSHSSQ (PREJ,
     1                      IDFISH, NYuse, REFNDF, SSQREF,
     2                      LPRINT, RRANGE)
C
C  FSHSSQ = the weighted sum of squared deviations
C    corresponding to the Fisher Varaiance-Ratio for PREJ (the
C    probability to reject), REFNDF (the degrees of freedom of the
C    Reference Solution), NYuse, and SSQREF (the SSQ for the Reference
C    Solution).
C  IDFISH .NE. 0 to always dump.
C         .LT. 0 to stop after first dump.
C
      CHARACTER CHSUBP*6
      DOUBLE PRECISION DGAMLN
      CHSUBP='FSHSSQ'
      mtry = 20
      ptol = 5.e-4
      IF (AMIN1(PREJ,SSQREF).LE.0. .OR. PREJ.GE.1. .OR. REFNDF.LE.1.
     1    .OR. REFNDF.GE.FLOAT(NYuse-1)) CALL ERRMES (1, 4, CHSUBP)
C
C     Get an approximate value for the inverse of the F-distribution
C     using formulas from Abramowitz and Stegun (A&S).
C
C     Put the inverse of the complementary normal distribution function
C     corresponding to PREJ into YP.  A&S 26.2.23 is used.
C
      EXMAX=ALOG(RRANGE)
      PCOMP=1.-PREJ
      IF (PCOMP .LE. .5) THEN
         PSMALL=PCOMP
      ELSE
         PSMALL=PREJ
      END IF
      TT=SQRT(-2.*ALOG(PSMALL))
      YP=TT-(2.515517+TT*(.802853+TT*.010328))/(1.+TT*(1.432788+TT*
     1   (.189269+TT*.001308)))
      IF (PCOMP .GT. .5) YP=-YP
C
C     Use A&S 26.6.16 and 26.5.22 to put approximate F value in FF.
C
      R2AM1=1./(FLOAT(NYuse-1)-REFNDF)
      R2BM1=1./(REFNDF-1.)
      SLAMBD=(YP*YP-3.)/6.
      HH=2./(R2AM1+R2BM1)
      DUM=HH+SLAMBD
      IF (DUM .LE. 0.) THEN
C
C        This is a special case where the degrees of freedom are so few
C        that this approximation is very poor.  DUM is arbitrarily set
C        to 0., and an error message warns that the result can be very
C        inaccurate.
C
         DUM=0.
         CALL ERRMES (2, 3, CHSUBP)
      END IF
      WW=YP*SQRT(DUM)/HH-(R2BM1-R2AM1)*(SLAMBD+(5.-4./HH)/6.)
      IF (2.*WW .GE. EXMAX) CALL ERRMES (3, 4, CHSUBP)
      FF=EXP(2.*WW)
C
C     Use Newton's method to refine value of FF so that it's PREJ is
C     within PTOL of the required value of PREJ.
C
      DPBEST=RRANGE
      DF1=REFNDF
      HDF1=.5*DF1
      DF2=FLOAT(NYuse)-REFNDF
      HDF2=.5*DF2
      FACLOG=DGAMLN(DBLE(HDF1+HDF2)) - DGAMLN(DBLE(HDF1)) -
     1       DGAMLN(DBLE(HDF2)) + HDF1*ALOG(DF1) + HDF2*ALOG(DF2)
      EXP1=HDF1-1.
      EXP2=-HDF1-HDF2
      fbest = rrange
      DO 210 JTRY=1,MTRY
         PRTRY=FISHNI(FF,DF1,DF2,LPRINT)
         DPRTRY=PRTRY-PREJ
         IF (ABS(DPRTRY) .LE. DPBEST) THEN
            DPBEST=ABS(DPRTRY)
            FBEST=FF
            IF (DPBEST .LE. PTOL) GO TO 300
         END IF
         DPDF=EXP(FACLOG+EXP1*ALOG(FF)+EXP2*ALOG(DF2+DF1*FF))
         IF (DPDF .LE. 0.) THEN
            CALL ERRMES (4, 2, CHSUBP)
            GO TO 300
         END IF
         FF=FF-DPRTRY/DPDF
         IF (FF .LE. 0.) THEN
            CALL ERRMES (5, 3, CHSUBP)
            GO TO 300
         END IF
  210 CONTINUE
      CALL ERRMES (6, 3, CHSUBP)
  300 FSHSSQ=SSQREF*(1.+FBEST*DF1/DF2)
C      IF (IABS(IDFISH) .GE. 1) THEN
C         SPREJ=PREJ
C         SFSH=FSHSSQ
C         SREFND=REFNDF
C         NAMELIST /NAM2/ SPREJ,PCOMP,R2AM1,R2BM1,SREFND,TT,YP,
C     1   SLAMBD,HH,WW,FF,SFSH,JTRY,DPRTRY,FBEST,FACLOG,DPDF
C         IF (LPRINT .GT. 0) WRITE (LPRINT,NML=NAM2)
C         IF (IDFISH .LE. 0) THEN
C            CALL ERRMES (7, 1, CHSUBP)
C            STOP
C         END IF
C      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PLINLS (ISTAGE, IERROR)
C
C  Nonlinear least-squares analysis with first NLIN parameters occurring
C    linearly and parameters NLIN+1,NLIN+NNONL=NPAR nonlinearly.
C  SOLUTN(JPAR) = linear parameters when JPAR=1,NLIN
C               = full steps for nonlinear parameters when
C                 JPAR=NLIN+1,NLIN+NNONL=NPAR
C  PARNLN(JNONL) = nonlinear parameters, JNONL=1,NNONL
C
C  ISTAGE = 1 for preliminary analysis with only 1 broadening and shift
C             parameter
C         = 2 for Reference Solution of full analysis or other difficult
C             solution, e.g., with one alpha = 0.
C         = 3 for regularized solution of full analysis
C  SDREF = rough estimate for SD of fit.  It is continually updated to the
C          minimum while LSTAGE=1 and during the Starting Solution when
C          LSTAGE=2 (i.e., if INISOL=T).
C          SDREF is used for weighting the priors, the regularizors, and
C            the Marquardt rows.
C          The continuous updating of SDREF will cause an extra decrease
C            in OBJECT during the Starting Solution, and this may require
C            extra iterations.
C
C  IERROR = 1 on normal return
C         = 2 on max. iterations.
C
      INCLUDE 'lcmodel.inc'
      DOUBLE PRECISION GRADST, GRAD2, OBJOLD, PMARQ, PMQACT, PMQOLD,
     1                 SDOLD, SOLOLD(MPAR), STEP2
      LOGICAL LERROR, TRYZER
      CHSUBP='PLINLS'
      IF (ISTAGE.LE.0 .OR. ISTAGE.GT.3) CALL ERRMES (1, 5, CHSUBP)
      IF (AMIN1(PMQST(ISTAGE),PMQSTL(ISTAGE),RSTPMX(ISTAGE),
     1          RCONVR(ISTAGE),RMQDEC(ISTAGE),RMQINC(1,ISTAGE),
     2          RMQINC(2,ISTAGE),COSMIN(ISTAGE),RINCRS(ISTAGE),
     3          RPMQMN(ISTAGE)) .LE. 0. .OR.
     3    AMAX1(RSTPMN(ISTAGE),RCONVR(ISTAGE),RMQDEC(ISTAGE),
     4          RINCRS(ISTAGE)) .GE. 1.) CALL ERRMES (2, 4, CHSUBP)
      LSTAGE=MIN0(2,ISTAGE)
      PMARQ=PMQST(ISTAGE)
      PMQACT=0.D0
      COSINE=-2.
      RSTEP=0.
      IERROR=1
      INTERP=0
      ACTRED=0.
      PRERED=0.
      TRYZER=.TRUE.
      NZBAD=0
      NZSKIP=0
C     -------------------------------------------------------------------------
C     Initialize nonnegativity constraints on delta(1/T2).
C     -------------------------------------------------------------------------
      DO 110 JNONL=LRT2ST,NNONL
         NONNEG(NLIN+JNONL)=PARNLN(JNONL) .LE. 0.D0
  110 CONTINUE
C     -------------------------------------------------------------------------
C     Get OBJECT = objective function.
C     -------------------------------------------------------------------------
      CALL SOLVE (LSTAGE, .FALSE., 0.D0, .FALSE., LERROR)
      IF (LERROR) THEN
         OBJECT=DRANGE
         STDDEV=DRANGE
      END IF
 5110 FORMAT (////' Iter', 3X, 'Obj. funct.', 3X, 'Rel. dec.', 3X,
     1         'Pred. dec.', 3X, 'NDF', 3X, 'Marquardt', 3X,
     2         'Rel. step', 6X, 'Cosine', 3X, 'Interps.', 3X,
     3         'Degrees', 3X, 'Deg/ppm')
      IF (LPRINT .GT. 0) WRITE (LPRINT,5110)
      DO 210 ITER=1,MITER(ISTAGE)+1
         PMQOLD=PMQACT
         OBJOLD=OBJECT
         SDOLD=STDDEV
         ITEROL=ITER-1
         DGZER=PARNLN(LPHAST)/RADIAN
         DGPPM=PARNLN(LPHAST+1)/RADIAN
 5210    FORMAT (1X, I4, 1PE14.6, E12.3, E13.3, I6, E12.2,
     1           0P2F12.6, I11, 2F10.2, 1PE15.4)
         IF (LPRINT .GT. 0) WRITE (LPRINT,5210) ITEROL, OBJECT, ACTRED,
     1         PRERED, NDF, PMQACT, RSTEP, COSINE, INTERP, DGZER, DGPPM
         IF (IDUMP(ISTAGE).GE.2 .AND. LPRINT.GT.0) CALL DUMP1 (LSTAGE)
         INTERP=0
         COSINE=-2.
         DO 215 JNONL=1,NNONL
            PAROLD(JNONL)=PARNLN(JNONL)
  215    CONTINUE
         DO 216 JPAR=1,NLIN
            SOLOLD(JPAR)=SOLUTN(JPAR)
  216    CONTINUE
         IF (TRYZER .AND. DOZERO(ISTAGE) .AND. NZSKIP.LE.0) THEN
C           -------------------------------------------------------------------
C           Compute new direction, first trying 0 Marquardt parameter.
C           -------------------------------------------------------------------
            CALL SOLVE (LSTAGE, .TRUE., 0.D0, .FALSE., LERROR)
            IF (LERROR) GO TO 222
            COSINE=-2.
            PMQACT=0.D0
            RSTEP=1.
            CALL PASTEP (RSTEP)
            CALL SOLVE (LSTAGE, .FALSE., 0.D0, .FALSE., LERROR)
            IF (OBJECT .GE. OBJOLD) THEN
               NZBAD=NZBAD+1
               NZSKIP=NZBAD
            ELSE
               NZBAD=0
            END IF
         ELSE
            NZSKIP=NZSKIP-1
         END IF
  222    IF (OBJECT.GE.OBJOLD .OR. .NOT.(TRYZER.AND.DOZERO(ISTAGE))
     1       .OR. LERROR) THEN
C           -------------------------------------------------------------------
C           Try again with PMARQ
C           -------------------------------------------------------------------
            DO 224 JNONL=1,NNONL
               PARNLN(JNONL)=PAROLD(JNONL)
  224       CONTINUE
            DO 225 JPAR=1,NLIN
               SOLUTN(JPAR)=SOLOLD(JPAR)
  225       CONTINUE
            CALL SOLVE (LSTAGE, .TRUE., PMARQ, .FALSE., LERROR)
            IF (LERROR) GO TO 280
            PMQACT=PMARQ
            RSTEP=1.
            CALL PASTEP (RSTEP)
            CALL SOLVE (LSTAGE, .FALSE., 0.D0, .FALSE., LERROR)
            IF (LERROR) GO TO 280
         END IF
         IF (OBJOLD .LE. 0.D0) THEN
            CALL ERRMES (3, 2, CHSUBP)
            GO TO 300
         ELSE
            PRERED=1.D0-OBJLIN/OBJOLD
            ACTRED=1.D0-OBJECT/OBJOLD
         END IF
C        ----------------------------------------------------------------------
C        Test objective function.
C        ----------------------------------------------------------------------
         PMQSAV=PMQACT
         TRYZER=ACTRED .GE. -RINCRS(ISTAGE)
         IF (TRYZER) THEN
C           -------------------------------------------------------------------
C           Objective function has not significantly increased with a full
C             step.  Decrease Marquardt parameter PMARQ.
C           -------------------------------------------------------------------
            PMARQ=DMAX1(RPMQMN(ISTAGE)*PRECIS,PMARQ*RMQDEC(ISTAGE))
C           -------------------------------------------------------------------
C           Convergence test.
C           -------------------------------------------------------------------
            IF (AMAX1(ABS(PRERED),ABS(ACTRED)) .LE. RCONVR(ISTAGE) .AND.
     1          ABS(ACTRED) .LT. 2.*ABS(PRERED) .AND.
     2          PMQACT .LT. PMQSTL(ISTAGE) .AND.
     3          RSTEP .GE. 1.) GO TO 300
            GO TO 210
         ELSE
C           -------------------------------------------------------------------
C           Objective function has significantly increased.
C           Compute dot products of gradient and step and test the cosine of
C             the angle between them.
C           -------------------------------------------------------------------
            GRAD2=0.D0
            STEP2=0.D0
            GRADST=0.D0
            DO 238 JNONL=1,NNONL
               GRAD2=GRAD2+GRAD(JNONL)**2
               STEP2=STEP2+SOLUTN(NLIN+JNONL)**2
               GRADST=GRADST+GRAD(JNONL)*SOLUTN(NLIN+JNONL)
  238       CONTINUE
            IF (GRADST .GT. 0.D0) CALL ERRMES (4, 1, CHSUBP)
            DTERM(1)=DSQRT(STEP2)*DSQRT(GRAD2)
            IF (DTERM(1) .LE. 0.D0) THEN
               CALL ERRMES (5, 2, CHSUBP)
               GO TO 280
            END IF
            COSINE=-GRADST/DTERM(1)
            IF (COSINE .LT. COSMIN(ISTAGE)) GO TO 280
C           -------------------------------------------------------------------
C           Gradient and step are close enough to being parallel.  Interpolate
C             to get RSTEP (fractional step size).
C           -------------------------------------------------------------------
            DO 240 INTERP=1,MINTER(ISTAGE)
               IF (PMQACT .GT. 0.D0) PMARQ=PMARQ*RMQINC(2,ISTAGE)
               DTERM(1)=OBJOLD-OBJECT+GRADST*RSTEP
               IF (DABS(DTERM(1)) .LE. 0.D0) THEN
C                 -------------------------------------------------------------
C                 In the unlikely case that the denominator is exactly zero,
C                   go to the next iteration with an increased Marquardt
C                   parameter.
C                 -------------------------------------------------------------
                  CALL ERRMES (6, 3, CHSUBP)
                  GO TO 280
               END IF
               RSTEP=.5D0*GRADST*RSTEP**2/DTERM(1)
C              ----------------------------------------------------------------
C              If the fractional step size is unreasonable, go to the next
C                iteration with an increased PMARQ.
C              ----------------------------------------------------------------
               IF (RSTEP.LT.RSTPMN(ISTAGE) .OR. RSTEP.GT.RSTPMX(ISTAGE))
     1               GO TO 280
C              ----------------------------------------------------------------
C              Compute the new parameter set.
C              ----------------------------------------------------------------
               CALL PASTEP (RSTEP)
C              ----------------------------------------------------------------
C              Compute the objective function.
C              ----------------------------------------------------------------
               CALL SOLVE (LSTAGE, .FALSE., 0.D0, .FALSE., LERROR)
               IF (LERROR) GO TO 280
               IF (OBJECT .LT. OBJOLD) GO TO 210
  240       CONTINUE
         END IF
C        ----------------------------------------------------------------------
C        The current solution was unsuccessful.  Restore PAROLD in PARNLN,
C          OBJOLD in OBJECT, and SDOLD in STDDEV.  Increase PMARQ and start
C          next iteration.
C        ----------------------------------------------------------------------
  280    DO 290 JNONL=1,NNONL
            PARNLN(JNONL)=PAROLD(JNONL)
  290    CONTINUE
         DO 295 JPAR=1,NLIN
            SOLUTN(JPAR)=SOLOLD(JPAR)
  295    CONTINUE
         OBJECT=OBJOLD
         STDDEV=SDOLD
         PMQSAV=PMQOLD
         IF (PMQACT .GT. 0.D0) PMARQ=dmax1(PMARQ*RMQINC(1,ISTAGE),
     1                                     dble(pmqstl(istage)))
  210 CONTINUE
      IERROR=2
      IF (ISTAGE .GE. 2) CALL ERRMES (7, 1, CHSUBP)
  300 ITEROL=ITER-1
      DGZER=PARNLN(LPHAST)/RADIAN
      DGPPM=PARNLN(LPHAST+1)/RADIAN
      IF (LPRINT .GT. 0) WRITE (LPRINT,5210) ITEROL, OBJECT, ACTRED,
     1      PRERED, NDF, PMQACT, RSTEP, COSINE, INTERP, DGZER, DGPPM
      IF (IDUMP(ISTAGE).GE.1 .AND. LPRINT.GT.0) CALL DUMP1 (LSTAGE)
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DUMP1 (LSTAGE)
C
      INCLUDE 'lcmodel.inc'
      IF (LPRINT .LE. 0) RETURN
 5110 FORMAT (' Conc =', 1P10E12.3/ (7X, 1P10E12.3))
      WRITE (LPRINT,5110) (SOLUTN(J),J=1,NMETAB)
 5120 FORMAT (' Backgr =', 1P10E12.3/ (9X, 1P10E12.3))
      if (nbackg .gt. 0)
     1   WRITE (LPRINT,5120) (SOLUTN(J),J=NMETAB+1,NMETAB+NBACKG)
 5140 FORMAT (' Lineshape =', 10F12.4/ (12X, 10F12.4))
      IF (NSIDE2 .GT. 0) WRITE (LPRINT,5140) (PARNLN(J),J=1,NSIDE2)
 5150 FORMAT (' 1000*Shift =', 3P10F11.2/ (13X, 10F11.2))
      WRITE (LPRINT,5150) (PARNLN(J)/(2.*PI*HZPPPM),J=LSHIST,
     1                       LSHIST+NEXPON-1)
      IF (LSTAGE .EQ. 1) THEN
 5160    FORMAT (' Gaussian FWHM (ppm) =', F8.4)
         WRITE (LPRINT,5160) TOFWHM*PARNLN(LRT2ST)
      ELSE
 5165    FORMAT (' delta(1/T2) =', 10F11.2/ (14X, 10F11.2))
         WRITE (LPRINT,5165) (PARNLN(J),J=LRT2ST,nnonl)
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASTEP (RSTEP)
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='PASTEP'
      DO 110 JNONL=1,NNONL
         IF(SOLUTN(NLIN+JNONL) .NE. 0.D0) RSTEP=AMIN1(RSTEP,
     1        ABS(FSTPMQ*SNGL(DPARMQ(JNONL)/SOLUTN(NLIN+JNONL))))
  110 CONTINUE
C     -------------------------------------------------------------------------
C     Enforce nonnegativity of delta(1/T2).
C     -------------------------------------------------------------------------
      LPAR=0
      DO 115 JNONL=LRT2ST,NNONL
         JPAR=NLIN+JNONL
         IF (SOLUTN(JPAR).LT.0.D0 .AND. .NOT.NONNEG(JPAR)) THEN
            TERM=DABS(PARNLN(JNONL)/SOLUTN(JPAR))
            IF (TERM .LT. RSTEP) THEN
               RSTEP=TERM
               LPAR=JPAR
            END IF
         END IF
         NONNEG(JPAR)=NONNEG(JPAR) .AND. SOLUTN(JPAR).LE.0.D0
  115 CONTINUE
      IF (LPAR .NE. 0) NONNEG(LPAR)=.TRUE.
      DO 120 JNONL=1,NNONL
         PARNLN(JNONL)=PAROLD(JNONL)+RSTEP*SOLUTN(NLIN+JNONL)
  120 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SOLVE (LSTAGE, DONONL, PMQACT, ONLYFT, LERROR)
C
C  Sets up and solves nonnegative linear least squares problem.
C
C  SOLVE must first be called with PMQACT=0 and DONONL=ONLYFT=F so that SDREF
C    can be initialized for computing the Marquardt rows.
C  SDREF must be initialized to DRANGE for the test of SDREF below to be
C    effective.
C
C  LSTAGE = 1 for preliminary unregularized analysis with only 1 Gaussian
C             peak-broadener, 1 shift, NSIDES=0, and no priors.
C         = 2 for full analysis.
C
C  ONLYFT = T to only compute for plotting:
C             YFITRE(JY,0) = the fit to the data,
C             YFITRE(JY,JMETAB) = the contribution of metabolite JMETAB +
C                                 background to the fit,
C           and return without changing DAMAT, so that error estimates can be
C           computed from DAMAT from previous call.
C         = F for normal computation of solution and variance from the real
C             data (from zero-filling).
C
C  DONONL = T for computation of SOLUTN and OBJLIN
C         = F for evaluation of OBJECT & initial linear part of SOLUTN only.
C
C  Needs (in COMMON) from previous iteration:
C    SOLUTN(JPAR), JPAR=1,NLIN (if DONONL=T)
C    PARNLN(JNONL), JNONL=1,NNONL (always)
C
C  Columns:
c    NLESS>0 will reduce indices above NLIN by NLESS (not shown here):
C    1, NMETAB : nonnegative concentration coefficients
C    NMETAB+1, NMETAB+NBACKG=NLIN : background coefficients
C     NLIN+1, NLIN+2*NSIDES+ : lineshape coefficients
C     LPHAST, LPHAST+1: 0- and 1st-order phase corrections
C      LSHIST, LSHIST-1+NMETAB : shift corrections
C      LRT2ST, LRT2ST-1+NMETAB : delta(1/T2) corrections
C
C  Rows:
C    1, NY : real parts of data from zero-filling
C    NY+1, NY+NRATIO_USED : priors for CONC ratios
C    NY+NRATIO_USED+1, NY+NRATIO_USED+NMETAB : priors for shift corrections
C    NY+NRATIO_USED+NMETAB+1, NY+NRATIO_USED+2*NMETAB : priors for
C                                                       delta(1/T2) corrections
C    NY+2*NMETAB+NRATIO_USED+1, NY+2*NMETAB+NRATIO_USED+[0--2] = NROWDA:
c                                             priors for phase corrections
C    NROWDA+1, NROWDA+NBACKG : regularizor (if ALPHAB>0)
C    NROWDA+(NBACKG)+1, NROWDA+(NBACKG)+(NSIDE2+3): regularizor
C                                                   (if ALPHAS>0 and DONONL=T).
C    NROWDA+(NBACKG)+(NSIDE2+3)+1, NROWDA+(NBACKG)+(NSIDE2+3)+NNONL:
C      Marquardt rows (if PMQACT>0 and DONONL=T).
C
      INCLUDE 'lcmodel.inc'
      external penlty
C      COMPLEX CYFIT(MY)!CT1!CT2
C      LOGICAL LTEST!CT1!CT2
C      COMMON /BLT1/ CYFIT, LTEST!CT1!CT2
      SAVE BASIYF, DCEXOL, LSTOLD
      COMPLEX*16 BASIYF(1-MINCSD*MSIDES:MY+MINCSD*MSIDES,mmetab,2),
     1           basiyf_power(MY, mcoeff_power),
     1           DCDATA(MDATA,5), DCEXOL(MMETAB,2), DCEXUS, DCFACT,
     2           DCPHAS, DCSUM(2), DCTERM, DCYFIT
C      COMPLEX*16 CDEXP!IRIX
      DOUBLE PRECISION DAMARQ(MNONL), DEXPRE,
     1                 DFZERO, DPPINC, DPPM, DRHS(MROW),
     2                 DDTIME, DTIME, DWORK(MROW+MPAR), PENLTY, OBJTOT,
     3                 PMQACT, SQRTWT
      INTEGER LSTOLD(2)
      LOGICAL DONONL, LERROR, NEWEXP, ONLYFT, PRIORP, PRIORZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON DAMARQ, DCDATA, DRHS, DWORK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      CDEXP(DCTERM)=ZEXP(DCTERM)!IRIX
      CHSUBP='SOLVE'
      dcsum(1) = (0.d0,0.d0)
      dcsum(2) = (0.d0,0.d0)
      if (initialize_solve) then
         initialize_solve = .false.
         do 50 k = 1, 2
            lstold(k) = 0
            do 55 j = 1, mmetab
               dcexol(j, k) = (0.D0,0.D0)
 55         continue
 50      continue
      end if
      IF (LSTAGE.LE.0 .OR. LSTAGE.GE.3) CALL ERRMES (1, 5, CHSUBP)
      PRIORZ=.NOT.(FXDEGZ .OR. SDDEGZ.GE.45.)
      PRIORP=.NOT.FXDEGP
C     -------------------------------------------------------------------------
C     Indices starting with K are only used instead of those with L for
C       temporarily contracting DAMAT, GRAD & NONNEG before the call to PNNLS
C       and right after restoring these.  This only is needed when FXDEGZ or
C       FXDEGP = T.
c     SOLUTN comes out of PNNLS contracted.
C     Outside of SOLVE these arrays are all the full size.
C     PARNLN is never compressed.
C     -------------------------------------------------------------------------
      IF (FXDEGZ) THEN
         NLESS=1
         KDEGP=LPHAST
      ELSE
         NLESS=0
         KDEGP=LPHAST+1
      END IF
      IF (FXDEGP) NLESS=NLESS+1
      KSHIST=LSHIST-NLESS
      KRT2ST=LRT2ST-NLESS
      KNONL=NNONL-NLESS
C     -------------------------------------------------------------------------
C     Compute broadened basis in frequency-domain.
C     -------------------------------------------------------------------------
      DDTIME=DBLE(DELTAT)
      IF (LSTAGE .EQ. 1) THEN
C        ----------------------------------------------------------------------
C        DCDATA(JDATA,1) = data for individual metabolites
C                     2  = sum over all metabolites of conc*(metabolite data)
C                     3  = broadening and shifting factors for time-domain
C                          data.
c                     4  = shifting factors for time-domain data (without
c                          broadening -- when LSHAPE=F).
C                     5  = sum over all broadened (LSHAPE=T) metabolites of
c                          conc*(metabolite data)
C        ----------------------------------------------------------------------
         IF (DONONL) THEN
            DO 101 JDATA=1,NDATA
               DCDATA(JDATA,2)=(0.D0,0.D0)
               DCDATA(JDATA,5)=(0.D0,0.D0)
  101       CONTINUE
         END IF
 
c         write (6, 9101) PARNLN(LRT2ST) * tofwhm
c 9101    format ('Broadening FWHM =', 1pe13.4)
 
         DEXPRE=PARNLN(LRT2ST)**2
         DCEXUS=DCMPLX(PARNLN(LRT2ST),PARNLN(LSHIST))
         NEWEXP=DCEXUS.NE.DCEXOL(1,1) .OR. LSTOLD(1).NE.1
         IF (NEWEXP .OR. DONONL) THEN
            DTIME=0.D0
            IF (NEWEXP) THEN
               DO 102 JDATA=1,NDATA
                  DCDATA(JDATA,4)=CDEXP(DCMPLX(0.d0,
     1                                         -dtime*PARNLN(LSHIST)))
                  dcdata(jdata,3) = dexp(-dexpre * dtime**2) *
     1                              dcdata(jdata,4)
                  DTIME=DTIME+DDTIME
  102          CONTINUE
            END IF
            DO 103 JMETAB=1,NMETAB
c              ----------------------------------------------------------------
c              LSHAPE = T to do Gaussian broadening (LSTAGE=1).
c              ----------------------------------------------------------------
               if (lshape(jmetab)) then
                  DO 104 JDATA=1,NDATA
                     DCDATA(JDATA,1)=DCDATA(JDATA,3)*
     1                               BASIST(JDATA,JMETAB)
 104              CONTINUE
                  if (dononl) then
                     do 1042 jdata = 1, ndata
                        dcterm = SOLUTN(JMETAB)*DCDATA(JDATA,1)
                        DCDATA(JDATA,2)=DCDATA(JDATA,2) + dcterm
                        DCDATA(JDATA,5)=DCDATA(JDATA,5) + dcterm
 1042                continue
                  end if
               else
c                 -------------------------------------------------------------
c                 No broadening (LSHAPE = F).
c                 -------------------------------------------------------------
                  DO 1045 JDATA=1,NDATA
                     DCDATA(JDATA,1)=DCDATA(JDATA,4)*
     1                               BASIST(JDATA,JMETAB)
 1045             CONTINUE
                  if (dononl) then
                     DO 1047 JDATA=1,NDATA
                        DCDATA(JDATA,2)=DCDATA(JDATA,2) +
     1                                  SOLUTN(JMETAB)*DCDATA(JDATA,1)
 1047                continue
                  end if
               end if
               IF (NEWEXP) THEN
                  LSTOLD(1)=1
                  DCEXOL(1,1)=DCEXUS
                  CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                        DWFFTC)
                  JDATA=LDATST
                  DO 105 JY=1,NY
                     BASIYF(JY,JMETAB,1)=DCDATA(JDATA,1)
                     JDATA=JDATA+1
  105             CONTINUE
               END IF
  103       CONTINUE
         END IF
         IF (DONONL) THEN
C           -------------------------------------------------------------------
C           i*exp(-phase)*BASIYF(JY,1,2) = derivative wrt shift parameter.
C             exp(-phase)*BASIYF(JY,2,2) = derivative wrt broadening parameter.
C           BASIYF(*,*,2) from previous calls cannot be used, because they
C             contain concentrations, which change in every call.
C           -------------------------------------------------------------------
            DTIME=0.D0
            DO 106 JDATA=1,NDATA
               DCDATA(JDATA,1)=DTIME*DCDATA(JDATA,2)
               DCDATA(JDATA,2)=DTIME**2 * DCDATA(JDATA,5)
               DTIME=DTIME-DDTIME
  106       CONTINUE
            CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                  DWFFTC)
            CALL DCFFT_R (DCDATA(1,2), DCDATA(1,2), NDATA, LDWFFT,
     1                  DWFFTC)
            JDATA=LDATST
            DO 107 JY=1,NY
               BASIYF(JY,1,2)=DCDATA(JDATA,1)
               BASIYF(JY,2,2)=-2.D0*PARNLN(LRT2ST)*DCDATA(JDATA,2)
               JDATA=JDATA+1
  107       CONTINUE
         END IF
      ELSE
C        ----------------------------------------------------------------------
C        LSTAGE = 2.
C        ----------------------------------------------------------------------
         DO 110 JMETAB=1,NMETAB
            if (imethd .eq. 3) then
               dtime = 0.d0
               do 112 jdata = 1, ndata
                  dterm(1) = 0.d0
                  do 114 jpower = 1, npower(jmetab)
                     dterm(1) = dterm(1) -
     1                          parnln(lpowen(jmetab - 1) + jpower) *
     2                          tpower(jpower, jdata)
 114              continue
                  dtime = dtime - ddtime
                  dcdata(jdata, 1) = cdexp(dcmplx(dterm(1),
     1               PARNLN(LSHIST-1+JMETAB) * dtime)) *
     2               BASIST(JDATA,JMETAB)
 112           continue
            else
               DCEXUS=DCMPLX(PARNLN(LRT2ST-1+JMETAB),PARNLN(LSHIST-1+
     1                                                      JMETAB))
               IF (DCEXUS.EQ.DCEXOL(JMETAB,1) .AND. LSTOLD(1).EQ.2)
     1            GO TO 110
               LSTOLD(1)=2
               DCEXOL(JMETAB,1)=DCEXUS
               DCFACT=CDEXP(-DCEXUS*DDTIME)
               DCTERM=(1.D0,0.D0)
               DO 120 JDATA=1,NDATA
                  DCDATA(JDATA,1)=DCTERM*BASIST(JDATA,JMETAB)
                  DCTERM=DCTERM*DCFACT
 120           CONTINUE
            end if
            CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                    DWFFTC)
C           -------------------------------------------------------------------
C           BASIYF(JYSIDE,JMETAB,1) = shifted frequency-domain basis spectra
C                                     in region of interest with extensions
C                                     of INCSID*NSIDES points at each end.
C           -------------------------------------------------------------------
            JDATA=LDATST-INCSID*NSIDES
            DO 130 JYSIDE=1-INCSID*NSIDES,NY+INCSID*NSIDES
               BASIYF(JYSIDE,JMETAB,1)=DCDATA(JDATA,1)
               JDATA=JDATA+1
  130       CONTINUE
  110    CONTINUE
         IF (DONONL) THEN
C           -------------------------------------------------------------------
C           Compute derivative wrt delta(1/T2) of broadened frequency-domain
C             basis spectra.
C           -------------------------------------------------------------------
            if (imethd .eq. 3) then
C              ----------------------------------------------------------------
C              BASIYF(JY,JMETAB,2) = (derivative wrt shift of exponential
c                                    term) / -sqrt(-1)
C              ----------------------------------------------------------------
               do 150 jmetab = 1, nmetab
                  dtime = 0.d0
                  do 151 jdata = 1, ndata
                     dterm(1) = 0.d0
                     do 152 jpower = 1, npower(jmetab)
                        dterm(1) = dterm(1) -
     1                       parnln(lpowen(jmetab - 1) + jpower) *
     2                       tpower(jpower, jdata)
 152                 continue
                     dtime = dtime - ddtime
                     dcdata(jdata, 1) = dtime *
     1                  cdexp(dcmplx(dterm(1),
     2                           PARNLN(LSHIST-1+JMETAB) * dtime)) *
     3                  BASIST(JDATA,JMETAB)
 151              continue
                  CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                          DWFFTC)
                  jdata = ldatst
                  DO 153 JY = 1, NY
                     BASIYF(JY, jmetab, 2)=DCDATA(JDATA,1)
                     JDATA=JDATA+1
 153              CONTINUE
 150           continue
C              ----------------------------------------------------------------
C              BASIYF_power(Jy,jcoeff_power) = derivative wrt coefficient of
c                                              POWER
C              ----------------------------------------------------------------
               jcoeff_power = 0
               do 155 jmetab = 1, nmetab
                  do 156 kpower = 1, npower(jmetab)
                     jcoeff_power = jcoeff_power + 1
                     dtime = 0.d0
                     do 157 jdata = 1, ndata
                        dterm(1) = 0.d0
                        do 159 jpower = 1, npower(jmetab)
                           dterm(1) = dterm(1) -
     1                          parnln(lpowen(jmetab - 1) + jpower) *
     2                          tpower(jpower, jdata)
 159                    continue
                        dtime = dtime - ddtime
                        dcdata(jdata, 1) = -tpower(kpower, jdata) *
     1                     cdexp(dcmplx(dterm(1),
     2                              PARNLN(LSHIST-1+JMETAB) * dtime)) *
     3                     BASIST(JDATA,JMETAB)
 157                 continue
                     CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                             DWFFTC)
                     jdata = ldatst
                     DO 158 JY = 1, NY
                        BASIYF_power(JY, jcoeff_power)=DCDATA(JDATA,1)
                        JDATA=JDATA+1
 158                 CONTINUE
 156              continue
 155           continue
            else
               DO 160 JMETAB=1,NMETAB
                  DCEXUS=DCEXOL(JMETAB,1)
                  IF (DCEXUS.EQ.DCEXOL(JMETAB,2) .AND. LSTOLD(2).EQ.2)
     1               go TO 160
                  LSTOLD(2)=2
                  DCEXOL(JMETAB,2)=DCEXUS
                  DCFACT=CDEXP(-DCEXUS*DDTIME)
                  DTIME=0.D0
                  DCTERM=(1.D0,0.D0)
                  DO 170 JDATA=1,NDATA
                     DCDATA(JDATA,1)=DTIME*DCTERM*BASIST(JDATA,JMETAB)
                     DCTERM=DCTERM*DCFACT
                     DTIME=DTIME-DDTIME
 170              CONTINUE
                  CALL DCFFT_R (DCDATA, DCDATA, NDATA, LDWFFT,
     1                          DWFFTC)
C                 -------------------------------------------------------------
C                 BASIYF(JYSIDE,JMETAB,2) = derivative wrt gamma (in MRM) of
C                                           shifted frequency-domain basis in
C                                           spectral region of interest with
C                                           extensions of INCSID*NSIDES points
C                                           at each end.
C                 -------------------------------------------------------------
                  JDATA=LDATST-INCSID*NSIDES
                  DO 180 JYSIDE=1-INCSID*NSIDES,NY+INCSID*NSIDES
                     BASIYF(JYSIDE,JMETAB,2)=DCDATA(JDATA,1)
                     JDATA=JDATA+1
 180              CONTINUE
 160           CONTINUE
            END IF
         END IF
      END IF
C     -------------------------------------------------------------------------
C     Load DAMAT and DRHS, or YFITRE.
C     -------------------------------------------------------------------------
      DSSQ=0.D0
      IF (.NOT.ONLYFT) THEN
         IF (DONONL) THEN
            DO 202 JNONL=1,KNONL
               GRAD(JNONL)=0.D0
  202       CONTINUE
         END IF
         DO 204 JPAR=1,NPAR
            DO 205 IROW=1,MROW
               DAMAT(IROW,JPAR)=0.D0
  205       CONTINUE
  204    CONTINUE
         DO 206 IROW=1,MROW
            DRHS(IROW)=0.D0
  206    CONTINUE
      END IF
      DFZERO=1.D0
      DO 208 JNONL=1,NSIDE2
         DFZERO=DFZERO-PARNLN(JNONL)
  208 CONTINUE
      DPPM=DBLE(DELPPM(1))
      DPPINC=DBLE(PPMINC)
      DCPHAS=CDEXP(DCMPLX(0.D0,
     1                    dble(-PARNLN(LPHAST)-PARNLN(LPHAST+1)*DPPM)))
      DCFACT=CDEXP(DCMPLX(0.D0,dble(PARNLN(LPHAST+1)*DPPINC)))
      nrow = 0
      DO 210 JY=1,NY
         if (lcy_skip(jy)) go to 210
         nrow = nrow + 1
         IF (ONLYFT) YFITRE(NROW,0)=0.
         DCYFIT=(0.D0,0.D0)
         jcoeff_power = 0
         DO 220 JMETAB=1,NMETAB
c           -------------------------------------------------------------------
c            LSHAPE(JMETAB) = F to eliminate convolution with Lineshape
c                               function for JMETAB.
c           -------------------------------------------------------------------
            if (lshape(jmetab)) then
               DCSUM(1)=DFZERO*BASIYF(JY,JMETAB,1)
               IF (DONONL .AND. LSTAGE.EQ.2) DCSUM(2)=DFZERO*
     1                                             BASIYF(JY,JMETAB,2)
               JNONL=1
               DO 230 JYSIDE=JY-INCSID*NSIDES,JY+INCSID*NSIDES,INCSID
                  IF (JYSIDE .EQ. JY) GO TO 230
                  DCSUM(1)=DCSUM(1)+PARNLN(JNONL)*
     1                              BASIYF(JYSIDE,JMETAB,1)
                  IF (DONONL) DCSUM(2)=DCSUM(2)+PARNLN(JNONL)*
     1                                          BASIYF(JYSIDE,JMETAB,2)
                  JNONL=JNONL+1
 230           CONTINUE
            else
               DCSUM(1)=BASIYF(JY,JMETAB,1)
               IF (DONONL .AND. LSTAGE.EQ.2) DCSUM(2)=
     1                                       BASIYF(JY,JMETAB,2)
            end if
            DCSUM(1)=DCSUM(1)*DCPHAS
            IF (ONLYFT) THEN
               YFITRE(NROW,JMETAB)=SOLUTN(JMETAB)*real(REAL(DCSUM(1)))
               YFITRE(NROW,0)=YFITRE(NROW,0)+YFITRE(NROW,JMETAB)
            ELSE
C              ----------------------------------------------------------------
C              Load derivatives wrt concentration coefficients.
C              ----------------------------------------------------------------
               DAMAT(NROW,JMETAB)=DREAL(DCSUM(1))
               IF (DONONL) THEN
                  DCYFIT=DCYFIT+SOLUTN(JMETAB)*DCSUM(1)
C                 -------------------------------------------------------------
C                 Load derivatives wrt shift and delta(1/T2)
C                 -------------------------------------------------------------
                  IF (LSTAGE .EQ. 2) THEN
                     DCSUM(2)=DCSUM(2)*DCPHAS*SOLUTN(JMETAB)
                     DAMAT(NROW,NLIN+KSHIST-1+JMETAB)=-DIMAG(DCSUM(2))
                     if (imethd .eq. 3) then
                        do 240 kpower = 1, npower(jmetab)
                           jcoeff_power = jcoeff_power + 1
                           DAMAT(NROW,NLIN+KRT2ST-1+jcoeff_power) =
     1                        dreal(DCPHAS * SOLUTN(JMETAB) *
     2                        BASIYF_power(JY, jcoeff_power))
 240                    continue
                     else
                        DAMAT(NROW,NLIN+KRT2ST-1+JMETAB)=DREAL(DCSUM(2))
                     end if
                  ELSE IF (JMETAB .EQ. 1) THEN
                     DAMAT(NROW,NLIN+KRT2ST)=
     1                  DREAL(DCPHAS*BASIYF(JY,2,2))
                     DAMAT(NROW,NLIN+KSHIST)=
     1                  -DIMAG(DCPHAS*BASIYF(JY,1,2))
                  END IF
               END IF
            END IF
  220    CONTINUE
C        ----------------------------------------------------------------------
C        Background is assumed real.
C        ----------------------------------------------------------------------
         BACFIT=0.
         JSOL=NMETAB
         DO 250 JBACKG=1,NBACKG
            JSOL=JSOL+1
            DCTERM=DCPHAS*DBLE(BACKGR(JY,JBACKG))
            IF (ONLYFT) THEN
               BACFIT=BACFIT+SNGL(SOLUTN(JSOL))*real(REAL(DCTERM))
            ELSE
               DAMAT(NROW,JSOL)=DREAL(DCTERM)
               DCYFIT=DCYFIT+SOLUTN(JSOL)*DCTERM
            END IF
  250    CONTINUE
         IF (ONLYFT) THEN
            BACKRE(NROW)=BACFIT
            IF (SUBBAS) THEN
               YREAL(NROW)=REAL(CY(JY))-BACFIT
            ELSE
               YREAL(NROW)=REAL(CY(JY))
               DO 260 JMETAB=0,NMETAB
                  YFITRE(NROW,JMETAB)=YFITRE(NROW,JMETAB)+BACFIT
  260          CONTINUE
            END IF
            GO TO 390
         END IF
C        ----------------------------------------------------------------------
C        Right-hand side.
C        ----------------------------------------------------------------------
         DRHS(NROW)=DBLE(REAL(CY(JY)))
         IF (DONONL) THEN
C           -------------------------------------------------------------------
C           Derivatives wrt phase corrections.
C           -------------------------------------------------------------------
            IF (.NOT.FXDEGZ) DAMAT(NROW,NLIN+LPHAST)=DIMAG(DCYFIT)
            IF (.NOT.FXDEGP) DAMAT(NROW,NLIN+KDEGP)=DPPM*DIMAG(DCYFIT)
            DPPM=DPPM-DPPINC
         END IF
         IF (DONONL) THEN
C           -------------------------------------------------------------------
C           Derivatives wrt to lineshape coefficients
C           -------------------------------------------------------------------
            ICOL=NLIN
            DO 310 JYSIDE=JY-INCSID*NSIDES,JY+INCSID*NSIDES,INCSID
               IF (JYSIDE .EQ. JY) GO TO 310
               DCSUM(1)=(0.D0,0.D0)
               DO 320 JMETAB=1,NMETAB
                  if (lshape(jmetab)) DCSUM(1)=DCSUM(1)+SOLUTN(JMETAB)*
     1                                         (BASIYF(JYSIDE,JMETAB,1)-
     2                                          BASIYF(JY,JMETAB,1))
  320          CONTINUE
               ICOL=ICOL+1
               DAMAT(NROW,ICOL)=DREAL(DCSUM(1)*DCPHAS)
  310       CONTINUE
C           -------------------------------------------------------------------
C           Gradients.
C           -------------------------------------------------------------------
            DTERM(1)=-2.D0*(DRHS(NROW)-DREAL(DCYFIT))
            DO 350 JNONL=1,KNONL
               GRAD(JNONL)=GRAD(JNONL)+DTERM(1)*DAMAT(NROW,NLIN+JNONL)
  350       CONTINUE
         END IF
  390    DCPHAS=DCPHAS*DCFACT
C      IF (LTEST) CYFIT(JY)=DCYFIT!CT1!CT2
  210 CONTINUE
      IF (ONLYFT) RETURN
C     -------------------------------------------------------------------------
c     Append priors for CONC ratios.
c     Weighting according to SDRATI*CSUM  has been done.  So, only need to
c       weight elements of CPRIOR by SDREF
c     For (linear) CONC, the row is simply the original constraint.
c       RHS=0 (not rhs_orig - fit).
C     -------------------------------------------------------------------------
      if (sdref .lt. drange) then
         do 401 jratio = 1, nratio_used
            nrow = nrow + 1
            do 402 jmetab = 1, nmetab
               damat(nrow, jmetab) = sdref *
     1                               dble(cprior(jratio, jmetab))
 402        continue
 401     continue
      end if
      IF (DONONL .AND. LSTAGE.EQ.2) THEN
C        ----------------------------------------------------------------------
C        Append rows for priors.
C        EXRT2(JMETAB) = prior mean of delta(1/T2)
C        SDRT2(JMETAB) = prior standard deviation of delta(1/T2)
C        SDREF = rough estimate for SD of fit.  It is continually updated
C                to the minimum while LSTAGE=1 and during the Starting
C                Solution when LSTAGE=2 (i.e., if INISOL=T).
C        First append rows for priors for shifts.
C        SDSHIF(JMETAB) = prior standard deviation of shift (prior mean is
C                         assumed 0).  No special provision for a larger
C                         SD in the overall shift of the data is (or
C                         should) be made, since this can be effected by
C                         a shift in the lineshape coefficients.
C        ----------------------------------------------------------------------
         DO 404 JMETAB=1,NMETAB
            NROW=NROW+1
            JNONL=KSHIST-1+JMETAB
            IF (SDSHIF(JMETAB) .LE. 0.) CALL ERRMES (2, 4, CHSUBP)
            SQRTWT=SDREF/SDSHIF(JMETAB)
            DAMAT(NROW,NLIN+JNONL)=SQRTWT
            DRHS(NROW)=-SQRTWT*PARNLN(JNONL+NLESS)
            GRAD(JNONL)=GRAD(JNONL)-2.D0*DRHS(NROW)*SQRTWT
 404     CONTINUE
C        ----------------------------------------------------------------------
C        Append rows for priors for group shifts.
C        ----------------------------------------------------------------------
         do 450 jrow_group_shift = 1, nrow_group_shift
            nrow = nrow + 1
            sqrtwt = sdref / sdgroup_shift_row(jrow_group_shift)
c           -------------------------------------------------------------------
c           DTERM(1) = weighted residual
c           -------------------------------------------------------------------
            dterm(1) = 0.d0
            do 460 jmetab = 1, nmetab
               jnonl = kshist + jmetab - 1
               dterm(2) = cgroup_shift(jrow_group_shift, jmetab) *
     1                    sqrtwt
               damat(nrow, nlin + jnonl) = dterm(2)
               dterm(1) = dterm(1) - dterm(2) * parnln(jnonl + nless)
 460        continue
            drhs(nrow) = dterm(1)
            do 470 jnonl = kshist, kshist + nmetab - 1
               grad(jnonl) = grad(jnonl) - 2.d0 * dterm(1) *
     1                                     damat(nrow, nlin + jnonl)
 470        continue
 450     continue
C        ----------------------------------------------------------------------
C        Append rows for priors for delta(1/T2).
C        ----------------------------------------------------------------------
         if (imethd.eq.3) then
            jcoeff_power = 0
            do 480 JMETAB=1,NMETAB
                do 485 kpower = 1, npower(jmetab)
                   NROW=NROW+1
                   jcoeff_power = jcoeff_power + 1
                   jnonl = krt2st - 1 + jcoeff_power
                   IF (coeff_power_sd(kpower, jmetab) .LE. 0.d0)
     1                CALL ERRMES (9, 4, CHSUBP)
                   SQRTWT = sdref / coeff_power_sd(kpower, jmetab)
                   DAMAT(NROW, nlin + jnonl)=SQRTWT
                   DRHS(NROW)=-SQRTWT*PARNLN(jnonl + nless)
                   GRAD(jnonl)=GRAD(jnonl)-2.D0*DRHS(NROW)*SQRTWT
 485            continue
 480         continue
         else
            DO 490 JMETAB=1,NMETAB
               NROW=NROW+1
               JNONL=KRT2ST-1+JMETAB
               IF (SDRT2(JMETAB) .LE. 0.) CALL ERRMES (5, 4, CHSUBP)
               SQRTWT=SDREF/SDRT2(JMETAB)
               DAMAT(NROW,NLIN+JNONL)=SQRTWT
               DRHS(NROW)=SQRTWT*(EXRT2(JMETAB)-PARNLN(JNONL+NLESS))
               GRAD(JNONL)=GRAD(JNONL)-2.D0*DRHS(NROW)*SQRTWT
 490        CONTINUE
         END IF
      END IF
      IF (DONONL) THEN
C        ----------------------------------------------------------------------
C        Append rows for priors for phases.
C        ----------------------------------------------------------------------
         IF (PRIORZ) THEN
            NROW=NROW+1
            SQRTWT=SDREF/(SDDEGZ*RADIAN)
            DAMAT(NROW,NLIN+LPHAST)=SQRTWT
            DRHS(NROW)=SQRTWT*(EXDEGZ*RADIAN-PARNLN(LPHAST)-PHITOT(1))
            GRAD(LPHAST)=GRAD(LPHAST)-2.D0*DRHS(NROW)*SQRTWT
         END IF
         IF (PRIORP) THEN
            NROW=NROW+1
            SQRTWT=SDREF/(SDDEGP*RADIAN)
            DAMAT(NROW,NLIN+KDEGP)=SQRTWT
            DRHS(NROW)=SQRTWT*(EXDEGP*RADIAN-PARNLN(LPHAST+1)
     1                                      -PHITOT(2))
            GRAD(KDEGP)=GRAD(KDEGP)-2.D0*DRHS(NROW)*SQRTWT
         END IF
      END IF
      NROWDA=NROW
      if (imethd .eq. 2   .and.   lstage .eq. 2   .and.
     1    alphas .gt. 0.d0) then
c        ---------------------------------------------------------------------
c        IMETHD = 2; Append rows for this special case with linewidths
c                    (rather than lineshape) being regularized.  It
c                    contains concentration and must therefore also be
c                    done with linear analysis (DONONL=F & T).
c        No idea why DRHS & GRAD should be commented out below, but it only
c           works this way.
c        ---------------------------------------------------------------------
         do 495 jmetab = 1, nmetab
            nrow = nrow + 1
            jnonl = krt2st - 1 + jmetab
            sqrtwt = sdref * alphas / conc_expect(jmetab)
            dterm(1) = rt2min(jmetab) + rlrntz * parnln(jnonl + nless)
            if (ipowrg .eq. 1) then
               damat(nrow, jmetab) = sqrtwt * dterm(1)
            else
               damat(nrow, jmetab) = sqrtwt * dterm(1)**2
            end if
            if (dononl) then
               if (ipowrg .eq. 1) then
                  DAMAT(NROW,NLIN+JNONL) = sqrtwt * solutn(jmetab) *
     1                                     rlrntz
               else
                  DAMAT(NROW,NLIN+JNONL) = 2.d0 * sqrtwt *
     1               solutn(jmetab) * dterm(1) * rlrntz
               end if
cc               drhs(nrow) = -solutn(jmetab) * damat(nrow, jmetab)
cc               grad(jnonl) = grad(jnonl) - 2.0 * drhs(nrow) * rlrntz *
cc     1                                     sqrtwt * solutn(jmetab)
            end if
 495     continue
      end if
      IF (ALPHAB.GT.0.D0 .AND. SDREF.LT.DRANGE) THEN
C        ----------------------------------------------------------------------
C        Append rows for background regularizor, REGB.
C        ----------------------------------------------------------------------
         DO 510 JBACKG=1,NBACKG
            NROW=NROW+1
            DO 520 JREG=1,NBACKG
               DAMAT(NROW,NMETAB+JREG)=SDREF*ALPHAB*REGB(JBACKG,JREG)
  520       CONTINUE
  510    CONTINUE
      END IF
      IF (ALPHAS.GT.0.D0  .AND. DONONL   .and.   imethd .ne. 2   .and.
     1    nside2 .gt. 0) then
C        ----------------------------------------------------------------------
C        Append rows for regularizor for lineshape coefficients.
C        REGF is assumed to be for the lineshape coefficients with the
C           center coefficient eliminated by the equality constraint.
C           Therefore, the 1st 3 rows are assumed to have the constants 2,
C           -1, -1 on the rhs, rather than 0.
c        NSIDE2 = 0 when LSTAGE = 1, but could conceivably have NSIDE2=0
c                   in other cases; so above test is safer.
C        ----------------------------------------------------------------------
         DO 530 IROWF=1,NREGF
            NROW=NROW+1
            IF (IROWF .EQ. 1) THEN
               DTERM(1)=2.*SDREF*ALPHAS
            ELSE IF (IROWF .LE. 3) THEN
               DTERM(1)=-SDREF*ALPHAS
            ELSE
               DTERM(1)=0.
            END IF
            ICOL=NLIN
            DO 540 JNONL=1,NSIDE2
               DTERM(2)=SDREF*ALPHAS*REGF(IROWF,JNONL)
               DTERM(1)=DTERM(1)-DTERM(2)*PARNLN(JNONL)
               ICOL=ICOL+1
               DAMAT(NROW,ICOL)=DTERM(2)
 540        CONTINUE
            DRHS(NROW)=DTERM(1)
            ICOL=NLIN
            DO 545 JNONL=1,NSIDE2
               ICOL=ICOL+1
               GRAD(JNONL)=GRAD(JNONL)-2.D0*DTERM(1)*DAMAT(NROW,ICOL)
 545        CONTINUE
 530     CONTINUE
      end if
      IF (PMQACT.GT.0.D0 .AND. DONONL) THEN
C        ----------------------------------------------------------------------
C        Append Marquardt rows.
C        ----------------------------------------------------------------------
         IF (SDREF .GE. DRANGE) CALL ERRMES (6, 4, CHSUBP)
         RTMARQ=SQRT(SNGL(PMQACT))
         INONL=0
         DO 550 JNONL=1,NNONL
            IF ((FXDEGZ .AND. JNONL.EQ.LPHAST) .OR.
     1          (FXDEGP .AND. JNONL.EQ.LPHAST+1)) GO TO 550
            INONL=INONL+1
            NROW=NROW+1
            IF (DPARMQ(JNONL) .LE. 0.) CALL ERRMES (7, 4, CHSUBP)
            DAMARQ(INONL)=RTMARQ*SDREF/DPARMQ(JNONL)
            DAMAT(NROW,NLIN+INONL)=DAMARQ(INONL)
  550    CONTINUE
      END IF
      IF (DONONL) THEN
         NCOL=NPAR-NLESS
         IF (NLESS .GT. 0) THEN
C           -------------------------------------------------------------------
C           Compress NONNEG
C           -------------------------------------------------------------------
            IF (FXDEGZ .AND. .NOT.FXDEGP) NONNEG(NLIN+LPHAST)=
     1                                    NONNEG(NLIN+LPHAST+1)
            DO 560 JNONL=KSHIST,KNONL
               NONNEG(NLIN+JNONL)=NONNEG(NLIN+JNONL+NLESS)
  560       CONTINUE
         END IF
      ELSE
         NCOL=NLIN
      END IF
C     -------------------------------------------------------------------------
C     Compute optimal SOLUTN with PNNLS.
C     -------------------------------------------------------------------------
      CALL PNNLS (DAMAT, MROW, NROW, NCOL, DRHS, SOLUTN, OBJTOT, DWORK,
     1            DWORK(NPAR+1), INDCOL, IERROR, DRANGE, NONNEG, 0.D0,
     2            NDF)
      LERROR=IERROR.NE.1
      IF (.NOT.LERROR) THEN
         IF (PMQACT.GT.0.D0 .AND. DONONL) THEN
C           -------------------------------------------------------------------
C           Subtract Marquardt penalty from OBJTOT.
C           -------------------------------------------------------------------
            INONL=0
            DO 610 JNONL=1,NNONL
               IF ((FXDEGZ .AND. JNONL.EQ.LPHAST) .OR.
     1             (FXDEGP .AND. JNONL.EQ.LPHAST+1)) GO TO 610
               INONL=INONL+1
               OBJTOT=OBJTOT-(DAMARQ(INONL)*SOLUTN(NLIN+INONL))**2
  610       CONTINUE
         END IF
         IF (NLESS .GT. 0   .and.   dononl) THEN
C           -------------------------------------------------------------------
C           Restore (expand) NONNEG, SOLUTION & GRAD.
C           -------------------------------------------------------------------
            DO 612 JNONL=NNONL,LSHIST,-1
               NONNEG(NLIN+JNONL)=NONNEG(NLIN+JNONL-NLESS)
               SOLUTN(NLIN+JNONL)=SOLUTN(NLIN+JNONL-NLESS)
               GRAD(JNONL)=GRAD(JNONL-NLESS)
 612        CONTINUE
            IF (FXDEGZ .AND. .NOT.FXDEGP) THEN
               NONNEG(NLIN+LPHAST+1)=NONNEG(NLIN+LPHAST)
               SOLUTN(NLIN+LPHAST+1)=SOLUTN(NLIN+LPHAST)
               GRAD(LPHAST+1)=GRAD(LPHAST)
            END IF
            IF (FXDEGZ) THEN
               NONNEG(NLIN+LPHAST)=.FALSE.
               SOLUTN(NLIN+LPHAST)=0.D0
               GRAD(LPHAST)=0.D0
            END IF
            IF (FXDEGP) THEN
               NONNEG(NLIN+LPHAST+1)=.FALSE.
               SOLUTN(NLIN+LPHAST+1)=0.D0
               GRAD(LPHAST+1)=0.D0
            END IF
         END IF
         IF (.NOT.DONONL .AND. SDREF.LT.DRANGE .AND. LSTAGE.EQ.2) THEN
C           -------------------------------------------------------------------
C           Add penalties for shift and 1/T2 priors to OBJTOT.
C           -------------------------------------------------------------------
            jcoeff_power = 0
            DO 615 JMETAB=1,NMETAB
               OBJTOT=OBJTOT+(PARNLN(LSHIST-1+JMETAB)*SDREF/
     1                        SDSHIF(JMETAB))**2
               if (imethd .eq. 3) then
                  do 616 kpower = 1, npower(jmetab)
                     jcoeff_power = jcoeff_power + 1
                     jnonl = LRT2ST-1+jcoeff_power
                     objtot = objtot + (parnln(jnonl) * sdref /
     1                        coeff_power_sd(kpower, jmetab))**2
 616              continue
               else
                  objtot = objtot + ((EXRT2(JMETAB)-
     3                         PARNLN(LRT2ST-1+JMETAB))*SDREF/
     4                        SDRT2(JMETAB))**2
               end if
 615        CONTINUE
            do 617 jrow_group_shift = 1, nrow_group_shift
               dterm(1) = 0.d0
               do 618 jmetab = 1, nmetab
                  jnonl = lshist + jmetab - 1
                  dterm(1) = dterm(1) + parnln(jnonl) *
     1                       cgroup_shift(jrow_group_shift, jmetab)
 618           continue
               objtot = objtot + (sdref * dterm(1) /
     1                           sdgroup_shift_row(jrow_group_shift))**2
 617        continue
         END IF
         IF (.NOT.DONONL .AND. SDREF.LT.DRANGE) THEN
C           -------------------------------------------------------------------
C           Add penalties for priors for phases to OBJTOT.
C           (The expressions below are simplified from the original
C             EXDEG?*RADIAN and SDDEG?*RADIAN.)
C           -------------------------------------------------------------------
            IF (PRIORZ) OBJTOT=OBJTOT+((EXDEGZ-(PARNLN(LPHAST)+
     1                         PHITOT(1))/RADIAN)*SDREF/SDDEGZ)**2
            IF (PRIORP) OBJTOT=OBJTOT+((EXDEGP-(PARNLN(LPHAST+1)+
     1                         PHITOT(2))/RADIAN)*SDREF/SDDEGP)**2
         END IF
C        ----------------------------------------------------------------------
C        DSSQ = OBJTOT - (regularizor penalty)
C        ----------------------------------------------------------------------
         DSSQ= OBJTOT-PENLTY(ALPHAB, 0.D0, SOLUTN, PARNLN)
         if (ALPHAS .GT. 0.D0   .and.   imethd .eq. 2   .and.
     1       lstage .eq. 2   .and.   sdref .lt. drange) then
            if (dononl) then
               do 630 jmetab = 1, nmetab
                  jnonl = lrt2st - 1 + jmetab
                  dwork(jnonl) = solutn(nlin + jnonl) + parnln(jnonl)
 630           continue
c              ---------------------------------------------------------------
c              Get identical final results if DPY is replaced with PARNLN in
c                 next statement (or if it is skipped).
c              Here IMETHD=2, and ALPHAS component is already in the
c                 solution when DONONL=F; so OBJTOT+=PENLTY(ALPHAS) should not
c                 be done.  Otherwise, much worse, often no convergence.
c              OBJ* affect the output of "Rel. dec." & "Pred. dec." in and
c                   the Marquardt criteria in PLINLS.  Negative "Rel. dec."
c                   occur frequently in Prel.  It also determines RSTEP.
c                   When IMETHD=2, PENLTY(ALPHAS) is already in OBJTOT from
c                   PNNLS and does not need to be added in below.
c              DSSQ determines the interpolations for ALPHA in RFALSI and
c                   the convergence criterion in PLINLS.  However, analysis
c                   is independent of following statement; PENLTY(...ALPS...)
c                   only affects output penB; i.e., different PENLTY function
c                   for ALPHAS yields identical solution.  So, DSSQ is
c                   apparently only used when DONONL=F.
c              ---------------------------------------------------------------
              DSSQ=DSSQ-PENLTY(0.D0, ALPHAS, SOLUTN, dwork)
 
c              dterm(1) = PENLTY(0.D0, ALPHAS, SOLUTN, dwork)
c              write (lprint, 9630) dssq,
c     1            dterm(1),
c     2            solutn(nlin + jnonl),
c     1            dwork(jnonl), parnln(jnonl), nlin, jnonl, jmetab
c 9630         format (1p5e13.3, 3i5)
c              if (nmetab .ne. 98) stop
 
            else
c               OBJTOT=OBJTOT+PENLTY(0.D0, ALPHAS, SOLUTN, PARNLN)
C              DSSQ=DSSQ-PENLTY(0.D0, ALPHAS, SOLUTN, parnln)!<<<<
            end if
         end if
         IF (NSIDE2.GT.0 .AND. ALPHAS.GT.0.D0 .and. imethd.ne.2) THEN
            IF (DONONL) THEN
               DO 640 JSIDE2=1,NSIDE2
                  DPY(JSIDE2)=SOLUTN(NLIN+JSIDE2)+PARNLN(JSIDE2)
  640          CONTINUE
               DSSQ=DSSQ-PENLTY(0.D0, ALPHAS, SOLUTN, DPY)
            ELSE
               OBJTOT=OBJTOT+PENLTY(0.D0, ALPHAS, SOLUTN, PARNLN)
            END IF
         END IF
         IF (DONONL) THEN
            OBJLIN=OBJTOT
         ELSE
            OBJECT=OBJTOT
C           -------------------------------------------------------------------
C           Correct NDF for KNONL nonlinear parameters.
C           -------------------------------------------------------------------
            NDF=NDF+KNONL
            IF (NROWDA.GT.NDF .AND. DSSQ.GT.0.D0) THEN
               STDDEV=DSQRT(DSSQ/DBLE(NROWDA-NDF))
C              ----------------------------------------------------------------
C              SDREF = rough estimate for SD of fit.  It is continually
C                      updated to the minimum while LSTAGE=1 and during the
C                      Starting Solution when LSTAGE=2 (i.e., if INISOL=T).
C              SDREF is used for weighting the priors, the regularizors, and
C                the Marquardt rows.
C              The continuous updating of SDREF will cause an extra decrease
C                in OBJECT during the Reference Solution, and this may require
C                extra iterations.
C              ----------------------------------------------------------------
               IF (INISOL) SDREF=DMIN1(SDREF,STDDEV)
            ELSE
               CALL ERRMES (8, 3, CHSUBP)
            END IF
         END IF
      END IF
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SAVBES (ILEVEL)
C
C  ILEVEL = 1 for saving the best solution from a Regula Falsi analysis for a
C             given alpha-ratio.
C         = 2 for saving the best of all solutions so far for the final output
C             and also loading PARNLN with PARBES.
c          =3 for saving the rephased solution with the max ALPHAB
c          =-3 for loading the *(3) values into *(2)
C
      INCLUDE 'lcmodel.inc'
      CHSUBP='SAVBES'
      IF (ILEVEL .EQ. 1) THEN
         PENBES(1)=OBJECT-DSSQ
         IF (SNGL(ALPHAB).LT.ALPBPN .AND. SDREF.LT.DRANGE) PENBES(1)=
     1        PENBES(1)+PNALPB*SDREF**2
         SSQBES(1)=DSSQ
         SDBEST(1)=STDDEV
         PMQBES(1)=PMQSAV
         ALPBBS(1)=ALPHAB
         ALPSBS(1)=ALPHAS
         DO 110 JNONL=1,NNONL
            PARBES(JNONL,1)=PARNLN(JNONL)
  110    CONTINUE
         DO 120 JPAR=1,NLIN
            SOLBES(JPAR,1)=SOLUTN(JPAR)
  120    CONTINUE
         phitot_sav(1, 1) = phitot(1)
         phitot_sav(1, 2) = phitot(2)
         do 130 j = 1, ny
            cy_sav(j, 1) = cy(j)
 130     continue
 
      ELSE IF (ILEVEL .EQ. 2) THEN
         PENBES(2)=PENBES(1)
         SSQBES(2)=SSQBES(1)
         SDBEST(2)=SDBEST(1)
         PMQBES(2)=PMQBES(1)
         ALPBBS(2)=ALPBBS(1)
         ALPSBS(2)=ALPSBS(1)
         DO 210 JNONL=1,NNONL
            PARBES(JNONL,2)=PARBES(JNONL,1)
            PARNLN(JNONL)=PARBES(JNONL,1)
  210    CONTINUE
         DO 220 JPAR=1,NLIN
            SOLBES(JPAR,2)=SOLBES(JPAR,1)
  220    CONTINUE
         phitot_sav(2, 1) = phitot_sav(1, 1)
         phitot_sav(2, 2) = phitot_sav(1, 2)
         do 230 j = 1, ny
            cy_sav(j, 2) = cy_sav(j, 1)
 230     continue
 
      else
         if (iabs(ilevel) .lt. 3   .or.
     1       iabs(ilevel) .gt. mmdegp3 + 7) call errmes (1, 5, chsubp)
         if (ilevel .gt. 0) then
            PENBES(ilevel)=PENBES(2)
            SSQBES(ilevel)=SSQBES(2)
            SDBEST(ilevel)=SDBEST(2)
            PMQBES(ilevel)=PMQBES(2)
            ALPBBS(ilevel)=ALPBBS(2)
            ALPSBS(ilevel)=ALPSBS(2)
            DO 310 JNONL=1,NNONL
               PARBES(JNONL,ilevel)=PARBES(JNONL,2)
               PARNLN(JNONL)=PARBES(JNONL,2)
 310        CONTINUE
            DO 320 JPAR=1,NLIN
               SOLBES(JPAR,ilevel)=SOLBES(JPAR,2)
 320        CONTINUE
            phitot_sav(ilevel, 1) = phitot_sav(2, 1)
            phitot_sav(ilevel, 2) = phitot_sav(2, 2)
            do 330 j = 1, ny
               cy_sav(j, ilevel) = cy_sav(j, 2)
 330        continue
         else
            PENBES(2)=PENBES(-ilevel)
            SSQBES(2)=SSQBES(-ilevel)
            SDBEST(2)=SDBEST(-ilevel)
            PMQBES(2)=PMQBES(-ilevel)
            ALPBBS(2)=ALPBBS(-ilevel)
            ALPSBS(2)=ALPSBS(-ilevel)
            DO 350 JNONL=1,NNONL
               PARBES(JNONL,2)=PARBES(JNONL,-ilevel)
               PARNLN(JNONL)=PARBES(JNONL,-ilevel)
 350        CONTINUE
            DO 360 JPAR=1,NLIN
               SOLBES(JPAR,2)=SOLBES(JPAR,-ilevel)
 360        CONTINUE
            phitot_sav(2, 1) = phitot_sav(-ilevel, 1)
            phitot_sav(2, 2) = phitot_sav(-ilevel, 2)
            do 370 j = 1, ny
               cy_sav(j, 2) = cy_sav(j, -ilevel)
 370        continue
         end if
      end if
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FINOUT ()
C
C  Output final summary.
C
      INCLUDE 'lcmodel.inc'
c     -----------------------------------------------------------------------
c     Excel only allows up to 256 columns; so must have MCSV <= 254/3 = 84 or
c        allow fewer than the three types (CONC, %SD, /NAMREL).
c     -----------------------------------------------------------------------
      external ilen, ldegmx
      parameter (mcsv=84)
      CHARACTER csv_element*56, csv_line(2)*(mcsv * 56 + 10),
     1          FMTC*8, fmtdat*(mchfmt), fmtpm*10, FMTR*6, id*(mchid)
      complex cfactor
      DOUBLE PRECISION DAPOSI(MPAR,MPAR), DSUM, dwork(mdwork_finout),
     1                 PHIOLD(2)
      integer nchar_metab(2)
      LOGICAL bruker, EFORM, INCLUD(MPAR), LDEGMX, LERROR, first_row,
     1        seqacq, skip_line
      REAL ERRCON(MCONC)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON DAPOSI, INCLUD, ERRCON, csv_line, dwork
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      save first_row
      data first_row/.true./
      namelist /seqpar/ hzpppm
      namelist /nmid/ bruker, fmtdat, id, seqacq, tramp, volume
      EFORM(TERM)=(abs(TERM).GT.999.999 .OR. ABS(TERM).lt..0995) .AND.
     1            TERM.NE.0.
      CHSUBP='FINOUT'
C     -------------------------------------------------------------------------
C     Rephase data and repeat analysis for plotting.
C     -------------------------------------------------------------------------
      IF (DOFULL) THEN
         ISTAGE=3
      ELSE
         ISTAGE=1
      END IF
      LSTAGE=MIN0(ISTAGE,2)
      DO 110 JREPHA=1,MREPHA(2)
         IF (.NOT.LDEGMX(2)) GO TO 140
         IF (LPRINT .GT. 0) WRITE (LPRINT,5110) JREPHA
 5110    FORMAT (////20X,
     1               'Rephase data for plotting and repeat analysis',
     2                5X, 'Iteration', I2)
         ALPHAB=ALPBBS(2)
         ALPHAS=ALPSBS(2)
         PHIOLD(1)=PARBES(LPHAST,2)
         PHIOLD(2)=PARBES(LPHAST+1,2)
         CALL REPHAS ()
         CALL PLINLS (ISTAGE, IERROR)
         IF (OBJECT .GE. DRANGE) THEN
            CALL ERRMES (1, 3, CHSUBP)
            PARBES(LPHAST,2)=-PHIOLD(1)
            PARBES(LPHAST+1,2)=-PHIOLD(2)
            CALL REPHAS()
            PARBES(LPHAST,2)=PHIOLD(1)
            PARBES(LPHAST+1,2)=PHIOLD(2)
            GO TO 140
         ELSE
            CALL SAVBES (1)
            CALL SAVBES (2)
         END IF
  110 CONTINUE
      IF (LDEGMX(2) .AND. MREPHA(2).GT.0) CALL ERRMES (2, 3, CHSUBP)
C     -------------------------------------------------------------------------
C     Reproduce best solution and matrix for covariance from values from
C       SAVBES.
C     Covariance matrix will actually correspond to the first stage of the
C       next iteration with DONONL=T.  Therefore, NONNEG is reset below to
C       remove the non-binding constraints.
C     The following are preserved throughout and do not have to be saved or
c       restored: EXDEGP & EXDEGZ
C     -------------------------------------------------------------------------
  140 DO 150 JPAR=1,NLIN
         SOLUTN(JPAR)=SOLBES(JPAR,2)
         NONNEG(JPAR)=SOLUTN(JPAR) .EQ. 0.D0
  150 CONTINUE
      DO 160 JNONL=1,NNONL
         PARNLN(JNONL)=PARBES(JNONL,2)
  160 CONTINUE
      do 161 j = 1, ny
         cy(j) = cy_sav(j, 2)
 161  continue
      phitot(1) = phitot_sav(2, 1)
      phitot(2) = phitot_sav(2, 2)
      ALPHAB=ALPBBS(2)
      ALPHAS=ALPSBS(2)
c     ------------------------------------------------------------------------
c     Save solution with SAVBES(3).  This overwrites old SAVBES(3) from
c        max ALPHAB solution, which is not needed anymore.
c     ------------------------------------------------------------------------
      CALL SOLVE (LSTAGE, .TRUE., PMQBES(2), .FALSE., LERROR)
      IF (LERROR) CALL ERRMES (3, 4, CHSUBP)
C     -------------------------------------------------------------------------
C     Now call SOLVE to get the following for plotting:
C     YREAL(JY) = real part of data, JY=1,NYuse
C     YFITRE(JY,0) =           fit
C     YFITRE(JY,JMETAB) = contribution of metabolite JMETAB and background to
C                         fit.
C     BACKRE(JY) = background
C     -------------------------------------------------------------------------
      DO 165 JPAR=1,NLIN
         SOLUTN(JPAR)=SOLBES(JPAR,2)
  165 CONTINUE
      CALL SOLVE (LSTAGE, .FALSE., 0.D0, .TRUE., LERROR)
      IF (LERROR) CALL ERRMES (4, 4, CHSUBP)
      ISTAGO=3
      IF (NSIDES .GT. 0) THEN
C        ----------------------------------------------------------------------
C        Print out lineshape smoothed with DGAUSS and use this to estimate
C          FWHMLS = approx. fwhm of lineshape.
C        DWORK(NSIDE2+3+1+JNONL) = smoothed lineshape at point JNONL.
C        The ordering of the lineshape is reversed for the plot.
C        ----------------------------------------------------------------------
         J=NEXTRE (PARNLN, NSIDE2, DWORK, DGAUSS, THRLIN, imethd)
         JNONL=NSIDE2+2
         DO 180 J=1,NSIDE2+1
            CPY(J)=INCSID*FLOAT(NSIDES+1-J)
            JNONL=JNONL-1
            CPY(NSIDE2+1+J)=SNGL(DWORK(NSIDE2+4+JNONL))
  180    CONTINUE
         HALFMX=-RRANGE
         DO 200 J=1,NSIDE2+1
            HALFMX=amax1(CPY(NSIDE2+1+J),halfmx)
  200    CONTINUE
         HALFMX=.5*HALFMX
         DO 210 jleft=1,NSIDE2+1
            IF (CPY(NSIDE2+1+jleft) .gE. HALFMX) GO TO 215
  210    CONTINUE
  215    DO 220 jright=NSIDE2+1,1,-1
            IF (CPY(NSIDE2+1+jright) .gE. HALFMX) GO TO 225
  220    CONTINUE
  225    FWHMLS=FLOAT(INCSID*(Jright-jleft+1))*PPMINC
         IF (LPRINT .GT. 0) THEN
            WRITE (LPRINT,5180) FWHMLS
 5180       FORMAT (////' Lineshape coefficients', 10X,
     1                  '(Approx. FWHM =', F6.3, ' ppm)')
            CALL PLPRIN (CPY, CPY(NSIDE2+1+1), CPY, NSIDE2+1, .TRUE.,
     1                   LPRINT, RRANGE, 0, 0, 0, CPY, .FALSE.)
         END IF
      END IF
C     -------------------------------------------------------------------------
C     Start compuation of error estimates and correlation coefficients.
C     First invert the NDF**2 upper triangular factor in DAMAT from PNNLS
C       onto itself and then square it using algorithm from Lawson & Hanson
C       pp. 68-69 to form the covariance matrix.
C     INDCOL(JDF) = the column no. of the JDFth parameter in the positive
C                     set; i.e., the column down to the diagonal has JDF
C                     elements.  INDCOL comes from PNNLS.
C     DAPOSI(JDF,KDF) = covariance between SOLUTN(INDCOL(JDF)) and
C                                          SOLUTN(INDCOL(KDF)).
C     -------------------------------------------------------------------------
      DO 310 I=1,NDF
         II=INDCOL(I)
         DAMAT(I,II)=1.D0/DAMAT(I,II)
  310 CONTINUE
      DO 320 I=1,NDF-1
         DO 330 J=I+1,NDF
            IJ=INDCOL(J)
            DSUM=0.D0
            DO 340 L=I,J-1
               IL=INDCOL(L)
               DSUM=DSUM+DAMAT(I,IL)*DAMAT(L,IJ)
  340       CONTINUE
            DAMAT(I,IJ)=-DAMAT(J,IJ)*DSUM
  330    CONTINUE
  320 CONTINUE
      DO 350 I=1,NDF
         DO 360 J=I,NDF
            DSUM=0.D0
            DO 370 L=J,NDF
               IL=INDCOL(L)
               DSUM=DSUM+DAMAT(I,IL)*DAMAT(J,IL)
  370       CONTINUE
            DAPOSI(I,J)=DSUM
            DAPOSI(J,I)=DSUM
  360    CONTINUE
  350 CONTINUE
C     -------------------------------------------------------------------------
C     CONCEN(NCONC-by-1) = BMAT(NCONC-by-NDF) * SOLUTN(NDF-by-1).
C     BMAT(JCONC,JDF) = 1 or 0, and is not expicitly formed, but rather
C                          determined by LCOMPO.
C     Cov_CONCEN = covariance between CONCEN(JCONC) and CONCEN(KCONC),
C                  JCONC.LE.KCONC
C                = BMAT * Cov_SOLUTN *BMAT**T, where Cov_SOLUTN is in DAPOSI.
C     DAMAT(JDF,KDF) = Cov_SOLUTN * BMAT**T
C     -------------------------------------------------------------------------
      DO 410 ICOL=1,NCONC
         DO 430 JDF=1,NDF
            DO 440 JCOMPO=1,NCOMPO(ICOL)
               INCLUD(JDF)=INDCOL(JDF) .EQ. LCOMPO(JCOMPO,ICOL)
               IF (INCLUD(JDF)) GO TO 430
  440       CONTINUE
  430    CONTINUE
         DO 450 IROW=1,NDF
            DSUM=0.D0
            DO 460 JDF=1,NDF
               IF (INCLUD(JDF)) DSUM=DSUM+DAPOSI(IROW,JDF)
  460       CONTINUE
            DAMAT(IROW,ICOL)=DSUM
  450    CONTINUE
  410 CONTINUE
C     -------------------------------------------------------------------------
C     DAPOSI(JCONC,KCONC) = covariance between CONCEN(JCONC) and
C                                              CONCEN(KCONC), JCONC.LE.KCONC
C                        = BMAT * DAMAT
C     -------------------------------------------------------------------------
      DO 470 IROW=1,NCONC
         DO 480 JDF=1,NDF
            DO 490 JCOMPO=1,NCOMPO(IROW)
               INCLUD(JDF)=INDCOL(JDF) .EQ. LCOMPO(JCOMPO,IROW)
               IF (INCLUD(JDF)) GO TO 480
  490       CONTINUE
  480    CONTINUE
         DO 510 ICOL=IROW,NCONC
            DSUM=0.D0
            DO 520 JDF=1,NDF
               IF (INCLUD(JDF)) DSUM=DSUM+DAMAT(JDF,ICOL)
  520       CONTINUE
            DAPOSI(IROW,ICOL)=DSUM
  510    CONTINUE
  470 CONTINUE
C     -------------------------------------------------------------------------
C     ERRCON(JCONC) = standard error estimate for CONCEN(JCONC) (temporarily
C                     not scaled by SDBEST(2) until after the
C                     correlation coefficients have been calculated).
C     -------------------------------------------------------------------------
      DO 530 JCONC=1,NCONC
         CONCEN(JCONC)=0.
         DO 540 JCOMPO=1,NCOMPO(JCONC)
            CONCEN(JCONC)=CONCEN(JCONC)+SOLBES(LCOMPO(JCOMPO,JCONC),2)
  540    CONTINUE
         TERM=DAPOSI(JCONC,JCONC)
         IF (TERM .LT. 0.) CALL ERRMES (5, 4, CHSUBP)
         ERRCON(JCONC)=SQRT(TERM)
  530 CONTINUE
C     -------------------------------------------------------------------------
C     DAPOSI(IROW,ICOL) = correlation coefficient between CONCEN(IROW) and
C                         CONCEN(ICOL), IROW<ICOL.
C     The output is in the form of the lower triangle.
C     -------------------------------------------------------------------------
 5560 FORMAT (//////20X,'Correlation coefficients')
      IF (LPRINT .GT. 0) WRITE (LPRINT,5560)
 5565 FORMAT (/(13X,17(1X,A6)))
      IF (LPRINT .GT. 0) WRITE (LPRINT,5565) (NACOMB(IROW),
     1                                        IROW=1,NCONC-1)
      DO 560 ICOL=1,NCONC
         DO 570 IROW=1,ICOL-1
            TERM=ERRCON(IROW)*ERRCON(ICOL)
            IF (TERM .LE. 0.) THEN
               DAPOSI(IROW,ICOL)=0.
            ELSE
               DAPOSI(IROW,ICOL)=DAPOSI(IROW,ICOL)/TERM
            END IF
  570    CONTINUE
 5570    FORMAT (1X,A11,17F7.3/(12X,17F7.3))
         IF (ICOL.GT.1 .AND. LPRINT .GT. 0) WRITE (LPRINT,5570)
     1         NACOMB(ICOL), (DAPOSI(IROW,ICOL),IROW=1,ICOL-1)
  560 CONTINUE
 5571 FORMAT ((13X,17(1X,A6)))
      IF (LPRINT .GT. 0) WRITE (LPRINT,5571) (NACOMB(IROW),
     1                                        IROW=1,NCONC-1)
C     -------------------------------------------------------------------------
c     Output number of SDs that CONC ratios deviate from prior expectation
c       values.
c     Priors are of form in 3rd line:
c       c_num / c_sum = exrati +- sdrati
c       c_num - c_sum * exrati = 0 +- c_sum * sdrati
c       [1 / (sdrati * csum)] [c_num - c_sum_true * exrati] = 0 +- 1
c       where (the approximate) CSUM is from the perliminary analysis with
c       ALP*ST.
c     FSD_CORR = sd_true / SD_USED = CSUM / CSUM_TRUE = factor to multiply
c               SD to correct for this approximation.
c     FSD_CORR < 1 means that the weighting was too strong.
C     -------------------------------------------------------------------------
      if (min0(nratio_used, lprint, ipdump - 2) .gt. 0) then
         write (lprint, 5572)
 5572    format (///5x, 'SDs', 3x, 'Corr Factor', 3x,
     1           'Concentration Ratio')
         do 572 jratio = 1, nratio_used
            csum = cprior(jratio, lmetab_prior(jratio)) *
     1             sdrati(jratio)
c           -------------------------------------------------------------------
c           Error -- neither of the above terms can be zero.
c           -------------------------------------------------------------------
            if (csum .le. 0.) call errmes (12, 5, chsubp)
            csum = 1. / csum
            prior_sum = sqrtwt_ratio_used(jratio) *
     1                  concen(lmetab_prior(jratio))
            sd_used = 0.
            do 573 jmetab = 1, nmetab
               term = cprior(jratio, jmetab) * concen(jmetab)
               sd_used = sd_used + term
               prior_sum = prior_sum - term
 573        continue
            if (prior_sum .le. 0.) then
               fsd_corr = 999.
            else
c              ----------------------------------------------------------------
c              csum_true = prior_sum * sdrati(jratio) * csum / exrati(jratio)
c              fsd_corr = csum / csum_true
c              ----------------------------------------------------------------
               fsd_corr = exrati(jratio) / (prior_sum * sdrati(jratio))
            end if
            len_chrato = ilen(chrato(jratio))
            write (lprint, 5573) sd_used, fsd_corr,
     1                           chrato(jratio)(:len_chrato)
 5573       format (f8.1, g14.2, 3x, a)
 572     continue
      end if
C     -------------------------------------------------------------------------
C     Output table of CONCEN and errors.  Load TABLE for plotting.
C     PCERR = standard percent error estimate for CONCEN(JCONC).
C     -------------------------------------------------------------------------
      NTABLE=NCONC+1
      IF (LCOORD .GT. 0) WRITE (LCOORD,5581) NTABLE
 5581 FORMAT (1X,I2,' lines in following concentration table = NCONC+1')
      IF (LTABLE .GT. 0) WRITE (LTABLE,6581) NTABLE
 6581 FORMAT (/'$$CONC',I3,
     1        ' lines in following concentration table = NCONC+1')
      LINTBL=1
      IPCERR(1)=999
      ratio_outlier(1) = .false.
      FACREL=0.
      DO 574 JCONC=1,NCONC
         IF ((NAMREL.EQ.NACOMB(JCONC)   .or.
     1       namrel .eq. nacom2(jconc))   .and.    CONCEN(JCONC).GT.0.)
     2      THEN
            FACREL=CONREL/CONCEN(JCONC)
            GO TO 575
         END IF
  574 CONTINUE
c     -------------------------------------------------------------------------
c     NAMREL not in NACOMB.  Check for synonyms
c     -------------------------------------------------------------------------
      do 5742 jconc = 1, nconc
         do 5744 jsyn = 1, mpmet
            IF (SYNUS1(1,JSYN).EQ.' ' .OR. SYNUS1(2,JSYN).EQ.' ')
     1           GO TO 5742
            IF (((SYNUS1(1,JSYN) .EQ. namrel   .AND.
     1            SYNUS1(2,JSYN) .EQ. nacomb(jconc))   .OR.
     2           (SYNUS1(2,JSYN) .EQ. nacomb(jconc)   .AND.
     3            SYNUS1(1,JSYN) .EQ. namrel))   .and.
     4          CONCEN(JCONC) .GT. 0.) then
              namrel = nacomb(jconc)
              FACREL=CONREL/CONCEN(JCONC)
              GO TO 575
           END IF
 5744   continue
 5742 continue
c     -----------------------------------------------------------------------
c     Use NAMREL='Cr+PCr' or 'Cre+PCr' instead of only Cr or Cre, if PCr is
c       also present.
c     This can be blocked by using CHCOMB(9)='PCr+Cr' instead of the default
c       'Cr+PCr'
c     -----------------------------------------------------------------------
  575 if (namrel .eq. 'Cr'  .or.  namrel .eq. 'Cre') then
         do 576 jconc = 1, nconc
            if ((nacomb(jconc) .eq. 'Cr+PCr'  .or.
     1           nacomb(jconc) .eq. 'Cre+PCr')  .and.
     2          concen(jconc) .gt. 0.) then
               namrel = nacomb(jconc)
               FACREL=CONREL/CONCEN(JCONC)
               go to 579
            end if
 576     continue
      end if
c     -----------------------------------------------------------------------
c     Use NAMREL='GPC+PCh+Cho' or as many as present, if fewer are in NAMREL.
c     This can be blocked by ordering the names differently in CHCOMB.
c     The ordering of the Metabolite Name pairs below must be that in CHCOMB.
c     -----------------------------------------------------------------------
      if (namrel .eq. 'PCh+GPC') namrel = 'GPC+PCh'
      if (namrel .eq. 'Cho+GPC') namrel = 'GPC+Cho'
      if (namrel .eq. 'Cho+PCh') namrel = 'PCh+Cho'
      if (namrel .eq. 'GPC+PCh'  .or.  namrel .eq. 'GPC+Cho'   .or.
     1    namrel .eq. 'PCh+Cho') then
          do 577 jconc = 1, nconc
            if (nacomb(jconc) .eq. 'Cho+GPC+PCh'  .and.
     1          concen(jconc) .gt. 0.) then
               namrel = nacomb(jconc)
               FACREL=CONREL/CONCEN(JCONC)
               go to 579
            end if
 577     continue
      end if
      if (namrel .eq. 'GPC'  .or.  namrel .eq. 'PCh'   .or.
     1    namrel .eq. 'Cho') then
         do 578 jconc = 1, nconc
            if ((nacomb(jconc) .eq. 'GPC+PCh'   .or.
     1           nacomb(jconc) .eq. 'GPC+Cho'   .or.
     2           nacomb(jconc) .eq. 'PCh+Cho'   .or.
     3           nacomb(jconc) .eq. 'Cho+GPC+PCh')   .and.
     4           concen(jconc) .gt. 0.) then
               namrel = nacomb(jconc)
               FACREL=CONREL/CONCEN(JCONC)
               go to 579
            end if
 578     continue
      end if
 579  do 5792 jpage = 1, 2
         nchar_metab(jpage) = nchlin(jpage) - 22
         if (nchar_metab(jpage) .lt. 12) then
            IF (NAMREL(5:6) .NE. '  ') THEN
               WRITE (TABLE(1,JPAGE),5582) NAMREL
 5582          FORMAT (3X, 'Conc.', 2X, '%SD', ' /', A6, 2X, 'Metab.')
            ELSE
               WRITE (TABLE(1,JPAGE),5583) NAMREL
 5583          FORMAT (3X, 'Conc.', 2X, '%SD', '   /', A4, 2X, 'Metab.')
            END IF
         else
            IF (NAMREL(5:6) .NE. '  ') THEN
               WRITE (TABLE(1,JPAGE),5584) NAMREL
 5584          FORMAT (3X, 'Conc.', 2X, '%SD', ' /', A6, 2X,
     1                 'Metabolite')
            ELSE
               WRITE (TABLE(1,JPAGE),5585) NAMREL
 5585          FORMAT (3X, 'Conc.', 2X, '%SD', '   /', A4, 2X,
     1                 'Metabolite')
            END IF
         end if
 5792 continue
      len_namrel = 0
      if (lcsv .gt. 0) then
         if (first_row   .and.   ioffset_current_in .le. 0) then
            csv_line(1) = 'Row, Col'
         end if
         len_namrel = ilen(namrel)
         write (csv_line(2), 5586) idrow, idcol
 5586    format (i3, ', ', i3)
      end if
      if (lprint .gt. 0) then
         if (.not.fxdegp   .and.   ndegppm3_used .gt. 0)
     1      write (lprint, 5603) (degp_degp(j), ssq_degp(j),
     2                            alpb_degp(j), dist_degp(j),
     3                           j = 0, ndegppm3_used)
 5603       format(///'  DEGPPM', 8x, 'SSQ', 5x,
     1             'ALPHAB   Spline Distance'/(0pf8.2, 1p2e11.2, e18.2))
 
         write (lprint, 5587)
 5587    format (///14x, '1000Shift', 3x, 'SDs', 2x,
     1              'delta(1/T2)', 3x, 'SDs', 8x, 'CONC',
     2              2x, '%SD')
      end if
      ncsv = 0
      DO 580 Jline = 1,NCONC
         lconc = iconc_line_table(jline)
         skip_line = onlyco   .and.   lconc .le. nmetab
         if (skip_line) then
            ipcerr(lintbl) = -99
            go to 581
         end if
c        --------------------------------------------------------------------
c        Output concentration table.
c        --------------------------------------------------------------------
         ncsv = ncsv + 1
         LINTBL=LINTBL+1
         ERRCON(LCONC)=ERRCON(LCONC)*SDBEST(2)
         RELCON=FACREL*CONCEN(LCONC)
         IF (CONCEN(LCONC) .LE. 0.) THEN
            concen(lconc) = 0.
            PCERR=9999.
         ELSE
            PCERR=100.*ERRCON(LCONC)/CONCEN(LCONC)
         END IF
         IPCERR(LINTBL)=NINT(AMIN1(PCERR,999.))
         IF (EFORM(CONCEN(LCONC))) THEN
            if (conc3f) then
               FMTC='SS1PE8.2'
            else
               FMTC='SS1PE8.1'
            end if
         ELSE
            FMTC='0PF8.3'
         END IF
         IF (EFORM(RELCON)) THEN
            FMTR='1PE8.1'
         ELSE
            FMTR='0PF8.3'
         END IF
         fmtpm = ', '' '', a)'
         do 582 jratio = 1, nratio_used
            if (lmetab_prior(jratio) .eq. lconc) go to 583
 582     continue
         go to 586
 583     sum = 0.
         do 584 jmetab = 1, nmetab
            sum = sum + cprior(jratio, jmetab) * concen(jmetab)
 584     continue
         ratipm = abs(ratipm)
         ratio_outlier(lintbl) = abs(sum) .gt. ratipm
         if (sum .gt. ratipm) then
            fmtpm = ', ''+'', a)'
         else if (sum .lt. -ratipm) then
            fmtpm = ', ''-'', a)'
         end if
c        ---------------------------------------------------------------------
c        NACOMB(J) is changed here, but only for J > NMETAB, and these are
c           needed again (and NACOMB is redefined with next CSI voxel)
c        ---------------------------------------------------------------------
 586     if (nacom2(lconc) .ne. ' ') nacomb(lconc) = nacom2(lconc)
         do 587 jpage = 1, 2
            WRITE (TABLE(LINTBL, jpage),
     1             '(' // FMTC // ',I4,''%'',' // FMTR // fmtpm )
     2            CONCEN(LCONC), IPCERR(LINTBL), RELCON,
     3            NACOMB(LCONC)(:nchar_metab(jpage))
 587     continue
c        -------------------------------------------------------------------
c        Print final table of shifts & broadenings for all lines, even lines
c           skipped in concentration table (where IPCERR was set to -99 above).
c        -------------------------------------------------------------------
 581     if (lprint .gt. 0   .and.   lconc .le. nmetab   .and.
     1       dofull) then
            if (imethd .eq. 3) then
               write (lprint, 5588) lconc, nacomb(lconc)(:6),
     1          parbes(lshist+lconc-1,2)/(2.*PI*HZPPPM),
     2          parbes(lshist+lconc-1,2)/sdshif(lconc),
     3          parbes(lpowen(lconc - 1) + 1, 2),
     4          parbes(lpowen(lconc - 1) + 1, 2) /
     5            coeff_power_sd(1, lconc),
     6          CONCEN(LCONC), IPCERR(LINTBL)
               do 588 jpower = 2, npower(lconc)
                  write (lprint, 5604)
     1               parbes(lpowen(lconc - 1) + jpower, 2),
     2               parbes(lpowen(lconc - 1) + jpower, 2) /
     3                  coeff_power_sd(jpower, lconc)
 5604             format (29x, f13.2, f6.1)
 588           continue
            else
               write (lprint, 5588) lconc, nacomb(lconc)(:6),
     1           parbes(lshist+lconc-1,2)/(2.*PI*HZPPPM),
     2           parbes(lshist+lconc-1,2)/sdshif(lconc),
     3           parbes(lrt2st+lconc-1,2),
     4           (parbes(lrt2st+lconc-1,2) - exrt2(lconc))/SDRT2(lconc),
     5           CONCEN(LCONC), IPCERR(LINTBL)
 5588          format(i3, 1x, 'met= ', a6, 3pf8.1, 0pf6.1, f13.2, f6.1,
     1                1pe12.2, i5)
            end if
         end if
c        -------------------------------------------------------------------
c        Output to CSV file.
c        -------------------------------------------------------------------
         if (skip_line) go to 580
         if (lcsv .gt. 0   .and.   ncsv .le. mcsv) then
            len_nacomb = min0(13, ilen(NACOMB(LCONC)))
            if (first_row   .and.   ioffset_current_in .le. 0) then
               len_line = ilen(csv_line(1))
               csv_line(1) = csv_line(1)(:len_line) // ', ' //
     1                       nacomb(lconc)(:len_nacomb) // ', ' //
     2                       nacomb(lconc)(:len_nacomb) // ' %SD, ' //
     3                       nacomb(lconc)(:len_nacomb) // '/' //
     4                       namrel(:len_namrel)
            end if
            len_line = ilen(csv_line(2))
            write (csv_element, '(' // fmtc // ''', '',i4, '', '',' //
     1                          fmtr // ')')
     2               concen(lconc), IPCERR(LINTBL), RELCON
            len_line = ilen(csv_line(2))
            csv_line(2) = csv_line(2)(:len_line) // ', ' // csv_element
         end if
  580 CONTINUE
 5590 FORMAT (1X,A)
      IF (LPRINT .GT. 0) THEN
 5003    FORMAT (//1X)
         WRITE (LPRINT,5003)
 5580    FORMAT (1X, A)
         WRITE (LPRINT,5580) VERSIO, (TITLE_line(j), j=1,nlines_title)
 5001    FORMAT (1X)
         WRITE (LPRINT,5001)
         WRITE (LPRINT,5590) (TABLE(J,2),J=1,LINTBL)
      END IF
      IF (LCOORD .GT. 0) WRITE (LCOORD,5590) (TABLE(J,2),J=1,LINTBL)
      IF (LTABLE .GT. 0) WRITE (LTABLE,5590) (TABLE(J,2),J=1,LINTBL)
      if (lcsv .gt. 0) then
         if (first_row   .and.   ioffset_current_in .le. 0) then
            first_row = .false.
            len_line = ilen(csv_line(1))
            write (lcsv, 5589) csv_line(1)(:len_line)
         end if
         len_line = ilen(csv_line(2))
         write (lcsv, 5589) csv_line(2)(:len_line)
 5589    format(a)
      end if
C     -------------------------------------------------------------------------
C     Load ETCOUT with further output.
C     LINETC = no. of lines in ETCOUT.
C     -------------------------------------------------------------------------
      IF (SDBEST(2) .GT. 0.D0) THEN
         SIGMAX=0.
         DO 590 JY=1,NYuse
            IF (SUBBAS) THEN
               SIGMAX=AMAX1(SIGMAX, abs(YFITRE(JY,0)))
            ELSE
               SIGMAX=AMAX1(SIGMAX, abs(YFITRE(JY,0)-BACKRE(JY)))
            END IF
  590    CONTINUE
         SIGTON=.5*SIGMAX/SDBEST(2)
      ELSE
         CALL ERRMES (6, 2, CHSUBP)
         SIGTON=-1.
      END IF
      if (nsides .gt. 0) then
         WRITE (ETCOUT,5591) fwhmls, NINT(SIGTON)
      else
         WRITE (ETCOUT,5591) SQRT(fwhmst_full**2+FWHMBA**2),
     1                       NINT(SIGTON)
      end if
 5591 FORMAT (' FWHM =', F6.3, ' ppm', 4X, 'S/N =', I4)
      SHIFD=PPMINC*FLOAT(ISHIFD)
      WRITE (ETCOUT(2),5593) SHIFD
 5593 FORMAT (' Data shift =', F6.3, ' ppm')
      WRITE (ETCOUT(4),5594) ALPBBS(2), ALPSBS(2)
 5594 FORMAT (' alphaB,S =',1PE8.1, ',', E10.1)
      IF (ALPBBS(2).LE.DBLE(ALPBPN) .AND. ALPBBS(2).gt.0.d0 .and.
     1    NBACKG.GE.10) CALL ERRMES (7, 1, CHSUBP)
c     -------------------------------------------------------------------------
c     Put PHITOT(1) in range +-180
c     -------------------------------------------------------------------------
      PHITOT(1)=DMOD((PHITOT(1)+PARBES(LPHAST,2))/RADIAN, 360.D0)
      IF (PHITOT(1) .LT. -180.d0) PHITOT(1)=PHITOT(1)+360.
      if (phitot(1) .gt.  180.d0) PHITOT(1)=PHITOT(1)-360.
      PHITOT(2)=(PHITOT(2)+PARBES(LPHAST+1,2))/RADIAN
      IF (PHITOT(2) .lt. dgppmn - sddegp  .or.
     1    PHITOT(2) .gt. dgppmx + sddegp) CALL ERRMES (8, 2, CHSUBP)
      WRITE (ETCOUT(3),5595) NINT(PHITOT(1)), PHITOT(2)
 5595 FORMAT (' Ph:', I4, ' deg', f10.1, ' deg/ppm')
      WRITE (ETCOUT(5),5596) NBACKG, NSIDES, INCSID
 5596 FORMAT (1X,I3, ' spline knots.', 3X, 'Ns =', I2, '(', I1, ')')
      LINETC=5
      IF (NSIDE2 .GT. 0) THEN
         IF (NINFL(2) .LT. 0) THEN
            IF (NEXTR(1) .GT. 1) CALL ERRMES (9, 1, CHSUBP)
            WRITE (ETCOUT(6),5597) NINFL(1), NEXTR(1)
 5597       FORMAT (1X, I2, ' inflections.',  I6, ' extrema')
         ELSE
            IF (NEXTR(2) .GT. 1) CALL ERRMES (9, 1, CHSUBP)
            WRITE (ETCOUT(6),5598) NINFL, NEXTR
 5598       FORMAT (' (', I2, ')', I2, ' infls.', 2X,
     1              ' (', I2, ')', I2, ' extrs.')
         END IF
         LINETC=6
      END IF
      if (.NOT.(FXDEGZ .OR. SDDEGZ.GE.45.)) then
c        ----------------------------------------------------------------------
c        Manual says that EXDEGZ (input as DEGZER) must be in range 0--360.
c          Bring PHITOT(1) back to this range.
c        ----------------------------------------------------------------------
         IF (PHITOT(1) .LT. 0.) PHITOT(1)=PHITOT(1)+360.
         TEST=AMOD(ABS(phitot(1)-EXDEGZ),360.)
         IF (AMIN1(TEST,360.-TEST).GT.4.*SDDEGZ)
     1        call errmes (10, 2, chsubp)
      end if
      if (.not.fxdegp) then
         if (abs(phitot(2) - exdegp) .gt. 6. * sddegp)
     1        call errmes (11, 2, chsubp)
      end if
      IF (LPRINT .GT. 0) THEN
         WRITE (LPRINT,5003)
         WRITE (LPRINT,5590) (ETCOUT(J),J=1,LINETC)
      END IF
      IF (LCOORD .GT. 0) THEN
         WRITE (LCOORD,5602) LINETC
 5602    FORMAT (1X,I3, ' lines in following misc. output table')
         WRITE (LCOORD,5590) (ETCOUT(J),J=1,LINETC)
C        ----------------------------------------------------------------------
C        Finish output onto unit LCOORD.
C        ----------------------------------------------------------------------
         WRITE (LCOORD,5610) NYuse, (PPM(JY),JY=1,NYuse)
 5610    FORMAT (1X,I4,' points on ppm-axis = NY'/(1X,10F13.6))
         WRITE (LCOORD,5620) (YREAL(JY),JY=1,NYuse)
 5620    FORMAT (' NY phased data points follow'/(1X,1P10E13.5))
         WRITE (LCOORD,5630) (YFITRE(JY,0),JY=1,NYuse)
 5630    FORMAT (' NY points of the fit to the data follow'/
     1           (1X,1P10E13.5))
         IF (NBACKG .GT. 0) THEN
C           -------------------------------------------------------------------
C           Output background.
C           -------------------------------------------------------------------
 5640       FORMAT (' NY background values follow'/(1X,1P10E13.3))
            WRITE (LCOORD,5640) (BACKRE(JY),JY=1,NYuse)
         END IF
         IF (NEACH .GE. 1) THEN
            DO 660 JMETAB=1,NMETAB
               IF (CONCEN(JMETAB) .le. 0.) go to 660
               do 670 JEACH=1,MIN0(MMETAB,NEACH)
                  IF (NAMEAC(JEACH).EQ.NACOMB(JMETAB)(1:6) .OR.
     1                NEACH.GT.NMETAB) THEN
                     WRITE (LCOORD,5660) NACOMB(JMETAB), CONCEN(JMETAB),
     1                                  (YFITRE(JY,JMETAB),JY=1,NYuse)
 5660                FORMAT (1X, A6, '   Conc. =', 1PE9.2/ (1X,10E13.5))
                     go to 660
                  end if
 670           continue
  660       CONTINUE
         END IF
      END IF
      IF (LTABLE .GT. 0) THEN
         WRITE (LTABLE,5662) LINETC
 5662    FORMAT (/'$$MISC',I3,
     1           ' lines in following misc. output table')
         WRITE (LTABLE,5590) (ETCOUT(J),J=1,LINETC)
      END IF
      if (lcoraw .gt. 0) then
C        ----------------------------------------------------------------------
C             Output corrected RAW file to FILCOR (which has already been
c        opened).
c             DATAT has already been scaled by FCALIB*TRAMP/VOLUME, ECC'd and
c        been corrected for BRUKER=T & SEQACQ=T.
c             Only need phase & shift corrections.  Use method of CHANGE-RAW.f
C        ----------------------------------------------------------------------
         cterm(1) = cexp(cmplx(0.,
     1              radian * (phitot(1) + sngl(parbes(lphast, 2)))))
         cfactor = cexp(cmplx(0., -2. * pi * hzpppm * deltat * shifd))
         do 710 junfil = 1, nunfil
            datat(junfil) = datat(junfil) * cterm(1)
            cterm(1) = cterm(1) * cfactor
 710     continue
         do 720 jdata = nunfil + 1, ndata
            datat(jdata) = (0., 0.)
 720     continue
         call cfft_r (datat, dataf, ndata, lwfft, wfftc)
         delta_ppm = float(nunfil) * ppminc
         do 730 jdata = 1, ndata
            dataf(jdata) = dataf(jdata) * cexp(cmplx(0., radian *
     1         delta_ppm * (phitot(2) + sngl(parbes(lphast + 1, 2)))))
            delta_ppm = delta_ppm - ppminc
 730     continue
         do 740 junfil = 1, nunfil
            cterm(1) = dataf(nunfil + junfil)
            dataf(nunfil + junfil) = dataf(junfil)
            dataf(junfil) = cterm(1)
 740     continue
         call cfftin (dataf, datat, ndata, lwfft, wfftc)
         write (lcoraw, nml=seqpar)
         bruker = .false.
         fmtdat = '(2e15.6)'
         id = 'FILCOR'
         seqacq = .false.
         tramp = 1.
         volume = 1.
         write (lcoraw, nml=nmid)
         write (lcoraw, fmt=fmtdat) (datat(junfil), junfil = 1, nunfil)
      END IF
      CALL EXITPS (.false.)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE EXITPS (lstop)
C
C  Output Postscript plot file.
C  Load error message table into ERRORS and output this.
C  LSTOP = T to exit with STOP.
C
      INCLUDE 'lcmodel.inc'
      logical lstop
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
      CHSUBP='EXITPS'
      IF (LPS .LE. 0) THEN
         CALL ERRTBL ()
      ELSE
         LLPS=LPS
         STRIP1=.NOT.CCNTRL
         CALL MAKEPS ()
         CLOSE (LPS)
      END IF
      IF (MAX0(LPRINT,LCOORD,LTABLE) .GT. 0) CALL ERRTBL ()
      IF (LPRINT .GT. 0) THEN
         WRITE (LPRINT, 5110) (ERRORS(J), J=1,LINERR)
 5110    FORMAT (////' Summary of diagnostics'/(A))
         WRITE (LPRINT, 5002)
 5002    FORMAT (/' ')
         CLOSE (LPRINT)
      END IF
      IF (LCOORD .GT. 0) THEN
         WRITE (LCOORD, 5120) LINERR, (ERRORS(J), J=1,LINERR)
 5120    FORMAT (1X,I3, ' lines in following diagnostic table:'/(A))
         WRITE (LCOORD, 5130) LINCHG(2), (CHANGE(J,2), J=1,LINCHG(2))
 5130    FORMAT (1X, I3, ' lines in following table of input changes:'/
     1           (A))
         CLOSE (LCOORD)
      END IF
      IF (LTABLE .GT. 0) THEN
         WRITE (LTABLE, 5122) LINERR, (ERRORS(J), J=1,LINERR)
 5122    FORMAT (/'$$DIAG',I3,
     1           ' lines in following diagnostic table:'/(A))
         WRITE (LTABLE, 5132) LINCHG(2), (CHANGE(J,2), J=1,LINCHG(2))
 5132    FORMAT (/'$$INPU', I3,
     1           ' lines in following table of input changes:'/ (A))
         CLOSE (LTABLE)
      END IF
      if (lcoraw .gt. 0) then
c        ---------------------------------------------------------------------
c        Intel Windows version outputs Namelists that cannot be read; so do
c           not bother to fix Cyg Namelist.
c        ---------------------------------------------------------------------
C        call fix_g77_namelist(lcoraw) !sun
         close (lcoraw)
      end if
      if (lstop) STOP
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ERRTBL ()
C
C  Load error message table into ERRORS.
C
      INCLUDE 'lcmodel.inc'
      CHARACTER*9 CHTYPE(5,2)
      LOGICAL ILLOG
      DATA CHTYPE/'info', 'warning', 'ERROR', 'FATAL', 'ILLOGICAL',
     1            'info''s', 'warnings', 'ERRORS', 'FATAL', 'ILLOGICAL'/
      CHSUBP='ERRTBL'
      DO 110 J=1,LINERR
         ILLOG=LEVERR(J).LT.1 .OR. LEVERR(J).GT.5
         IF (ILLOG) THEN
            NERROR(J)=1
            LEVERR(J)=5
            CHERR(J)=CHSUBP
            IERRNO(J)=1
         END IF
         WRITE (ERRORS(J),5110) NERROR(J),
     1     CHTYPE(LEVERR(J),MIN0(2,NERROR(J))), CHERR(J), IERRNO(J)
 5110    FORMAT (i4, 1X, A9, 1X, A6, I3)
         IF (ILLOG) GO TO 120
  110 CONTINUE
  120 LINERR=MIN0(J,LINERR)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MAKEPS ()
C
C  Compute and output Postscript plot file onto unit LPS.
C
C  ISTAGO tells what output is available for plotting as follows:
C  ISTAGO = 0 when nothing yet is available in DATAT for plotting;
C         = 1 when DATAT is loaded with time-domain data;
C         = 2 when CY is loaded with frequency-domain data (not necessarily
C             properly referenced or phased);
C         = 3 when YFITRE, BACKRE, and PPM are available
C         < 3 means that the run was aborted, and only the absolute value of
C             CY (if ISTAGO > 0) can be plotted.
C
      INCLUDE 'lcmodel.inc'
      CHARACTER SUBTIT*25
      CHSUBP='MAKEPS'
      call toupper_lower (.true., pgnorm)
      IF (PGNORM .EQ. 'US') THEN
         PAGEWD=8.5*2.54
         PAGEHT=11.*2.54
      ELSE IF (PGNORM.EQ.'EU' .OR. PGNORM.EQ.'A4') THEN
         PAGEWD=21.
         PAGEHT=29.7
      ELSE
         IF (AMIN1(PAGEWD,PAGEHT) .LE. 0.) CALL ERRMES (1, -4, CHSUBP)
      END IF
      IF (LANDSC) THEN
         PAGEX=PAGEHT
         PAGEY=PAGEWD
      ELSE
         PAGEX=PAGEWD
         PAGEY=PAGEHT
      END IF
C     -------------------------------------------------------------------------
C     Initial setup of PostScript file.
C     -------------------------------------------------------------------------
      IF (FILPS .NE. ' ') THEN
C         OPEN (LPS, FILE=FILPS, STATUS='UNKNOWN', !OSF
C     1         CARRIAGECONTROL='LIST', err=802)!OSF
          OPEN (LPS, FILE=FILPS, STATUS='UNKNOWN', err=802)!Cyg
C        OPEN (LPS, FILE=FILPS, STATUS='UNKNOWN', err=802)!sun
C         OPEN (LPS, FILE=FILPS, STATUS='UNKNOWN', err=802)!IRIX
         CALL PSETUP (.true., PAGEWD, PAGEHT, LANDSC)
      END IF
      IF (ISTAGO. EQ. 0) THEN
         CALL ONEPAG (0, CPY, CPY, NCHLIN(1), 0,
     1                PAGEX, PAGEY, SUBTIT)
      ELSE IF (ISTAGO .LE. 2) THEN
C        ----------------------------------------------------------------------
C        Put absolute value of (possibly unreferenced) spectrum into
C          YREAL(JY)
C        ----------------------------------------------------------------------
         IF (ISTAGO .EQ. 1) CALL FTDATA (0)
         DO 110 JY=1,NY
            YREAL(JY)=CABS(CY(JY))
  110    CONTINUE
C        ----------------------------------------------------------------------
C        Reverse order of arrays for easier plotting.
C        ----------------------------------------------------------------------
         CALL REVERS (PPM, NYUSE)
         CALL REVERS (YREAL, NYUSE)
         CALL ONEPAG (1, YREAL, CPY, NCHLIN(1), 0,
     1                PAGEX, PAGEY, SUBTIT)
      ELSE
         CALL REVERS (PPM, NYUSE)
         CALL REVERS (YREAL, NYUSE)
         CALL REVERS (BACKRE, NYUSE)
         CALL REVERS (YFITRE(1,0), NYUSE)
         CALL ONEPAG (4, YFITRE(1,0), YREAL, NCHLIN(1), 0,
     1                PAGEX, PAGEY, SUBTIT)
         DO 210 JMETAB=1,NMETAB
            IF (CONCEN(JMETAB) .LE. 0.) GO TO 210
            DO 220 JEACH=1,MIN0(MMETAB,NEACH)
               IF (NAMEAC(JEACH).EQ.NACOMB(JMETAB)(1:6) .OR.
     1             NEACH.GT.NMETAB) THEN
                  if (ldump(4) .and. subbas) then
                     sumy = 0.
                     do 230 jy = 1, nyuse
                        sumy = sumy + yfitre(jy, jmetab)
 230                 continue
 5230                format ('metab = ', a6, 3x, 'Conc =', 1pe10.3, 3x,
     1                       'Area/Conc =', 1pe10.3)
                     write (lprint, 5230) nacomb(JMETAB)(1:6),
     1                     concen(jmetab), sumy / CONCEN(JMETAB)
                  end if
                  WRITE (SUBTIT,5220) NACOMB(JMETAB)(1:6),
     1                                CONCEN(JMETAB)
 5220             FORMAT (A6, '   Conc. =', 1PE9.2)
                  CALL REVERS (YFITRE(1,JMETAB), NYUSE)
                  CALL ONEPAG (3, YFITRE(1,JMETAB), YREAL, 0, 1,
     1                         PAGEX, PAGEY, SUBTIT)
                  GO TO 210
               END IF
  220       CONTINUE
  210    CONTINUE
      END IF
      CALL ENDPS ()
      return
 802  call errmes (2, -4, chsubp)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE REVERS (X, N)
C
C  Reverse order of array X.
C
      REAL X(N)
      K=N
      DO 110 J=1,N/2
         TERM=X(K)
         X(K)=X(J)
         X(J)=TERM
         K=K-1
  110 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ONEPAG (NCURVE, YFIT, YDATA, LCHLIN, NSUBTI,
     1                   PAGEX, PAGEY, SUBTIT)
C
C  Make one PostScript page.
C
C  NCURVE = no. of curves to be plotted;
C           0 if only the table is to be output.
C           1 if only the table (with diagnostics and input) and the absolute
C             value of the data (in YFIT in this routine) are to be plotted;
C           3 if the contribution of one metabolite (in YFIT), the data (in
C             YDATA), the background (in BACKRE) are to be plotted.  No table
C             is plotted.
C           4 if the fit (in FIT), the data (in YDATA), the background (in
C             BACKRE), the residuals, and the table are to be plotted.
C
c  LCHLIN is only used to compute axis length for plot on page 1.
c         So, normally the actual argument in the calling program is NCHLIN(1).
c
C  All distances are in cm, except variables beginning with PT, which are in
C    points.
C                      YFITRE(*,JMETAB)] and the data are to be plotted.
C
      INCLUDE 'lcmodel.inc'
      CHARACTER CHREAL*8, CHSTR*132, SUBTIT*(*), XLABEL*20
      logical outside
      REAL YDATA(*), YFIT(*)
      DATA RBOXSP/.3/, XLABEL/'Chemical Shift (ppm)'/
      CHSUBP='ONEPAG'
      IF (NCURVE.LE.-1 .OR. NCURVE.EQ.2 .OR. NCURVE.GE.5) CALL ERRMES
     1      (1, -5, CHSUBP)
      IF (AMIN1(PTLABL,PTTITL,PTOUTP,PTVERS).LE.0. .OR.
     1    AMIN1(RHLABL,RHTITL,RHOUTP,RHVERS).LT.1.) CALL ERRMES
     1      (2, -4, CHSUBP)
      CALL PSETUP (.false., PAGEWD, PAGEHT, LANDSC)
C     -------------------------------------------------------------------------
C     SZ* = the character height in cm.
C     SP* = the total height of a line, including the space between the
C           current line and the line toward the center of the plot.
C         = (SZ*)*(RH*)
C     -------------------------------------------------------------------------
      PTTOCM=2.54/72.
      SZLABL=PTLABL*PTTOCM
      SZTITL=PTTITL*PTTOCM
      SZOUTP=PTOUTP*PTTOCM
      SZVERS=PTVERS*PTTOCM
      SPLABL=SZLABL*RHLABL
      SPTITL=SZTITL*RHTITL
      SPOUTP=SZOUTP*RHOUTP
      SPVERS=SZVERS*RHVERS
      TCKLNG=.6*SZLABL
      IF (MIN0(LCHLIN,NSUBTI).LT.0 .OR. NSUBTI.GT.1 .OR.
     1    LCHLIN.GT.NCHLIN(1)) CALL ERRMES (3, -5, CHSUBP)
      XORIG=XLEFT +SPLABL +TCKLNG
      XAXIS=PAGEX -XORIG -RWFONT*FLOAT(LCHLIN)*SZOUTP -XRIGHT
      YORIG=YBOTT +2.*SPLABL +TCKLNG
      YAXIS=PAGEY -YORIG -FLOAT(NSUBTI+nlines_title)*SPTITL -2.*SPVERS
     1      -YTOP
      IF (AMIN1(XORIG,XAXIS,YORIG,YAXIS) .LT. 0.) CALL ERRMES (4, -4,
     1                                                         CHSUBP)
C     -------------------------------------------------------------------------
C     Outer box for plots, titles, and version.
C     DSHPAT(J,1) = dash pattern for fine grid lines (J=1,2);
C              2  =                  data.
C       CALL DASH (0,DSHPAT) produces a solid curve.
C       CALL DASH (J,DSHPAT) produces DSHPAT(*,J) if J=1,2
C     WDLINE(1) = linewidth of fit to data (except for absolute value of the
C                                           data when NCURVE=1);
C            2  =              baseline;
C            3  =              data;
C            4  =              residuals;
C            5  =              axes and frame;
C            6  =              fine grid;
C     -------------------------------------------------------------------------
      CALL DASH (0, DSHPAT)
      CALL LINEWD (WDLINE(5))
      CALL RGB (RGBLIN(1,5))
      CALL BOX (XORIG, YORIG, XAXIS, YAXIS)
      call rgb (black)
      CALL FONT (PTTITL, 'Helvetica')
      XCENTR=.5*(PAGEX -XLEFT -XRIGHT) +XLEFT
      YCURR=PAGEY -YTOP -SZTITL +sptitl
      do 110 jline = 1, nlines_title
         CALL STRCHK (TITLE_line(jline), CHSTR)
         YCURR=ycurr - sptitl
         if (nlines_title .eq. 1) then
            CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
         else
            CALL STRING ('l', 0, XORIG+.3, YCURR, CHSTR)
         end if
 110  continue
      IF (NSUBTI .EQ. 1) THEN
         CALL STRCHK (SUBTIT, CHSTR)
         YCURR=YCURR-SPTITL
         CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
      END IF
      CALL FONT (PTVERS, 'Helvetica-Oblique')
      CALL STRCHK (OWNOUT, CHSTR)
      YCURR=YCURR-SPTITL+SZTITL-SZVERS
      CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
      CALL STRCHK (VERSIO, CHSTR)
      YCURR=YCURR-SPVERS
      CALL STRING ('l', 0, XORIG+.3, YCURR, CHSTR)
      CALL STRCHK (CHDATE, CHSTR)
      CALL STRING ('r', 0, PAGEX-XRIGHT, YCURR, CHSTR)
      IF (NCURVE .GT. 0) THEN
C        ----------------------------------------------------------------------
C        X-axis (ppm).
C        ----------------------------------------------------------------------
         CALL FONT (PTLABL, 'Helvetica')
         CALL STRING ('c', 0, XORIG+.5*XAXIS, YBOTT, XLABEL)
         IF (XSTEP .LE. 0.) THEN
            XSTEP=.2
            CALL ERRMES (5, 3, CHSUBP)
         END IF
C        ----------------------------------------------------------------------
C        PPM is now in increasing order (since call to REVERS in MAKEPS).
C        X-axis frame has ends that are integer multiples of XSTEP.
C        Allow the curve to extend out of the x-axis frame by up to PPMINC.
C        ----------------------------------------------------------------------
         CALL ENDRND (PPM(1), PPM(NYUSE), XSTEP, PPMINC, XBEG, XEND)
         IF (XEND .LE. XBEG) CALL ERRMES (6, -4, CHSUBP)
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         CALL AXIS (0, XORIG+XAXIS, YORIG, -XAXIS, XBEG, XEND, XBEG,
     1              XSTEP, 1., .TRUE.)
         IF (NSUBTK .LE. 0) THEN
            NSUBTK=1
            CALL ERRMES (7, 3, CHSUBP)
         END IF
         CALL TICK (0, XORIG+XAXIS, YORIG, -XAXIS, XBEG, XEND, XBEG,
     1              XSTEP/FLOAT(NSUBTK), -.3*SZLABL)
         CALL LINEWD (WDLINE(6))
         CALL RGB (RGBLIN(1,6))
         IF (WDLINE(6) .GT. 0.) then
            CALL LINEWD (WDLINE(6))
            CALL RGB (RGBLIN(1,6))
            CALL DASH (1, DSHPAT)
            CALL TICK (0, XORIG+XAXIS, YORIG, -XAXIS, XBEG, XEND, XBEG,
     1                 XSTEP/FLOAT(NSUBTK), YAXIS)
         end if
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         YMIN=0.
         YMAX=0.
         IF (NCURVE .EQ. 1) THEN
            DO 210 J=1,NYUSE
               YMIN=AMIN1(YMIN, YFIT(J))
               YMAX=AMAX1(YMAX, YFIT(J))
  210       CONTINUE
         ELSE IF (SUBBAS) THEN
            DO 215 J=1,NYUSE
               YMIN=AMIN1(YMIN, YFIT(J), YDATA(J))
               YMAX=AMAX1(YMAX, YFIT(J), YDATA(J))
  215       CONTINUE
         ELSE
            DO 220 J=1,NYUSE
               YMIN=AMIN1(YMIN, YFIT(J), YDATA(J), BACKRE(J))
               YMAX=AMAX1(YMAX, YFIT(J), YDATA(J), BACKRE(J))
  220       CONTINUE
         END IF
         XSPACE=SPLABL-SZLABL+TCKLNG
         IF (NCURVE .LT. 4) THEN
            YAXIS1=YAXIS
            YTOTAL=YMAX-YMIN
            IF (YTOTAL .LE. 0.) CALL ERRMES (8, -4, CHSUBP)
            YTOCM=YAXIS/YTOTAL
         ELSE
C           -------------------------------------------------------------------
C           Plot residuals above main plot.
C           CPY = residuals.
C           -------------------------------------------------------------------
            RESMAX=0.
            RESMIN=0.
            DO 230 J=1,NYUSE
               CPY(J)=YDATA(J)-YFIT(J)
               RESMIN=AMIN1(RESMIN,CPY(J))
               RESMAX=AMAX1(RESMAX,CPY(J))
  230       CONTINUE
            YTOTAL=RESMAX-RESMIN+YMAX-YMIN
            IF (YTOTAL .LE. 0.) CALL ERRMES (9, -4, CHSUBP)
            YAXIS1=YAXIS*(YMAX-YMIN)/YTOTAL
            YAXIS2=YAXIS*(RESMAX-RESMIN)/YTOTAL
            YTOCM=YAXIS/YTOTAL
            IF (RESMAX .LE. RESMIN) THEN
               CALL ERRMES (10, 3, CHSUBP)
            ELSE
               YZERO=YORIG+YAXIS-RESMAX*YTOCM
               CALL LINEWD (WDLINE(6))
               CALL RGB (RGBLIN(1,6))
               CALL LINE (XORIG, YZERO, XAXIS, 0.)
               IF (AMIN1(RESMAX,-RESMIN)*YTOCM .GT. 1.) CALL STRING
     1               ('c', 90, XORIG-XSPACE, YZERO, '0')
               if (resmax * ytocm .gt. .6) then
                  CALL STRING ('c', 90, XORIG-XSPACE, YORIG+YAXIS,
     1                         CHREAL(RESMAX, RESMAX, .falsE.))
               else
                  CALL STRING ('l', 90, XORIG-XSPACE, YORIG+YAXIS,
     1                         CHREAL(RESMAX, RESMAX, .truE.))
               end if
               CALL DASH (0, DSHPAT)
               CALL LINEWD (WDLINE(5))
               CALL RGB (RGBLIN(1,5))
               CALL LINE (XORIG-TCKLNG, YZERO, TCKLNG, 0.)
               CALL LINE (XORIG-TCKLNG, YORIG+YAXIS, TCKLNG, 0.)
               CALL LINEWD (WDLINE(4))
               CALL RGB (RGBLIN(1,4))
               CALL PLOT_gap (NYUSE, PPM, CPY, XBEG, XEND, RESMIN,
     1            RESMAX, XORIG+XAXIS, YORIG+YAXIS1, -XAXIS, YAXIS2)
            END IF
         END IF
         if (solgrd   .and.   wdline(6) .gt. 0.) then
            CALL LINEWD (WDLINE(6))
            CALL RGB (RGBLIN(1,6))
            CALL DASH (0, DSHPAT)
            CALL TICK (0, XORIG+XAXIS, YORIG, -XAXIS, XBEG, XEND,
     1                 XBEG, XSTEP, YAXIS1)
         end if
C        ----------------------------------------------------------------------
C        Label main Y-axis, and plot first curve (in YFIT).
C        ----------------------------------------------------------------------
         CALL DASH (1, DSHPAT)
         CALL LINEWD (WDLINE(6))
         CALL RGB (RGBLIN(1,6))
         YZERO=YORIG+YAXIS1-YMAX*YTOCM
         CALL LINE (XORIG, YZERO, XAXIS, 0.)
         CALL DASH (0, DSHPAT)
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         IF (YMAX*YTOCM .GT. 1.) CALL STRING ('c', 90, XORIG-XSPACE,
     1                                        YZERO, '0')
         if (-resmin * ytocm .gt. .6) then
            CALL STRING ('c', 90, XORIG-XSPACE, YORIG+YAXIS1,
     1                   CHREAL(YMAX, YMAX, .falsE.))
         else
            CALL STRING ('r', 90, XORIG-XSPACE, YORIG+YAXIS1,
     1                   CHREAL(YMAX, YMAX, .falsE.))
         end if
         CALL LINE (XORIG-TCKLNG, YZERO, TCKLNG, 0.)
         CALL LINE (XORIG-TCKLNG, YORIG+YAXIS1, TCKLNG+XAXIS, 0.)
         CALL LINEWD (WDLINE(1))
         CALL RGB (RGBLIN(1,1))
         IF (YMAX .LE. YMIN) CALL ERRMES (11, -4, CHSUBP)
         CALL PLOT_gap (NYUSE, PPM, YFIT, XBEG, XEND, YMIN, YMAX,
     1                  XORIG+XAXIS, YORIG, -XAXIS, YAXIS1)
         IF (NCURVE .GE. 3) THEN
C           -------------------------------------------------------------------
C           Plot curves 2 and 3 (in BACKRE and YDATA).
C           -------------------------------------------------------------------
            IF (WDLINE(2).GT.0. .AND. .NOT.SUBBAS   .and.
     1          nbackg .gt. 0) THEN
               CALL LINEWD (WDLINE(2))
               CALL RGB (RGBLIN(1,2))
               CALL PLOT_gap (NYUSE, PPM, BACKRE, XBEG, XEND, YMIN,
     1                        YMAX, XORIG+XAXIS, YORIG, -XAXIS, YAXIS1)
            END IF
            CALL DASH (2, DSHPAT)
            CALL LINEWD (WDLINE(3))
            CALL RGB (RGBLIN(1,3))
            CALL PLOT_gap (NYUSE, PPM, YDATA, XBEG, XEND, YMIN, YMAX,
     1                     XORIG+XAXIS, YORIG, -XAXIS, YAXIS1)
         END IF
      END IF
      IF (LCHLIN .le. 0) go to 800
C     ----------------------------------------------------------------------
C     Output data for licensing
c            table of concentrations (in TABLE),
C            misc. output (in ETCOUT),
C            diagnostics (in ERRORS), and
C            input changes to Control Variables (in CHANGES).
C     Goes right to the bottom of the page, if necessary, ignoring YBOTT.
C     ----------------------------------------------------------------------
      CALL DASH (0, DSHPAT)
      CALL LINEWD (WDLINE(5))
      CALL RGB (RGBLIN(1,5))
      XBOXLO=XORIG+XAXIS
      XBOXHI=PAGEX-XRIGHT
      XCENTR=.5*(XBOXHI+XBOXLO)
      YBOXHI=PAGEY -YTOP - FLOAT(NSUBTI+nlines_title)*SPTITL -2.*SPVERS
      YCURR=YBOXHI-SPOUTP
      NLINE=0
      YBOXSP=RBOXSP*SPOUTP
      LOUTMX=INT(YBOXHI/SPOUTP)
      NETCOU=MAX0(0,MIN0(LINETC,IETCOU))
C     ------------------------------------------------------------------
C     Output data that users must return for licensing.
c     NOLINE = T means that there is no license.
c     LETT = 0 if this is not Linux, which means that DATA FOR LICENSING table
c              must be output.
C     ------------------------------------------------------------------
      IF (NOLINE) THEN
         nline=1
         CALL FONT (PTOUTP, 'Courier-Bold')
         if (lett .eq. 0) then
            CALL STRING ('c', 0, XCENTR, YCURR, 'DATA FOR LICENSING')
         else
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, ycurr,
     1                   'Send the following file:')
         end if
         CALL FONT (PTOUTP, 'Courier-Bold')
         call rgb (rgberr)
         if (lett .eq. 0) then
            DO 305 JDEV=1,4
               IF (LDEV1(JDEV) .LT. 0.) GO TO 306
               nline=nline+1
               WRITE (CHSTR, '(I10)') LDEV1(JDEV)
               YCURR=YCURR-SPOUTP
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                      CHSTR)
 305        CONTINUE
         else
            nline = nline + 1
            YCURR=YCURR-SPOUTP
            CALL STRING ('c', 0, XCENTR, YCURR,
     1                   '~/.lcmodel/data-for-licensing')
         end if
 306     CALL FONT (PTOUTP, 'Courier')
         call rgb (black)
         YBOXLO=YCURR-YBOXSP
         IF (YBOXLO .LE. 0.) GO TO 700
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         CALL BOX (XBOXLO, YBOXLO, XBOXHI-XBOXLO, YBOXHI-YBOXLO)
         CALL RGB (BLACK)
         YCURR=YCURR-SPOUTP+YBOXSP
         nline=nline+1
      END IF
c
      IF (LINTBL+nline .GE. LOUTMX   .and.   ipage2 .le. 0) THEN
         LINTBL=LOUTMX-1-nline
         CALL ERRMES (12, 3, CHSUBP)
      END IF
      IF (LINTBL .GE. 1) THEN
C     -------------------------------------------------------------------
C     Concentration table.
C     -------------------------------------------------------------------
         call rgb (black)
         CALL FONT (PTOUTP, 'Courier-Bold')
         if (noline) YCURR=YCURR-SPOUTP
         CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                TABLE(1,1))
         DO 310 J=2,LINTBL
            IF (IPCERR(J) .LT. ISDBOL) THEN
               if (ratio_outlier(j)) then
                  call rgb(rgbrat)
                  CALL FONT (PTOUTP, 'Courier-BoldOblique')
               else
                  CALL RGB (RGBBOL)
                  CALL FONT (PTOUTP, 'Courier-Bold')
               end if
            ELSE
               if (ratio_outlier(j)) then
                  call rgb(rgbrat)
                  CALL FONT (PTOUTP, 'Courier-BoldOblique')
               else
                  CALL RGB (BLACK)
                  CALL FONT (PTOUTP, 'Courier')
               end if
            END IF
            YCURR=YCURR-SPOUTP
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   TABLE(J,1))
            if (j - 1 .eq. ntable_top   .and.   j .lt. lintbl) then
               ycurr = ycurr - yboxsp
               CALL RGB (BLACK)
               call line (xboxlo, ycurr, xboxhi - xboxlo, 0.)
            end if
 310     CONTINUE
         YBOXLO=YCURR-YBOXSP
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         CALL BOX (XBOXLO, YBOXLO, XBOXHI-XBOXLO, YBOXHI-YBOXLO)
         CALL RGB (BLACK)
         YBOXHI=YBOXLO
         YCURR=YCURR-SPOUTP+YBOXSP
         NLINE=nline+LINTBL+1
      END IF
      CALL FONT (PTOUTP, 'Courier-Bold')
      NLINE=NLINE+1
      IF (NLINE .GE. LOUTMX) GO TO 700
      YCURR=YCURR-SPOUTP
      IF (LINERR.GE.1 .OR. lline .or. .NOT.DOFULL .or.  wsdone) THEN
C        -------------------------------------------------------------------
C        Diagnostics.
C        -------------------------------------------------------------------
         CALL STRING ('c', 0, XCENTR, YCURR, 'DIAGNOSTICS')
         CALL ERRTBL ()
         DO 320 J=1,LINERR
            NLINE=NLINE+1
            IF (NLINE .GE. LOUTMX) GO TO 700
            IF (INDEX(ERRORS(J),'info') .gt. 0) THEN
               CALL FONT (PTOUTP, 'Courier')
               CALL RGB (black)
            else if (INDEX(ERRORS(J),'warning') .gt. 0) then
               CALL FONT (PTOUTP, 'Courier-Bold')
               CALL RGB (rgbbol)
            ELSE
               CALL FONT (PTOUTP, 'Courier-Bold')
               CALL RGB (rgberr)
            END IF
            YCURR=YCURR-SPOUTP
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   ERRORS(J))
 320     CONTINUE
         call rgb(black)
         IF (LLINE) THEN
            NLINE=NLINE+1
            IF (NLINE .GE. LOUTMX) GO TO 700
            YCURR=YCURR-SPOUTP
            CALL FONT (PTOUTP, 'Courier-Bold')
            call rgb (rgberr)
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1           '  LICENSE REQUIRED FOR THIS.')
            CALL FONT (PTOUTP, 'Courier')
            call rgb (black)
         else
            IF (.NOT.DOFULL) THEN
               NLINE=NLINE+1
               IF (NLINE .GE. LOUTMX) GO TO 700
               YCURR=YCURR-SPOUTP
               CALL FONT (PTOUTP, 'Courier-Bold')
               call rgb (rgberr)
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1              '  CRUDE MODEL')
               CALL FONT (PTOUTP, 'Courier')
               call rgb (black)
            END IF
            if (wsdone) then
               NLINE=NLINE+1
               IF (NLINE .GE. LOUTMX) GO TO 700
               YCURR=YCURR-SPOUTP
               CALL FONT (PTOUTP, 'Courier')
               call rgb (black)
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1              '  Doing Water-Scaling')
            END IF
         END IF
      ELSE
         CALL STRING ('c', 0, XCENTR, YCURR, 'NO DIAGNOSTICS')
      END IF
      YBOXLO=YCURR-YBOXSP
      IF (YBOXLO .LE. 0.) GO TO 700
      CALL LINEWD (WDLINE(5))
      CALL RGB (RGBLIN(1,5))
      CALL BOX (XBOXLO, YBOXLO, XBOXHI-XBOXLO, YBOXHI-YBOXLO)
      CALL RGB (BLACK)
      YBOXHI=YBOXLO
      YCURR=YCURR-SPOUTP+YBOXSP
      NLINE=NLINE+1
      IF (NETCOU .GE. 1) THEN
C     -------------------------------------------------------------------
C     Misc. output.
C     -------------------------------------------------------------------
         CALL FONT (PTOUTP, 'Courier-Bold')
         NLINE=NLINE+1
         IF (NLINE .GE. LOUTMX) GO TO 700
         YCURR=YCURR-SPOUTP
         CALL STRING ('c', 0, XCENTR, YCURR, 'MISCELLANEOUS OUTPUT')
         CALL FONT (PTOUTP, 'Courier')
         DO 330 J=1,NETCOU
            NLINE=NLINE+1
            IF (NLINE .GE. LOUTMX) GO TO 700
            YCURR=YCURR-SPOUTP
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   ETCOUT(J))
 330     CONTINUE
         YBOXLO=YCURR-YBOXSP
         IF (YBOXLO .LE. 0.) GO TO 700
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         CALL BOX (XBOXLO, YBOXLO, XBOXHI-XBOXLO, YBOXHI-YBOXLO)
         YBOXHI=YBOXLO
         CALL RGB (BLACK)
         YCURR=YCURR-SPOUTP+YBOXSP
         NLINE=NLINE+1
      END IF
      IF (LINCHG(1) .GE. 1) THEN
C     -------------------------------------------------------------------
C     Input changes.
C     -------------------------------------------------------------------
         CALL FONT (PTOUTP, 'Courier-Bold')
         NLINE=NLINE+1
         IF (NLINE .GE. LOUTMX) GO TO 700
         YCURR=YCURR-SPOUTP
         CALL STRING ('c', 0, XCENTR, YCURR, 'INPUT CHANGES')
         CALL FONT (PTOUTP, 'Courier')
         DO 340 J=1,LINCHG(1)
            NLINE=NLINE+1
            IF (NLINE .GE. LOUTMX) GO TO 700
            YCURR=YCURR-SPOUTP
            CALL STRCHK (CHANGE(J,1), CHSTR)
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR, CHSTR)
 340     CONTINUE
         YBOXLO=YCURR-YBOXSP
         IF (YBOXLO .LE. 0.) GO TO 700
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
         CALL BOX (XBOXLO, YBOXLO, XBOXHI-XBOXLO, YBOXHI-YBOXLO)
         CALL RGB (BLACK)
         YBOXHI=YBOXLO
         YCURR=YCURR-SPOUTP+YBOXSP
      END IF
  700 CALL SHOWPG()
      IF (ipage2 .LE. 0) return
      IF (IPAGE2.EQ.1 .AND. NLINE.LT.LOUTMX) return
C     -----------------------------------------------------------------------
C     Plot tables on page 2.
C     -----------------------------------------------------------------------
      CALL PSETUP (.false., PAGEWD, PAGEHT, LANDSC)
      XORIG=XLEFT
      XAXIS=PAGEX -XORIG -XRIGHT
C     -----------------------------------------------------------------------
C     Title, owner, version & date (as page 1)
C     -----------------------------------------------------------------------
      CALL DASH (0, DSHPAT)
      CALL FONT (PTTITL, 'Helvetica')
      XCENTR=.5*(PAGEX -XLEFT -XRIGHT) +XLEFT
      YCURR=PAGEY -YTOP -SZTITL +sptitl
      do 410 jline = 1, nlines_title
         CALL STRCHK (TITLE_line(jline), CHSTR)
         YCURR=ycurr - sptitl
         if (nlines_title .eq. 1) then
            CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
         else
            CALL STRING ('l', 0, XORIG+.3, YCURR, CHSTR)
         end if
 410  continue
      IF (NSUBTI .EQ. 1) THEN
         CALL STRCHK (SUBTIT, CHSTR)
         YCURR=YCURR-SPTITL
         CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
      END IF
      CALL FONT (PTVERS, 'Helvetica-Oblique')
      CALL STRCHK (OWNOUT, CHSTR)
      YCURR=YCURR-SPTITL+SZTITL-SZVERS
      CALL STRING ('c', 0, XCENTR, YCURR, CHSTR)
      CALL STRCHK (VERSIO, CHSTR)
      YCURR=YCURR-SPVERS
      CALL STRING ('l', 0, XORIG+.3, YCURR, CHSTR)
      CALL STRCHK (CHDATE, CHSTR)
      CALL STRING ('r', 0, PAGEX-XRIGHT, YCURR, CHSTR)
C     -----------------------------------------------------------------------
C     Set column dimensions.
C     -----------------------------------------------------------------------
      column_width=rwfont*float(nchlin(2))*szoutp
      xboxlo = xorig
      xboxlo_max = pagex - xright - column_width
      ytop_column = PAGEY -YTOP - FLOAT(NSUBTI+nlines_title)*SPTITL
     1              -2.*SPVERS
      ycurr = ytop_column
C     -----------------------------------------------------------------------
C     For each table:
C       Draw top line.
C       If at bottom, draw side lines (and update positions)
C       If at end of table, draw side lines & bottom line.
C     -----------------------------------------------------------------------
      CALL LINEWD (WDLINE(5))
      CALL RGB (RGBLIN(1,5))
      call line (xboxlo, ycurr, column_width, 0.)
C     ------------------------------------------------------------------
C     Output data that users must return for licensing.
c     NOLINE = T means that there is no license.
c     LETT = 0 if this is not Linux, which means that DATA FOR LICENSING table
c              must be output.
C     ------------------------------------------------------------------
      IF (NOLINE) THEN
         ycurr = ycurr - spoutp
         CALL FONT (PTOUTP, 'Courier-Bold')
         xcentr = xboxlo + .5*column_width
         if (lett .eq. 0) then
            CALL STRING ('c', 0, XCENTR, YCURR, 'DATA FOR LICENSING')
         else
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, ycurr,
     1                   'Send the following file:')
         end if
         call rgb (rgberr)
         if (lett .eq. 0) then
            DO 705 JDEV=1,4
               IF (LDEV1(JDEV) .LT. 0.) GO TO 706
               WRITE (CHSTR, '(I10)') LDEV1(JDEV)
               call check_bottom (ycurr, spoutp, xboxlo, outside,
     1              ytop_column, column_width, xboxlo_max)
               if (outside) go to 800
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                      CHSTR)
 705        CONTINUE
         else
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1                         ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            CALL STRING ('c', 0, XCENTR, YCURR,
     1                   '~/.lcmodel/data-for-licensing')
         end if
 706     CALL FONT (PTOUTP, 'Courier')
         call rgb (black)
         call end_table (ycurr-yboxsp, xboxlo, ytop_column,
     1        column_width)
      END IF
      IF (LINTBL .GE. 1) THEN
C     -------------------------------------------------------------------
C     Concentration table.
C     -------------------------------------------------------------------
         if (noline) then
            call check_bottom (ycurr, spoutp-yboxsp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
         end if
         CALL FONT (PTOUTP, 'Courier-Bold')
         call check_bottom (ycurr, spoutp, xboxlo, outside,
     1        ytop_column, column_width, xboxlo_max)
         if (outside) go to 800
         xcentr = xboxlo + .5*column_width
         call rgb (black)
         CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                TABLE(1,2))
         DO 710 J=2,LINTBL
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            IF (IPCERR(J) .LT. ISDBOL) THEN
               if (ratio_outlier(j)) then
                  call rgb(rgbrat)
                  CALL FONT (PTOUTP, 'Courier-BoldOblique')
               else
                  CALL RGB (RGBBOL)
                  CALL FONT (PTOUTP, 'Courier-Bold')
               end if
            ELSE
               if (ratio_outlier(j)) then
                  call rgb(rgbrat)
                  CALL FONT (PTOUTP, 'Courier-Oblique')
               else
                  CALL RGB (BLACK)
                  CALL FONT (PTOUTP, 'Courier')
               end if
            END IF
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   TABLE(J,2))
            if (j - 1 .eq. ntable_top   .and.   j .lt. lintbl   .and.
     1          ycurr .lt. ytop_column - 1.01 * spoutp) then
               CALL RGB (BLACK)
               ycurr = ycurr - yboxsp
               call line (xboxlo, ycurr, column_width, 0.)
            end if
 710     CONTINUE
         call end_table (ycurr-yboxsp, xboxlo, ytop_column,
     1        column_width)
      end if
C     -------------------------------------------------------------------
C     Diagnostics.
C     -------------------------------------------------------------------
      call check_bottom (ycurr, spoutp-yboxsp, xboxlo, outside,
     1     ytop_column, column_width, xboxlo_max)
      if (outside) go to 800
      CALL RGB (BLACK)
      CALL FONT (PTOUTP, 'Courier-Bold')
      call check_bottom (ycurr, spoutp, xboxlo, outside,
     1     ytop_column, column_width, xboxlo_max)
      if (outside) go to 800
      xcentr = xboxlo + .5*column_width
      IF (LINERR.GE.1 .OR. lline .or..NOT.DOFULL .or. wsdone) THEN
         CALL STRING ('c', 0, XCENTR, YCURR, 'DIAGNOSTICS')
         CALL ERRTBL ()
         DO 720 J=1,LINERR
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            IF (INDEX(ERRORS(J),'info') .gt. 0) THEN
               CALL FONT (PTOUTP, 'Courier')
               CALL RGB (black)
            else if (INDEX(ERRORS(J),'warning') .gt. 0) then
               CALL FONT (PTOUTP, 'Courier-Bold')
               CALL RGB (rgbbol)
            ELSE
               CALL FONT (PTOUTP, 'Courier-Bold')
               CALL RGB (rgberr)
            END IF
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   ERRORS(J))
 720     CONTINUE
         call rgb(black)
         IF (LLINE) THEN
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            CALL FONT (PTOUTP, 'Courier-Bold')
            call rgb (rgberr)
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1           '  LICENSE REQUIRED FOR THIS.')
            CALL FONT (PTOUTP, 'Courier')
            call rgb (black)
         else
            IF (.NOT.DOFULL) THEN
               call check_bottom (ycurr, spoutp, xboxlo, outside,
     1              ytop_column, column_width, xboxlo_max)
               if (outside) go to 800
               CALL FONT (PTOUTP, 'Courier-Bold')
               call rgb (rgberr)
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1              '  CRUDE MODEL')
               CALL FONT (PTOUTP, 'Courier')
               call rgb (black)
            END IF
            if (wsdone) then
               call check_bottom (ycurr, spoutp, xboxlo, outside,
     1              ytop_column, column_width, xboxlo_max)
               if (outside) go to 800
               CALL FONT (PTOUTP, 'Courier')
               call rgb (black)
               CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1              '  Doing Water-Scaling')
            END IF
         END IF
      ELSE
         CALL STRING ('c', 0, XCENTR, YCURR, 'NO DIAGNOSTICS')
      END IF
      call end_table (ycurr-yboxsp, xboxlo, ytop_column, column_width)
      IF (NETCOU .GE. 1) THEN
C        -------------------------------------------------------------------
C        Misc. output.
C        -------------------------------------------------------------------
         call check_bottom (ycurr, spoutp-yboxsp, xboxlo, outside,
     1        ytop_column, column_width, xboxlo_max)
         if (outside) go to 800
         CALL FONT (PTOUTP, 'Courier-Bold')
         call check_bottom (ycurr, spoutp, xboxlo, outside,
     1        ytop_column, column_width, xboxlo_max)
         if (outside) go to 800
         xcentr = xboxlo + .5*column_width
         CALL STRING ('c', 0, XCENTR, YCURR, 'MISCELLANEOUS OUTPUT')
         CALL FONT (PTOUTP, 'Courier')
         DO 730 J=1,NETCOU
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR,
     1                   ETCOUT(J))
 730     CONTINUE
         call end_table (ycurr-yboxsp, xboxlo, ytop_column,
     1        column_width)
      END IF
      IF (LINCHG(2) .GE. 1) THEN
C     -------------------------------------------------------------------
C     Input changes.
C     -------------------------------------------------------------------
         call check_bottom (ycurr, spoutp-yboxsp, xboxlo, outside,
     1        ytop_column, column_width, xboxlo_max)
         if (outside) go to 800
         CALL FONT (PTOUTP, 'Courier-Bold')
         call check_bottom (ycurr, spoutp, xboxlo, outside,
     1        ytop_column, column_width, xboxlo_max)
         if (outside) go to 800
         xcentr = xboxlo + .5*column_width
         CALL STRING ('c', 0, XCENTR, YCURR, 'INPUT CHANGES')
         CALL FONT (PTOUTP, 'Courier')
         DO 740 J=1,LINCHG(2)
            call check_bottom (ycurr, spoutp, xboxlo, outside,
     1           ytop_column, column_width, xboxlo_max)
            if (outside) go to 800
            CALL STRCHK (CHANGE(J,2), CHSTR)
            CALL STRING ('l', 0, XBOXLO+.5*RWFONT*SZOUTP, YCURR, CHSTR)
 740     CONTINUE
         call end_table (ycurr-yboxsp, xboxlo, ytop_column,
     1        column_width)
      end if
 800  call showpg ()
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine end_table (ycurr, xboxlo, ytop_column, column_width)
c
c  Draw side lines and bottom line at end of table.
c
      include 'lcmodel.inc'
      CALL LINEWD (WDLINE(5))
      CALL RGB (RGBLIN(1,5))
      call line (xboxlo, ybott, 0., ytop_column - ybott)
      call line (xboxlo + column_width, ycurr, 0.,
     1        ytop_column - ycurr)
      call line (xboxlo, ycurr, column_width, 0.)
      call rgb (black)
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine check_bottom (ycurr, decrement, xboxlo, outside,
     1     ytop_column, column_width, xboxlo_max)
c
c  If YCURR is below YBOTT, then draw side lines and shift to the top of next
c    column.
c    OUTSIDE=T if this new column is out of bounds (i.e., if the updated
c              xboxlo > xboxlo_max
c
      include 'lcmodel.inc'
      logical outside
      outside = .false.
      ycurr = ycurr - decrement
      if (ycurr .lt. ybott) then
         CALL LINEWD (WDLINE(5))
         CALL RGB (RGBLIN(1,5))
c        ---------------------------------------------------------------------
c        Draw sidelines for full column.
c        ---------------------------------------------------------------------
         call line (xboxlo, ybott, 0., ytop_column - ybott)
         xboxlo = xboxlo + column_width
         call line (xboxlo, ybott, 0., ytop_column - ybott)
         ycurr = ytop_column - decrement
         outside = xboxlo .gt. xboxlo_max
         if (.not.outside) then
c           ------------------------------------------------------------------
c           Draw top line for new column.
c           ------------------------------------------------------------------
            call line (xboxlo, ytop_column, column_width, 0.)
         end if
         call rgb (black)
      end if
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE HEX(VAL,INUM,FLUSH)
      SAVE CHBUF
      CHARACTER CHBUF(19)*4, CHDIGT(0:15)*1, CHSUBP*6
      INTEGER INUM
      REAL VAL, POWER(4)
        LOGICAL FLUSH
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
      DATA CHDIGT/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
     1            'A', 'B', 'C', 'D', 'E', 'F'/,
     2     POWER/4096., 256., 16., 1./
      CHSUBP='HEX'
        IF (.NOT.FLUSH) THEN
         INUM=INUM+1
         REMAIN=INT(MAX(-1.0,MIN(1.0,VAL))*32767.)+32768
         DO 110 J=1,4
            IDIGIT=INT(REMAIN/POWER(J))
            IF (IDIGIT.LE.-1 .OR. IDIGIT.GE.16) CALL ERRMES (1, -5,
     1                                                       CHSUBP)
            CHBUF(INUM)(J:J)=CHDIGT(IDIGIT)
            REMAIN=REMAIN-FLOAT(IDIGIT)*POWER(J)
  110    CONTINUE
         ENDIF
        IF (INUM.EQ.19.OR.(FLUSH.AND.INUM.GT.0))THEN
10       FORMAT(1X, 19A4)
         WRITE(BUFOUT,10)(CHBUF(I),I=1,INUM)
         CALL STRPOU ()
         IF(.NOT.FLUSH)INUM=0
         ENDIF
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine remove_blank_start (str)
c
      implicit none
      character str*(*)
      integer ilen, j, length, nblanks
      external ilen
c
      length = ilen(str)
      do 110 nblanks = 1, length
         if (str(nblanks:nblanks) .ne. ' ') go to 200
 110  continue
 200  nblanks = nblanks - 1
      if (nblanks .gt. 0) then
         do 210 j = 1, length - nblanks
            str(j:j) = str(j + nblanks : j + nblanks)
 210     continue
      end if
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine toupper_lower (lupper_out, str)
c
c LUPPER_OUT = T to put STR in all upper case
c              F                   lower
c
      implicit none
      character ch(26,2)*(1), str*(*)
      integer iin, iout, jalpha, jstr
      logical lupper_out
      data ch/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
     1        'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
     2        'u', 'v', 'w', 'x', 'y', 'z',
     3        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
     1        'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
     2        'U', 'V', 'W', 'X', 'Y', 'Z'/
c
      if (lupper_out) then
         iin = 1
         iout = 2
      else
         iin = 2
         iout = 1
      end if
      do 110 jstr = 1, len(str)
         do 120 jalpha = 1, 26
            if (str(jstr:jstr) .eq. ch(jalpha, iin)) then
               str(jstr:jstr) = ch(jalpha, iout)
               go to 110
            end if
 120     continue
 110  continue
      end
c
c
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Converts g77 Namelist format to that readable by Sun, SGI, DEC.
c
      subroutine fix_g77_namelist (lunit)
      parameter (lscratch=31)
      character line*526
      external ilen
      open (lscratch, status='SCRATCH')
      rewind lunit
c
c     1st copy lunit to lscratch.
c
      do 110 jline = 1, 9999999
         read (lunit, 5110, end=200) line
 5110    format (a)
         llen = ilen(line)
         write (lscratch, 5110) line(1:llen)
 110  continue
c
c     Read lscratch, correct the format, and overwrite lunit.
c
 200  rewind lunit
      rewind lscratch
      do 210 jline = 1, 9999999
         read (lscratch, 5110, end=310) line
         llen = ilen(line)
         if (llen .gt. 520) then
            write (*, 5210) line
 5210       format ('***** The following line exceeds 520 charcters:'/a)
         endif
c        -----------------------------------------------------------------
c        Replace gfortran's quote delimiters with apostrophes.
c        -----------------------------------------------------------------
         do 215 jchar = 1, len(line)
            if (line(jchar:jchar) .eq. '"') line(jchar:jchar) = ''''
 215     continue
         if (index(line, '&') .eq. 1) then
            llen = llen + 1
            do 220 jchar = llen, 3, -1
               line(jchar:jchar) = line(jchar-1:jchar-1)
 220        continue
            line(1:2) = ' $'
         endif
c
         if (line(llen:llen) .eq. '/') then
            write (lunit, 5110) line(1:llen-1)
            line = ' $END'
            llen = 5
         endif
         write (lunit, 5110) line(1:llen)
 210  continue
 310  close (lscratch)
      return
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ENDRND (XMN, XMX, XSTEP, XINC, XMNRND, XMXRND)
C
C  Finds endpoints enclosing the interval [XMN,XMX] (to within a grace margin
C    of XINC) that are integer multiples of XSTEP.
C
C  XSTEP should be positive, and cannot be zero.
C
      XMNRND=XSTEP*FLOAT(NINT(XMN/XSTEP))
      XGRACE=XMNRND-ABS(XINC)
      IF (XGRACE .GT. XMN) XMNRND=XMNRND-ABS(XSTEP)
      XMXRND=XSTEP*FLOAT(NINT(XMX/XSTEP))
      XGRACE=XMXRND+ABS(XINC)
      IF (XGRACE .LT. XMX) XMXRND=XMXRND+ABS(XSTEP)
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CHARACTER*8 FUNCTION CHREAL (X, XSTEP, LJUST)
C
C  Outputs X into CHREAL in the shortest FORMAT, assuming that the E FORMAT
C    requires 7 characters for a positive number
C    and that at least 2 significant figures are required
C    and that X+XSTEP must be distinguishable from X.
C
      CHARACTER CH*8, FMTR*10
      LOGICAL LJUST
c     -------------------------------------------------------------------------
c     LJUST = T to left-justify numbers.
c               (Previously, it just returned a blank (commented out below) to
c                avoid rightmost axis number running against the tables in the
c                One_page output.)
c     -------------------------------------------------------------------------
c      if (ljust) then
c         chreal = ' '
c         return
c      end if
      AX=ABS(X)
      ASTEP=ABS(XSTEP)
      IF ((AX.LT.1.E-7 .AND. ASTEP.GT.1.E-5) .OR. AX.LE.0.) THEN
         IF (LJUST) THEN
            CHREAL='0.0'
         ELSE
            CHREAL='  0.0'
         END IF
         RETURN
      END IF
      CH=' '
      IF (AX.LT..0001 .OR. AX.GE.999999.) THEN
         IF (LJUST .AND. X.GE.0.) THEN
            FMTR='(1PE7.1E2)'
         ELSE
            FMTR='(1PE8.1E2)'
         END IF
         WRITE (CH,FMTR) X
         CHREAL=CH
         RETURN
      END IF
      IF (AX .LT. .001) THEN
         RFMT=8.5
      ELSE IF (AX .LT. .01) THEN
         RFMT=7.4
      ELSE IF (AX .LT. .1) THEN
         RFMT=6.3
      ELSE IF (AX .LT. 1.) THEN
         RFMT=5.2
      ELSE IF (AX .LT. 9.95) THEN
         RFMT=4.1
      ELSE
C        ----------------------------------------------------------------------
C        Here provision is made for Fn.0 so that it can be expanded below if
C          necessary to distinguish between neihboring numbers.  Otherwise,
C          I(n-1) will be used.
C        ----------------------------------------------------------------------
         RFMT=FLOAT(INT(ALOG10(FLOAT(NINT(AX))))+3)+.001
      END IF
      DO 210 J=1,4
         FACT=10.**NINT(10.*AMOD(RFMT,1.))
         IF (INT(FACT*(AX+.9999*ASTEP)).NE.INT(FACT*AX) .OR.
     1       RFMT.GE.8.) GO TO 220
         RFMT=RFMT+1.1
  210 CONTINUE
  220 IF (AMOD(RFMT,1.) .LT. .01) THEN
         L=INT(RFMT)
         IF (LJUST) THEN
            L=L-1
            IF (X .GE. 0.) L=L-1
         END IF
 5100    FORMAT ('(I', I1, ')')
         WRITE (FMTR,5100) L
         WRITE (CH,FMTR) NINT(X)
         CHREAL=CH
         RETURN
      END IF
      IF (RFMT.LT.8. .AND. .NOT.LJUST) RFMT=RFMT+1.
      IF (LJUST .AND. X.GE.0.) RFMT=RFMT-1.
 5110 FORMAT ('(F', F3.1, ')')
      WRITE (FMTR,5110) RFMT
      WRITE (CH,FMTR) X
      CHREAL=CH
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE STRCHK(ST,PS)
C
C  PS will be truncated to a max. of 123 characters to prevent BUFOUT from
C    being longer than 132 characters.
C
      EXTERNAL ILEN
      CHARACTER ST*(*),PS*(*)
      PS=' '
      IPS=0
      DO 110 I=1,ILEN(ST)
         IF (IPS .GE. 122) RETURN
         IF (ST(I:I).EQ.'('
     &     .OR.ST(I:I).EQ.')'
     &     .OR.ST(I:I).EQ.'%'
     &     .OR.ST(I:I).EQ.'\') THEN
            PS(IPS+1:IPS+2)='\'//ST(I:I)
            IPS=IPS+2
         ELSE
            IF (IPS .GE. 123) RETURN
               PS(IPS+1:IPS+1)=ST(I:I)
            IPS=IPS+1
         END IF
  110 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:TICK
C
C
C
        SUBROUTINE TICK(ANG,OX,OY,LENGTH,GMN,GMX,TCKBEG,TCKINC,GRID)
        INTEGER ANG
        REAL OX,OY,LENGTH,GMN,GMX,TCKBEG,TCKINC,GRID
        INTEGER I,J
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
        I=1
        VAL=TCKBEG
        DO 110 K=1,199
         IF (VAL-GMN .GE. -.001*ABS(TCKINC)) GO TO 115
         VAL=TCKBEG+FLOAT(I)*TCKINC
         I=I+1
  110   CONTINUE
  115   J=0
        DO 120 K=1,199
         IF (VAL-GMX .GT. .001*ABS(TCKINC)) GO TO 125
         VAL=TCKBEG+FLOAT(I)*TCKINC
         I=I+1
         J=J+1
  120   CONTINUE
10      FORMAT(1X, i3,3(1x,f8.4),i5,2(1x,f8.4),' ticks')
  125   WRITE(BUFOUT,10)J,LENGTH*TCKINC/ABS(GMX-GMN),
     &  LENGTH*(TCKBEG-GMN)/ABS(GMX-GMN),GRID,ANG,OX,OY
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:AXIS
C
C
C  LJUST1 = T to skip call to CHREAL and output of rightmost number to avoid
C             it running into the tables at the right.
C
        SUBROUTINE AXIS(ANG,OX,OY,LENGTH,GMN,GMX,
     &    TCKBEG,TCKINC,POS,LJUST1)
      EXTERNAL CHREAL
      INTEGER ANG
      CHARACTER CHREAL*8, chstr*8
      LOGICAL LJUST1
      REAL POS,OX,OY,LENGTH,GMN,GMX,TCKBEG,TCKINC
      REAL VAL
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10    FORMAT(1X, '(',a,')')
      I=1
      VAL=TCKBEG
      DO 110 K=1,199
       IF (VAL-GMN .GE. -.001*ABS(TCKINC)) GO TO 115
       VAL=TCKBEG+FLOAT(I)*TCKINC
       I=I+1
  110 CONTINUE
  115 J=0
      DO 120 K=1,199
       IF (VAL-GMX .GT. .001*ABS(TCKINC)) GO TO 125
       if (k .eq. 1   .and.   ljust1) then
          chstr=' '
       else
          chstr = CHREAL(VAL, TCKINC, .false.)
       end if
       WRITE(BUFOUT,10) chstr
       CALL STRPOU ()
       VAL=TCKBEG+FLOAT(I)*TCKINC
       I=I+1
       J=J+1
  120 CONTINUE
11    FORMAT(1X, i2,4(1x,f8.4),i5,2(1x,f8.4),' axis')
  125 WRITE(BUFOUT,11)J,POS,LENGTH*TCKINC/ABS(GMX-GMN),
     &LENGTH*(TCKBEG-GMN)/ABS(GMX-GMN),LENGTH,ANG,OX,OY
      CALL STRPOU ()
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:LINEWD
C
C LINEWD sets the line width using the current units (default cm).
C
        SUBROUTINE LINEWD(WIDTH)
        REAL WIDTH
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
 5110 FORMAT (1X, F8.4, A)
        WRITE(BUFOUT,5110)WIDTH,' setlinewidth'
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:RGB
C
C RGB sets the current color.
C
        SUBROUTINE RGB (RGBV)
        REAL RGBV(3)
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
 5110 FORMAT (1X, 3F8.4, A)
        WRITE(BUFOUT,5110) RGBV, ' setrgbcolor'
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:DASH
C
        SUBROUTINE DASH (IDSHPT, DSHPAT)
        REAL DSHPAT(2,*)
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, '[] 0 setdash')
12      FORMAT(1X, '[',f8.4,1x,f8.4,'] 0 setdash')
      IF (IDSHPT.LE.0 .OR. IDSHPT.GE.3) THEN
         WRITE(BUFOUT,10)
      ELSE IF (AMIN1(DSHPAT(1,IDSHPT),DSHPAT(2,IDSHPT)) .GT. 0.) THEN
         WRITE(BUFOUT,12) (DSHPAT(J,IDSHPT),J=1,2)
      ELSE
         WRITE(BUFOUT,10)
      END IF
      CALL STRPOU ()
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:PLOT
C
C
C
        SUBROUTINE PLOT(N,X,Y,XMN,XMX,YMN,YMX,OX,OY,WD,HT)
        INTEGER N,I,INUM
        REAL X(N),Y(N),XMN,XMX,YMN,YMX,OX,OY,WD,HT
        REAL DX,DY,XN,YN,XX,YX,XPAND,YPAND,XOFF,YOFF
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
        DX=ARBBOX(N,X,XN,XX)
        DY=ARBBOX(N,Y,YN,YX)
        ABSX=ABS(XMX-XMN)
        ABSY=ABS(YMX-YMN)
        IF (AMIN1(ABSX,ABSY,ABS(DX),ABS(DY)) .LE. 0.) RETURN
        XPAND=WD*DX/ABSX
        YPAND=HT*DY/ABSY
        XOFF=WD*(XN-MIN(XMN,XMX))/ABSX
        YOFF=HT*(YN-MIN(YMN,YMX))/ABSY
10      FORMAT(1X, i5,6(1x,f9.4),' plot')
        WRITE(BUFOUT,10)N,XOFF,YOFF,XPAND,YPAND,OX,OY
        CALL STRPOU ()
        INUM=0
        CALL HEX((X(1)-XN)/DX,INUM,.FALSE.)
        CALL HEX((Y(1)-YN)/DY,INUM,.FALSE.)
        DO 110 I=1,N
         CALL HEX((X(I)-XN)/DX,INUM,.FALSE.)
         CALL HEX((Y(I)-YN)/DY,INUM,.FALSE.)
  110   CONTINUE
        CALL HEX(0.,INUM,.TRUE.)
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:PLOT_gap
C
c Skips connecting lines thru gaps [detected by changes in X(J)-X(J-1)]
C
      SUBROUTINE PLOT_gap(N,X,Y,XMN,XMX,YMN,YMX,OX,OY,WD,HT)
      INTEGER N,I,INUM
      REAL X(N),Y(N),XMN,XMX,YMN,YMX,OX,OY,WD,HT
      REAL DX,DY,XN,YN,XX,YX,XPAND,YPAND,XOFF,YOFF
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
      DX=ARBBOX(N,X,XN,XX)
      DY=ARBBOX(N,Y,YN,YX)
      ABSX=ABS(XMX-XMN)
      ABSY=ABS(YMX-YMN)
      IF (AMIN1(ABSX,ABSY,ABS(DX),ABS(DY)) .LE. 0.) RETURN
      XPAND=WD*DX/ABSX
      YPAND=HT*DY/ABSY
      XOFF=WD*(XN-MIN(XMN,XMX))/ABSX
      YOFF=HT*(YN-MIN(YMN,YMX))/ABSY
      nend = 0
      do 110 jgap = 1, 10
         nstart = nend + 1
         if (nstart .ge. n) return
         xtol = .5 * abs(x(nstart + 1) - x(nstart))
         do 120 jx = nstart + 2, n
            if (abs(x(jx) - 2. * x(jx-1) + x(jx-2)) .gt. xtol)
     1         go to 130
 120     continue
 130     nend = jx - 1
 5110    FORMAT(1X, i5,6(1x,f9.4),' plot')
         WRITE(BUFOUT,5110) nend-nstart+1, XOFF,YOFF,XPAND,YPAND,OX,OY
         CALL STRPOU ()
         INUM=0
         CALL HEX((X(nstart)-XN)/DX,INUM,.FALSE.)
         CALL HEX((Y(nstart)-YN)/DY,INUM,.FALSE.)
         DO 150 I=nstart, nend
            CALL HEX((X(I)-XN)/DX,INUM,.FALSE.)
            CALL HEX((Y(I)-YN)/DY,INUM,.FALSE.)
 150     CONTINUE
         CALL HEX(0.,INUM,.TRUE.)
 110  continue
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        REAL FUNCTION ARBBOX(N,A,AMN,AMX)
        INTEGER N
        REAL A(N),AMN,AMX
        AMN=A(1)
        AMX=A(N)
        DO 110 I=1,N
         AMN=MIN(AMN,A(I))
         AMX=MAX(AMX,A(I))
  110   CONTINUE
        ARBBOX=ABS(AMX-AMN)
        IF(ARBBOX.EQ.0.)ARBBOX=1.0
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:BOX
C
C BOX draws box defined by the origin, the width and the height.
C
        SUBROUTINE BOX(OX,OY,WD,HT)
        REAL OX,OY,WD,HT
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, 4(f8.4,1x),'box')
        WRITE(BUFOUT,10)WD,HT,OX,OY
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:LINE
C
C LINE draws a line from the origin to the specified width and height.
C Intensive use of LINE could produce very large files.
C
        SUBROUTINE LINE(OX,OY,WD,HT)
        REAL OX,OY,WD,HT
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, 4(f8.4,1x),'li')
        WRITE(BUFOUT,10)WD,HT,OX,OY
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:PSETUP
C
      SUBROUTINE PSETUP (top_of_file, WD, HT, LANDSC)
      save kpage
      LOGICAL LANDSC, top_of_file
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
      data kpage/0/
C     -------------------------------------------------------------------------
C     Initialize output buffer BUFOUT.
C     -------------------------------------------------------------------------
      DO 110 J=1,20
         BUFOUT(J)=' '
  110 CONTINUE
      if (top_of_file) then
 10      FORMAT(1X, '%!PS-Adobe-2.0',/,1X,
     &          '%%Creator: LCModel',/,1X,
     &          '%%BoundingBox: 0 0 ',i4,' ',i4,/,1X,
     &          '%%EndComments')
         WRITE (BUFOUT,10) NINT(WD*72.0/2.54),NINT(HT*72.0/2.54)
         CALL STRPOU ()
         kpage = 0
         return
      end if
11      FORMAT(1X, '/MPGdict 100 dict def MPGdict begin',/,1X,
     &  '/bd{bind def}def',/,1X,
     &  '/set{exch def}bd',/,1X,
     &  '/cm{initmatrix 72 2.54 div dup scale 0 ',
     &  'setgray .015 setlinewidth',/,1X,
     &  '  1 selectplot /Helvetica 11 font}bd',/,1X,
     &  '/font{/pitch exch 72 div 2.54 mul def ',
     &  'findfont pitch scalefont setfont}bd',/,1X,
     &  '/s{show}bd',/,1X,
     &  '/lst{exec}bd')
12      FORMAT(1X, '/rst{stlen neg 0 rmoveto exec}bd',/,1X,
     &  '/cst{stlen -2 div 0 rmoveto exec}bd',/,1X,
     &  '/r{rotate}bd',/,1X,
     &  '/m{moveto}bd',/,1X,
     &  '/li{newpath moveto rlineto stroke}bd',/,1X,
     &  '/boxpath{newpath moveto dup 3 1 roll 0 exch',
     &  ' rlineto 0 rlineto 0 exch neg ',/,1X,
     &  ' rlineto closepath} bd',/,1X,
     &  '/box{boxpath stroke}def',/,1X,
     &  '/go{save exch 1 add 1 roll translate}bd')
13      FORMAT(1X, '/axis{',/,1X,
     &  '  8 go rotate',/,1X,
     &  '  /wd set',/,1X,
     &  '  /beg set',/,1X,
     &  '  /inc set',/,1X,
     &  '  /tick exch neg pitch mul def',/,1X,
     &  '  /n set',/,1X,
     &  '  newpath 0 0 moveto wd 0 lineto stroke',/,1X,
     &  '  newpath beg n 1 sub inc mul add 0 moveto',/,1X,
     &  '  n{0 tick 0.6 mul rlineto',/,1X,
     &  '    0 tick 0 gt{0.25}{1}ifelse tick mul rmoveto',/,1X,
     &  '    exch currentpoint pop inc sub exch',/,1X,
     &  '    dup stringwidth pop -0.5 mul 0 rmoveto show',/,1X,
     &  '    stroke 0 moveto}repeat',/,1X,
     &  '  restore}bd')
14      FORMAT(1X, '/readflt{currentfile 2 string readhexstring pop ',
     &  'dup 1 get exch 0 get',/,1X,
     &  '  256 mul add 32768 sub 32767 div mul}bd',/,1X,
     &  '/selectplot{/n set /plotproc',/,1X,
     &  '  n 1 eq{{lineto}}if',/,1X,
     &  '  def}bd',/,1X,
     &  '/plot{',/,1X,
     &  '  7 go /ypand set /xpand set /yoff set /xoff set /n set')
15      FORMAT(1X, '  /xy{xpand readflt xoff add ypand ',
     &  'readflt yoff add}bd',/,1X,
     &  '  xy moveto',/,1X,
     &  '  n{xy lineto currentpoint stroke moveto}repeat',/,1X,
     &  '  restore',/,1X,
     &  '}bd')
17      FORMAT(1X, '/ticks{',/,1X,
     &  '  7 go rotate',/,1X,
     &  '  /wd set',/,1X,
     &  '  /beg set',/,1X,
     &  '  /inc set',/,1X,
     &  '  /n set',/,1X,
     &  '  newpath beg 0 moveto',/,1X,
     &  '  n{0 wd rlineto currentpoint pop inc add 0',/,1X,
     &  '  stroke moveto}repeat',/,1X,
     &  ' restore}bd')
20      FORMAT(1X, '/stlen{',/,1X,
     &  '  dup currentpoint 3 -1 roll',/,1X,
     &  '  /len 0 def',/,1X,
     &  '  /dlen 0 store',/,1X,
     &  '  /s{',/,1X,
     &  '   stringwidth pop len add dlen add /len exch def',/,1X,
     &  '   /dlen 0 store}bd',/,1X,
     &  '  exec ',/,1X,
     &  '  moveto /s{show}bd len',/,1X,
     &  '  }bd')
 
        WRITE(BUFOUT,11)
        CALL STRPOU ()
        WRITE(BUFOUT,12)
        CALL STRPOU ()
        WRITE(BUFOUT,13)
        CALL STRPOU ()
        WRITE(BUFOUT,14)
        CALL STRPOU ()
        WRITE(BUFOUT,15)
        CALL STRPOU ()
        WRITE(BUFOUT,17)
        CALL STRPOU ()
        WRITE(BUFOUT,20)
        CALL STRPOU ()
 
        kpage = kpage + 1
30      FORMAT(1X, 'end',/,1X, '%%EndProlog',/, ' %%Page:', 2i4/
     1       1X, 'MPGdict begin cm')
        WRITE(BUFOUT,30) kpage, kpage
        CALL STRPOU ()
31      FORMAT(1X, '0 ',f7.4,' translate -90 rotate')
        IF (LANDSC) WRITE(BUFOUT,31)HT
        CALL STRPOU ()
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:FONT
C
C For convenience the FONT is called using the conventional
C pitch size of popular word-processing softwares; thus a
C 11 pitch size font will be 0.388 cm high (i.e. 11./72.*2.54).
C The name of the fonts can be any fonts available on your
C PostScript printer; typically;
C
C       CALL FONT(PITCH,'Helvetica-Bold')
C       CALL FONT(PITCH,'Times-Roman')
C       CALL FONT(PITCH,'Times-Italic')
C       CALL FONT(PITCH,'Times-Bold')
C       CALL FONT(PITCH,'Courier')
C       CALL FONT(PITCH,'Courier-Bold')
C       CALL FONT(PITCH,'Palatino-Bold')
C       CALL FONT(PITCH,'Palatino-Italic')
C       CALL FONT(PITCH,'Symbol')
C
        SUBROUTINE FONT(PITCH,POLICE)
      EXTERNAL ILEN
        CHARACTER POLICE*(*)
        REAL PITCH
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, '/',a,' ',f8.4,' font')
        WRITE(BUFOUT,10)POLICE(1:ILEN(POLICE)),PITCH
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:STRING
C
C STRING prints a character chain.  The string can be left
C justified ('l'), right justified ('r') or centered ('c').
C It can be written at an angle (integer in degrees) from
C the specified origin.
C Some editing commands have been incorporated for convenience, but a plain
C text string will always be easier to print.  A special care
C with '(', ')' and '%' signs should be made since they are
C interprated by PostScript.  They should be replaced by
C respectively '/('  '/)' and '/%'.
C
        SUBROUTINE STRING(FLUSH,ANG,OX,OY,ST)
        INTEGER ANG
        REAL OX,OY
        CHARACTER ST*(*),FLUSH*1
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, i4,' r')
11      FORMAT(1X, 2(f7.3,1x), 'm ',i4,' r')
12      FORMAT(1X, 2(f7.3,1x),'m')
        IF (ABS(ANG) .GT. 1.E-4) THEN
         WRITE(BUFOUT,11)OX,OY,ANG
         CALL STRPOU ()
         CALL SHOW(FLUSH,ST)
         WRITE(BUFOUT,10)-ANG
         CALL STRPOU ()
        ELSE
         WRITE(BUFOUT,12)OX,OY
         CALL STRPOU ()
         CALL SHOW(FLUSH,ST)
        ENDIF
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SHOW(FLUSH,ST)
      EXTERNAL ILEN
        CHARACTER ST*(*),FLUSH*1
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, '{(',a,')s}',a1, 'st')
        WRITE(BUFOUT,10)ST(1:ILEN(ST)),FLUSH
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C:SHOWPG
C
        SUBROUTINE SHOWPG()
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
10      FORMAT(1X, 'showpage')
        WRITE(BUFOUT,10)
        CALL STRPOU ()
        END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ENDPS ()
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
 5110 FORMAT (1X, '%%Trailer',/,1X, '%%EOF')
      WRITE (BUFOUT, 5110)
      CALL STRPOU ()
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        INTEGER FUNCTION ILEN(ST)
        CHARACTER ST*(*)
        ILEN=LEN(ST)
        DO 110 J=1,IABS(ILEN)
           IF (ST(ILEN:ILEN).NE.' ') RETURN
           ILEN=ILEN-1
  110   CONTINUE
        ILEN=1
        END
C
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE STRPOU ()
C
C  Strips the first (blank) character from each record in BUFOUT if STRIP1=T.
C  Writes BUFOUT onto LLPS.
C  Fills first characters of BUFOUT with blanks.
C
      EXTERNAL ILEN
      CHARACTER BUFOUT(20)*132
      LOGICAL STRIP1
      COMMON /BLPS/ BUFOUT, LLPS, STRIP1
      IF (STRIP1) THEN
         ICOL1=2
      ELSE
         ICOL1=1
      END IF
      DO 160 J=1,20
         L=ILEN(BUFOUT(J))
         IF (L.LE.ICOL1 .AND. BUFOUT(J)(ICOL1:ICOL1).EQ.' ') GO TO 200
         WRITE (LLPS,5160) BUFOUT(J)(ICOL1:L)
 5160    FORMAT (A)
  160 CONTINUE
  200 DO 210 J=1,20
         BUFOUT(J)=' '
  210 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE csft_r (DATAT, ft, ncap)
C
C  Input:
C    DATAT(J), J=1,Ncap
C  Output:
C    FT = FT of DATAT WITH rearrangement.  The FT is normalized; i.e.,
C         the sum is divided by sqrt(Ncap).  This is compatible with DCFFT
C         and CFFTIN and Siemens FFTs.
c
c Uses formula and notation of IMSL Math, p 716.
c
      complex cfact, cterm, datat(ncap), ft(ncap)
      TwoPI = 6.28318530717959d0
      do 150 m = 1, ncap
         cterm = (1., 0.)
         ft(m) = (0., 0.)
         cfact = cexp(cmplx(0., -twopi * float(m - 1) / float(ncap)))
         cterm = (1., 0.)
         do 180 n = 1, ncap
            ft(m) = ft(m) + datat(n) * cterm
            cterm = cterm * cfact
 180     continue
 150  continue
C     ------------------------------------------------------------------------
C     Scale & rearrange
C     ------------------------------------------------------------------------
      FACT=1./SQRT(FLOAT(Ncap))
      nunfil = ncap / 2
      DO 210 J=1,nunfil
         CTERM=FT(NUNFIL+J)
         FT(NUNFIL+J)=FT(J) * fact
         FT(J)=CTERM * fact
  210 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE csftin_r (ft, ftwork, ftinv, ncap)
C
C  Input:
C    Rearranged FT(J), J=1,N.  This will be unrearranged (in FTWORK) for the
c       inverse FFT.
c    FTWORK is an extra array used for the unrearranged FT.
C  Output:
C    FTINV = inverse FFT of FT without rearrangement.  FTINV is normalized;
C             i.e., the sum is divided by sqrt(N).  This is compatible with
C             DCFFT and CFFT and Siemens FFTs.
c Uses formula and notation of IMSL Math, p 716.
c
c     -------------------------------------------------------------------------
c     Unrearrange.
c     -------------------------------------------------------------------------
      complex cfact, cterm, ft(ncap), ftinv(ncap), ftwork(ncap)
      TwoPI = 6.28318530717959d0
      nunfil = ncap / 2
      do 110 J=1,nunfil
         ftwork(j)=FT(NUNFIL+J)
         FTwork(NUNFIL+J)=FT(J)
  110 CONTINUE
      do 150 n = 1, ncap
         cterm = (1., 0.)
         ftinv(n) = (0., 0.)
         cfact = cexp(cmplx(0., twopi * float(n - 1) / float(ncap)))
         cterm = (1., 0.)
         do 180 m = 1, ncap
            ftinv(n) = ftinv(n) + ftwork(m) * cterm
            cterm = cterm * cfact
 180     continue
 150  continue
      FACT=1./SQRT(FLOAT(Ncap))
      DO 210 J=1,Ncap
         ftinv(J)=ftinv(J)*FACT
  210 CONTINUE
      END
C
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SEQTOT (DATAT, DATAF, NUNFIL, LWFFT, WFFTC)
C
C  DATAT = NUNFIL sequentially acquired Bruker time-domain data on input.
C        = NUNFIL time-domain data after the correction (which includes
C          complex-conjugation of the time-domain data, since the Bruker
C          quadrature data are apparently at -90, rather than 90 deg)
C          on output.
C  DATAF : must be at least 2*NUNFIL long.
C  There cannot be any peaks in the last RENDS*NUNFIL/2 spectrum points at
C    each end.
C
      COMPLEX ENDAVG, DATAT(NUNFIL), DATAF(2*NUNFIL)
      REAL WFFTC(*)
      DATA NZERO/20/, RENDS/.05/, RZERO/.01/
C     -------------------------------------------------------------------------
C     Scale by sqrt(2) because the inverse FFT will be twice as long and
C       scaled by 1/sqrt(2*NUNFIL) rather than 1/sqrt(NUNFIL).
C     -------------------------------------------------------------------------
      DO 110 J=1,NUNFIL
         DATAF(J)=CMPLX(0., -AIMAG(DATAT(J))*1.414214)
  110 CONTINUE
      CALL CFFT (DATAF, DATAF, NUNFIL, LWFFT, WFFTC)
C     -------------------------------------------------------------------------
C     Fill FFT of imaginary part of DATAT with the mean of the high-frequency
C       ends.
C     -------------------------------------------------------------------------
      NHALF=NUNFIL/2
      NENDS=NINT(RENDS*FLOAT(NHALF))
      ENDAVG=(0.,0.)
      DO 115 J=NHALF-NENDS+1, NHALF+NENDS
         ENDAVG=ENDAVG+DATAF(J)
  115 CONTINUE
      ENDAVG=ENDAVG/FLOAT(2*NENDS)
      DO 120 J=NHALF+1,NUNFIL
         DATAF(J+NUNFIL)=DATAF(J)
         DATAF(J)=ENDAVG
         DATAF(J+NHALF)=ENDAVG
  120 CONTINUE
      NDATA=2*NUNFIL
      CALL CFFTIN (DATAF, DATAF, NDATA, LWFFT, WFFTC)
C     -------------------------------------------------------------------------
C     Insert interpolated imag. part into DATAT [except for first point, which
C       will be left (extrapolated) as the first original imag. point; this
C       will be corrected below anyway].
C     -------------------------------------------------------------------------
      DATAT(1)=CONJG(DATAT(1))
      JF=0
      DO 130 J=2,NUNFIL
         JF=JF+2
         DATAT(J)=CMPLX(REAL(DATAT(J)), AIMAG(DATAF(JF)))
  130 CONTINUE
C     -------------------------------------------------------------------------
C     Adjust DATAT(1) so that the high-frequency ends of the spectrum are about
C       zero.
C     -------------------------------------------------------------------------
      CALL CFFT (DATAT, DATAF, NUNFIL, LWFFT, WFFTC)
      ENDAVG=(0.,0.)
      DO 150 J=NHALF-NENDS+1, NHALF+NENDS
         ENDAVG=ENDAVG+DATAF(J)
  150 CONTINUE
      ENDAVG=ENDAVG/FLOAT(2*NENDS)
      DATAT(1)=DATAT(1)-ENDAVG
C     -------------------------------------------------------------------------
C     Zero up to RZERO*NUNFIL points if CABS of the preceding NZERO points are
C       all smaller.
C     -------------------------------------------------------------------------
      DO 210 JZERO=MAX0(1,NINT((1.-RZERO)*FLOAT(NUNFIL))), NUNFIL
         TERM=REAL(DATAT(JZERO))**2+AIMAG(DATAT(JZERO))**2
         IF (NZERO .GE. JZERO) RETURN
         DO 220 J=JZERO-1, JZERO-NZERO, -1
            IF (TERM .LE. REAL(DATAT(J))**2+AIMAG(DATAT(J))**2) GO TO
     1                                                          210
  220    CONTINUE
         DO 230 J=JZERO,NUNFIL
            DATAT(J)=(0.,0.)
  230    CONTINUE
         RETURN
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFTIN (FT, FTINV, N, LWFFT,
     1                   WFFTC)
C
C  Input:
C    FT(J), J=1,N
C    LWFFT = argument in last call to CALL FFTCI (LWFFT, WFFTC)
C  Output:
C    FTINV = inverse FFT of FT without rearrangement.  FFTINV is normalized;
C            i.e., the sum is divided by sqrt(N).  This is compatible with
C            DCFFT and CFFT and Siemens FFTs.
C
      COMPLEX FT(N), FTINV(N)
      REAL WFFTC(4*N+15)
      IF (N. NE. LWFFT) THEN
         CALL FFTCI (N, WFFTC)
         LWFFT=N
      END IF
      CALL F2TCB (N, FT, FTINV, WFFTC)
      FACT=1./SQRT(FLOAT(N))
      DO 210 J=1,N
         FTINV(J)=FTINV(J)*FACT
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFTIN_r (FT, ftwork, FTINV, N, LWFFT,
     1                     WFFTC)
C
C  Input:
C    Rearranged FT(J), J=1,N.  This will be unrearranged (in FTWORK) for the
c       inverse FFT.
c    FTWORK is an extra array used for the unrearranged FT.  So, FT is
c           unchanged (unless FT & FTINV are identical in calling program).
C    LWFFT = argument in last call to CALL FFTCI (LWFFT, WFFTC)
C  Output:
C    FTINV = inverse FFT of FT without rearrangement.  FFTINV is normalized;
C            i.e., the sum is divided by sqrt(N).  This is compatible with
C            DCFFT and CFFT and Siemens FFTs.
C
      COMPLEX FT(N), FTINV(N), ftwork(n)
      REAL WFFTC(4*N+15)
c     -------------------------------------------------------------------------
c     Unrearrange.
c     -------------------------------------------------------------------------
      nunfil = n / 2
      do 110 J=1,nunfil
         ftwork(j)=FT(NUNFIL+J)
         FTwork(NUNFIL+J)=FT(J)
  110 CONTINUE
      IF (N .NE. LWFFT) THEN
         CALL FFTCI (N, WFFTC)
         LWFFT=N
      END IF
      CALL F2TCB (N, FTwork, FTINV, WFFTC)
      FACT=1./SQRT(FLOAT(N))
      DO 210 J=1,N
         FTINV(J)=FTINV(J)*FACT
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFT (DATAT, FT, N, LWFFT,
     1                 WFFTC)
C
C  Input:
C    DATAT(J), J=1,N
C    LWFFT = argument in last call to CALL FFTCI (LWFFT, WFFTC)
C  Output:
C    FT = FFT of DATAT without rearrangement.  The FFT is normalized; i.e.,
C             the sum is divided by sqrt(N).  This is compatible with DCFFT
C             and CFFTIN and Siemens FFTs.
C
      COMPLEX DATAT(N), FT(N)
      REAL WFFTC(4*N+15)
      IF (N .NE. LWFFT) THEN
         CALL FFTCI (N, WFFTC)
         LWFFT=N
      END IF
      CALL F2TCF (N, DATAT, FT, WFFTC)
      FACT=1./SQRT(FLOAT(N))
      DO 210 J=1,N
         FT(J)=FT(J)*FACT
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFT_r (DATAT, FT, N, LWFFT,
     1                   WFFTC)
C
C  Input:
C    DATAT(J), J=1,N
C    LWFFT = argument in last call to CALL FFTCI (LWFFT, WFFTC)
C  Output:
C    FT = FFT of DATAT WITH rearrangement.  The FFT is normalized; i.e.,
C             the sum is divided by sqrt(N).  This is compatible with DCFFT
C             and CFFTIN and Siemens FFTs.
C
      COMPLEX cterm, DATAT(N), FT(N)
      REAL WFFTC(4*N+15)
      IF (N .NE. LWFFT) THEN
         CALL FFTCI (N, WFFTC)
         LWFFT=N
      END IF
      CALL F2TCF (N, DATAT, FT, WFFTC)
C     ------------------------------------------------------------------------
C     Scale & rearrange
C     ------------------------------------------------------------------------
      FACT=1./SQRT(FLOAT(N))
      nunfil = n / 2
      DO 210 J=1,nunfil
         CTERM=FT(NUNFIL+J)
         FT(NUNFIL+J)=FT(J) * fact
         FT(J)=CTERM * fact
  210 CONTINUE
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DF2TCF (N, C, YOUT, WSAVE)
C
C  From Paul Swarztrauber's FFTPACK
C
      DOUBLE PRECISION C(*), YOUT(*), WSAVE(*)
      DO 110 J=1,2*N
         YOUT(J)=C(J)
  110 CONTINUE
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTF1 (N,YOUT,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DCFTF1 (N,C,CH,WA,IFAC)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPASF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPASF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,-.951056516295154D0,
     1-.809016994374947D0,-.587785252292473D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPASF3 (IDO,L1,CC,CH,WA1,WA2)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5D0,-.866025403784439D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPASF2 (IDO,L1,CC,CH,WA1)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DPASF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
C
C  From Paul Swarztrauber's FFTPACK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DFFTCI (N,WSAVE)
C
C  From Paul Swarztrauber's FFTPACK.
C
      DOUBLE PRECISION WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE DCFTI1 (N,WA,IFAC)
C
C  From Paul Swarztrauber's FFTPACK.
C
      DOUBLE PRECISION ARG, ARGH, ARGLD, FI, TPI, WA(*)
      INTEGER IFAC(*), NTRYH(4)
      DATA NTRYH/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959D0
      ARGH = TPI/DBLE(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.D0
            WA(I) = 0.D0
            LD = LD+L1
            FI = 0.D0
            ARGLD = DBLE(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = DCOS(ARG)
               WA(I) = DSIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE F2TCF (N, C, YOUT, WSAVE)
C
C  From Paul Swarztrauber's FFTPACK
C
      REAL C(*), YOUT(*), WSAVE(*)
      DO 110 J=1,2*N
         YOUT(J)=C(J)
  110 CONTINUE
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTF1 (N,YOUT,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,-.951056516295154,
     1-.809016994374947,-.587785252292473/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSF3 (IDO,L1,CC,CH,WA1,WA2)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5,-.866025403784439/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSF2 (IDO,L1,CC,CH,WA1)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE F2TCB (N, C, YOUT, WSAVE)
C
C  From Paul Swarztrauber's FFTPACK
C
      REAL C(*), YOUT(*), WSAVE(*)
      DO 110 J=1,2*N
         YOUT(J)=C(J)
  110 CONTINUE
      IF (N .LE. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTB1 (N,YOUT,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFTB1 (N,C,CH,WA,IFAC)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSB3 (IDO,L1,CC,CH,WA1,WA2)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5,.866025403784439/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSB2 (IDO,L1,CC,CH,WA1)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
C
C  From Paul Swarztrauber's FFTPACK
C
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FFTCI (N,WSAVE)
C
C  From Paul Swarztrauber's FFTPACK.
C
      DIMENSION       WSAVE(*)
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CFFTI1 (N,WA,IFAC)
C
C  From Paul Swarztrauber's FFTPACK.
C
      DIMENSION WA(*), IFAC(*), NTRYH(4)
      DATA NTRYH/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959D0
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD+L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++    4222
C  FUNCTION RANDOM.  PRODUCES A PSEUDORANDOM REAL ON THE OPEN INTERVAL      4223
C      (0.,1.).                                                             4224
C  DIX (IN DOUBLE PRECISION) MUST BE INITIALIZED TO A WHOLE NUMBER          4225
C      BETWEEN 1.D0 AND 2147483646.D0 BEFORE THE FIRST CALL TO RANDOM       4226
C      AND NOT CHANGED BETWEEN SUCCESSIVE CALLS TO RANDOM.                  4227
C  BASED ON L. SCHRAGE, ACM TRANS. ON MATH. SOFTWARE 5, 132 (1979).         4228
C-----------------------------------------------------------------------    4229
      FUNCTION RANDOM(DIX)                                                  4230
C                                                                           4231
C  PORTABLE RANDOM NUMBER GENERATOR                                         4232
C   USING THE RECURSION                                                     4233
C    DIX = DIX*A MOD P                                                      4234
C                                                                           4235
      DOUBLE PRECISION A,P,DIX,B15,B16,XHI,XALO,LEFTLO,FHI,K                4236
C                                                                           4237
C  7**5, 2**15, 2**16, 2**31-1                                              4238
      DATA A/16807.D0/,B15/32768.D0/,B16/65536.D0/,P/2147483647.D0/         4239
C                                                                           4240
C  GET 15 HI ORDER BITS OF DIX                                              4241
      XHI = DIX / B16                                                       4242
      XHI = XHI - DMOD(XHI,1.D0)                                            4243
C  GET 16 LO BITS IF DIX AND FORM LO PRODUCT                                4244
      XALO=(DIX-XHI*B16)*A                                                  4245
C  GET 15 HI ORDER BITS OF LO PRODUCT                                       4246
      LEFTLO = XALO/B16                                                     4247
      LEFTLO = LEFTLO - DMOD(LEFTLO,1.D0)                                   4248
C  FORM THE 31 HIGHEST BITS OF FULL PRODUCT                                 4249
      FHI = XHI*A + LEFTLO                                                  4250
C  GET OVERFLO PAST 31ST BIT OF FULL PRODUCT                                4251
      K = FHI/B15                                                           4252
      K = K - DMOD(K,1.D0)                                                  4253
C  ASSEMBLE ALL THE PARTS AND PRESUBTRACT P                                 4254
C   THE PARENTHESES ARE ESSENTIAL                                           4255
      DIX = (((XALO-LEFTLO*B16) - P) + (FHI-K*B15)*B16) + K                 4256
C  ADD P BACK IN IF NECESSARY                                               4257
      IF (DIX .LT. 0.D0) DIX = DIX + P                                      4258
C  MULTIPLY BY 1/(2**31-1)                                                  4259
      RANDOM=DIX*4.656612875D-10                                            4260
      RETURN                                                                4261
      END                                                                   4262
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  FUNCTION FISHNI.  CALCULATES FISHER F-DISTRIBUTION FOR ARGUMENT F
C      WITH DF1 AND DF2 (NOT NECESSARILY INTEGER) DEGREES OF FREEDOM.
C      (SEE ABRAMOWITZ AND STEGUN, EQS. 26.6.2 AND 26.5.2.)
C-------------------------------------------------------------------------------
C  CALLS SUBPROGRAMS - ERRMES, BETAIN
C  WHICH IN TURN CALL - GAMLN
C-------------------------------------------------------------------------------
      FUNCTION FISHNI (F,DF1,DF2,NOUT)
      CHARACTER CHSUBP*6
      CHSUBP='FISHNI'
      IF (AMIN1(DF1,DF2) .LE. 0.) CALL ERRMES (1, 4, CHSUBP)
      HDF1=.5*DF1
      HDF2=.5*DF2
      DUM=DF1*F
      FISHNI=BETAIN(DUM/(DF2+DUM),HDF1,HDF2,NOUT)
      RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  FUNCTION DGAMLN.  COMPUTES LN OF THE GAMMA FUNCTION FOR POSITIVE
C      XARG USING A VERSION OF CACM ALGORITHM NO. 291.
      DOUBLE PRECISION FUNCTION DGAMLN (XARG)
      DOUBLE PRECISION XARG, X, P, Z
      X=XARG
      P=1.D0
      DGAMLN=0.D0
  110 IF (X .GE. 30.D0) GO TO 150
      P=P*X
      X=X+1.D0
      GO TO 110
  150 IF (XARG .LT. 30.D0) DGAMLN=-DLOG(P)
      Z=1.D0/X**2
      DGAMLN=DGAMLN+(X-.5D0)*DLOG(X)-X+.918938533204672742D0-
     1       (((Z/1680.D0 - 1.D0/1260.D0)*Z + 1.D0/360.D0)*Z -
     2        1.D0/12.D0)/X
      RETURN
      END
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  FUNCTION BETAIN_LCMODL.  APPROXIMATES THE INCOMPLETE BETA FUNCTION RATIO
C      I(SUB X)(A,B) USING ABRAMOWITZ AND STEGUN (26.5.5).
C  Gamma functions evaluated with DGAMLN.
C  TOL made smaller.
C  GOOD TO ABOUT 1 PART IN 1/[max(A,B)*TOL] (TOL SET IN DATA STATEMENT BELOW).
C  FOR VERY LARGE (.GT. 1.E+4) A OR B, OF THE ORDER OF MAX(A,B) TERMS
C      ARE NEEDED, AND AN ASYMPTOTIC FORMULA WOULD BE BETTER.
C      TAKES AN ERROR EXIT IF A OR B .GE. 2.E+4.
C  ERRMES is the LCMODL (not CONTIN) version.
C-------------------------------------------------------------------------------
C  CALLS SUBPROGRAMS - ERRMEs, dGAMLN
C-------------------------------------------------------------------------------
      FUNCTION BETAIN (X,A,B,NOUT)
      CHARACTER CHSUBP*6
      DOUBLE PRECISION DAA, DBB, DGAMLN, DRI
      LOGICAL SWAP
      CHSUBP='BETAIN'
      tol = 1.e-8
      IF (X.LT.0. .OR. X.GT.1. .OR. AMIN1(A,B).LE.0. .OR.
     1    AMAX1(A,B).GE.2.E+4) CALL ERRMES (1, 4, CHSUBP)
      BETAIN=X
      SWAP=X .GT. .5
      IF (SWAP) GO TO 150
      XX=X
      AA=A
      BB=B
      GO TO 200
C-------------------------------------------------------------------------------
C      WHEN SWAP=.TRUE., I(SUB 1-X)(B,A)=1-I(SUB X)(A,B) IS EVALUATED
C      FIRST.
C-------------------------------------------------------------------------------
  150 XX=1.D0-X
      AA=B
      BB=A
  200 CX=1.D0-XX
      IF (AMIN1(XX,X).LE.0. .OR. AMAX1(CX,X).GE.1.) RETURN
      R=XX/CX
C-------------------------------------------------------------------------------
C  TERM IMAX IS APPROXIMATELY THE MAXIMUM TERM IN THE SUM.
C  0 < R < 1 implies IMAX < BB.
C-------------------------------------------------------------------------------
      IMAX=MAX0(0,INT((R*BB-AA-1.)/(R+1.)))
      RI=FLOAT(IMAX)
      SUM=0.
      DAA=AA
      DRI=DBLE(IMAX)
      DBB=BB
      TERMAX=(DAA+DRI)*DLOG(DBLE(XX))+(DBB-DRI-1.D0)*DLOG(DBLE(CX))+
     1 DGAMLN(DAA+DBB)-DGAMLN(DAA+DRI+1.D0)-DGAMLN(DBB-DRI)
      IF (TERMAX .LT. -50.) GO TO 700
      TERMAX=EXP(TERMAX)
      TERM=TERMAX
      SUM=TERM
C-------------------------------------------------------------------------------
C  SUM TERMS FOR I=IMAX+1,IMAX+2,... UNTIL CONVERGENCE.
C-------------------------------------------------------------------------------
      I1=IMAX+1
      DO 250 I=I1,40000
        TNUMER=BB-FLOAT(I)
        TERM=TERM*R*TNUMER/(AA+FLOAT(I))
        SUM=SUM+TERM
C       -----------------------------------------------------------------------
C       TNUMER = 0 implies that the denominators of all following terms contain
C                  beta functions with zero arguments, which make the terms 0.
C                  This occurs when BB is an integer.
C       -----------------------------------------------------------------------
        IF (ABS(TERM).LE.TOL*SUM .OR. ABS(TNUMER).LE.1.E-3) GO TO 300
  250 CONTINUE
      CALL ERRMES (2, 3, CHSUBP)
  300 IF (IMAX .EQ. 0) GO TO 700
C-------------------------------------------------------------------------------
C  SUM TERMS FOR I=IMAX-1,IMAX-2,... UNTIL CONVERGENCE.
C-------------------------------------------------------------------------------
      TERM=TERMAX
      DO 320 I=IMAX,1,-1
        RI=FLOAT(I)
        TERM=TERM*(AA+RI)/(R*(BB-RI))
        SUM=SUM+TERM
        IF (ABS(TERM) .LE. TOL*SUM) GO TO 700
  320 CONTINUE
  700 BETAIN=SUM
      IF (SWAP) BETAIN=1.D0-BETAIN
      RETURN
      END
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE PNNLS  (A,MDA,M,N,B,X,VAR,W,ZZ,INDEX,MODE,RANGE,
C                         NONNEG,DVARAC,NSETP)
C  NSETP = no. of positive solution points
C        = no. of degrees of freedom
C  NONNEG(J) = T if parameter J is to be constrained to be nonnegative
C  DVAR = DVARAC + RNORM**2 is returned rather than RNORM.
C  FACTOR is decreased by a factor of RFACTR up to MFACTR times when there are
C    more than 2*N iterations in the inner loop.  See the DATA statement below.
C     BASED ON C.L.LAWSON AND R.J.HANSON,
C     'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C
C         **********   NONNEGATIVE LEAST SQUARES   **********
C
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, WHICH SOLVES THE LEAST SQUARES PROBLEM
C
C                      A * X = B  SUBJECT TO X .GE. 0
C
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N
C                     MATRIX, A.           ON EXIT A() CONTAINS
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY
C                     THIS SUBROUTINE.
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON-
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL
C             CONTAIN THE SOLUTION VECTOR.
C     DVAR    ON EXIT DVAR CONTAINS THE squared EUCLIDEAN NORM OF THE
C             RESIDUAL VECTOR.  DVAR is DOUBLE PRECISION.
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z
C     ZZ()     AN M-ARRAY OF WORKING SPACE.
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
C                 P AND Z AS FOLLOWS..
C
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING
C             MEANINGS.
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 2*N*MFACTR ITERATIONS.
C
C  RANGE IS 2 OR 3 ORDERS OF MAGNITUDE SMALLER THAN BIG, WHERE BIG IS
C      THE LARGEST NUMBER THAT DOES NOT OVERFLOW AND 1/BIG DOES NOT
C      UNDERFLOW.  FOR THE DOUBLE PRECISION VERSION, BIG AND RANGE
C      ARE IN DOUBLE PRECISION.  FOR THE SINGLE PRECISION VERSION,
C      THEY ARE IN SINGLE PRECISION (AND THEREFORE RANGE=SRANGE).
C-------------------------------------------------------------------------------
C  CALLS SUBPROGRAMS - DIFF, H12, G1, G2
C-------------------------------------------------------------------------------
      SUBROUTINE PNNLS (A,MDA,M,N,B,X,DVAR,W,ZZ,INDEX,MODE,RANGE,
     1                   NONNEG,DVARAC,NSETP)
      DOUBLE PRECISION A, ABS, ALPHA, ASAVE, B, CC, DIFF, DUMMY,
     1 DVAR, DVARAC, FACTOR, RANGE, SM, SQRT, SS, T, TWO, UNORM,
     2 UP, W, WMAX, X, ZERO, ZTEST, ZZ
      DIMENSION A(MDA,N), B(*), X(*), W(*), ZZ(*)
      INTEGER INDEX(N)
      CHARACTER CHSUBP*6
      LOGICAL NONNEG(N)
      DATA MFACTR/4/, RFACTR/.01/
      ABS(T)=DABS(T)
      SQRT(T)=DSQRT(T)
      ZERO=0.D0
      TWO=2.D0
      FACTOR=1.D-2
C
      CHSUBP='PNNLS'
      NFACTR=0
      MODE=1
      IF (M.GT.0.AND.N.GT.0) GO TO 10
      MODE=2
      RETURN
   10 ITER=0
      ITMAX=2*N
C
C                    INITIALIZE THE ARRAYS INDEX() AND X().
C
          DO 20 I=1,N
          X(I)=ZERO
   20     INDEX(I)=I
C
      IZ2=N
      IZ1=1
      NSETP=0
      NPP1=1
C                             ******  MAIN LOOP BEGINS HERE  ******
   30 CONTINUE
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.
C
      IF (IZ1.GT.IZ2.OR.NSETP.GE.M) GO TO 350
C
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C
          DO 50 IZ=IZ1,IZ2
          J=INDEX(IZ)
          SM=ZERO
              DO 40 L=NPP1,M
   40         SM=SM+A(L,J)*B(L)
          IF (NONNEG(J)) THEN
             W(J)=SM
          ELSE
             W(J)=ABS(SM)
          END IF
   50     CONTINUE
 
C                                   FIND LARGEST POSITIVE W(J).
   60 WMAX=ZERO
          DO 70 IZ=IZ1,IZ2
          J=INDEX(IZ)
          IF (W(J).LE.WMAX) GO TO 70
          WMAX=W(J)
          IZMAX=IZ
   70     CONTINUE
C
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C
      IF (WMAX) 350,350,80
   80 IZ=IZMAX
      J=INDEX(IZ)
C
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
C     NEAR LINEAR DEPENDENCE.
C
      ASAVE=A(NPP1,J)
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0,RANGE)
      UNORM=ZERO
      IF (NSETP.EQ.0) GO TO 100
          DO 90 L=1,NSETP
   90     UNORM=UNORM+A(L,J)**2
      UNORM=SQRT(UNORM)
  100 IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM)) 130,130,110
C
C     COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ AND
C     SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
C
  110     DO 120 L=1,M
  120     ZZ(L)=B(L)
      CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1,RANGE)
      ZTEST=ZZ(NPP1)/A(NPP1,J)
      IF (.NOT.NONNEG(J)) ZTEST=ABS(ZTEST)
C
C                                     SEE IF ZTEST IS POSITIVE
C
      IF (ZTEST) 130,130,140
C
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.
C
  130 A(NPP1,J)=ASAVE
      W(J)=ZERO
      GO TO 60
C
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
C     COL J,  SET W(J)=0.
C
  140     DO 150 L=1,M
  150     B(L)=ZZ(L)
C
      INDEX(IZ)=INDEX(IZ1)
      INDEX(IZ1)=J
      IZ1=IZ1+1
      NSETP=NPP1
      NPP1=NPP1+1
C
      IF (IZ1.GT.IZ2) GO TO 170
          DO 160 JZ=IZ1,IZ2
          JJ=INDEX(JZ)
  160     CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1,RANGE)
  170 CONTINUE
C
      IF (NSETP.EQ.M) GO TO 190
          DO 180 L=NPP1,M
  180     A(L,J)=ZERO
  190 CONTINUE
C
      W(J)=ZERO
C                                SOLVE THE TRIANGULAR SYSTEM.
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      next = 200
      GO TO 400
  200 CONTINUE
C
C                       ******  SECONDARY LOOP BEGINS HERE ******
C
C                          ITERATION COUNTER.
C
  210 ITER=ITER+1
      IF (ITER .GT. ITMAX) THEN
         IF (NFACTR .GE. MFACTR) THEN
            MODE=3
            GO TO 350
         ELSE
            CALL ERRMES (1, 1, CHSUBP)
            NFACTR=NFACTR+1
            FACTOR=FACTOR*RFACTR
            ITER=0
         END IF
      END IF
C
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
C                                  IF NOT COMPUTE ALPHA.
C
      ALPHA=TWO
          DO 240 IP=1,NSETP
          L=INDEX(IP)
          IF (NONNEG(L) .AND. ZZ(IP).LE.ZERO) THEN
             T=-X(L)/(ZZ(IP)-X(L))
             IF (ALPHA.LE.T) GO TO 240
             ALPHA=T
             JJ=IP
          END IF
  240     CONTINUE
C
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.
C
      IF (ALPHA.EQ.TWO) GO TO 330
C
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.
C
          DO 250 IP=1,NSETP
          L=INDEX(IP)
  250     X(L)=X(L)+ALPHA*(ZZ(IP)-X(L))
C
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
C        FROM SET P TO SET Z.
C
      I=INDEX(JJ)
  260 X(I)=ZERO
C
      IF (JJ.EQ.NSETP) GO TO 290
      JJ=JJ+1
          DO 280 J=JJ,NSETP
          II=INDEX(J)
          INDEX(J-1)=II
          CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))
          A(J,II)=ZERO
              DO 270 L=1,N
              IF (L.NE.II) CALL G2 (CC,SS,A(J-1,L),A(J,L))
  270         CONTINUE
  280     CALL G2 (CC,SS,B(J-1),B(J))
  290 NPP1=NSETP
      NSETP=NSETP-1
      IZ1=IZ1-1
      INDEX(IZ1)=I
C
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO
C        AND MOVED FROM SET P TO SET Z.
C
          DO 300 JJ=1,NSETP
          I=INDEX(JJ)
          IF (NONNEG(I) .AND. X(I).LE.ZERO) GO TO 260
  300     CONTINUE
C
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C
          DO 310 I=1,M
  310     ZZ(I)=B(I)
      next = 320
      GO TO 400
  320 CONTINUE
      GO TO 210
C                      ******  END OF SECONDARY LOOP  ******
C
  330     DO 340 IP=1,NSETP
          I=INDEX(IP)
  340     X(I)=ZZ(IP)
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
      GO TO 30
C
C                        ******  END OF MAIN LOOP  ******
C
C                        COME TO HERE FOR TERMINATION.
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.
C
  350 SM=ZERO
      IF (NPP1.GT.M) GO TO 370
          DO 360 I=NPP1,M
  360     SM=SM+B(I)**2
      GO TO 390
  370     DO 380 J=1,N
  380     W(J)=ZERO
  390 DVAR=SM+DVARAC
      RETURN
C
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().
C
  400     DO 430 L=1,NSETP
          IP=NSETP+1-L
          IF (L.EQ.1) GO TO 420
              DO 410 II=1,IP
  410         ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)
  420     JJ=INDEX(IP)
  430     ZZ(IP)=ZZ(IP)/A(IP,JJ)
      if (next .eq. 200) go to 200
      if (next .eq. 320) go to 320
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  FUNCTION DIFF.
C     BASED ON C.L.LAWSON AND R.J.HANSON,
C     'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C     FUNCTION DIFF(X,Y)!SP
      DOUBLE PRECISION FUNCTION DIFF(X,Y)
      DOUBLE PRECISION X, Y
      DIFF=X-Y
      RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  SUBROUTINE G1.
      SUBROUTINE G1 (A,B,COS,SIN,SIG)
C     BASED ON C.L.LAWSON AND R.J.HANSON,
C     'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C
C
C     COMPUTE ORTHOGONAL ROTATION MATRIX..
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
C                        (-S,C)         (-S,C)(B)   (   0          )
C     COMPUTE SIG = SQRT(A**2+B**2)
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C
      DOUBLE PRECISION A, ABS, B, COS, ONE, SIG, SIGN, SIN, SQRT,
     1 XR, YR, ZERO
      ABS(A)=DABS(A)
      SQRT(A)=DSQRT(A)
      SIGN(A,B)=DSIGN(A,B)
C     ZERO=0.E0!SP
      ZERO=0.D0
C     ONE=1.E0!SP
      ONE=1.D0
      IF (ABS(A).LE.ABS(B)) GO TO 10
      XR=B/A
      YR=SQRT(ONE+XR**2)
      COS=SIGN(ONE/YR,A)
      SIN=COS*XR
      SIG=ABS(A)*YR
      RETURN
   10 IF (B) 20,30,20
   20 XR=A/B
      YR=SQRT(ONE+XR**2)
      SIN=SIGN(ONE/YR,B)
      COS=SIN*XR
      SIG=ABS(B)*YR
      RETURN
   30 SIG=ZERO
      COS=ZERO
      SIN=ONE
      RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  SUBROUTINE G2.
      SUBROUTINE G2    (COS,SIN,X,Y)
C     BASED ON C.L.LAWSON AND R.J.HANSON,
C     'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C          APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).
      DOUBLE PRECISION COS, SIN, X, XR, Y
      XR=COS*X+SIN*Y
      Y=-SIN*X+COS*Y
      X=XR
      RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV,RANGE)
C     BASED ON C.L.LAWSON AND R.J.HANSON,
C     'SOLVING LEAST SQUARES PROBLEMS', PRENTICE-HALL, 1974
C
C     CONSTRUCTION AND/OR APPLICATION OF A SINGLE
C     HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B
C
C     MODE    = 1 OR 2   TO SELECT ALGORITHM  H1  OR  H2 .
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    ON ENTRY TO H1 U() CONTAINS THE PIVOT VECTOR.
C                   IUE IS THE STORAGE INCREMENT BETWEEN ELEMENTS.
C                                       ON EXIT FROM H1 U() AND UP
C                   CONTAIN QUANTITIES DEFINING THE VECTOR U OF THE
C                   HOUSEHOLDER TRANSFORMATION.   ON ENTRY TO H2 U()
C                   AND UP SHOULD CONTAIN QUANTITIES PREVIOUSLY COMPUTED
C                   BY H1.  THESE WILL NOT BE MODIFIED BY H2.
C     C()    ON ENTRY TO H1 OR H2 C() CONTAINS A MATRIX WHICH WILL BE
C            REGARDED AS A SET OF VECTORS TO WHICH THE HOUSEHOLDER
C            TRANSFORMATION IS TO BE APPLIED.  ON EXIT C() CONTAINS THE
C            SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0
C            NO OPERATIONS WILL BE DONE ON C().
C  RANGE IS 2 OR 3 ORDERS OF MAGNITUDE SMALLER THAN BIG, WHERE BIG IS
C      THE LARGEST NUMBER THAT DOES NOT OVERFLOW AND 1/BIG DOES NOT
C      UNDERFLOW.  FOR THE DOUBLE PRECISION VERSION, BIG AND RANGE
C      ARE IN DOUBLE PRECISION.  FOR THE SINGLE PRECISION VERSION,
C      THEY ARE IN SINGLE PRECISION (AND THEREFORE RANGE=SRANGE).
C
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV,RANGE)
      DIMENSION U(IUE,*), C(*)
      DOUBLE PRECISION SM,B
      DOUBLE PRECISION ABS, AMAX1, C, CL, CLINV, DOUBLE, ONE, RANGE,
     1 RANGIN, SIGN, SM1, SQRT, U, UP
      ABS(SM)=DABS(SM)
      AMAX1(SM,ONE)=DMAX1(SM,ONE)
C     DOUBLE(ONE)=DBLE(ONE)!SP
      DOUBLE(SM)=SM
      SQRT(SM)=DSQRT(SM)
      SIGN(SM,ONE)=DSIGN(SM,ONE)
C     ONE=1.E0!SP
      ONE=1.D0
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      RANGIN=ONE/RANGE
      CL=ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GO TO 60
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M
   10     CL=AMAX1(ABS(U(1,J)),CL)
      IF (CL .LE. RANGIN) GO TO 130
      CLINV=ONE/CL
      SM=(DOUBLE(U(1,LPIVOT))*CLINV)**2
          DO 30 J=L1,M
   30     SM=SM+(DOUBLE(U(1,J))*CLINV)**2
C                          CONVERT DOUBLE PREC. SM TO SNGL. PREC. SM1
      SM1=SM
      CL=-SIGN(CL*SQRT(SM1),U(1,LPIVOT))
      UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GO TO 70
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
   60 IF (CL .LE. RANGIN) GO TO 130
   70 IF (NCV.LE.0) RETURN
      B=DOUBLE(UP)*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
      IF (B .GE. -RANGIN) GO TO 130
      B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
          DO 120 J=1,NCV
          I2=I2+ICV
          I3=I2+INCR
          I4=I3
          SM=C(I2)*DOUBLE(UP)
              DO 90 I=L1,M
              SM=SM+C(I3)*DOUBLE(U(1,I))
   90         I3=I3+ICE
          IF (SM) 100,120,100
  100     SM=SM*B
          C(I2)=C(I2)+SM*DOUBLE(UP)
              DO 110 I=L1,M
              C(I4)=C(I4)+SM*DOUBLE(U(1,I))
  110         I4=I4+ICE
  120     CONTINUE
  130 RETURN
      END
C++++++++++++++++ DOUBLE PRECISION VERSION 2DP (MAR 1984) ++++++++++++++
C  SUBROUTINE PLPRIN.  PLOTS Y1 (AND Y2, IF ONLY1=.FALSE.) VS. X ON
C      LINE PRINTER, AND PRINTS Y1 AND X.
C  IF PLTERR=T, AN ERROR BAND IS ALSO PLOTTED USING THE ERRORS
C      SUPPLIED IN YERR.
C  IF NLINF.GT.0, THEN THE NLINF LINEAR COEFFICIENTS Y(J), J=NG+1,MY1
C      ARE PRINTED OUT.  IN THIS CASE MY1=NGL NORMALLY.
      SUBROUTINE PLPRIN (X,Y1,Y2,N,ONLY1,NOUT,SRANGE,NLINF,NG,MY1,
     1 YERR, PLTERR)
      character*1 char(5), ih(109)
      DOUBLE PRECISION YERR, DUB
      LOGICAL ONLY1, PLTERR
      DIMENSION X(*), Y1(*), Y2(*), YERR(*)
      DATA CHAR/' ', 'X', 'O', '*', '.'/
C     SINGLE(DUB)=DUB!SP
      SINGLE(DUB)=SNGL(DUB)
      DUB=1.D0
      YMIN=SRANGE
      YMAX=-SRANGE
      DO 120 J=1,N
        YMIN=AMIN1(YMIN,Y1(J))
        YMAX=AMAX1(YMAX,Y1(J))
        IF (ONLY1) GO TO 120
        YMIN=AMIN1(YMIN,Y2(J))
        YMAX=AMAX1(YMAX,Y2(J))
  120 CONTINUE
      DUM=YMAX-YMIN
      NCHAR=109
      IF (DUM .LE. FLOAT(NCHAR)/SRANGE) DUM=1.
      IF (PLTERR) GO TO 130
      WRITE (NOUT,5120)
 5120 FORMAT (/4X,8HORDINATE,2X,8HABSCISSA)
      GO TO 140
  130 NCHAR=100
 5130 FORMAT (/4X,8HORDINATE,4X,5HERROR,2X,8HABSCISSA)
      WRITE (NOUT,5130)
  140 R=(FLOAT(NCHAR)-.001)/DUM
      DO 150 J=1,N
        DO 155 L1=1,NCHAR
          IH(L1)=CHAR(1)
  155   CONTINUE
        IF (.NOT.PLTERR) GO TO 158
        LMIN=INT((AMAX1(YMIN,Y1(J)-ABS(SINGLE(YERR(J))))-YMIN)*R)+1
        LMAX=INT((AMIN1(YMAX,Y1(J)+ABS(SINGLE(YERR(J))))-YMIN)*R)+1
        IF (LMIN .GE. LMAX) GO TO 158
        DO 156 L1=LMIN,LMAX
          IH(L1)=CHAR(5)
  156   CONTINUE
  158   L1=INT((Y1(J)-YMIN)*R)+1
        IH(L1)=CHAR(2)
        IF (ONLY1) GO TO 160
        L2=INT((Y2(J)-YMIN)*R)+1
        IH(L2)=CHAR(3)
        IF (L1 .EQ. L2) IH(L2)=CHAR(4)
  160   IF (.NOT.PLTERR) WRITE (NOUT,5160) Y1(J),X(J),IH
 5160   FORMAT (1X,1PE11.3,E10.2,109A1)
        IF (PLTERR) WRITE (NOUT,5161) Y1(J),YERR(J),X(J),
     1  (IH(L1),L1=1,NCHAR)
C5161   FORMAT (1X,1PE11.3,E9.1,E10.2,100A1)!SP
 5161   FORMAT (1X,1PE11.3,D9.1,E10.2,100A1)
  150 CONTINUE
      IF (NLINF .LE. 0) GO TO 800
      L2=NG+1
 5200 FORMAT (22H0LINEAR COEFFICIENTS =,1P8E13.4/(22X,8E13.4))
      IF (.NOT.PLTERR) WRITE (NOUT,5200) (Y1(J),J=L2,MY1)
      IF (PLTERR) WRITE(NOUT,5201) (Y1(J),YERR(J),J=L2,MY1)
C5201 FORMAT (22H0LINEAR COEFFICIENTS =,!SP
C    1 1PE13.4,3H +-,E9.1,E20.4,3H +-,E9.1,E20.4,3H +-,E9.1/!SP
C    2 (22X,1PE13.4,3H +-,E9.1,E20.4,3H +-,E9.1,E20.4,3H +-,E9.1))!SP
 5201 FORMAT (22H0LINEAR COEFFICIENTS =,
     1 1PE13.4,3H +-,D9.1,E20.4,3H +-,D9.1,E20.4,3H +-,D9.1/
     2 (22X,1PE13.4,3H +-,D9.1,E20.4,3H +-,D9.1,E20.4,3H +-,D9.1))
  800 RETURN
      END
      subroutine EIGVrs(nm,n,a,w,z,fv1,fv2,ierr)
c
      integer n,nm,ierr
      real a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
c     .......... find both eigenvalues and eigenvectors ..........
   10 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real d(n),e(n),z(nm,n)
      real c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0e0
      tst1 = 0.0e0
      e(n) = 0.0e0
c
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0e0 * e(l))
         r = pythag(p,1.0e0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0e0
         c2 = c
         el1 = e(l1)
         s = 0.0e0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      real a(nm,n),d(n),e(n),z(nm,n)
      real f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0e0
         scale = 0.0e0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(d(k))
c
         if (scale .ne. 0.0e0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0e0
            z(j,i) = 0.0e0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0e0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0e0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0e0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0e0
         h = d(i)
         if (h .eq. 0.0e0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0e0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0e0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0e0
  520 continue
c
      z(n,n) = 1.0e0
      e(1) = 0.0e0
      return
      end
C
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function pythag(a,b)
      real a,b
c
c     finds sqrt(a**2+b**2) without overflow or destructive underflow
c
      real p,r,s,t,u
      p = amax1(abs(a),abs(b))
      if (p .eq. 0.0e0) go to 20
      r = (amin1(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         u = 1.0e0 + 2.0e0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
