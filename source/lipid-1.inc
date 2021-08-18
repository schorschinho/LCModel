c     For breast-1  lipid-1  old-lipid-1  breast-2  only-cho-1  lipid-2  
c     prostate-a-e
c     111111111111111111111111111111111111111111111111111111111111111111111111
c     Modify liver-1 settings (from liver-1.inc) for breast-1, and for initial 
c        settings for lipid-1, which will be modified below.
c     *** Note *** The subscripts below must agree with liver-1.inc.
c     111111111111111111111111111111111111111111111111111111111111111111111111
      if (sptype(:9) .eq. 'breast-1 '  .or.
     1    sptype(:8) .eq. 'lipid-1 ') then

         chnot1(27)= 'cho4'
         chnot1(28)= 'cho5'
         chnot1(29)= 'cho6'

         chsimu(41)= 'cho1 @ 3.20+-.015 FWHM=.06<.08+-.01 AMP=1.'
         chsimu(42)= 'cho2 @ 3.20+-.015 FWHM=.08<.10+-.01 AMP=1.'
         chsimu(43)= 'cho3 @ 3.20+-.015 FWHM=.10<.12+-.01 AMP=1.'
         chsimu(44)= 'cho4 @ 3.20+-.015 FWHM=.12<.14+-.01 AMP=1.'
         chsimu(45)= 'cho5 @ 3.20+-.015 FWHM=.14<.16+-.01 AMP=1.'
         chsimu(46)= 'cho6 @ 3.20+-.015 FWHM=.16<.18+-.01 AMP=1.'

         nnot1= 29
         nsimul= 46
         ppmst= 4.0
      end if

c     222222222222222222222222222222222222222222222222222222222222222222222222
c     Modify liver-2 settings (from liver-1.inc) for breast-2 & only-cho-1.
c     *** Note *** The subscripts below must agree with liver-1.inc.
c     222222222222222222222222222222222222222222222222222222222222222222222222
      if (sptype(:8) .eq. 'breast-2'   .or.   
     1    sptype(:10) .eq. 'only-cho-1'   .or.
     2    sptype(:8) .eq. 'breast-3'   .or.   
     3    sptype(:10) .eq. 'only-cho-2') then
         chsimu(38)= 'cho1 @ 3.185+-.015 FWHM=.030<.050+-.012 AMP=1.'
         chsimu(39)= 'cho2 @ 3.185+-.015 FWHM=.074<.086+-.012 AMP=1.'
         chsimu(40)= 'cho3 @ 3.185+-.015 FWHM=.110<.122+-.012 AMP=1.'
         chsimu(41)= 'cho4 @ 3.185+-.015 FWHM=.146<.158+-.012 AMP=1.'
         chsimu(42)= 'cho5 @ 3.185+-.015 FWHM=.182<.194+-.012 AMP=1.'
         chsimu(43)= 'cho6 @ 3.185+-.015 FWHM=.218<.230+-.012 AMP=1.'
         chsimu(44)= 'cho7 @ 3.185+-.015 FWHM=.254<.266+-.012 AMP=1.'
c        ---------------------------------------------------------------------
c        Vienna email in 080312 claims there are spectra with only Cho.
c        Also, Cho is often prominent in breast.
c        ---------------------------------------------------------------------
         if (ppmst .lt. 5.) then
            nrefpk(2)= 2
            ppmref(1,2)= 3.185
            ppmref(2,2)= 1.28
c           ------------------------------------------------------------------
c           SHIFMN(2): do not hit water instead of Cho
c           SHIFMX(2): do not hit Lip13 instead of Cho
c           Only go halfway, since these limits assume that there is no 
c              referencing shift.
c           ------------------------------------------------------------------
            shifmn(2)= -.7
            shifmx(2)= 1.
         end if
         nsimul= 44
         ppmst= 3.8
         if (sptype(:10) .eq. 'only-cho-1'   .or.
     1       sptype(:10) .eq. 'only-cho-2') then
            chomit(1)= 'L20a'
            chomit(2)= 'L20b'
            chomit(3)= 'L20c'
            chomit(4)= 'L20d'
            chomit(5)= 'L20e'
            chomit(6)= 'L20f'
            dorefs(1)= .true.
            dorefs(2)= .false.
            fh2omx= .8
            namrel= 'Cho'
            nomit= 6
            nrefpk(2)= 1
            ppmend= 2.7
            rsdgp3= 1.05
         end if
      end if
c     111111111111111111111111111111111111111111111111111111111111111111111111
c     Modify breast-1 settings above for lipid-1.
c     *** Note *** The subscripts below must agree with liver-1.inc.
c     111111111111111111111111111111111111111111111111111111111111111111111111
      if (sptype(:7) .eq. 'lipid-1') then
         nsimul= 40
         ppmst= 3.4
      end if

c     222222222222222222222222222222222222222222222222222222222222222222222222
c     Modify liver-2 settings (from liver-1.inc) for lipid-2.
c     *** Note *** The subscripts below must agree with liver-1.inc.
c     222222222222222222222222222222222222222222222222222222222222222222222222
      if (sptype(:7) .eq. 'lipid-2'   .or.
     1    sptype(:7) .eq. 'lipid-3') then
         if (ppmst .lt. 5.) idgppm= 0
         nsimul= 37
         ppmst= 3.6
      end if

c     444444444444444444444444444444444444444444444444444444444444444444444444
c     Modify liver-4 settings for breast-4
c     444444444444444444444444444444444444444444444444444444444444444444444444
      if (sptype(:8) .eq. 'breast-4') then
         chsimu(42)= 'cho1 @ 3.185+-.015 FWHM=.030<.050+-.012 AMP=1.'
         chsimu(43)= 'cho2 @ 3.185+-.015 FWHM=.074<.086+-.012 AMP=1.'
         chsimu(44)= 'cho3 @ 3.185+-.015 FWHM=.110<.122+-.012 AMP=1.'
         chsimu(45)= 'cho4 @ 3.185+-.015 FWHM=.146<.158+-.012 AMP=1.'
         chsimu(46)= 'cho5 @ 3.185+-.015 FWHM=.182<.194+-.012 AMP=1.'
         chsimu(47)= 'cho6 @ 3.185+-.015 FWHM=.218<.230+-.012 AMP=1.'
         chsimu(48)= 'cho7 @ 3.185+-.015 FWHM=.254<.266+-.012 AMP=1.'
c        ---------------------------------------------------------------------
c        Vienna email in 080312 claims there are spectra with only Cho.
c        Also, Cho is often prominent in breast.
c        ---------------------------------------------------------------------
         if (ppmst .lt. 5.) then
            nrefpk(2)= 2
            ppmref(1,2)= 3.185
            ppmref(2,2)= 1.28
c           ------------------------------------------------------------------
c           SHIFMN(2): do not hit water instead of Cho
c           SHIFMX(2): do not hit Lip13 instead of Cho
c           Only go halfway, since these limits assume that there is no 
c              referencing shift.
c           ------------------------------------------------------------------
            shifmn(2)= -.7
            shifmx(2)= 1.
         end if
         nsimul= 48
         ppmst= 3.8
      end if

c     444444444444444444444444444444444444444444444444444444444444444444444444
c     Modify liver-4 settings for lipid-4
c     444444444444444444444444444444444444444444444444444444444444444444444444
      if (sptype(:7) .eq. 'lipid-4') then
         nsimul=41
      end if

c     -----------------------------------------------------------------------
c     Modify liver-5-10 settings for breast-5-10
c     -----------------------------------------------------------------------
      if (sptype(:8) .eq. 'breast-5'   .or.
     1    sptype(:8) .eq. 'breast-6'   .or.
     1    sptype(:8) .eq. 'breast-7'   .or.
     1    sptype(:8) .eq. 'breast-8'   .or.
     1    sptype(:8) .eq. 'breast-9'   .or.
     2    sptype(:9) .eq. 'breast-10') then
         if (sptype(:8) .eq. 'breast-5') then
            chsimu(48)= 'cho1 @ 3.185+-.015 FWHM=.03<.05+-.01 AMP=1.'
            chsimu(49)= 'cho2 @ 3.185+-.015 FWHM=.07<.09+-.01 AMP=1.'
            chsimu(50)= 'cho3 @ 3.185+-.015 FWHM=.11<.13+-.01 AMP=1.'
            chsimu(51)= 'cho4 @ 3.185+-.015 FWHM=.15<.17+-.01 AMP=1.'
            chsimu(52)= 'cho5 @ 3.185+-.015 FWHM=.19<.21+-.01 AMP=1.'
            nsimul= 52
            ppmref(1,2)= 3.185
         else if (sptype(:8) .eq. 'breast-6') then
            chsimu(52)= 'cho1 @ 3.185+-.015 FWHM=.03<.05+-.01 AMP=1.'
            chsimu(53)= 'cho2 @ 3.185+-.015 FWHM=.07<.09+-.01 AMP=1.'
            chsimu(54)= 'cho3 @ 3.185+-.015 FWHM=.11<.13+-.01 AMP=1.'
            chsimu(55)= 'cho4 @ 3.185+-.015 FWHM=.15<.17+-.01 AMP=1.'
            chsimu(56)= 'cho5 @ 3.185+-.015 FWHM=.19<.21+-.01 AMP=1.'
            nsimul= 56
            ppmref(1,2)= 3.185
         else if (sptype(:8) .eq. 'breast-7'   .or.
     1            sptype(:8) .eq. 'breast-8'   .or.
     2            sptype(:9) .eq. 'breast-10') then
            nsimul= 42
            ppmref(1,2)= 3.2
         else if (sptype(:8) .eq. 'breast-9') then
            nsimul= 15
            ppmref(1,2)= 3.2
         end if
c        ---------------------------------------------------------------------
c        Vienna email in 080312 claims there are spectra with only Cho.
c        Also, Cho is often prominent in breast.
c        ---------------------------------------------------------------------
         nrefpk(2)= 2
         ppmref(2,2)= 1.28
c        ---------------------------------------------------------------------
c        SHIFMN(2): do not hit water instead of Cho
c        SHIFMX(2): do not hit Lip13 instead of Cho
c        Only go halfway, since these limits assume that there is no 
c           referencing shift.
c        ---------------------------------------------------------------------
         shifmn(2)= -.7
         shifmx(2)= 1.
         if (.not.dowatr) ppmgap(2,1)= 3.8
         ppmst= 3.8
      end if

c     555555555555555555555555555555555555555555555555555555555555555555555555
c     Modify liver-5 settings for lipid-5
c     555555555555555555555555555555555555555555555555555555555555555555555555
      if (sptype(:7) .eq. 'lipid-5') nsimul= 47

c     666666666666666666666666666666666666666666666666666666666666666666666666
c     Modify liver-6 settings for lipid-6
c     666666666666666666666666666666666666666666666666666666666666666666666666
      if (sptype(:7) .eq. 'lipid-6') then
         nsimul= 51
         if (.not.dowatr   .or.   ppmst .lt. 5.) nuse1= 11
      end if

c     ------------------------------------------------------------------------
c     Modify liver-7 8 & 10 settings for lipid-7 8 & 10
c     ------------------------------------------------------------------------
      if (sptype(:7) .eq. 'lipid-7'   .or.
     1    sptype(:7) .eq. 'lipid-8'  .or.
     1    sptype(:8) .eq. 'lipid-10') nsimul= 39

c     999999999999999999999999999999999999999999999999999999999999999999999999
c     Modify liver-9 settings for lipid-9
c     999999999999999999999999999999999999999999999999999999999999999999999999
      if (sptype(:7) .eq. 'lipid-9') nsimul= 14

c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-a' thru 'prostate-e'
c     prostate-a: 4 peaks (for 3T) with no constraints.  Used initially for 
c                 Dijon; superceded (for Dijon) by prostate-d.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-a'   .or.   
     1    sptype(:10) .eq. 'prostate-b'   .or.
     2    sptype(:10) .eq. 'prostate-c'   .or.
     3    sptype(:10) .eq. 'prostate-d'   .or.
     4    sptype(:10) .eq. 'prostate-e') then
         is_sptype= .true.

         atth2o= 1.
         badref= .true.

         chcom2(1)= 'Cho'
         chcom2(2)= 'Cr'
         chcom2(3)= 'PA'
         chcom2(4)= 'Cit'
         chcom2(5)= 'PA+Cr'
         chcom2(6)= 'Cho+Cr'
         chcom2(7)= 'Cho+PA+Cr'

         chcomb(1)= 'Cho1+Cho2'
         chcomb(2)= 'Cr1+Cr2'
         chcomb(3)= 'PA1+PA2+PA3+PA4'
         chcomb(4)= 'Cit1a+Cit2a+Cit3a+Cit4a+Cit1b+Cit2b+Cit3b+Cit4b+Cit
     11c+Cit2c+Cit3c+Cit4c'
         chcomb(5)= 'PA1+PA2+PA3+PA4+Cr1+Cr2'
         chcomb(6)= 'Cho1+Cho2+Cr1+Cr2'
         chcomb(7)= 'Cho1+Cho2+PA1+PA2+PA3+PA4+Cr1+Cr2'

         chgrsh(1)= 'Cit1'
         chgrsh(2)= 'Cit2'
         chgrsh(3)= 'Cit3'
         chgrsh(4)= 'Cit4'
         chgrsh(5)= 'PA'
         chgrsh(6)= 'Cho'
         chgrsh(7)= 'Cr'

         chnot1(1)= 'Cit1b'
         chnot1(2)= 'Cit2b'
         chnot1(3)= 'Cit3b'
         chnot1(4)= 'Cit4b'
         chnot1(5)= 'Cit1c'
         chnot1(6)= 'Cit2c'
         chnot1(7)= 'Cit3c'
         chnot1(8)= 'Cit4c'
         chnot1(9)= 'PA1'
         chnot1(10)= 'PA2'
         chnot1(11)= 'PA3'
         chnot1(12)= 'PA4'
         chnot1(13)= 'Cho2'
         chnot1(14)= 'Cr2'

         chnot2(1)= 'PA0'

      chsimu(1)= 'Cho1 @ 3.2 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(2)= 'Cho2 @ 3.2 +- 0.015 FWHM= 0.06 < 0.07 +- 0.01 AMP=1.0'
      chsimu(3)= 'Cr1 @ 3.03 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(4)= 'Cr2 @ 3.03 +- 0.015 FWHM= 0.06 < 0.07 +- 0.01 AMP=1.0'
      chsimu(5)= 'PA0 @ 3.11 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(6)= 'PA1 @ 3.11 +- 0.015 FWHM= 0.04 < 0.06 +- 0.01 AMP=1.0'
      chsimu(7)= 'PA2 @ 3.11 +- 0.015 FWHM= 0.08 < 0.10 +- 0.01 AMP=1.0'
      chsimu(8)= 'PA3 @ 3.11 +- 0.015 FWHM= 0.12 < 0.14 +- 0.01 AMP=1.0'
      chsimu(9)= 'PA4 @ 3.11 +- 0.015 FWHM= 0.16 < 0.18 +- 0.01 AMP=1.0'
      chsimu(10)= 'Cit1a @ 2.81 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(11)= 'Cit1b @ 2.81 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(12)= 'Cit1c @ 2.81 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(13)= 'Cit2a @ 2.66 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(14)= 'Cit2b @ 2.66 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(15)= 'Cit2c @ 2.66 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(16)= 'Cit3a @ 2.62 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(17)= 'Cit3b @ 2.62 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(18)= 'Cit3c @ 2.62 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(19)= 'Cit4a @ 2.47 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(20)= 'Cit4b @ 2.47 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(21)= 'Cit4c @ 2.47 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
         dorefs(1)= .true.
         dorefs(2)= .false.
         endpha= .true.
         gauss_rt2= .false.
         incsmx= 1
         namrel= 'Cit'
         ncombi= 7
         ngrsh= 7
         nnot1= 14
         nnot2= 1
         nobasi= .true.
         nsidmn= 1
         nsidmx= 1
         nsimul= 21
         onlyco= .true.
         ppmend= 2.1
         ppmst= 4.0
         sddegp= 6.0
         sdgrsh(1)= 0.002
         sdgrsh(2)= 0.002
         sdgrsh(3)= 0.002
         sdgrsh(4)= 0.002
         sdgrsh(5)= 0.004
         sdgrsh(6)= 0.002
         sdgrsh(7)= 0.002
         vitro= .true.
         wconc= .5
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-b'
c     2 peaks (for 1.5T) with weak symmetry constraints.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-b') then
         chcomb(4)= 'Cit1a+Cit1b+Cit1c+Cit2a+Cit2b+Cit2c'
         chgrsh(1)= 'Cit1'  
         chgrsh(2)= 'Cit2'  
         chgrsh(3)= 'PA'  
         chgrsh(4)= 'Cho'  
         chgrsh(5)= 'Cr'  
         chrato(1)= 'Cit1a/Cit2a = 1. +- .1'
         chrato(2)= 'Cit1b/Cit2b = 1. +- .1'
         chrato(3)= 'Cit1c/Cit2c = 1. +- .1'
         chsimu(10)= 'Cit1a @ 2.67 +- .03 FWHM=.02 < .04 +- .01 AMP=1.'
         chsimu(11)= 'Cit1b @ 2.67 +- .03 FWHM=.06 < .08 +- .01 AMP=1.'
         chsimu(12)= 'Cit1c @ 2.67 +- .03 FWHM=.10 < .12 +- .01 AMP=1.'
         chsimu(13)= 'Cit2a @ 2.61 +- .03 FWHM=.02 < .04 +- .01 AMP=1.'
         chsimu(14)= 'Cit2b @ 2.61 +- .03 FWHM=.06 < .08 +- .01 AMP=1.'
         chsimu(15)= 'Cit2c @ 2.61 +- .03 FWHM=.10 < .12 +- .01 AMP=1.'
         idgppm = 2
         ngrsh= 5
         nratio= 3
         nsimul= 15
         ppmend= 2.1
         ppmst= 3.85
         sddegp= 6.
         sdgrsh(1)= .002
         sdgrsh(2)= .002
         sdgrsh(3)= .004
         sdgrsh(4)= .002
         sdgrsh(5)= .002
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-c'
c     4 peaks (for 3T) with exact symmetry with wide DEGPPM range.  
c     For Uppsala.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-c') then
         chcomb(4)= 'Cir1a+Cir2a+Cir3a+Cir4a+Cir5a+Cir6a+Cit7a+Cit8a+Cir
     11b+Cir2b+Cir3b+Cir4b+Cir5b+Cir6b+Cit7b+Cit8b+Cir1c+Cir2c+Cir3c+Cir
     24c+Cir5c+Cir6c+Cit7c+Cit8c'
         chgrsh(1)= 'PA'
         chgrsh(2)= 'Cho'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'Cir'
         chgrsh(5)= 'Cit7'
         chgrsh(6)= 'Cit8'
         chgrsh(7)= 'Cit'
         chgrsh(8)= 'Ci'
         chsimu(10)= 'Cir1a @ 2.830 +- 0.03 FWHM= 0.03 < 0.05+-0.01
     1             AMP=0.5  @ 2.450 FWHM=0.03 AMP=0.5'
         chsimu(11)= 'Cir1b @ 2.830 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.450 FWHM=0.07 AMP=0.5'
         chsimu(12)= 'Cir1c @ 2.830 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.450 FWHM=0.11 AMP=0.5'
         chsimu(13)= 'Cir2a @ 2.815 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.465 FWHM=0.03 AMP=0.5'
         chsimu(14)= 'Cir2b @ 2.815 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.465 FWHM=0.07 AMP=0.5'
         chsimu(15)= 'Cir2c @ 2.815 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1            AMP=0.5  @ 2.465 FWHM=0.11 AMP=0.5'
         chsimu(16)= 'Cir3a @ 2.800 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.480 FWHM=0.03 AMP=0.5'
         chsimu(17)= 'Cir3b @ 2.800 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.480 FWHM=0.07 AMP=0.5'
         chsimu(18)= 'Cir3c @ 2.800 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1            AMP=0.5  @ 2.480 FWHM=0.11 AMP=0.5'
         chsimu(19)= 'Cir4a @ 2.785 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.03 AMP=0.5'
         chsimu(20)= 'Cir4b @ 2.785 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.07 AMP=0.5'
         chsimu(21)= 'Cir4c @ 2.785 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.11 AMP=0.5'
         chsimu(22)= 'Cir5a @ 2.770 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.03 AMP=0.5'
         chsimu(23)= 'Cir5b @ 2.770 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.07 AMP=0.5'
         chsimu(24)= 'Cir5c @ 2.770 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.11 AMP=0.5'
         chsimu(25)= 'Cir6a @ 2.755 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.03 AMP=0.5'
         chsimu(26)= 'Cir6b @ 2.755 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.07 AMP=0.5'
         chsimu(27)= 'Cir6c @ 2.755 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.11 AMP=0.5'
         chsimu(28)= 'Cit7a @ 2.660 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1                AMP=1.0'
         chsimu(29)= 'Cit7b @ 2.660 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1                AMP=1.0'
         chsimu(30)= 'Cit7c @ 2.660 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1                AMP=1.0'
         chsimu(31)= 'Cit8a @ 2.620 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1                AMP=1.0'
         chsimu(32)= 'Cit8b @ 2.620 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1                AMP=1.0'
         chsimu(33)= 'Cit8c @ 2.620 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1                AMP=1.0'
         ddegp3= 20.
         degppm= 90.
         dgppmn= 30.
         dgppmx= 150.
         fh2omx= .8
         idgppm= 1
         ndgppm(2)= 5
         ngrsh= 8
         nsimul= 33
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.125
         rsdgp3= 1.05
         sddegp= 20.
         sdgrsh(1)= 0.004
         sdgrsh(2)= 0.002
         sdgrsh(3)= 0.002
         sdgrsh(4)= 0.0005
         sdgrsh(5)= 0.002
         sdgrsh(6)= 0.002
         sdgrsh(7)= 0.008
         sdgrsh(8)= 0.010
         sptype= 'prostate-a'
         useany= .true.
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-d' 
c     (Probably) superceded by prostate-g
c     4 peaks (strictly 6, since satellites have 2 each) with wide DEGPPM 
c     range.  For Dijon.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-d') then
         badref= .false.
         chcomb(4)= 'Lef1a+Lef2a+Cit3a+Cit4a+Rgt5a+Rgt6a+Lef1b+Lef2b+Cit
     13b+Cit4b+Rgt5b+Rgt6b+Lef1c+Lef2c+Cit3c+Cit4c+Rgt5c+Rgt6c'
         chgrsh(1)= 'PA'
         chgrsh(2)= 'Cho'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'Lef1'
         chgrsh(5)= 'Lef2'
         chgrsh(6)= 'Cit3'
         chgrsh(7)= 'Cit4'
         chgrsh(8)= 'Rgt5'
         chgrsh(9)= 'Rgt6'
         chgrsh(10)= 'Lef'
         chgrsh(11)= 'Cit'
         chgrsh(12)= 'Rgt'
         chnot2(2)= 'Cit0'
         chnot2(3)= 'Cit7'
      chsimu(10)= 'Lef1a @ 2.82 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(11)= 'Lef1b @ 2.82 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(12)= 'Lef1c @ 2.82 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(13)= 'Lef2a @ 2.77 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(14)= 'Lef2b @ 2.77 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(15)= 'Lef2c @ 2.77 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(16)= 'Cit3a @ 2.66 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(17)= 'Cit3b @ 2.66 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(18)= 'Cit3c @ 2.66 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(19)= 'Cit4a @ 2.62 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(20)= 'Cit4b @ 2.62 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(21)= 'Cit4c @ 2.62 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(22)= 'Rgt5a @ 2.51 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(23)= 'Rgt5b @ 2.51 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(24)= 'Rgt5c @ 2.51 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(25)= 'Rgt6a @ 2.46 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(26)= 'Rgt6b @ 2.46 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(27)= 'Rgt6c @ 2.46 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(28)= 'Cit0 @ 2.79 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(29)= 'Cit7 @ 2.49 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
         chuse1(1)= 'Cho1'
         chuse1(2)= 'PA0'
         chuse1(3)= 'Cr1'
         chuse1(4)= 'Cit0'
         chuse1(5)= 'Cit3a'
         chuse1(6)= 'Cit4a'
         chuse1(7)= 'Cit7'
         ddegp3= 20.
         degppm= 90.
         dgppmn= 30.
         dgppmx= 150.
         fh2omx= .8
         idgppm= 1
         ndgppm(2)= 5
         ngrsh= 12
         nnot1= 0
         nnot2= 3
         nsimul= 29
         nuse1= 7
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.125
         rsdgp3= 1.05
         sddegp= 20.
         sdgrsh(1)= 0.004
         sdgrsh(2)= 0.002
         sdgrsh(3)= 0.002
         sdgrsh(4)= 0.002
         sdgrsh(5)= 0.002
         sdgrsh(6)= 0.002
         sdgrsh(7)= 0.002
         sdgrsh(8)= 0.002
         sdgrsh(9)= 0.002
         sdgrsh(10)= 0.008
         sdgrsh(11)= 0.007
         sdgrsh(12)= 0.008
         sptype= 'prostate-a'
         useany= .true.
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-e'
c     Cit & Lip24 as singlets with DGPPM*=0,90.
c     For Dijon.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-e') then
         badref= .false.
         chcom2(8)= 'Lip24'
         chcomb(4)= 'Cit1a+Cit1b+Cit1c+Cit1d'
         chcomb(8)= 'L24a+L24b+L24c+L24d'
         chgrsh(1)= 'Cit1'  
         chgrsh(2)= 'L24'
         chgrsh(3)= 'PA'  
         chgrsh(4)= 'Cho'  
         chgrsh(5)= 'Cr'  
         chsimu(10)= 'Cit1a @ 2.64 +- .03 FWHM=.03 < .05 +- .01 AMP=1.'
         chsimu(11)= 'Cit1b @ 2.64 +- .03 FWHM=.07 < .09 +- .01 AMP=1.'
         chsimu(12)= 'Cit1c @ 2.64 +- .03 FWHM=.11 < .13 +- .01 AMP=1.'
         chsimu(13)= 'Cit1d @ 2.64 +- .03 FWHM=.15 < .17 +- .01 AMP=1.'
         chsimu(14)= 'L24a @ 2.35 +- .03 FWHM=.03 < .05 +- .01 AMP=1.'
         chsimu(15)= 'L24b @ 2.35 +- .03 FWHM=.07 < .09 +- .01 AMP=1.'
         chsimu(16)= 'L24c @ 2.35 +- .03 FWHM=.11 < .13 +- .01 AMP=1.'
         chsimu(17)= 'L24d @ 2.35 +- .03 FWHM=.15 < .17 +- .01 AMP=1.'
         chuse1(1)= 'Cho1'
         chuse1(2)= 'PA0'
         chuse1(3)= 'Cr1'
         chuse1(4)= 'Cit1a'
         chuse1(5)= 'L24a'
         ddegp3= 15.
         degppm= 45.
         dgppmn= 0.
         dgppmx= 90.
         fh2omx= 1.
         idgppm = 2
         ncombi= 8
         ndgppm(2)= 4
         ngrsh= 5
         nsimul= 17
         nuse1= 5
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.
         rsdgp3= 1.05
         sddegp= 20.
         sdgrsh(1)= .004
         sdgrsh(2)= .004
         sdgrsh(3)= .004
         sdgrsh(4)= .002
         sdgrsh(5)= .002
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-f'
c     Supercedes prostate-b.  For 1.5T.
c     All Cho PA Cr, plus symmetric Cit in Prel;
c     Broader Cho3-5, Cit*d-f; starts only at water peak; (self-contained).
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-f') then
         is_sptype= .true.

         atth2o= 1.
         badref= .true.

         chcom2(1)= 'Cho'
         chcom2(2)= 'Cr'
         chcom2(3)= 'PA'
         chcom2(4)= 'Cit'
         chcom2(5)= 'PA+Cr'
         chcom2(6)= 'Cho+Cr'
         chcom2(7)= 'Cho+PA+Cr'

         chcomb(1)= 'Cho1+Cho2+Cho3+Cho4+Cho5'
         chcomb(2)= 'Cr1+Cr2'
         chcomb(3)= 'PA1+PA2+PA3+PA4'
         chcomb(4)= 'Cit1a+Cit1b+Cit1c+Cit1d+Cit1e+Cit1f+Cit2a+Cit2b+Cit
     12c+Cit2d+Cit2e+Cit2f'
         chcomb(5)= 'PA1+PA2+PA3+PA4+Cr1+Cr2'
         chcomb(6)= 'Cho1+Cho2+Cho3+Cho4+Cho5+Cr1+Cr2'
         chcomb(7)= 'Cho1+Cho2+Cho3+Cho4+Cho5+PA1+PA2+PA3+PA4+Cr1+Cr2'

         chgrsh(1)= 'Cit1'
         chgrsh(2)= 'Cit2'
         chgrsh(3)= 'PA'
         chgrsh(4)= 'Cho'
         chgrsh(5)= 'Cr'

         chnot1(1)= 'Cit1a'
         chnot1(2)= 'Cit1b'
         chnot1(3)= 'Cit1c'
         chnot1(4)= 'Cit1d'
         chnot1(5)= 'Cit1e'
         chnot1(6)= 'Cit1f'
         chnot1(7)= 'Cit2a'
         chnot1(8)= 'Cit2b'
         chnot1(9)= 'Cit2c'
         chnot1(10)= 'Cit2d'
         chnot1(11)= 'Cit2e'
         chnot1(12)= 'Cit2f'

         chnot2(1)= 'PA0'
         chnot2(2)= 'Cita'
         chnot2(3)= 'Citb'
         chnot2(4)= 'Citc'
         chnot2(5)= 'Citd'
         chnot2(6)= 'Cite'
         chnot2(7)= 'Citf'

         chrato(1)= 'Cit1a/Cit2a = 1. +- .1'
         chrato(2)= 'Cit1b/Cit2b = 1. +- .1'
         chrato(3)= 'Cit1c/Cit2c = 1. +- .1'
         chrato(4)= 'Cit1d/Cit2d = 1. +- .1'
         chrato(5)= 'Cit1e/Cit2e = 1. +- .1'
         chrato(6)= 'Cit1f/Cit2f = 1. +- .1'
         chsimu(1)= 'Cho1 @ 3.21 +- 0.02 FWHM= 0.03 < 0.04 +- 0.01 
     1               AMP=1.0'
         chsimu(2)= 'Cho2 @ 3.21 +- 0.02 FWHM= 0.06 < 0.07 +- 0.01 
     1               AMP=1.0'
         chsimu(3)= 'Cho3 @ 3.21 +- 0.02 FWHM= 0.09 < 0.10 +- 0.01 
     1               AMP=1.0'
         chsimu(4)= 'Cho4 @ 3.21 +- 0.02 FWHM= 0.12 < 0.14 +- 0.01 
     1               AMP=1.0'
         chsimu(5)= 'Cho5 @ 3.21 +- 0.02 FWHM= 0.16 < 0.18 +- 0.01 
     1               AMP=1.0'

         chsimu(6)= 'Cr1 @ 3.03 +- 0.02 FWHM= 0.03 < 0.04 +- 0.01 
     1               AMP=1.0'
         chsimu(7)= 'Cr2 @ 3.03 +- 0.02 FWHM= 0.06 < 0.07 +- 0.01 
     1               AMP=1.0'
         chsimu(8)= 'PA0 @ 3.11 +- 0.02 FWHM= 0.03 < 0.04 +- 0.01 
     1               AMP=1.0'
         chsimu(9)= 'PA1 @ 3.11 +- 0.02 FWHM= 0.04 < 0.06 +- 0.01 
     1               AMP=1.0'
         chsimu(10)= 'PA2 @ 3.11 +- 0.02 FWHM= 0.08 < 0.10 +- 0.01 
     1                AMP=1.0'
         chsimu(11)= 'PA3 @ 3.11 +- 0.02 FWHM= 0.12 < 0.14 +- 0.01 
     1                AMP=1.0'
         chsimu(12)= 'PA4 @ 3.11 +- 0.02 FWHM= 0.16 < 0.18 +- 0.01 
     1                AMP=1.0'
         chsimu(13)= 'Cit1a @ 2.67 +- .02 FWHM=.02 < .04 +- .01 AMP=1.'
         chsimu(14)= 'Cit1b @ 2.67 +- .02 FWHM=.06 < .08 +- .01 AMP=1.'
         chsimu(15)= 'Cit1c @ 2.67 +- .02 FWHM=.10 < .12 +- .01 AMP=1.'
         chsimu(16)= 'Cit1d @ 2.67 +- .02 FWHM=.14 < .16 +- .01 AMP=1.'
         chsimu(17)= 'Cit1e @ 2.67 +- .02 FWHM=.18 < .20 +- .01 AMP=1.'
         chsimu(18)= 'Cit1f @ 2.67 +- .02 FWHM=.22 < .24 +- .01 AMP=1.'

         chsimu(19)= 'Cit2a @ 2.61 +- .02 FWHM=.02 < .04 +- .01 AMP=1.'
         chsimu(20)= 'Cit2b @ 2.61 +- .02 FWHM=.06 < .08 +- .01 AMP=1.'
         chsimu(21)= 'Cit2c @ 2.61 +- .02 FWHM=.10 < .12 +- .01 AMP=1.'
         chsimu(22)= 'Cit2d @ 2.61 +- .02 FWHM=.14 < .16 +- .01 AMP=1.'
         chsimu(23)= 'Cit2e @ 2.61 +- .02 FWHM=.18 < .20 +- .01 AMP=1.'
         chsimu(24)= 'Cit2f @ 2.61 +- .02 FWHM=.22 < .24 +- .01 AMP=1.'

         chsimu(25)= 'Cita@2.67+-.02 FWHM=.02<.04+-.01 AMP=.5 
     1                    @2.61 FWHM=.02 AMP=.5'
         chsimu(26)= 'Citb@2.67+-.02 FWHM=.06<.08+-.01 AMP=.5 
     1                    @2.61 FWHM=.06 AMP=.5'
         chsimu(27)= 'Citc@2.67+-.02 FWHM=.10<.12+-.01 AMP=.5 
     1                    @2.61 FWHM=.10 AMP=.5'
         chsimu(28)= 'Citd@2.67+-.02 FWHM=.14<.16+-.01 AMP=.5 
     1                   @2.61 FWHM=.14 AMP=.5'
         chsimu(29)= 'Cite@2.67+-.02 FWHM=.18<.20+-.01 AMP=.5 
     1                    @2.61 FWHM=.18 AMP=.5'
         chsimu(30)= 'Citf@2.67+-.02 FWHM=.22<.24+-.01 AMP=.5 
     1                    @2.61 FWHM=.22 AMP=.5'

         dorefs(1)= .true.
         dorefs(2)= .false.
         endpha= .true.
         fh2omx= 1.
         fwhmst= .1
         gauss_rt2= .false.
         idgppm= 2
         incsmx= 1
         namrel= 'Cit'
         ncombi= 7
         ngrsh= 5
         nnot1= 12
         nnot2= 7
         nratio= 6
         nobasi= .true.
         nsidmn= 1
         nsidmx= 1
         nsimul= 30
         onlyco= .true.
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.
         sddegp= 6.
         sdgrsh(1)= .005
         sdgrsh(2)= .005
         sdgrsh(3)= .010
         sdgrsh(4)= .005
         sdgrsh(5)= .005
         vitro= .true.
         wconc= .5
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-g' (self-contained)
c     (Probably) supercedes prostate-d.
c     For Mannheim 3T, which has highly variable Cit structures.
c     Prel uses single sharpest Cho1 PA0 Cr1 Cit3a(@2.66) Cit4a(@2.62).
c     Final has full variety of Cit (as in prostate-d) to empirically fit 
c        highly variable structure.
c     3 components for Cho1 PA0 Cr1 (up to FWHM=.11) (i.e., 1 less PA and 1
c        more Cho & Cr than prostate-d).
c     By default, starts only at water peak (but Mannheim has no water peak).
c     Missing ATTH2O=1 (error) added 100820.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-g') then
         is_sptype= .true.

         atth2o= 1.

         chcom2(1)= 'Cho'
         chcom2(2)= 'Cr'
         chcom2(3)= 'PA'
         chcom2(4)= 'Cit'
         chcom2(5)= 'PA+Cr'
         chcom2(6)= 'Cho+Cr'
         chcom2(7)= 'Cho+PA+Cr'

         chcomb(1)= 'Cho1+Cho2+Cho3'
         chcomb(2)= 'Cr1+Cr2+Cr3'
         chcomb(3)= 'PA1+PA2+PA3'
         chcomb(4)= 'Lef1a+Lef2a+Cit3a+Cit4a+Rgt5a+Rgt6a+Lef1b+Lef2b+Cit
     13b+Cit4b+Rgt5b+Rgt6b+Lef1c+Lef2c+Cit3c+Cit4c+Rgt5c+Rgt6c'
         chcomb(5)= 'PA1+PA2+PA3+Cr1+Cr2+Cr3'
         chcomb(6)= 'Cho1+Cho2+Cho3+Cr1+Cr2+Cr3'
         chcomb(7)= 'Cho1+Cho2+Cho3+PA1+PA2+PA3+Cr1+Cr2+Cr3'

         chgrsh(1)= 'PA'
         chgrsh(2)= 'Cho'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'Lef1'
         chgrsh(5)= 'Lef2'
         chgrsh(6)= 'Cit3'
         chgrsh(7)= 'Cit4'
         chgrsh(8)= 'Rgt5'
         chgrsh(9)= 'Rgt6'
         chgrsh(10)= 'Lef'
         chgrsh(11)= 'Cit'
         chgrsh(12)= 'Rgt'

         chnot2(1)= 'PA0'
      chsimu(1)= 'Cho1 @ 3.2 +- 0.015 FWHM= 0.03 < 0.04+-0.01 AMP=1.0'
      chsimu(2)= 'Cho2 @ 3.2 +- 0.015 FWHM= 0.06 < 0.07+-0.01 AMP=1.0'
      chsimu(3)= 'Cho3 @ 3.2 +- 0.015 FWHM= 0.09 < 0.11+-0.01 AMP=1.0'
      chsimu(4)= 'Cr1 @ 3.03 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(5)= 'Cr2 @ 3.03 +- 0.015 FWHM= 0.06 < 0.07 +- 0.01 AMP=1.0'
      chsimu(6)= 'Cr3 @ 3.03 +- 0.015 FWHM= 0.09 < 0.11 +- 0.01 AMP=1.0'
      chsimu(7)= 'PA0 @ 3.11 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(8)= 'PA1 @ 3.11 +- 0.015 FWHM= 0.04 < 0.06 +- 0.01 AMP=1.0'
      chsimu(9)= 'PA2 @ 3.11 +- 0.015 FWHM= 0.08 < 0.10 +- 0.01 AMP=1.0'
      chsimu(10)= 'PA3 @ 3.11 +- 0.015 FWHM= 0.12 < 0.14+-0.01 AMP=1.0'
      chsimu(11)= 'Lef1a @ 2.82 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(12)= 'Lef1b @ 2.82 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(13)= 'Lef1c @ 2.82 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(14)= 'Lef2a @ 2.77 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(15)= 'Lef2b @ 2.77 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(16)= 'Lef2c @ 2.77 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(17)= 'Cit3a @ 2.66 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(18)= 'Cit3b @ 2.66 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(19)= 'Cit3c @ 2.66 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(20)= 'Cit4a @ 2.62 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(21)= 'Cit4b @ 2.62 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(22)= 'Cit4c @ 2.62 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(23)= 'Rgt5a @ 2.51 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(24)= 'Rgt5b @ 2.51 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(25)= 'Rgt5c @ 2.51 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'
      chsimu(26)= 'Rgt6a @ 2.46 +- 0.03 FWHM= 0.03 < 0.05+-0.01 AMP=1.0'
      chsimu(27)= 'Rgt6b @ 2.46 +- 0.03 FWHM= 0.07 < 0.09+-0.01 AMP=1.0'
      chsimu(28)= 'Rgt6c @ 2.46 +- 0.03 FWHM= 0.11 < 0.13+-0.01 AMP=1.0'

         chuse1(1)= 'Cho1'
         chuse1(2)= 'PA0'
         chuse1(3)= 'Cr1'
         chuse1(4)= 'Cit3a'
         chuse1(5)= 'Cit4a'

         ddegp3= 10.
         dorefs(1)= .true.
         dorefs(2)= .false.
         fh2omx= .8
         gauss_rt2= .false.
         idgppm= 2
         incsmx= 1
         namrel= 'Cit'
         ncombi= 7
         ndgppm(2)= 5
         ngrsh= 12
         nnot2= 1
         nobasi= .true.
         nsidmn= 1
         nsidmx= 1
         nsimul= 28
         nuse1= 5
         onlyco= .true.
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.125
         rsdgp3= 1.05
         sdgrsh(1)= 0.008
         sdgrsh(2)= 0.004
         sdgrsh(3)= 0.004
         sdgrsh(4)= 0.004
         sdgrsh(5)= 0.004
         sdgrsh(6)= 0.004
         sdgrsh(7)= 0.004
         sdgrsh(8)= 0.004
         sdgrsh(9)= 0.004
         sdgrsh(10)= 0.010
         sdgrsh(11)= 0.008
         sdgrsh(12)= 0.010
         useany= .true.
         vitro= .true.
         wconc= .5
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     sptype = 'prostate-h' (self-contained); start of prostate-i
c     Symmetric Cit satellites and wide DEGPPM range [30, 150] (of
c        prostate-c). 
c     More balanced simulations of Cho, PA & Cr (of prostate-g) (to avoid PA 
c        overestimation); rest mostly prostate-g.
c     For Leipzig 3T of 1008, which has fairly uniform and symmetric Cit (but 
c        only of fair quality, with large negative side peaks).  Leipzig has 
c        no useable water peaks; so, PPMSHF=0 was added externally.
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-h'   .or.
     1    sptype(:10) .eq. 'prostate-i') then
         is_sptype= .true.
         atth2o= 1.

         chcom2(1)= 'Cho'
         chcom2(2)= 'Cr'
         chcom2(3)= 'PA'
         chcom2(4)= 'Cit'
         chcom2(5)= 'PA+Cr'
         chcom2(6)= 'Cho+Cr'
         chcom2(7)= 'Cho+PA+Cr'

         chcomb(1)= 'Cho1+Cho2+Cho3'
         chcomb(2)= 'Cr1+Cr2+Cr3'
         chcomb(3)= 'PA1+PA2+PA3'
         chcomb(4)= 'Cir1a+Cir2a+Cir3a+Cir4a+Cir5a+Cir6a+Cit7a+Cit8a+Cir
     11b+Cir2b+Cir3b+Cir4b+Cir5b+Cir6b+Cit7b+Cit8b+Cir1c+Cir2c+Cir3c+Cir
     24c+Cir5c+Cir6c+Cit7c+Cit8c'
         chcomb(5)= 'PA1+PA2+PA3+Cr1+Cr2+Cr3'
         chcomb(6)= 'Cho1+Cho2+Cho3+Cr1+Cr2+Cr3'
         chcomb(7)= 'Cho1+Cho2+Cho3+PA1+PA2+PA3+Cr1+Cr2+Cr3'

         chgrsh(1)= 'PA'
         chgrsh(2)= 'Cho'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'Cir'
         chgrsh(5)= 'Cit7'
         chgrsh(6)= 'Cit8'
         chgrsh(7)= 'Cit'
         chgrsh(8)= 'Ci'

         chnot2(1)= 'PA0'
         chnot2(2)= 'Cit0'

         chrato(1)= 'Cit7a/Cit8a = 1. +- .1'
         chrato(2)= 'Cit7b/Cit8b = 1. +- .1'
         chrato(3)= 'Cit7c/Cit8c = 1. +- .1'

      chsimu(1)= 'Cho1 @ 3.2 +- 0.015 FWHM= 0.03 < 0.04+-0.01 AMP=1.0'
      chsimu(2)= 'Cho2 @ 3.2 +- 0.015 FWHM= 0.06 < 0.07+-0.01 AMP=1.0'
      chsimu(3)= 'Cho3 @ 3.2 +- 0.015 FWHM= 0.09 < 0.11+-0.01 AMP=1.0'
      chsimu(4)= 'Cr1 @ 3.03 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(5)= 'Cr2 @ 3.03 +- 0.015 FWHM= 0.06 < 0.07 +- 0.01 AMP=1.0'
      chsimu(6)= 'Cr3 @ 3.03 +- 0.015 FWHM= 0.09 < 0.11 +- 0.01 AMP=1.0'
      chsimu(7)= 'PA0 @ 3.11 +- 0.015 FWHM= 0.03 < 0.04 +- 0.01 AMP=1.0'
      chsimu(8)= 'PA1 @ 3.11 +- 0.015 FWHM= 0.04 < 0.06 +- 0.01 AMP=1.0'
      chsimu(9)= 'PA2 @ 3.11 +- 0.015 FWHM= 0.08 < 0.10 +- 0.01 AMP=1.0'
      chsimu(10)= 'PA3 @ 3.11 +- 0.015 FWHM= 0.12 < 0.14+-0.01 AMP=1.0'
         chsimu(11)= 'Cir1a @ 2.830 +- 0.03 FWHM= 0.03 < 0.05+-0.01
     1             AMP=0.5  @ 2.450 FWHM=0.03 AMP=0.5'
         chsimu(12)= 'Cir1b @ 2.830 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.450 FWHM=0.07 AMP=0.5'
         chsimu(13)= 'Cir1c @ 2.830 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.450 FWHM=0.11 AMP=0.5'
         chsimu(14)= 'Cir2a @ 2.815 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.465 FWHM=0.03 AMP=0.5'
         chsimu(15)= 'Cir2b @ 2.815 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.465 FWHM=0.07 AMP=0.5'
         chsimu(16)= 'Cir2c @ 2.815 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1            AMP=0.5  @ 2.465 FWHM=0.11 AMP=0.5'
         chsimu(17)= 'Cir3a @ 2.800 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.480 FWHM=0.03 AMP=0.5'
         chsimu(18)= 'Cir3b @ 2.800 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.480 FWHM=0.07 AMP=0.5'
         chsimu(19)= 'Cir3c @ 2.800 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1            AMP=0.5  @ 2.480 FWHM=0.11 AMP=0.5'
         chsimu(20)= 'Cir4a @ 2.785 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.03 AMP=0.5'
         chsimu(21)= 'Cir4b @ 2.785 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.07 AMP=0.5'
         chsimu(22)= 'Cir4c @ 2.785 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.495 FWHM=0.11 AMP=0.5'
         chsimu(23)= 'Cir5a @ 2.770 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.03 AMP=0.5'
         chsimu(24)= 'Cir5b @ 2.770 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.07 AMP=0.5'
         chsimu(25)= 'Cir5c @ 2.770 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.510 FWHM=0.11 AMP=0.5'
         chsimu(26)= 'Cir6a @ 2.755 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.03 AMP=0.5'
         chsimu(27)= 'Cir6b @ 2.755 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.07 AMP=0.5'
         chsimu(28)= 'Cir6c @ 2.755 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1             AMP=0.5  @ 2.525 FWHM=0.11 AMP=0.5'
         chsimu(29)= 'Cit7a @ 2.660 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1                AMP=1.0'
         chsimu(30)= 'Cit7b @ 2.660 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1                AMP=1.0'
         chsimu(31)= 'Cit7c @ 2.660 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1                AMP=1.0'
         chsimu(32)= 'Cit8a @ 2.620 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1                AMP=1.0'
         chsimu(33)= 'Cit8b @ 2.620 +- 0.03 FWHM= 0.07 < 0.09+-0.01 
     1                AMP=1.0'
         chsimu(34)= 'Cit8c @ 2.620 +- 0.03 FWHM= 0.11 < 0.13+-0.01 
     1                AMP=1.0'
         chsimu(35)= 'Cit0 @ 2.640 +- 0.03 FWHM= 0.03 < 0.05+-0.01 
     1                AMP=1.0'
         chuse1(1)= 'Cho1'
         chuse1(2)= 'PA0'
         chuse1(3)= 'Cr1'
         chuse1(4)= 'Cit0'

         ddegp3= 20.
         degppm= 90.
         dgppmn= 30.
         dgppmx= 150.
         dorefs(1)= .true.
         dorefs(2)= .false.
         fh2omx= .8
         gauss_rt2= .false.
         idgppm= 1
         incsmx= 1
         namrel= 'Cit'
         ncombi= 7
         ndgppm(2)= 5
         ngrsh= 8
         nnot2= 2
         nobasi= .true.
         nratio= 3
         nsidmn= 1
         nsidmx= 1
         nsimul= 35
         nuse1= 4
         onlyco= .true.
         ppmend= 2.1
         ppmst= 3.85
         rfwhcc= 0.125
         rsdgp3= 1.05
         sdgrsh(1)= 0.006
         sdgrsh(2)= 0.003
         sdgrsh(3)= 0.003
         sdgrsh(4)= 0.00075
         sdgrsh(5)= 0.003
         sdgrsh(6)= 0.003
         sdgrsh(7)= 0.012
         sdgrsh(8)= 0.015
         useany= .true.
         vitro= .true.
         wconc= .5
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     SPTYPE = 'prostate-i'
c     Modify prostate-h to be more general-purpose, with better phase range & 
c        baseline.
c     Sent to Karolinska (1012)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:10) .eq. 'prostate-i') then
         ddegp3= 10.
         degmax(1)= 5.
         degppm= 0.
         dgppmn= -60.
         dgppmx= 60.
         isdbol= 5
         mrepha(1)= 1
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     SPTYPE = 'mega-press-1' or 'mega=-press-2' or 'mega=-press-3'
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:12) .eq. 'mega-press-1'   .or.   
     1    sptype(:12) .eq. 'mega-press-2'   .or.
     1    sptype(:12) .eq. 'mega-press-3') then
         is_sptype= .true.
         alsdsh(1)= .004
         alsdsh(2)= .008
         atth2o = .4
         chrato(1)= 'NAAG/NAA = .15 +- .15'
         chuse1(1)= 'NAA'
         chuse1(2)= 'Glu'
         isdbol= 5
         namrel= 'NAA+NAAG'
         nobase= .true.
         nratio= 1
         nrefpk(2)= 1
         nsimul= 11
         nuse1= 2
         ppmend= 1.9
         ppmst= 4.2
         ratipm= 10.
         sddegp= 1.
         wsmet= 'NAA'
         wsppm= 2.01
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     SPTYPE = 'mega-press-2' -3' (minor changes)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:12) .eq. 'mega-press-2'   .or.
     1    sptype(:12) .eq. 'mega-press-3') then
         alsdsh(2)= .01
         atth2o= .43
         chcomb(17)= 'GSH+Glu+Gln'
         chsimu(1)= 'MM09 @ .91 +- .02 FWHM=.12 < .15 +- .015 AMP=3.'
         dgppmn= -15.
         dgppmx= 15.
         ncombi= 17
         nsimul= 1
         sddegp= 4.
         shifmn(2)= -1.3
         shifmx(2)= 1.3
      end if
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     SPTYPE = 'mega-press-3' more freedom for GABA; 
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (sptype(:12) .eq. 'mega-press-3') then
         chsdsh(3)= 'GABA'
         ppmend= 1.95
      end if
