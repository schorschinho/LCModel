c     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
c     Define liver-1 and set initial parts for breast-* & lipid-* (which will 
c        be modified in lipid.inc).
c     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      if (sptype(:8) .eq. 'liver-1 '   .or.
     1    sptype(:8) .eq. 'lipid-1 '   .or.
     2    sptype(:9) .eq. 'breast-1 '  .or.
     3    sptype(:7) .eq. 'liver-2'   .or.
     4    sptype(:7) .eq. 'lipid-2'   .or.
     5    sptype(:8) .eq. 'breast-2'   .or.
     6    sptype(:10) .eq. 'only-cho-1'   .or.
     7    sptype(:7) .eq. 'lipid-3'   .or.
     8    sptype(:7) .eq. 'liver-3'   .or.
     9    sptype(:8) .eq. 'breast-3'   .or.
     1    sptype(:10) .eq. 'only-cho-2') then
         is_sptype= .true.

         if (sptype(:8) .eq. 'liver-1 '   .or.
     1       sptype(:8) .eq. 'lipid-1 '   .or.
     2       sptype(:9) .eq. 'breast-1 ') then
            if (ppmst .lt. -9998.) ppmst= 4.2
            if (ppmend .gt. 9998.) ppmend= -2.
            if (ppmst .ge. 6.) then
               namrel= 'Water'
               nrefpk(2)= 2
               ppmref(2,2)= 4.65
            else
               namrel= 'Lip09+Lip13'
               nrefpk(2)= 1
            end if
            if (ppmst - ppmend .gt. 6.3) then
               nsubtk= 5
               xstep= 1.
            end if
         end if

         atth2o= 1.
         badref= .true.
         chcomb(1)= 'L20a+L20b+L20c+L20d+L20e+L20f+L13a+L13b+L13c+L13d+L
     113e+L13f+L13g+L13h+L13i+L13j+L09a+L09b+L09c+L09d+L09e+L09f+L09g'
         chcomb(2)= 'L13j+L13i+L13h+L09g+L13g+L09f+L13f+L09e+L13e+L09d+L
     113d+L09c+L13c+L09b+L13b+L09a+L13a'
         chcomb(3)= 'L13j+L13i+L13h+L13g+L13f+L13e+L13d+L13c+L13b+L13a'
         chcomb(4)= 'L09g+L09f+L09e+L09d+L09c+L09b+L09a'
         chcomb(5)= 'L20f+L20e+L20d+L20c+L20b+L20a'
         chcomb(6)= 'L53f+L53e+L53d+L53c+L53b+L53a'
         chcomb(7)= 'W11+W10+W9+W8+L53f+W7+L53e+W6+L53d+W5+L53c+W4+L53b+
     1W3+L53a+W2+W1'
         chcomb(8)= 'W11+W10+W9+W8+W7+W6+W5+W4+W3+W2+W1'
         chcomb(9)= 'cho1+cho2+cho3+cho4+cho5+cho6'
         chcomb(10)= 'glycg1+glycg2+glycg3+glycg4+glycg5+glycg6+glycg7'

         chcom2(1)= 'L20+L13+L09'
         chcom2(2)= 'Lip09+Lip13'
         chcom2(3)= 'Lip13'
         chcom2(4)= 'Lip09'
         chcom2(5)= 'Lip20'
         chcom2(6)= 'Lip53'
         chcom2(7)= 'Lip53+Water'
         chcom2(8)= 'Water'
         chcom2(9)= 'Cho'
         chcom2(10)= 'Glycg'

         chgrsh(1)= 'L53'
         chgrsh(2)= 'L20'
         chgrsh(3)= 'L13'
         chgrsh(4)= 'L09'
         chgrsh(5)= 'cho'

         chnot1(1)= 'L53a'
         chnot1(2)= 'L53b'
         chnot1(3)= 'L53c'
         chnot1(4)= 'L53d'
         chnot1(5)= 'L53e'
         chnot1(6)= 'L53f'

         chnot1(7)= 'L20a'
         chnot1(8)= 'L20b'
         chnot1(9)= 'L20c'
         chnot1(10)= 'L20d'
         chnot1(11)= 'L20e'
         chnot1(12)= 'L20f'

         chnot1(13)= 'L13d'
         chnot1(14)= 'L13e'
         chnot1(15)= 'L13f'
         chnot1(16)= 'L13g'
         chnot1(17)= 'L13h'
         chnot1(18)= 'L13i'
         chnot1(19)= 'L13j'

         chnot1(20)= 'L09a'
         chnot1(21)= 'L09b'
         chnot1(22)= 'L09c'
         chnot1(23)= 'L09d'
         chnot1(24)= 'L09e'
         chnot1(25)= 'L09f'
         chnot1(26)= 'L09g'

         chnot1(27)= 'cho1'
         chnot1(28)= 'cho2'
         chnot1(29)= 'cho3'
         chnot1(30)= 'cho4'
         chnot1(31)= 'cho5'
         chnot1(32)= 'cho6'

         chnot1(33)= 'glycg1'
         chnot1(34)= 'glycg2'
         chnot1(35)= 'glycg3'
         chnot1(36)= 'glycg4'
         chnot1(37)= 'glycg5'
         chnot1(38)= 'glycg6'
         chnot1(39)= 'glycg7'

         chsimu(1)= 'L53a @ 5.33+-.02 FWHM=.2<.23+-.03 AMP=1.'
         chsimu(2)= 'L53b @ 5.33+-.02 FWHM=.3<.33+-.03 AMP=1.'
         chsimu(3)= 'L53c @ 5.33+-.02 FWHM=.4<.43+-.03 AMP=1.'
         chsimu(4)= 'L53d @ 5.33+-.02 FWHM=.5<.53+-.03 AMP=1.'
         chsimu(5)= 'L53e @ 5.33+-.02 FWHM=.6<.63+-.03 AMP=1.'
         chsimu(6)= 'L53f @ 5.33+-.02 FWHM=.7<.73+-.03 AMP=1.'

         chsimu(7)= 'W1 @ 4.65+-.05 FWHM=.01<.04+-.03 AMP=1.'
         chsimu(8)= 'W2 @ 4.65+-.05 FWHM=.05<.08+-.03 AMP=1.'
         chsimu(9)= 'W3 @ 4.65+-.05 FWHM=.15<.18+-.03 AMP=1.'
         chsimu(10)= 'W4 @ 4.65+-.05 FWHM=.25<.28+-.03 AMP=1.'
         chsimu(11)= 'W5 @ 4.65+-.05 FWHM=.35<.38+-.03 AMP=1.'
         chsimu(12)= 'W6 @ 4.65+-.05 FWHM=.45<.48+-.03 AMP=1.'
         chsimu(13)= 'W7 @ 4.65+-.05 FWHM=.55<.58+-.03 AMP=1.'
         chsimu(14)= 'W8 @ 4.65+-.05 FWHM=.65<.68+-.03 AMP=1.'
         chsimu(15)= 'W9 @ 4.65+-.05 FWHM=.75<.78+-.03 AMP=1.'
         chsimu(16)= 'W10 @ 4.65+-.05 FWHM=.85<.88+-.03 AMP=1.'
         chsimu(17)= 'W11 @ 4.65+-.05 FWHM=.95<.98+-.03 AMP=1.'

         chsimu(18)= 'L20a @ 2.06+-.02 FWHM=.1<.13+-.03 AMP=.58 
     1                @ 2.25 FWHM=.1 AMP=.29  @ 2.7 FWHM=.15 AMP=.13'
         chsimu(19)= 'L20b @ 2.06+-.02 FWHM=.2<.23+-.03 AMP=.58 
     1                @ 2.25 FWHM=.2 AMP=.29  @ 2.7 FWHM=.25 AMP=.13'
         chsimu(20)= 'L20c @ 2.06+-.02 FWHM=.3<.33+-.03 AMP=.58  
     1                @ 2.25 FWHM=.3 AMP=.29  @ 2.7 FWHM=.35 AMP=.13'
         chsimu(21)= 'L20d @ 2.06+-.02 FWHM=.4<.43+-.03 AMP=.58  
     1                @ 2.25 FWHM=.4 AMP=.29  @ 2.7 FWHM=.45 AMP=.13'
         chsimu(22)= 'L20e @ 2.06+-.02 FWHM=.5<.53+-.03 AMP=.58  
     1                @ 2.25 FWHM=.5 AMP=.29  @ 2.7 FWHM=.55 AMP=.13'
         chsimu(23)= 'L20f @ 2.06+-.02 FWHM=.6<.63+-.03 AMP=.58  
     1                @ 2.25 FWHM=.6 AMP=.29  @ 2.7 FWHM=.65 AMP=.13'

         chsimu(24)= 'L13a @ 1.28+-.02 FWHM=.05<.08+-.03 AMP=1.'
         chsimu(25)= 'L13b @ 1.28+-.02 FWHM=.15<.18+-.03 AMP=1.'
         chsimu(26)= 'L13c @ 1.28+-.02 FWHM=.25<.28+-.03 AMP=1.'
         chsimu(27)= 'L13d @ 1.28+-.02 FWHM=.35<.38+-.03 AMP=1.'
         chsimu(28)= 'L13e @ 1.28+-.02 FWHM=.45<.48+-.03 AMP=1.'
         chsimu(29)= 'L13f @ 1.28+-.02 FWHM=.55<.58+-.03 AMP=1.'
         chsimu(30)= 'L13g @ 1.28+-.02 FWHM=.65<.68+-.03 AMP=1.'
         chsimu(31)= 'L13h @ 1.28+-.02 FWHM=.75<.78+-.03 AMP=1.'
         chsimu(32)= 'L13i @ 1.28+-.02 FWHM=.85<.88+-.03 AMP=1.'
         chsimu(33)= 'L13j @ 1.28+-.02 FWHM=.95<.98+-.03 AMP=1.'

         chsimu(34)= 'L09a @ .87+-.02 FWHM=.15<.18+-.03 AMP=1.'
         chsimu(35)= 'L09b @ .87+-.02 FWHM=.25<.28+-.03 AMP=1.'
         chsimu(36)= 'L09c @ .87+-.02 FWHM=.35<.38+-.03 AMP=1.'
         chsimu(37)= 'L09d @ .87+-.02 FWHM=.45<.48+-.03 AMP=1.'
         chsimu(38)= 'L09e @ .87+-.02 FWHM=.55<.58+-.03 AMP=1.'
         chsimu(39)= 'L09f @ .87+-.02 FWHM=.65<.68+-.03 AMP=1.'
         chsimu(40)= 'L09g @ .87+-.02 FWHM=.75<.78+-.03 AMP=1.'

         chsimu(41)= 'cho1 @ 3.18+-.01 FWHM=.06<.08+-.01 AMP=1.'
         chsimu(42)= 'cho2 @ 3.18+-.01 FWHM=.08<.10+-.01 AMP=1.'
         chsimu(43)= 'cho3 @ 3.18+-.01 FWHM=.10<.12+-.01 AMP=1.'
         chsimu(44)= 'cho4 @ 3.18+-.01 FWHM=.12<.14+-.01 AMP=1.'
         chsimu(45)= 'cho5 @ 3.18+-.01 FWHM=.14<.16+-.01 AMP=1.'
         chsimu(46)= 'cho6 @ 3.18+-.01 FWHM=.16<.18+-.01 AMP=1.'

         chsimu(47)= 'glycg1 @ 3.82+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(48)= 'glycg2 @ 3.75+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(49)= 'glycg3 @ 3.68+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(50)= 'glycg4 @ 3.61+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(51)= 'glycg5 @ 3.54+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(52)= 'glycg6 @ 3.47+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(53)= 'glycg7 @ 3.40+-.01 FWHM=.09<.11+-.01 AMP=1.'

         dkntmn(2)= 2.4
         endpha= .true.
         fwhmmx= .4
         fwhmst= .3
         gauss_rt2= .false.
         incsmx= 1
         ncombi= 10
         ngrsh= 5
         nnot1= 39
         nobasi= .true.
         nratio= 0
         nsidmn= 1
         nsidmx= 1
         nsimul= 53
         nuse1= 0
         onlyco= .true.
         ppmref(1,2)= 1.28
         sdgrsh(1)= .003
         sdgrsh(2)= .004
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .002
         shifmn(2)= -1.
         shifmx(2)= 1.
         vitro= .true.
         wconc= .5
      end if
c     222222222222222222222222222222222222222222222222222222222222222222222222
c     Modify *-1 settings for SPTYPE = '*-2' & 'only-cho-1, which have (only) 
c        L13j L13i glycg7 glycg6 & glycg5 omitted
c     222222222222222222222222222222222222222222222222222222222222222222222222
      if (sptype(:7) .eq. 'liver-2'   .or.
     1    sptype(:7) .eq. 'lipid-2'   .or.
     2    sptype(:8) .eq. 'breast-2'   .or. 
     3    sptype(:10) .eq. 'only-cho-1'   .or.
     4    sptype(:7) .eq. 'lipid-3'   .or.
     5    sptype(:7) .eq. 'liver-3'   .or.
     6    sptype(:8) .eq. 'breast-3'   .or. 
     7    sptype(:10) .eq. 'only-cho-2') then

         chcomb(1)= 'L20a+L20b+L20c+L20d+L20e+L20f+L13a+L13b+L13c+L13d+L
     113e+L13f+L13g+L13h+L09a+L09b+L09c+L09d+L09e+L09f'
         chcomb(2)= 'L13h+L13g+L09f+L13f+L09e+L13e+L09d+L13d+L09c+L13c+L
     109b+L13b+L09a+L13a'
         chcomb(3)= 'L13h+L13g+L13f+L13e+L13d+L13c+L13b+L13a'
         chcomb(4)= 'L09f+L09e+L09d+L09c+L09b+L09a'
         chcomb(9)= 'cho1+cho2+cho3+cho4+cho5+cho6+cho7'
         chcomb(10)= 'glycg1+glycg2+glycg3+glycg4'

         do 7210 j = 13, 17
            chnot1(j)= 'xqlytc'
7210     continue
         chnot1(27)='cho7'
         chnot1(28)='xqlytc'
         chnot1(29)='xqlytc'

         chsimu(18)= 'L20a @ 2.02+-.03 FWHM=.1<.13+-.03 AMP=.58 
     1                @ 2.23 FWHM=.1 AMP=.29  @ 2.76 FWHM=.15 AMP=.13'
         chsimu(19)= 'L20b @ 2.02+-.03 FWHM=.2<.23+-.03 AMP=.58 
     1                @ 2.23 FWHM=.2 AMP=.29  @ 2.76 FWHM=.25 AMP=.13'
         chsimu(20)= 'L20c @ 2.02+-.03 FWHM=.3<.33+-.03 AMP=.58  
     1                @ 2.23 FWHM=.3 AMP=.29  @ 2.76 FWHM=.35 AMP=.13'
         chsimu(21)= 'L20d @ 2.02+-.03 FWHM=.4<.43+-.03 AMP=.58  
     1                @ 2.23 FWHM=.4 AMP=.29  @ 2.76 FWHM=.45 AMP=.13'
         chsimu(22)= 'L20e @ 2.02+-.03 FWHM=.5<.53+-.03 AMP=.58  
     1                @ 2.23 FWHM=.5 AMP=.29  @ 2.76 FWHM=.55 AMP=.13'
         chsimu(23)= 'L20f @ 2.02+-.03 FWHM=.6<.63+-.03 AMP=.58  
     1                @ 2.23 FWHM=.6 AMP=.29  @ 2.76 FWHM=.65 AMP=.13'

         chsimu(24)= 'L13a @ 1.28+-.02 FWHM=.04<.08+-.03 AMP=1.'

         chsimu(32)= 'L09a @ .87+-.02 FWHM=.15<.18+-.03 AMP=1.'
         chsimu(33)= 'L09b @ .87+-.02 FWHM=.25<.28+-.03 AMP=1.'
         chsimu(34)= 'L09c @ .87+-.02 FWHM=.35<.38+-.03 AMP=1.'
         chsimu(35)= 'L09d @ .87+-.02 FWHM=.45<.48+-.03 AMP=1.'
         chsimu(36)= 'L09e @ .87+-.02 FWHM=.55<.58+-.03 AMP=1.'
         chsimu(37)= 'L09f @ .87+-.02 FWHM=.65<.68+-.03 AMP=1.'

         chsimu(38)= 'cho1 @ 3.185+-.01 FWHM=.030<.050+-.012 AMP=1.'
         chsimu(39)= 'cho2 @ 3.185+-.01 FWHM=.074<.086+-.012 AMP=1.'
         chsimu(40)= 'cho3 @ 3.185+-.01 FWHM=.110<.122+-.012 AMP=1.'
         chsimu(41)= 'cho4 @ 3.185+-.01 FWHM=.146<.158+-.012 AMP=1.'
         chsimu(42)= 'cho5 @ 3.185+-.01 FWHM=.182<.194+-.012 AMP=1.'
         chsimu(43)= 'cho6 @ 3.185+-.01 FWHM=.218<.230+-.012 AMP=1.'
         chsimu(44)= 'cho7 @ 3.185+-.01 FWHM=.254<.266+-.012 AMP=1.'

         chsimu(45)= 'glycg1 @ 3.82+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(46)= 'glycg2 @ 3.75+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(47)= 'glycg3 @ 3.68+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(48)= 'glycg4 @ 3.61+-.01 FWHM=.09<.11+-.01 AMP=1.'

         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         if (ppmst .ge. 5.) then
            namrel= 'Water'
         else
            if (sptype(:7) .eq. 'liver-2'   .or.
     1          sptype(:8) .eq. 'breast-2'   .or. 
     2          sptype(:10) .eq. 'only-cho-1') idgppm = 1
            namrel= 'Lip09+Lip13'
            nchles= 2
            nnot1= 0
            ralinc= 1.4
            rfwhcc= 0.
            rrt2mq= .25
         end if
         if (sptype(:7) .eq. 'lipid-3'   .or.
     1       sptype(:7) .eq. 'liver-3'   .or.
     2       sptype(:8) .eq. 'breast-3'   .or. 
     3       sptype(:10) .eq. 'only-cho-2') idgppm = 2

         nrefpk(2)= 1
         nsimul= 48
         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         sdgrsh(2)= .006
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
      end if

c     444444444444444444444444444444444444444444444444444444444444444444444444
c     Original settings for *-4
c     444444444444444444444444444444444444444444444444444444444444444444444444
      if (sptype(:7) .eq. 'lipid-4'   .or.
     1    sptype(:7) .eq. 'liver-4'   .or.
     2    sptype(:8) .eq. 'breast-4')  then
         is_sptype= .true.
         atth2o= 1.
         chcomb(1)= 'L20a+L20b+L20c+L20d+L20e+L20f+L16a+L16b+L16c+L13a+L
     113b+L13c+L13d+L13e+L13f+L13g+L13h+L13i+L09a+L09b+L09c+L09d+L09e+L0
     29f'
         chcomb(2)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e+L13f+L13g+L
     113h+L13i+L09a+L09b+L09c+L09d+L09e+L09f'
         chcomb(3)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e+L13f+L13g+L
     113h+L13i'
         chcomb(4)= 'L13i+L13h+L13g+L13f+L13e+L13d+L13c+L13b+L13a'
         chcomb(5)= 'L09f+L09e+L09d+L09c+L09b+L09a'
         chcomb(6)= 'L16a+L16b+L16c'
         chcomb(7)= 'L20f+L20e+L20d+L20c+L20b+L20a'
         chcomb(8)= 'L53f+L53e+L53d+L53c+L53b+L53a'
         chcomb(9)= 'L53a+L53b+L53c+L53d+L53e+L53f+W1+W2+W3+W4+W5+W6+W7+
     1W8+W9+W10+W11'
         chcomb(10)= 'W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11'
         chcomb(11)= 'cho1+cho2+cho3+cho4+cho5+cho6+cho7'
         chcomb(12)= 'glycg1+glycg2+glycg3+glycg4'

         chcom2(1)= 'L20+L16+L09+L13'
         chcom2(2)= 'L16+L09+L13'
         chcom2(3)= 'Lip16+Lip13'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip20'
         chcom2(8)= 'Lip53'
         chcom2(9)= 'Lip53+Water'
         chcom2(10)= 'Water'
         chcom2(11)= 'Cho'
         chcom2(12)= 'Glycg'

         chgrsh(1)= 'L53'
         chgrsh(2)= 'L20'
         chgrsh(3)= 'L13'
         chgrsh(4)= 'L09'
         chgrsh(5)= 'cho'
         chgrsh(6)= 'L16'

         chsimu(1)= 'L53a @ 5.33+-.02 FWHM=.2<.23+-.03 AMP=1.'
         chsimu(2)= 'L53b @ 5.33+-.02 FWHM=.3<.33+-.03 AMP=1.'
         chsimu(3)= 'L53c @ 5.33+-.02 FWHM=.4<.43+-.03 AMP=1.'
         chsimu(4)= 'L53d @ 5.33+-.02 FWHM=.5<.53+-.03 AMP=1.'
         chsimu(5)= 'L53e @ 5.33+-.02 FWHM=.6<.63+-.03 AMP=1.'
         chsimu(6)= 'L53f @ 5.33+-.02 FWHM=.7<.73+-.03 AMP=1.'

         chsimu(7)= 'W1 @ 4.65+-.05 FWHM=.01<.04+-.03 AMP=1.'
         chsimu(8)= 'W2 @ 4.65+-.05 FWHM=.05<.08+-.03 AMP=1.'
         chsimu(9)= 'W3 @ 4.65+-.05 FWHM=.15<.18+-.03 AMP=1.'
         chsimu(10)= 'W4 @ 4.65+-.05 FWHM=.25<.28+-.03 AMP=1.'
         chsimu(11)= 'W5 @ 4.65+-.05 FWHM=.35<.38+-.03 AMP=1.'
         chsimu(12)= 'W6 @ 4.65+-.05 FWHM=.45<.48+-.03 AMP=1.'
         chsimu(13)= 'W7 @ 4.65+-.05 FWHM=.55<.58+-.03 AMP=1.'
         chsimu(14)= 'W8 @ 4.65+-.05 FWHM=.65<.68+-.03 AMP=1.'
         chsimu(15)= 'W9 @ 4.65+-.05 FWHM=.75<.78+-.03 AMP=1.'
         chsimu(16)= 'W10 @ 4.65+-.05 FWHM=.85<.88+-.03 AMP=1.'
         chsimu(17)= 'W11 @ 4.65+-.05 FWHM=.95<.98+-.03 AMP=1.'

         chsimu(18)= 'L20a @ 2.02+-.03 FWHM=.07<.1+-.03 AMP=.58   
     1               @ 2.23 FWHM=.07 AMP=.29  @ 2.75 FWHM=.060 AMP=.13'
         chsimu(19)= 'L20b @ 2.02+-.03 FWHM=.17<.2+-.03 AMP=.58  
     1               @ 2.23 FWHM=.17 AMP=.29  @ 2.75 FWHM=.146 AMP=.13'
         chsimu(20)= 'L20c @ 2.02+-.03 FWHM=.27<.3+-.03 AMP=.58 
     1               @ 2.23 FWHM=.27 AMP=.29  @ 2.75 FWHM=.231 AMP=.13'
         chsimu(21)= 'L20d @ 2.02+-.03 FWHM=.37<.4+-.03 AMP=.58  
     1               @ 2.23 FWHM=.37 AMP=.29  @ 2.75 FWHM=.317 AMP=.13'
         chsimu(22)= 'L20e @ 2.02+-.03 FWHM=.47<.5+-.03 AMP=.58  
     1               @ 2.23 FWHM=.47 AMP=.29  @ 2.75 FWHM=.403 AMP=.13'
         chsimu(23)= 'L20f @ 2.02+-.03 FWHM=.57<.6+-.03 AMP=.58  
     1               @ 2.23 FWHM=.57 AMP=.29  @ 2.75 FWHM=.488 AMP=.13'

         chsimu(24)= 'L16a @ 1.58+-.01 FWHM=.06<.09+-.03 AMP=1.'
         chsimu(25)= 'L16b @ 1.58+-.01 FWHM=.16<.19+-.03 AMP=1.'
         chsimu(26)= 'L16c @ 1.58+-.01 FWHM=.26<.29+-.03 AMP=1.'

         chsimu(27)= 'L13a @ 1.28+-.02 FWHM=.04<.05+-.01 AMP=1.'
         chsimu(28)= 'L13b @ 1.28+-.02 FWHM=.07<.09+-.02 AMP=1.'
         chsimu(29)= 'L13c @ 1.28+-.02 FWHM=.13<.16+-.03 AMP=1.'
         chsimu(30)= 'L13d @ 1.28+-.02 FWHM=.23<.26+-.03 AMP=1.'
         chsimu(31)= 'L13e @ 1.28+-.02 FWHM=.33<.36+-.03 AMP=1.'
         chsimu(32)= 'L13f @ 1.28+-.02 FWHM=.43<.46+-.03 AMP=1.'
         chsimu(33)= 'L13g @ 1.28+-.02 FWHM=.53<.56+-.03 AMP=1.'
         chsimu(34)= 'L13h @ 1.28+-.02 FWHM=.63<.66+-.03 AMP=1.'
         chsimu(35)= 'L13i @ 1.28+-.02 FWHM=.73<.76+-.03 AMP=1.'

         chsimu(36)= 'L09a @ .87+-.02 FWHM=.1<.13+-.03 AMP=1.'
         chsimu(37)= 'L09b @ .87+-.02 FWHM=.2<.23+-.03 AMP=1.'
         chsimu(38)= 'L09c @ .87+-.02 FWHM=.3<.33+-.03 AMP=1.'
         chsimu(39)= 'L09d @ .87+-.02 FWHM=.4<.43+-.03 AMP=1.'
         chsimu(40)= 'L09e @ .87+-.02 FWHM=.5<.53+-.03 AMP=1.'
         chsimu(41)= 'L09f @ .87+-.02 FWHM=.6<.63+-.03 AMP=1.'

         chsimu(42)= 'cho1 @ 3.185+-.01 FWHM=.030<.050+-.012 AMP=1.'
         chsimu(43)= 'cho2 @ 3.185+-.01 FWHM=.074<.086+-.012 AMP=1.'
         chsimu(44)= 'cho3 @ 3.185+-.01 FWHM=.110<.122+-.012 AMP=1.'
         chsimu(45)= 'cho4 @ 3.185+-.01 FWHM=.146<.158+-.012 AMP=1.'
         chsimu(46)= 'cho5 @ 3.185+-.01 FWHM=.182<.194+-.012 AMP=1.'
         chsimu(47)= 'cho6 @ 3.185+-.01 FWHM=.218<.230+-.012 AMP=1.'
         chsimu(48)= 'cho7 @ 3.185+-.01 FWHM=.254<.266+-.012 AMP=1.'

         chsimu(49)= 'glycg1 @ 3.82+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(50)= 'glycg2 @ 3.75+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(51)= 'glycg3 @ 3.68+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(52)= 'glycg4 @ 3.61+-.01 FWHM=.09<.11+-.01 AMP=1.'

         dkntmn(2)= 2.4
         endpha= .true.
         fwhmmx= .4
         fwhmst= .3
         gauss_rt2= .false.
         incsmx= 1
         idgppm= 2
         ncombi= 12
         ngrsh= 6
         nobasi= .true.
         nratio= 0
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 52
         onlyco= .true.
         ppmref(1,2)= 1.28

         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -1.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         if (ppmst .lt. 5.) then
            badref= .true.
            chless(1)= 'L09'
            chless(2)= 'L20'
            chnot1(1)= 'L16a'
            chnot1(2)= 'L16b'
            chnot1(3)= 'L16c'
            namrel= 'L16+L09+L13'
            nchles= 2
            nnot1= 3
            ralinc= 1.4
            rfwhcc= 0.
            rrt2mq= .25
         else
            chuse1(1)= 'L13a'
            chuse1(2)= 'L13b'
            chuse1(3)= 'L13c'
            chuse1(4)= 'L13d'
            chuse1(5)= 'L13e'
            chuse1(6)= 'L13f'
            chuse1(7)= 'L13g'
            chuse1(8)= 'L13h'
            chuse1(9)= 'L13i'
            chuse1(10)= 'W1'
            chuse1(11)= 'W2'
            chuse1(12)= 'W3'
            chuse1(13)= 'W4'
            chuse1(14)= 'W5'
            chuse1(15)= 'W6'
            chuse1(16)= 'W7'
            chuse1(17)= 'W8'
            chuse1(18)= 'W9'
            chuse1(19)= 'W10'
            chuse1(20)= 'W11'
            namrel= 'Water'
            nuse1= 20
         end if

         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         sdgrsh(1)= .003
         sdgrsh(2)= .006
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .002
         sdgrsh(6)= .003
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     555555555555555555555555555555555555555555555555555555555555555555555555
c     Complete settings for *-5 (resolved L21, L23 & L28) and start of *-6.
c     DOWATR=T to include water; else, there is a PPMGAP in [5.0,3.6]
c     555555555555555555555555555555555555555555555555555555555555555555555555
      if (sptype(:7) .eq. 'lipid-5'   .or.
     1    sptype(:8) .eq. 'breast-5'   .or.
     2    sptype(:7) .eq. 'liver-5') then
         is_sptype= .true.
         atth2o= 1.

         chcomb(1)= 'L16a+L16b+L16c+L16d+L13a+L13b+L13c+L13d+L13e+L13f+L
     109a+L09b+L09c+L09d+L09e'
         chcomb(2)= 'L16a+L16b+L16c+L16d+L13a+L13b+L13c+L13d+L13e+L13f'
         chcomb(3)= 'L28a+L28b+L28c+L28d+L23a+L23b+L23c+L23d+L21a+L21b+L
     121c+L21d'
         chcomb(4)= 'L13f+L13e+L13d+L13c+L13b+L13a'
         chcomb(5)= 'L09e+L09d+L09c+L09b+L09a'
         chcomb(6)= 'L16a+L16b+L16c+L16d'
         chcomb(7)= 'L21a+L21b+L21c+L21d'
         chcomb(8)= 'L23a+L23b+L23c+L23d'
         chcomb(9)= 'L28a+L28b+L28c+L28d'
         chcomb(10)= 'L53a+L53b+L53c+L53d+L52a+L52b+L52c+L52d'
         chcomb(11)= 'L53a+L53b+L53c+L53d+L52a+L52b+L52c+L52d+W1+W2+W3+W
     14+W5+W6+W7+W8+W9'
         chcomb(12)= 'W1+W2+W3+W4+W5+W6+W7+W8+W9'
         chcomb(13)= 'cho1+cho2+cho3+cho4+cho5'
         chcomb(14)= 'glycg1+glycg2+glycg3+glycg4'

         chcom2(1)= 'L16+L09+L13'
         chcom2(2)= 'Lip16+Lip13'
         chcom2(3)= 'L28+L23+L21'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip21'
         chcom2(8)= 'Lip23'
         chcom2(9)= 'Lip28'
         chcom2(10)= 'Lip53+Lip52'
         chcom2(11)= 'L53+L52+Water'
         chcom2(12)= 'Water'
         chcom2(13)= 'Cho'
         chcom2(14)= 'Glycg'

         chgrsh(1)= 'L28'
         chgrsh(2)= 'L23'
         chgrsh(3)= 'L21'
         chgrsh(4)= 'L13'
         chgrsh(5)= 'L09'
         chgrsh(6)= 'cho'
         chgrsh(7)= 'L16'

         chless(1)= 'L09'
         chless(2)= 'L20'

         chnot2(1)= 'L20a'
         chnot2(2)= 'L20b'
         chnot2(3)= 'L20c'

         chsimu(1)= 'L53a @ 5.31+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(2)= 'L53b @ 5.31+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(3)= 'L53c @ 5.31+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(4)= 'L53d @ 5.31+-.02 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(5)= 'L52a @ 5.21+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(6)= 'L52b @ 5.21+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(7)= 'L52c @ 5.21+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(8)= 'L52d @ 5.21+-.02 FWHM=.27<.290+-.0100 AMP=1.'

         if (roomt) then
            chsimu(9)= 'W1 @ 4.83+-.05 FWHM=.01<.04+-.03 AMP=1.'
            chsimu(10)= 'W2 @ 4.83+-.05 FWHM=.05<.08+-.03 AMP=1.'
            chsimu(11)= 'W3 @ 4.83+-.05 FWHM=.15<.18+-.03 AMP=1.'
            chsimu(12)= 'W4 @ 4.83+-.05 FWHM=.25<.28+-.03 AMP=1.'
            chsimu(13)= 'W5 @ 4.83+-.05 FWHM=.35<.38+-.03 AMP=1.'
            chsimu(14)= 'W6 @ 4.83+-.05 FWHM=.45<.48+-.03 AMP=1.'
            chsimu(15)= 'W7 @ 4.83+-.05 FWHM=.55<.58+-.03 AMP=1.'
            chsimu(16)= 'W8 @ 4.83+-.05 FWHM=.65<.68+-.03 AMP=1.'
            chsimu(17)= 'W9 @ 4.83+-.05 FWHM=.75<.78+-.03 AMP=1.'
         else
            chsimu(9)= 'W1 @ 4.65+-.05 FWHM=.01<.04+-.03 AMP=1.'
            chsimu(10)= 'W2 @ 4.65+-.05 FWHM=.05<.08+-.03 AMP=1.'
            chsimu(11)= 'W3 @ 4.65+-.05 FWHM=.15<.18+-.03 AMP=1.'
            chsimu(12)= 'W4 @ 4.65+-.05 FWHM=.25<.28+-.03 AMP=1.'
            chsimu(13)= 'W5 @ 4.65+-.05 FWHM=.35<.38+-.03 AMP=1.'
            chsimu(14)= 'W6 @ 4.65+-.05 FWHM=.45<.48+-.03 AMP=1.'
            chsimu(15)= 'W7 @ 4.65+-.05 FWHM=.55<.58+-.03 AMP=1.'
            chsimu(16)= 'W8 @ 4.65+-.05 FWHM=.65<.68+-.03 AMP=1.'
            chsimu(17)= 'W9 @ 4.65+-.05 FWHM=.75<.78+-.03 AMP=1.'
         end if

         chsimu(18)= 'L28a @ 2.75+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(19)= 'L28b @ 2.75+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(20)= 'L28c @ 2.75+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(21)= 'L28d @ 2.75+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(22)= 'L23a @ 2.23+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(23)= 'L23b @ 2.23+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(24)= 'L23c @ 2.23+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(25)= 'L23d @ 2.23+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(26)= 'L21a @ 2.02+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(27)= 'L21b @ 2.02+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(28)= 'L21c @ 2.02+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(29)= 'L21d @ 2.02+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(30)= 'L20a @ 2.02+-.03 FWHM=.07<.1+-.03 AMP=.58   
     1               @ 2.23 FWHM=.07 AMP=.29  @ 2.75 FWHM=.060 AMP=.13'
         chsimu(31)= 'L20b @ 2.02+-.03 FWHM=.17<.2+-.03 AMP=.58  
     1               @ 2.23 FWHM=.17 AMP=.29  @ 2.75 FWHM=.146 AMP=.13'
         chsimu(32)= 'L20c @ 2.02+-.03 FWHM=.27<.3+-.03 AMP=.58 
     1               @ 2.23 FWHM=.27 AMP=.29  @ 2.75 FWHM=.231 AMP=.13'

         chsimu(33)= 'L16a @ 1.58+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(34)= 'L16b @ 1.58+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(35)= 'L16c @ 1.58+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(36)= 'L16d @ 1.58+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(37)= 'L13a @ 1.28+-.02 FWHM=.03<.050+-.0100 AMP=1.'
         chsimu(38)= 'L13b @ 1.28+-.02 FWHM=.07<.105+-.0175 AMP=1.'
         chsimu(39)= 'L13c @ 1.28+-.02 FWHM=.14<.175+-.0175 AMP=1.'
         chsimu(40)= 'L13d @ 1.28+-.02 FWHM=.21<.245+-.0175 AMP=1.'
         chsimu(41)= 'L13e @ 1.28+-.02 FWHM=.28<.315+-.0175 AMP=1.'
         chsimu(42)= 'L13f @ 1.28+-.02 FWHM=.35<.370+-.0100 AMP=1.'

         chsimu(43)= 'L09a @ .87+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(44)= 'L09b @ .87+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(45)= 'L09c @ .87+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(46)= 'L09d @ .87+-.02 FWHM=.27<.305+-.0175 AMP=1.'
         chsimu(47)= 'L09e @ .87+-.02 FWHM=.34<.360+-.0100 AMP=1.'

         chsimu(48)= 'cho1 @ 3.185+-.01 FWHM=.03<.05+-.01 AMP=1.'
         chsimu(49)= 'cho2 @ 3.185+-.01 FWHM=.07<.09+-.01 AMP=1.'
         chsimu(50)= 'cho3 @ 3.185+-.01 FWHM=.11<.13+-.01 AMP=1.'
         chsimu(51)= 'cho4 @ 3.185+-.01 FWHM=.15<.17+-.01 AMP=1.'
         chsimu(52)= 'cho5 @ 3.185+-.01 FWHM=.19<.21+-.01 AMP=1.'

         chsimu(53)= 'glycg1 @ 3.82+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(54)= 'glycg2 @ 3.75+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(55)= 'glycg3 @ 3.68+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(56)= 'glycg4 @ 3.61+-.01 FWHM=.09<.11+-.01 AMP=1.'

         dkntmn(2)= 99.
         endpha= .true.
         fwhmmx= .4
         fwhmst= .3
         gauss_rt2= .false.
         incsmx= 1
         idgppm= 2
         ncombi= 14
         ngrsh= 7
         nnot2= 3
         nobasi= .true.
         nratio= 0
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 56
         onlyco= .true.
         ppmref(1,2)= 1.28

         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -1.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         if (dowatr   .and.   ppmst .ge. 5.) then
            chuse1(1)= 'L13a'
            chuse1(2)= 'L13b'
            chuse1(3)= 'L13c'
            chuse1(4)= 'L13d'
            chuse1(5)= 'L13e'
            chuse1(6)= 'L13f'
            chuse1(7)= 'W1'
            chuse1(8)= 'W2'
            chuse1(9)= 'W3'
            chuse1(10)= 'W4'
            chuse1(11)= 'W5'
            chuse1(12)= 'W6'
            chuse1(13)= 'W7'
            chuse1(14)= 'W8'
            chuse1(15)= 'W9'
            namrel= 'Water'
            nuse1= 15
         else
            chuse1(1)= 'L20a'
            chuse1(2)= 'L20b'
            chuse1(3)= 'L20c'
            chuse1(4)= 'L13a'
            chuse1(5)= 'L13b'
            chuse1(6)= 'L13c'
            chuse1(7)= 'L13d'
            chuse1(8)= 'L13e'
            chuse1(9)= 'L13f'
            chuse1(10)= 'L09a'
            chuse1(11)= 'L09b'
            chuse1(12)= 'L09c'
            chuse1(13)= 'L09d'
            chuse1(14)= 'L09e'
            nchles= 2
            nuse1= 14
            if (ppmst .ge. 5.) then
               namrel= 'Lip09'
               ppmgap(1,1)= 5.
               ppmgap(2,1)= 3.6
            else
               namrel= 'Lip23'
               ralinc= 1.4
               rfwhcc= 0.
               rrt2mq= .25
            end if
         end if

         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         sdgrsh(1)= .003
         sdgrsh(2)= .003
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .003
         sdgrsh(6)= .002
         sdgrsh(7)= .003
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     66666666666666666666666666666666666666666666666666666666666666666666666
c     Complete settings for *-6.  L43 & L41; L20 out of & cho into Prel; 
c        RLESMO=1.5
c     66666666666666666666666666666666666666666666666666666666666666666666666
      if (sptype(:7) .eq. 'lipid-6'   .or.
     1    sptype(:7) .eq. 'liver-6'   .or.
     2    sptype(:8) .eq. 'breast-6')  then
         is_sptype= .true.
         atth2o= 1.

         chcomb(1)= 'L16a+L16b+L16c+L16d+L13a+L13b+L13c+L13d+L13e+L13f+L
     109a+L09b+L09c+L09d+L09e'
         chcomb(2)= 'L16a+L16b+L16c+L16d+L13a+L13b+L13c+L13d+L13e+L13f'
         chcomb(3)= 'L28a+L28b+L28c+L28d+L23a+L23b+L23c+L23d+L21a+L21b+L
     121c+L21d'
         chcomb(4)= 'L13f+L13e+L13d+L13c+L13b+L13a'
         chcomb(5)= 'L09e+L09d+L09c+L09b+L09a'
         chcomb(6)= 'L16a+L16b+L16c+L16d'
         chcomb(7)= 'L21a+L21b+L21c+L21d'
         chcomb(8)= 'L23a+L23b+L23c+L23d'
         chcomb(9)= 'L28a+L28b+L28c+L28d'
         chcomb(10)= 'L53a+L53b+L53c+L53d+L52a+L52b+L52c+L52d'
         chcomb(11)= 'W1+W2+W3+W4+W5+W6+W7+W8+W9'
         chcomb(12)= 'cho1+cho2+cho3+cho4+cho5'
         chcomb(13)= 'L43a+L43b'
         chcomb(14)= 'L41a+L41b'
         chcomb(15)= 'glycg1+glycg2+glycg3+glycg4'
         chcomb(16)= 'L53a+L53b+L53c+L53d'
         chcomb(17)= 'L52a+L52b+L52c+L52d'

         chcom2(1)= 'L16+L09+L13'
         chcom2(2)= 'Lip16+Lip13'
         chcom2(3)= 'L28+L23+L21'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip21'
         chcom2(8)= 'Lip23'
         chcom2(9)= 'Lip28'
         chcom2(10)= 'Lip53+Lip52'
         chcom2(11)= 'Water'
         chcom2(12)= 'Cho'
         chcom2(13)= 'Lip43'
         chcom2(14)= 'Lip41'
         chcom2(15)= 'Glycg'
         chcom2(16)= 'Lip53'
         chcom2(17)= 'Lip52'

         chgrsh(1)= 'L28'
         chgrsh(2)= 'L23'
         chgrsh(3)= 'L21'
         chgrsh(4)= 'L13'
         chgrsh(5)= 'L09'
         chgrsh(6)= 'cho'
         chgrsh(7)= 'L16'
         chgrsh(8)= 'L43'
         chgrsh(9)= 'L41'
         chgrsh(10)= 'W'

         chless(1)= 'L09'
         chless(2)= 'L20'

         chnot2(1)= 'L20a'
         chnot2(2)= 'L20b'
         chnot2(3)= 'L20c'

         chrato(1)= 'L41a/L43a = 1. +- .02'
         chrato(2)= 'L41b/L43b = 1. +- .02'

         chsimu(1)= 'L53a @ 5.31+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(2)= 'L53b @ 5.31+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(3)= 'L53c @ 5.31+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(4)= 'L53d @ 5.31+-.02 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(5)= 'L52a @ 5.21+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(6)= 'L52b @ 5.21+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(7)= 'L52c @ 5.21+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(8)= 'L52d @ 5.21+-.02 FWHM=.27<.290+-.0100 AMP=1.'

         if (roomt) then
            chsimu(9)= 'W1 @ 4.83+-.03 FWHM=.01<.04+-.015 AMP=1.'
            chsimu(10)= 'W2 @ 4.83+-.03 FWHM=.05<.1+-.025 AMP=1.'
            chsimu(11)= 'W3 @ 4.83+-.03 FWHM=.15<.2+-.025 AMP=1.'
            chsimu(12)= 'W4 @ 4.83+-.03 FWHM=.25<.3+-.025 AMP=1.'
            chsimu(13)= 'W5 @ 4.83+-.03 FWHM=.35<.4+-.025 AMP=1.'
            chsimu(14)= 'W6 @ 4.83+-.03 FWHM=.45<.5+-.025 AMP=1.'
            chsimu(15)= 'W7 @ 4.83+-.03 FWHM=.55<.6+-.025 AMP=1.'
            chsimu(16)= 'W8 @ 4.83+-.03 FWHM=.65<.7+-.025 AMP=1.'
            chsimu(17)= 'W9 @ 4.83+-.03 FWHM=.75<.8+-.025 AMP=1.'
         else
            chsimu(9)= 'W1 @ 4.65+-.03 FWHM=.01<.04+-.015 AMP=1.'
            chsimu(10)= 'W2 @ 4.65+-.03 FWHM=.05<.1+-.025 AMP=1.'
            chsimu(11)= 'W3 @ 4.65+-.03 FWHM=.15<.2+-.025 AMP=1.'
            chsimu(12)= 'W4 @ 4.65+-.03 FWHM=.25<.3+-.025 AMP=1.'
            chsimu(13)= 'W5 @ 4.65+-.03 FWHM=.35<.4+-.025 AMP=1.'
            chsimu(14)= 'W6 @ 4.65+-.03 FWHM=.45<.5+-.025 AMP=1.'
            chsimu(15)= 'W7 @ 4.65+-.03 FWHM=.55<.6+-.025 AMP=1.'
            chsimu(16)= 'W8 @ 4.65+-.03 FWHM=.65<.7+-.025 AMP=1.'
            chsimu(17)= 'W9 @ 4.65+-.03 FWHM=.75<.8+-.025 AMP=1.'
         end if

         chsimu(18)= 'L43a@4.28+-.004 FWHM=.04<.08+-.02 AMP=1.'
         chsimu(19)= 'L43b@4.28+-.004 FWHM=.12<.16+-.02 AMP=1.'
         chsimu(20)= 'L41a@4.08+-.004 FWHM=.04<.08+-.02 AMP=1.'
         chsimu(21)= 'L41b@4.08+-.004 FWHM=.12<.16+-.02 AMP=1.'

         chsimu(22)= 'L28a @ 2.75+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(23)= 'L28b @ 2.75+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(24)= 'L28c @ 2.75+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(25)= 'L28d @ 2.75+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(26)= 'L23a @ 2.23+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(27)= 'L23b @ 2.23+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(28)= 'L23c @ 2.23+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(29)= 'L23d @ 2.23+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(30)= 'L21a @ 2.02+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(31)= 'L21b @ 2.02+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(32)= 'L21c @ 2.02+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(33)= 'L21d @ 2.02+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(34)= 'L20a @ 2.02+-.03 FWHM=.07<.1+-.03 AMP=.58   
     1               @ 2.23 FWHM=.07 AMP=.29  @ 2.75 FWHM=.060 AMP=.13'
         chsimu(35)= 'L20b @ 2.02+-.03 FWHM=.17<.2+-.03 AMP=.58  
     1               @ 2.23 FWHM=.17 AMP=.29  @ 2.75 FWHM=.146 AMP=.13'
         chsimu(36)= 'L20c @ 2.02+-.03 FWHM=.27<.3+-.03 AMP=.58 
     1               @ 2.23 FWHM=.27 AMP=.29  @ 2.75 FWHM=.231 AMP=.13'

         chsimu(37)= 'L16a @ 1.58+-.01 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(38)= 'L16b @ 1.58+-.01 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(39)= 'L16c @ 1.58+-.01 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(40)= 'L16d @ 1.58+-.01 FWHM=.27<.290+-.0100 AMP=1.'

         chsimu(41)= 'L13a @ 1.28+-.02 FWHM=.03<.050+-.0100 AMP=1.'
         chsimu(42)= 'L13b @ 1.28+-.02 FWHM=.07<.105+-.0175 AMP=1.'
         chsimu(43)= 'L13c @ 1.28+-.02 FWHM=.14<.175+-.0175 AMP=1.'
         chsimu(44)= 'L13d @ 1.28+-.02 FWHM=.21<.245+-.0175 AMP=1.'
         chsimu(45)= 'L13e @ 1.28+-.02 FWHM=.28<.315+-.0175 AMP=1.'
         chsimu(46)= 'L13f @ 1.28+-.02 FWHM=.35<.370+-.0100 AMP=1.'

         chsimu(47)= 'L09a @ .87+-.02 FWHM=.06<.095+-.0175 AMP=1.'
         chsimu(48)= 'L09b @ .87+-.02 FWHM=.13<.165+-.0175 AMP=1.'
         chsimu(49)= 'L09c @ .87+-.02 FWHM=.20<.235+-.0175 AMP=1.'
         chsimu(50)= 'L09d @ .87+-.02 FWHM=.27<.305+-.0175 AMP=1.'
         chsimu(51)= 'L09e @ .87+-.02 FWHM=.34<.360+-.0100 AMP=1.'

         chsimu(52)= 'cho1 @ 3.185+-.01 FWHM=.03<.05+-.01 AMP=1.'
         chsimu(53)= 'cho2 @ 3.185+-.01 FWHM=.07<.09+-.01 AMP=1.'
         chsimu(54)= 'cho3 @ 3.185+-.01 FWHM=.11<.13+-.01 AMP=1.'
         chsimu(55)= 'cho4 @ 3.185+-.01 FWHM=.15<.17+-.01 AMP=1.'
         chsimu(56)= 'cho5 @ 3.185+-.01 FWHM=.19<.21+-.01 AMP=1.'

         chsimu(57)= 'glycg1 @ 3.82+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(58)= 'glycg2 @ 3.75+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(59)= 'glycg3 @ 3.68+-.01 FWHM=.09<.11+-.01 AMP=1.'
         chsimu(60)= 'glycg4 @ 3.61+-.01 FWHM=.09<.11+-.01 AMP=1.'

         dkntmn(2)= 1.867
         endpha= .true.
         fwhmmx= .4
         fwhmst= .3
         gauss_rt2= .false.
         incsmx= 1
         isdbol= 5
         ncombi= 17
         ngrsh= 10
         nnot2= 3
         nobasi= .true.
         nratio= 2
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 60
         onlyco= .true.
         ppmref(1,2)= 1.28

         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         if (dowatr   .and.   ppmst .ge. 5.) then
            chuse1(1)= 'L13a'
            chuse1(2)= 'L13b'
            chuse1(3)= 'L13c'
            chuse1(4)= 'L13d'
            chuse1(5)= 'L13e'
            chuse1(6)= 'L13f'
            chuse1(7)= 'W1'
            chuse1(8)= 'W2'
            chuse1(9)= 'W3'
            chuse1(10)= 'W4'
            chuse1(11)= 'W5'
            chuse1(12)= 'W6'
            chuse1(13)= 'W7'
            chuse1(14)= 'W8'
            chuse1(15)= 'W9'
            namrel= 'Water'
            nuse1= 15
         else
            chuse1(1)= 'L13a'
            chuse1(2)= 'L13b'
            chuse1(3)= 'L13c'
            chuse1(4)= 'L13d'
            chuse1(5)= 'L13e'
            chuse1(6)= 'L13f'
            chuse1(7)= 'L09a'
            chuse1(8)= 'L09b'
            chuse1(9)= 'L09c'
            chuse1(10)= 'L09d'
            chuse1(11)= 'L09e'
            chuse1(12)= 'cho1'
            chuse1(13)= 'cho2'
            chuse1(14)= 'cho3'
            chuse1(15)= 'cho4'
            chuse1(16)= 'cho5'
            chuse1(17)= 'L20a'
            chuse1(18)= 'L20b'
            chuse1(19)= 'L20c'
            namrel= 'Lip13'
            nchles= 1
            nuse1= 16
            rlesmo= 1.5
            if (ppmst .ge. 5.) then
               ppmgap(1,1)= 5.
               ppmgap(2,1)= 3.6
               ppmsep(1)= 4.3
            else
               ralinc= 1.4
               rfwhcc= 0.
               rrt2mq= .25
            end if
         end if

         if (ppmst .lt. 5.) idgppm= 2

         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         sdgrsh(1)= .003
         sdgrsh(2)= .003
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .003
         sdgrsh(6)= .002
         sdgrsh(7)= .003
         sdgrsh(8)= .003
         sdgrsh(9)= .003
         sdgrsh(10)= .01
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     77777777777777777777777777777777777777777777777777777777777777777777777
c     Complete settings for *-7 and start of *-8 & liver-11.
c     Regularize peak widths.
c     88888888888888888888888888888888888888888888888888888888888888888888888
      if (sptype(:7) .eq. 'lipid-7'   .or.
     1    sptype(:7) .eq. 'liver-7'   .or.
     2    sptype(:8) .eq. 'breast-7'  .or.
     3    sptype(:7) .eq. 'lipid-8'   .or.
     4    sptype(:7) .eq. 'liver-8'   .or.
     5    sptype(:8) .eq. 'breast-8'   .or.
     6    sptype(:8) .eq. 'liver-11' )  then
         is_sptype= .true.
         accept_step2= .true.
         alpsmn= 1.d-6
         alpsmx= 2.d0
         alpsst= 1.d-4
         atth2o= 1.

         chcomb(1)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e+L09a+L09b+L
     109c+L09d'
         chcomb(2)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e'
         chcomb(3)= 'L28a+L28b+L28c+L23a+L23b+L23c+L21a+L21b+L21c'
         chcomb(4)= 'L13e+L13d+L13c+L13b+L13a'
         chcomb(5)= 'L09d+L09c+L09b+L09a'
         chcomb(6)= 'L16a+L16b+L16c'
         chcomb(7)= 'L21a+L21b+L21c'
         chcomb(8)= 'L23a+L23b+L23c'
         chcomb(9)= 'L28a+L28b+L28c'
         chcomb(10)= 'L53a+L53b+L53c+L52a+L52b+L52c'
         chcomb(11)= 'W1+W2+W3+W4+W5+W6'
         chcomb(12)= 'cho1+cho2+cho3'
         chcomb(13)= 'L43a+L43b'
         chcomb(14)= 'L41a+L41b'
         chcomb(15)= 'glycg1+glycg2+glycg3+glycg4'
         chcomb(16)= 'L53a+L53b+L53c'
         chcomb(17)= 'L52a+L52b+L52c'

         chcom2(1)= 'L16+L09+L13'
         chcom2(2)= 'Lip16+Lip13'
         chcom2(3)= 'L28+L23+L21'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip21'
         chcom2(8)= 'Lip23'
         chcom2(9)= 'Lip28'
         chcom2(10)= 'Lip53+Lip52'
         chcom2(11)= 'Water'
         chcom2(12)= 'Cho'
         chcom2(13)= 'Lip43'
         chcom2(14)= 'Lip41'
         chcom2(15)= 'Glycg'
         chcom2(16)= 'Lip53'
         chcom2(17)= 'Lip52'

         chgrsh(1)= 'L28'
         chgrsh(2)= 'L23'
         chgrsh(3)= 'L21'
         chgrsh(4)= 'L13'
         chgrsh(5)= 'L09'
         chgrsh(6)= 'cho'
         chgrsh(7)= 'L16'
         chgrsh(8)= 'L43'
         chgrsh(9)= 'L41'
         chgrsh(10)= 'W'

         chnot2(1)= 'W1pre'
         chnot2(2)= 'L13pre'

         chrato(1)= 'L41a/L43a = 1. +- .02'
         chrato(2)= 'L41b/L43b = 1. +- .02'

         chsimu(1)= 'L53a @ 5.31+-.02 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(2)= 'L53b @ 5.31+-.02 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(3)= 'L53c @ 5.31+-.02 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'

         chsimu(4)= 'L52a @ 5.21+-.02 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(5)= 'L52b @ 5.21+-.02 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(6)= 'L52c @ 5.21+-.02 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'

         if (roomt) then
            chsimu(7)= 'W1pre @ 4.83+-.03 FWHM=.03<.04+-.005 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(8)= 'W1 @ 4.83+-.03 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(9)= 'W2 @ 4.83+-.03 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(10)= 'W3 @ 4.83+-.03 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(11)= 'W4 @ 4.83+-.03 FWHM=1.6<1.8+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(12)= 'W5 @ 4.83+-.03 FWHM=2.2<2.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(13)= 'W6 @ 4.83+-.03 FWHM=2.8<3.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(7)= 'W1pre @ 4.65+-.03 FWHM=.03<.04+-.005 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(8)= 'W1 @ 4.65+-.03 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(9)= 'W2 @ 4.65+-.03 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(10)= 'W3 @ 4.65+-.03 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(11)= 'W4 @ 4.65+-.03 FWHM=1.6<1.8+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(12)= 'W5 @ 4.65+-.03 FWHM=2.2<2.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(13)= 'W6 @ 4.65+-.03 FWHM=2.8<3.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
         end if

         chsimu(14)= 'L43a @ 4.28+-.004 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(15)= 'L43b @ 4.28+-.004 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(16)= 'L41a @ 4.08+-.004 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(17)= 'L41b @ 4.08+-.004 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'

         chsimu(18)= 'L28a @ 2.75+-.01 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(19)= 'L28b @ 2.75+-.01 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(20)= 'L28c @ 2.75+-.01 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'

         chsimu(21)= 'L23a @ 2.23+-.01 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(22)= 'L23b @ 2.23+-.01 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(23)= 'L23c @ 2.23+-.01 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'

         chsimu(24)= 'L21a @ 2.02+-.01 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.1'
         chsimu(25)= 'L21b @ 2.02+-.01 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.1'
         chsimu(26)= 'L21c @ 2.02+-.01 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.1'

         chsimu(27)= 'L16a @ 1.58+-.01 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(28)= 'L16b @ 1.58+-.01 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'
         chsimu(29)= 'L16c @ 1.58+-.01 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.07'

         chsimu(30)= 'L13pre @ 1.28+-.02 FWHM=.03<.04+-.005 AMP=1.
     1               @999. FWHM=1. AMP=1.'
         chsimu(31)= 'L13a @ 1.28+-.02 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=1.'
         chsimu(32)= 'L13b @ 1.28+-.02 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=1.'
         chsimu(33)= 'L13c @ 1.28+-.02 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=1.'
         chsimu(34)= 'L13d @ 1.28+-.02 FWHM=1.6<1.8+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=1.'
         chsimu(35)= 'L13e @ 1.28+-.02 FWHM=2.2<2.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=1.'

         chsimu(36)= 'L09a @ .87+-.02 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.13'
         chsimu(37)= 'L09b @ .87+-.02 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.13'
         chsimu(38)= 'L09c @ .87+-.02 FWHM=1.4<1.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.13'
         chsimu(39)= 'L09d @ .87+-.02 FWHM=2.<2.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.13'
c        -------------------------------------------------------------------
c        Cho can be much sharper than Lip13.  Use only the sharpest models.
c        -------------------------------------------------------------------
         chsimu(40)= 'cho1 @ 3.2+-.01 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=.2'
         chsimu(41)= 'cho2 @ 3.2+-.01 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.2'
         chsimu(42)= 'cho3 @ 3.2+-.01 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.2'
c        -------------------------------------------------------------------
c        Only one model at each ppm.
c        -------------------------------------------------------------------
         chsimu(43)= 'glycg1 @ 3.82+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1               @999. FWHM=1. AMP=.1'
         chsimu(44)= 'glycg2 @ 3.75+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1               @999. FWHM=1. AMP=.1'
         chsimu(45)= 'glycg3 @ 3.68+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1               @999. FWHM=1. AMP=.1'
         chsimu(46)= 'glycg4 @ 3.61+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1               @999. FWHM=1. AMP=.1'

         chuse1(1)= 'L13pre'
         degmax(1)= 2.
         dkntmn(1)= .45
         fwhmmn= .1
         fwhmmx= .3
         fwhmst= .15
         gauss_rt2= .false.
         idgppm= 2
         imethd= 2
         ipowrg= 2
         isdbol= 5
         mrepha(1)= 1
         nbas_ccf= 0
         ncombi= 17
         ngrsh= 10
         nnot2= 2
         nobasi= .true.
         nratio = 0
         if (ppmst .ge. 4.28) nratio = 2
         nrefpk(2)= 1
         nsimul= 46
         onlyco= .true.
         ppmref(1,2)= 1.28
         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         namrel= 'L16+L09+L13'
         if (dowatr   .and.   ppmst .ge. 5.) then
            chuse1(2)= 'W1pre'
            if (.not.(dows   .or.   doecc)) namrel= 'Water'
            nuse1= 2
         else
            nuse1= 1
            if (ppmst .ge. 5.) then
               ppmgap(1,1)= 5.
               ppmgap(2,1)= 3.6
               ppmsep(1)= 4.3
            else
               rfwhcc= 0.
               rrt2mq= .25
            end if
         end if
         prmnmx(1,1)= .2
         prmnmx(2,1)= .5
         prmnmx(1,2)= .7
         prmnmx(2,2)= .9
         prmnmx(1,3)= .3
         prmnmx(2,3)= .69
         rbackg(1)= 22.2
         rbackg(2)= 11.1
         rlrntz= 2.
         scafwh= .true.
         sdgrsh(1)= .003
         sdgrsh(2)= .003
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .003
         sdgrsh(6)= .002
         sdgrsh(7)= .003
         sdgrsh(8)= .003
         sdgrsh(9)= .003
         sdgrsh(10)= .01
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
         wconc= .5
      end if

c     88888888888888888888888888888888888888888888888888888888888888888888888
c     Modify *-7 settings for *-8 and start of liver-11.
c     As liver-7, but with weaker regularization of peak widths.
c     88888888888888888888888888888888888888888888888888888888888888888888888
      if (sptype(:7) .eq. 'lipid-8'   .or.
     1    sptype(:7) .eq. 'liver-8'   .or.
     2    sptype(:8) .eq. 'breast-8'   .or.
     3    sptype(:8) .eq. 'liver-11')  then
         alpsmn= 2.d-4
         alpsmx= 1.d1
         alpsst= 1.d-2
         chcomb(15)= 'L53a+L53b+L53c'
         chcomb(16)= 'L52a+L52b+L52c'
         chcom2(15)= 'Lip53'
         chcom2(16)= 'Lip52'
         if (roomt) then
            chsimu(7)= 'W1pre @ 4.83+-.02 FWHM=.03<.04+-.005 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(8)= 'W1 @ 4.83+-.02 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(9)= 'W2 @ 4.83+-.02 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(10)= 'W3 @ 4.83+-.02 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(11)= 'W4 @ 4.83+-.02 FWHM=1.6<1.8+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(12)= 'W5 @ 4.83+-.02 FWHM=2.2<2.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(13)= 'W6 @ 4.83+-.02 FWHM=2.8<3.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(7)= 'W1pre @ 4.65+-.02 FWHM=.03<.04+-.005 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(8)= 'W1 @ 4.65+-.02 FWHM=.2<.3+-.025 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(9)= 'W2 @ 4.65+-.02 FWHM=.4<.6+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(10)= 'W3 @ 4.65+-.02 FWHM=1.<1.2+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(11)= 'W4 @ 4.65+-.02 FWHM=1.6<1.8+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(12)= 'W5 @ 4.65+-.02 FWHM=2.2<2.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
            chsimu(13)= 'W6 @ 4.65+-.02 FWHM=2.8<3.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=10.'
         end if

         chsimu(14)= 'L43a @ 4.27+-.005 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(15)= 'L43b @ 4.27+-.005 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(16)= 'L41a @ 4.07+-.005 FWHM=.2<.4+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         chsimu(17)= 'L41b @ 4.07+-.005 FWHM=.8<1.+-.0667 AMP=1.
     1               @999. FWHM=1. AMP=.035'
         degmax(1)= 12.5
         dkntmn(1)= 1.
         endpha= .true.
         fwhmmx= .4
         ipowrg= 1
         ncombi= 16
         nsimul= 42
         prmnmx(1,1)= .02
         prmnmx(2,1)= .08
         prmnmx(1,2)= .081
         prmnmx(2,2)= .25
         prmnmx(1,3)= .02
         prmnmx(2,3)= .08
         rlrntz= 1.
         skip_step3= .true.
      end if

c     99999999999999999999999999999999999999999999999999999999999999999999999
c     Complete settings for *-9.  Generalized Voigt functions.
c     NPOWER(JMETAB) = SIFWSD
C     FRACT_POWER_SD(JPOWER,JMETAB) = SIFWMN(JGAU,JMETAB) (when SIPPM=999)
c     99999999999999999999999999999999999999999999999999999999999999999999999
      if (sptype(:7) .eq. 'lipid-9'   .or.
     1    sptype(:7) .eq. 'liver-9'   .or.
     2    sptype(:8) .eq. 'breast-9')  then
         is_sptype= .true.
         atth2o= 1.

         chcom2(1)= 'L16+L09+L13'
         chcom2(2)= 'Lip16+Lip13'
         chcom2(3)= 'L28+L23+L21'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip21'
         chcom2(8)= 'Lip23'
         chcom2(9)= 'Lip28'
         chcom2(10)= 'Lip53+Lip52'
         chcom2(11)= 'Water'
         chcom2(12)= 'Cho'
         chcom2(13)= 'Lip43'
         chcom2(14)= 'Lip41'
         chcom2(15)= 'Glycg'
         chcom2(16)= 'Lip53'
         chcom2(17)= 'Lip52'

         chcomb(1)= 'L16+L13a+L13b+L09'
         chcomb(2)= 'L16+L13a+L13b'
         chcomb(3)= 'L28+L23+L21'
         chcomb(4)= 'L13a+L13b'
         chcomb(5)= 'L09'
         chcomb(6)= 'L16'
         chcomb(7)= 'L21'
         chcomb(8)= 'L23'
         chcomb(9)= 'L28'
         chcomb(10)= 'L53+L52'
         chcomb(11)= 'W1+W2+W3'
         chcomb(12)= 'Cho'
         chcomb(13)= 'L43'
         chcomb(14)= 'L41'
         chcomb(15)= 'glycg1+glycg2+glycg3+glycg4'
         chcomb(16)= 'L53'
         chcomb(17)= 'L52'

         chgrsh(1)= 'L13'
         chgrsh(2)= 'L5'
         chgrsh(3)= 'L4'

         chrato(1)= 'L41/L43 = 1. +- .02'

         chsimu(1)= 'L53 @ 5.31+-.02 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.4 AMP=-1.
     2                @999. FWHM=.08 AMP=-1.'
         chsimu(2)= 'L52 @ 5.21+-.02 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.4 AMP=-1.
     2                @999. FWHM=.08 AMP=-1.'

         if (roomt) then
            chsimu(3)= 'W1 @4.83+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
            chsimu(4)= 'W2 @4.82+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
            chsimu(5)= 'W3 @4.81+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(3)= 'W1 @4.65+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
            chsimu(4)= 'W2 @4.66+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
            chsimu(5)= 'W3 @4.67+-.02 FWHM=.02<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.'
         end if

         chsimu(6)= 'L43 @ 4.25+-.01 FWHM=.05<9.+-1. AMP=1.
     1                @999. FWHM=.4 AMP=-1.'
         chsimu(7)= 'L41 @ 4.05+-.01 FWHM=.05<9.+-1. AMP=1.
     1                @999. FWHM=.4 AMP=-1.'

         chsimu(8)= 'L28 @ 2.75+-.015 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'
         chsimu(9)= 'L23 @ 2.23+-.015 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'
         chsimu(10)= 'L21 @ 2.02+-.015 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'
         chsimu(11)= 'L16 @ 1.58+-.015 FWHM=.05<9.+-1. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'

         chsimu(12)= 'L13a @ 1.28+-.02 FWHM=.02<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'
         chsimu(13)= 'L13b @ 1.27+-.02 FWHM=.02<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'

         chsimu(14)= 'L09 @ .87+-.02 FWHM=.05<9.+-2. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.1 AMP=-1.'

         chsimu(15)= 'Cho @ 3.2+-.01 FWHM=.02<9.+-2. AMP=1.
     1                @999. FWHM=.25 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'

         chsimu(16)= 'glycg1 @ 3.82+-.01 FWHM=.1<9.+-1. AMP=1.
     1                @999. FWHM=.25 AMP=-1.'
         chsimu(17)= 'glycg2 @ 3.75+-.01 FWHM=.1<9.+-1. AMP=1.
     1                @999. FWHM=.25 AMP=-1.'
         chsimu(18)= 'glycg3 @ 3.68+-.01 FWHM=.1<9.+-1. AMP=1.
     1                @999. FWHM=.25 AMP=-1.'
         chsimu(19)= 'glycg4 @ 3.61+-.01 FWHM=.1<9.+-1. AMP=1.
     1                @999. FWHM=.25 AMP=-1.'

         chuse1(1)= 'L13a'
         dkntmn(2)= 1.867
         endpha= .true.
         fwhmmn= .02
         fwhmmx= .36
         fwhmst= .2
         gauss_rt2= .false.
         idgppm= 2
         imethd= 3
         incsmx= 1
         isdbol= 5
         mpower= 2
         mrepha(1)= 1
         nbas_ccf= 0
         ncombi= 17
         ngrsh= 3
         nobasi= .true.
         nratio = 0
         if (ppmst .ge. 4.25) nratio = 1
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 19
         onlyco= .true.
         ppmref(1,2)= 1.28
         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if

         namrel= 'L16+L09+L13'
         if (dowatr   .and.   ppmst .ge. 5.) then
            chuse1(2)= 'W1'
            if (.not.(dows   .or.   doecc)) namrel= 'Water'
            nuse1= 2
         else
            nuse1= 1
            if (ppmst .ge. 5.) then
               ppmgap(1,1)= 5.
               ppmgap(2,1)= 3.6
               ppmsep(1)= 4.3
            else
               rfwhcc= 0.
               rrt2mq= .25
            end if
         end if

         prmnmx(1,1)= .05
         prmnmx(2,1)= .9
         ralinc= 1.4
         rpowmq= .1
         sdgrsh(1)= .003
         sdgrsh(2)= .008
         sdgrsh(3)= .008
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     ======================================================================
c     Complete settings for *-10.  Similar to *-6, but scales SIFW* 
c        using FWHMST from Prel (as *-7).  As *-7, but with no 
c        linewidth regularization.
c     ======================================================================
      if (sptype(:8) .eq. 'lipid-10'   .or.
     1    sptype(:8) .eq. 'liver-10'   .or.
     2    sptype(:9) .eq. 'breast-10')  then
         is_sptype= .true.
         atth2o= 1.

         chcomb(1)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e+L09a+L09b+L
     109c+L09d'
         chcomb(2)= 'L16a+L16b+L16c+L13a+L13b+L13c+L13d+L13e'
         chcomb(3)= 'L28a+L28b+L28c+L23a+L23b+L23c+L21a+L21b+L21c'
         chcomb(4)= 'L13e+L13d+L13c+L13b+L13a'
         chcomb(5)= 'L09d+L09c+L09b+L09a'
         chcomb(6)= 'L16a+L16b+L16c'
         chcomb(7)= 'L21a+L21b+L21c'
         chcomb(8)= 'L23a+L23b+L23c'
         chcomb(9)= 'L28a+L28b+L28c'
         chcomb(10)= 'L53a+L53b+L53c+L52a+L52b+L52c'
         chcomb(11)= 'W1+W2+W3+W4+W5+W6'
         chcomb(12)= 'cho1+cho2+cho3'
         chcomb(13)= 'L43a+L43b'
         chcomb(14)= 'L41a+L41b'
         chcomb(15)= 'glycg1+glycg2+glycg3+glycg4'
         chcomb(16)= 'L53a+L53b+L53c'
         chcomb(17)= 'L52a+L52b+L52c'

         chcom2(1)= 'L16+L09+L13'
         chcom2(2)= 'Lip16+Lip13'
         chcom2(3)= 'L28+L23+L21'
         chcom2(4)= 'Lip13'
         chcom2(5)= 'Lip09'
         chcom2(6)= 'Lip16'
         chcom2(7)= 'Lip21'
         chcom2(8)= 'Lip23'
         chcom2(9)= 'Lip28'
         chcom2(10)= 'Lip53+Lip52'
         chcom2(11)= 'Water'
         chcom2(12)= 'Cho'
         chcom2(13)= 'Lip43'
         chcom2(14)= 'Lip41'
         chcom2(15)= 'Glycg'
         chcom2(16)= 'Lip53'
         chcom2(17)= 'Lip52'

         chgrsh(1)= 'L28'
         chgrsh(2)= 'L23'
         chgrsh(3)= 'L21'
         chgrsh(4)= 'L13'
         chgrsh(5)= 'L09'
         chgrsh(6)= 'cho'
         chgrsh(7)= 'L16'
         chgrsh(8)= 'L43'
         chgrsh(9)= 'L41'
         chgrsh(10)= 'L4'
         chgrsh(11)= 'W'

         chnot2(1)= 'W1pre'
         chnot2(2)= 'L13pre'

         chrato(1)= 'L41a/L43a = 1. +- .02'
         chrato(2)= 'L41b/L43b = 1. +- .02'

         chsimu(1)= 'L53a @ 5.31+-.02 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(2)= 'L53b @ 5.31+-.02 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(3)= 'L53c @ 5.31+-.02 FWHM=1.4<1.6+-.0381 AMP=1.'

         chsimu(4)= 'L52a @ 5.21+-.02 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(5)= 'L52b @ 5.21+-.02 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(6)= 'L52c @ 5.21+-.02 FWHM=1.4<1.6+-.0381 AMP=1.'

         if (roomt) then
            chsimu(7)= 'W1pre @ 4.83+-.03 FWHM=.03<.04+-.005 AMP=1.'
            chsimu(8)= 'W1 @ 4.83+-.03 FWHM=.2<.3+-.025 AMP=1.'
            chsimu(9)= 'W2 @ 4.83+-.03 FWHM=.4<.6+-.0667 AMP=1.'
            chsimu(10)= 'W3 @ 4.83+-.03 FWHM=1.<1.2+-.0667 AMP=1.'
            chsimu(11)= 'W4 @ 4.83+-.03 FWHM=1.6<1.8+-.0667 AMP=1.'
            chsimu(12)= 'W5 @ 4.83+-.03 FWHM=2.2<2.4+-.0667 AMP=1.'
            chsimu(13)= 'W6 @ 4.83+-.03 FWHM=2.8<3.+-.0381 AMP=1.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(7)= 'W1pre @ 4.65+-.03 FWHM=.03<.04+-.005 AMP=1.'
            chsimu(8)= 'W1 @ 4.65+-.03 FWHM=.2<.3+-.025 AMP=1.'
            chsimu(9)= 'W2 @ 4.65+-.03 FWHM=.4<.6+-.0667 AMP=1.'
            chsimu(10)= 'W3 @ 4.65+-.03 FWHM=1.<1.2+-.0667 AMP=1.'
            chsimu(11)= 'W4 @ 4.65+-.03 FWHM=1.6<1.8+-.0667 AMP=1.'
            chsimu(12)= 'W5 @ 4.65+-.03 FWHM=2.2<2.4+-.0667 AMP=1.'
            chsimu(13)= 'W6 @ 4.65+-.03 FWHM=2.8<3.+-.0381 AMP=1.'
         end if

         chsimu(14)= 'L43a @ 4.26+-.006 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(15)= 'L43b @ 4.26+-.006 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(16)= 'L41a @ 4.06+-.006 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(17)= 'L41b @ 4.06+-.006 FWHM=.8<1.+-.0381 AMP=1.'

         chsimu(18)= 'L28a @ 2.75+-.015 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(19)= 'L28b @ 2.75+-.015 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(20)= 'L28c @ 2.75+-.015 FWHM=1.4<1.6+-.0381 AMP=1.'

         chsimu(21)= 'L23a @ 2.23+-.015 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(22)= 'L23b @ 2.23+-.015 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(23)= 'L23c @ 2.23+-.015 FWHM=1.4<1.6+-.0381 AMP=1.'

         chsimu(24)= 'L21a @ 2.02+-.015 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(25)= 'L21b @ 2.02+-.015 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(26)= 'L21c @ 2.02+-.015 FWHM=1.4<1.6+-.0381 AMP=1.'

         chsimu(27)= 'L16a @ 1.58+-.015 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(28)= 'L16b @ 1.58+-.015 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(29)= 'L16c @ 1.58+-.015 FWHM=1.4<1.6+-.0667 AMP=1.'

         chsimu(30)= 'L13pre @ 1.28+-.02 FWHM=.03<.04+-.005 AMP=1.'
         chsimu(31)= 'L13a @ 1.28+-.02 FWHM=.2<.3+-.025 AMP=1.'
         chsimu(32)= 'L13b @ 1.28+-.02 FWHM=.4<.6+-.0667 AMP=1.'
         chsimu(33)= 'L13c @ 1.28+-.02 FWHM=1.<1.2+-.0667 AMP=1.'
         chsimu(34)= 'L13d @ 1.28+-.02 FWHM=1.6<1.8+-.0667 AMP=1.'
         chsimu(35)= 'L13e @ 1.28+-.02 FWHM=2.2<2.4+-.0381 AMP=1.'

         chsimu(36)= 'L09a @ .87+-.02 FWHM=.2<.4+-.0667 AMP=1.'
         chsimu(37)= 'L09b @ .87+-.02 FWHM=.8<1.+-.0667 AMP=1.'
         chsimu(38)= 'L09c @ .87+-.02 FWHM=1.4<1.6+-.0667 AMP=1.'
         chsimu(39)= 'L09d @ .87+-.02 FWHM=2.<2.2+-.0381 AMP=1.'
c        -------------------------------------------------------------------
c        Cho can be much sharper than Lip13.  Use only the sharpest models.
c        -------------------------------------------------------------------
         chsimu(40)= 'cho1 @ 3.2+-.01 FWHM=.2<.3+-.025 AMP=1.'
         chsimu(41)= 'cho2 @ 3.2+-.01 FWHM=.4<.6+-.0667 AMP=1.'
         chsimu(42)= 'cho3 @ 3.2+-.01 FWHM=1.<1.2+-.0381 AMP=1.'
c        -------------------------------------------------------------------
c        Only one model at each ppm.
c        -------------------------------------------------------------------
         chsimu(43)= 'glycg1 @ 3.82+-.01 FWHM=.5<.75+-.0625 AMP=1.'
         chsimu(44)= 'glycg2 @ 3.75+-.01 FWHM=.5<.75+-.0625 AMP=1.'
         chsimu(45)= 'glycg3 @ 3.68+-.01 FWHM=.5<.75+-.0625 AMP=1.'
         chsimu(46)= 'glycg4 @ 3.61+-.01 FWHM=.5<.75+-.0625 AMP=1.'

         chuse1(1)= 'L13pre'
         dkntmn(2)= 1.867
         endpha= .true.
         fwhmmn= .1
         fwhmmx= .4
         fwhmst= .15
         gauss_rt2= .false.
         idgppm= 2
         isdbol= 5
         mrepha(1)= 1
         nbas_ccf= 0
         ncombi= 17
         ngrsh= 11
         nnot2= 2
         nobasi= .true.
         nratio = 0
         if (ppmst .ge. 4.28) nratio = 2
         nrefpk(2)= 1
         nsimul= 46
         onlyco= .true.
         ppmref(1,2)= 1.28
         if (ppmst .lt. -9998.) ppmst = 3.6
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         namrel= 'L16+L09+L13'
         if (dowatr   .and.   ppmst .ge. 5.) then
            chuse1(2)= 'W1pre'
            if (.not.(dows   .or.   doecc)) namrel= 'Water'
            nuse1= 2
         else
            nuse1= 1
            if (ppmst .ge. 5.) then
               ppmgap(1,1)= 5.
               ppmgap(2,1)= 3.6
               ppmsep(1)= 4.3
            else
               rfwhcc= 0.
               rrt2mq= .25
            end if
         end if
         prmnmx(1,1)= .05
         prmnmx(2,1)= .9
         rbackg(1)= 22.2
         scafwh= .true.
         sdgrsh(1)= .003
         sdgrsh(2)= .003
         sdgrsh(3)= .003
         sdgrsh(4)= .003
         sdgrsh(5)= .003
         sdgrsh(6)= .002
         sdgrsh(7)= .003
         sdgrsh(8)= .003
         sdgrsh(9)= .003
         sdgrsh(10)= .005
         sdgrsh(11)= .01
         sdsmoo(3)=.03
c        ---------------------------------------------------------------------
c        SHIFMN(1) (PPMREF=4.65): searches above 4.65; practically no limit.
c        SHIFMN(2) (PPMREF=1.28): do not hit Cho instead of Lip13.
c        SHIFMX(1): do not hit Cho instead of water (but usually water>Cho)
c        SHIFMX(2): searches below 1.28; practically no limit.
c        Only go halfway, since these limits implicitly assume that there is 
c           no referencing shift.
c        ---------------------------------------------------------------------
         shifmn(1)= -4.
         shifmn(2)= -1.
         shifmx(1)= 1.
         shifmx(2)= 2.5
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     ======================================================================
c     liver-11.  Modify liver-8 to again include (more extensive) Glycg.
c     ======================================================================
      if (sptype(:8) .eq. 'liver-11')  then
         chcom2(17)= 'Glycg'
         chcomb(17)= 'glycg1+glycg2+glycg3+glycg4+glycg5+glycg6+glycg7'
         chsimu(43)= 'glycg1 @ 3.89+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(44)= 'glycg2 @ 3.82+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(45)= 'glycg3 @ 3.75+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(46)= 'glycg4 @ 3.68+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(47)= 'glycg5 @ 3.61+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(48)= 'glycg6 @ 3.54+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         chsimu(49)= 'glycg7 @ 3.47+-.01 FWHM=.5<.75+-.0625 AMP=1.
     1                @999.0 FWHM=1.0 AMP=0.1'
         ncombi= 17
         nsimul= 49
         if (.not.dowatr   .and.   ppmst .ge. 5.) ppmgap(2,1) = 4.
         ppmst = 4.
      end if
