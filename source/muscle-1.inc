      if (sptype(:8) .eq. 'muscle-1') then
         is_sptype= .true.
         badref = .true.
         chcomb(1) = 'I13d+I13c+I13b+I13a'
         chcomb(2) = 'E15d+E15c+E15b+E15a'
         chcomb(3) = 'Cr2+Cr1'
         chcomb(4) = 'tau5+tau4+tau3+tau2+tau1'
         chcomb(5) = 'cho4+cho3+cho2+cho1'
         chcomb(6) = 'cr28e+cr28d+cr28c+cr28b+cr28a'
         chcomb(7) = 'I09'
         chcomb(8) = 'E11'
         chcomb(9) = 'I21'
         chcomb(10) = 'E23'

         chcom2(1) = 'IMCL13'
         chcom2(2) = 'EMCL15'
         chcom2(3) = 'Cr'
         chcom2(4) = 'tau'
         chcom2(5) = 'cho'
         chcom2(6) = 'cr28'
         chcom2(7) = 'I09'
         chcom2(8) = 'E11'
         chcom2(9) = 'I21'
         chcom2(10) = 'E23'

         chgrsh(1) = 'E15'
         chgrsh(2) = 'I13'
         chgrsh(3) = 'Cr'
         chnot1(1) = 'tau1'
         chnot1(2) = 'tau2'
         chnot1(3) = 'tau3'
         chnot1(4) = 'tau4'
         chnot1(5) = 'tau5'
         chnot2(1) = 'E15e'
         chnot2(2) = 'E15f'
         chnot2(3) = 'I13e'
         chnot2(4) = 'I13f'
         chsimu(1)= 'E15a @1.52+-.01 FWHM=.19<.2+-.006 AMP=2.'
         chsimu(2)= 'E15b @1.52+-.01 FWHM=.16<.17+-.006 AMP=2.'
         chsimu(3)= 'E15c @1.52+-.01 FWHM=.13<.14+-.006 AMP=2.'
         chsimu(4)= 'E15d @1.52+-.01 FWHM=.10<.11+-.006 AMP=2.'
         chsimu(5)= 'E15e @1.54+-.01 FWHM=.16<.17+-.03 AMP=2.'
         chsimu(6)= 'E15f @1.50+-.01 FWHM=.16<.17+-.03 AMP=2.'
         chsimu(7)= 'I13a @1.27+-.01 FWHM=.13<.14+-.004 AMP=2.'
         chsimu(8)= 'I13b @1.27+-.01 FWHM=.10<.11+-.004 AMP=2.'
         chsimu(9)= 'I13c @1.27+-.01 FWHM=.07<.08+-.004 AMP=2.'
         chsimu(10)= 'I13d @1.27+-.01 FWHM=.04<.05+-.004AMP=2.'
         chsimu(11)= 'I13e @1.29+-.01 FWHM=.10<.11+-.02 AMP=2.'
         chsimu(12)= 'I13f @1.25+-.01 FWHM=.10<.11+-.02 AMP=2.'
         chsimu(13)= 'E11 @1.07+-.02 FWHM=.07<.08+-.02 AMP=3.'
         chsimu(14)= 'I09 @.89+-.02 FWHM=.07<.08+-.02 AMP=3.'
         chsimu(15)= 'E23 @2.33+-.03 FWHM=.25<.3+-.03 AMP=1. 
     1                 @1.8 FWHM=.74 AMP=1.'
         chsimu(16)= 'I21 @2.15+-.03 FWHM=.23<.28+-.03 AMP=1.1111  
     1                @1.6 FWHM=.74 AMP=.8889'
         chsimu(17)= 'cr28a @2.90+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(18)= 'cr28b @2.86+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(19)= 'cr28c @2.82+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(20)= 'cr28d @2.78+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(21)= 'cr28e @2.74+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(22)= 'cho1@3.24+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(23)= 'cho2@3.21+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(24)= 'cho3@3.18+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(25)= 'cho4@3.15+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(26)= 'Cr1@3.03+-.004 FWHM=.025<.03+-.003  AMP=3.'
         chsimu(27)= 'Cr2@3.03+-.004 FWHM=.03<.04+-.005  AMP=3.'
         chsimu(28)= 'tau1 @3.58+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(29)= 'tau2 @3.52+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(30)= 'tau3 @3.46+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(31)= 'tau4 @3.40+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(32)= 'tau5 @3.34+-.006 FWHM=.04<.06+-.01 AMP=3.'
         endpha= .true.
         incsmx = 1
         namrel= 'Cr'
         ncombi= 10
         ngrsh= 3
         nnot1= 5
         nnot2= 4
         nobasi = .true.
         nratio= 0
         nrefpk(2)= 1
         nsidmn = 1
         nsidmx = 1
         nsimul= 32
         onlyco= .true.
         ppmend= -1.0
         ppmref(1,2)= 3.03
         ppmst= 3.8
         sdgrsh(1) = .002
         sdgrsh(2) = .002
         sdgrsh(3) = .001
         vitro = .true.
      end if

c     22222222222222222222222222222222222222222222222222222222222222222222222
c     muscle-2
c     22222222222222222222222222222222222222222222222222222222222222222222222
      if (sptype(:8) .eq. 'muscle-2') then
         is_sptype= .true.
         badref = .true.
         chcom2(1)= 'IMCL13'
         chcom2(2)= 'EMCL15'
         chcom2(3)= 'Cr'
         chcom2(4)= 'tau'
         chcom2(5)= 'cho'
         chcom2(6)= 'cr28'
         chcom2(7)= 'I09'
         chcom2(8)= 'E11'
         chcom2(9)= 'I21'
         chcom2(10)= 'E23'

         chcomb(1)= 'I13e+I13d+I13c+I13b+I13a'
         chcomb(2)= 'E15g+E15f+E15e+E15d+E15c+E15b+E15a'
         chcomb(3)= 'Cr2+Cr1'
         chcomb(4)= 'tau5+tau4+tau3+tau2+tau1'
         chcomb(5)= 'cho4+cho3+cho2+cho1'
         chcomb(6)= 'cr28e+cr28d+cr28c+cr28b+cr28a'
         chcomb(7)= 'I09a+I09b'
         chcomb(8)= 'E11'
         chcomb(9)= 'I21a+I21b+I21c+I21d+I21e'
         chcomb(10)= 'E23a+E23b+E23c+E23d+E23e'

         chgrsh(1)= 'E15'
         chgrsh(2)= 'I13'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'E23'
         chgrsh(5)= 'I21'
         chgrsh(6)= 'I09'

         chnot1(1)= 'tau1'
         chnot1(2)= 'tau2'
         chnot1(3)= 'tau3'
         chnot1(4)= 'tau4'
         chnot1(5)= 'tau5'

         chnot2(1)= 'E15h'
         chnot2(2)= 'E15i'
         chnot2(3)= 'I13f'
         chnot2(4)= 'I13g'

         chsimu(1)= 'E15a @1.52+-.03 FWHM=.03<.04+-.006 AMP=2.'
         chsimu(2)= 'E15b @1.52+-.03 FWHM=.06<.07+-.006 AMP=2.'
         chsimu(3)= 'E15c @1.52+-.03 FWHM=.09<.10+-.006 AMP=2.'
         chsimu(4)= 'E15d @1.52+-.03 FWHM=.12<.13+-.006 AMP=2.'
         chsimu(5)= 'E15e @1.52+-.03 FWHM=.15<.16+-.006 AMP=2.'
         chsimu(6)= 'E15f @1.52+-.03 FWHM=.18<.19+-.006 AMP=2.'
         chsimu(7)= 'E15g @1.52+-.03 FWHM=.21<.22+-.006 AMP=2.'
         chsimu(8)= 'E15h @1.54+-.03 FWHM=.16<.18+-.03 AMP=2.'
         chsimu(9)= 'E15i @1.50+-.03 FWHM=.16<.18+-.03 AMP=2.'
         chsimu(10)= 'I13a @1.27+-.03 FWHM=.03<.04+-.006 AMP=2.'
         chsimu(11)= 'I13b @1.27+-.03 FWHM=.06<.07+-.006 AMP=2.'
         chsimu(12)= 'I13c @1.27+-.03 FWHM=.09<.10+-.006 AMP=2.'
         chsimu(13)= 'I13d @1.27+-.03 FWHM=.12<.13+-.006 AMP=2.'
         chsimu(14)= 'I13e @1.27+-.03 FWHM=.15<.16+-.006 AMP=2.'
         chsimu(15)= 'I13f @1.29+-.03 FWHM=.09<.11+-.02 AMP=2.'
         chsimu(16)= 'I13g @1.25+-.03 FWHM=.09<.11+-.02 AMP=2.'
         chsimu(17)= 'E11 @1.07+-.02 FWHM=.07<.08+-.02 AMP=3.'
         chsimu(18)= 'I09a @.89+-.02 FWHM=.05<.06+-.01 AMP=3.'
         chsimu(19)= 'I09b @.89+-.02 FWHM=.08<.09+-.01 AMP=3.'
         chsimu(20)= 'E23a @2.33+-.03 FWHM=.04<.06+-.02 AMP=1. 
     1                @1.8 FWHM=.59 AMP=1.'
         chsimu(21)= 'E23b @2.33+-.03 FWHM=.10<.12+-.02 AMP=1. 
     1                @1.8 FWHM=.59 AMP=1.'
         chsimu(22)= 'E23c @2.33+-.03 FWHM=.16<.18+-.02 AMP=1. 
     1                @1.8 FWHM=.69 AMP=1.'
         chsimu(23)= 'E23d @2.33+-.03 FWHM=.22<.24+-.02 AMP=1. 
     1                @1.8 FWHM=.79 AMP=1.'
         chsimu(24)= 'E23e @2.33+-.03 FWHM=.28<.30+-.02 AMP=1. 
     1                @1.8 FWHM=.79 AMP=1.'
         chsimu(25)= 'I21a @2.15+-.03 FWHM=.04<.06+-.02 AMP=1.1111 
     1                @1.6 FWHM=.61 AMP=.8889'
         chsimu(26)= 'I21b @2.15+-.03 FWHM=.10<.12+-.02 AMP=1.1111 
     1                @1.6 FWHM=.61 AMP=.8889'
         chsimu(27)= 'I21c @2.15+-.03 FWHM=.16<.18+-.02 AMP=1.1111 
     1                @1.6 FWHM=.71 AMP=.8889'
         chsimu(28)= 'I21d @2.15+-.03 FWHM=.22<.24+-.02 AMP=1.1111 
     1                @1.6 FWHM=.71 AMP=.8889'
         chsimu(29)= 'I21e @2.15+-.03 FWHM=.28<.30+-.02 AMP=1.1111 
     1                @1.6 FWHM=.71 AMP=.8889'
         chsimu(30)= 'cr28a @2.90+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(31)= 'cr28b @2.86+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(32)= 'cr28c @2.82+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(33)= 'cr28d @2.78+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(34)= 'cr28e @2.74+-.004 FWHM=.03<.04+-.004 AMP=3.'
         chsimu(35)= 'cho1 @3.24+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(36)= 'cho2 @3.21+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(37)= 'cho3 @3.18+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(38)= 'cho4 @3.15+-.003 FWHM=.025<.03+-.003  AMP=9.'
         chsimu(39)= 'Cr1 @3.03+-.02 FWHM=.025<.035+-.01  AMP=3.'
         chsimu(40)= 'Cr2 @3.03+-.02 FWHM=.045<.055+-.01  AMP=3.'
         chsimu(41)= 'tau1 @3.58+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(42)= 'tau2 @3.52+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(43)= 'tau3 @3.46+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(44)= 'tau4 @3.40+-.006 FWHM=.04<.06+-.01 AMP=3.'
         chsimu(45)= 'tau5 @3.34+-.006 FWHM=.04<.06+-.01 AMP=3.'
         dkntmn(2)= 2.4
         idgppm= 2
         incsmx= 1
         namrel= 'Cr'
         ncombi= 10
         ngrsh= 6
         nnot1= 5
         nnot2= 4
         nobasi= .true.
         nratio= 0
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 45
         onlyco= .true.
         ppmend= -1.0
         ppmref(1,2)= 3.03
         ppmst= 3.8
         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         ralinc= 1.4
         rrt2mq= .25
         sdgrsh(1)= .002
         sdgrsh(2)= .002
         sdgrsh(3)= .001
         sdgrsh(4)= .004
         sdgrsh(5)= .004
         sdgrsh(6)= .002
         useany= .true.
         vitro= .true.
      end if

c     33333333333333333333333333333333333333333333333333333333333333333333333
c     muscle-3.  Estimate E55, I53 & Water.
c     33333333333333333333333333333333333333333333333333333333333333333333333
      if (sptype(:8) .eq. 'muscle-3') then
         is_sptype= .true.
         atth2o= 1.
         badref = .true.
         chcom2(1)= 'IMCL13'
         chcom2(2)= 'EMCL15'
         chcom2(3)= 'Cr'
         chcom2(4)= 'E55+I53'
         chcom2(5)= 'I09'
         chcom2(6)= 'E11'
         chcom2(7)= 'cho'
         chcom2(8)= 'I21'
         chcom2(9)= 'E23'
         chcom2(10)= 'tau'
         chcom2(11)= 'cr39'
         chcom2(12)= 'cr28'
         chcom2(13)= 'E55'
         chcom2(14)= 'I53'
         chcom2(15)= 'Water'

         chcomb(1)= 'I13d+I13c+I13b+I13a'
         chcomb(2)= 'E15e+E15d+E15c+E15b+E15a'
         chcomb(3)= 'Cr2+Cr1'
         chcomb(4)= 'E55a+E55b+E55c+E55d+I53a+I53b+I53c+I53d'
         chcomb(5)= 'I09a+I09b+I09c'
         chcomb(6)= 'E11a+E11b+E11c'
         chcomb(7)= 'cho4+cho3+cho2+cho1'
         chcomb(8)= 'I21a+I21b+I21c+I21d+I21e'
         chcomb(9)= 'E23a+E23b+E23c+E23d+E23e'
         chcomb(10)= 'tau5+tau4+tau3+tau2+tau1'
         chcomb(11)= 'cr39a+cr39b'
         chcomb(12)= 'cr28e+cr28d+cr28c+cr28b+cr28a'
         chcomb(13)= 'E55a+E55b+E55c+E55d'
         chcomb(14)= 'I53a+I53b+I53c+I53d'
         chcomb(15)= 'W1+W2+W3+W4+W5+W6+W7'

         chgrsh(1)= 'E15'
         chgrsh(2)= 'I13'
         chgrsh(3)= 'Cr'
         chgrsh(4)= 'I09'
         chgrsh(5)= 'E11'
         chgrsh(6)= 'I21'
         chgrsh(7)= 'E23'
         chgrsh(8)= 'cr39'
         chgrsh(9)= 'E55'
         chgrsh(10)= 'I53'
         chgrsh(11)= 'W'

         chnot1(1)= 'tau1'
         chnot1(2)= 'tau2'
         chnot1(3)= 'tau3'
         chnot1(4)= 'tau4'
         chnot1(5)= 'tau5'
         chnot1(6)= 'cho1'
         chnot1(7)= 'cho2'
         chnot1(8)= 'cho3'
         chnot1(9)= 'cho4'
         chnot1(10)= 'cr28a'
         chnot1(11)= 'cr28b'
         chnot1(12)= 'cr28c'
         chnot1(13)= 'cr28d'
         chnot1(14)= 'cr28e'
         chnot1(15)= 'E11a'
         chnot1(16)= 'E11b'
         chnot1(17)= 'E11c'
         chnot1(18)= 'I09a'
         chnot1(19)= 'I09b'
         chnot1(20)= 'I09c'

         chnot2(1)= 'E15f'
         chnot2(2)= 'E15g'
         chnot2(3)= 'I13e'
         chnot2(4)= 'I13f'

         chsimu(1)= 'E15a @1.52+-.01 FWHM=.03<.04+-.005 AMP=1.'
         chsimu(2)= 'E15b @1.52+-.01 FWHM=.06<.07+-.005 AMP=1.'
         chsimu(3)= 'E15c @1.52+-.01 FWHM=.13<.14+-.005 AMP=1.'
         chsimu(4)= 'E15d @1.52+-.01 FWHM=.20<.21+-.005 AMP=1.'
         chsimu(5)= 'E15e @1.52+-.01 FWHM=.27<.28+-.005 AMP=1.'
         chsimu(6)= 'E15f @1.54+-.01 FWHM=.16<.18+-.02 AMP=1.'
         chsimu(7)= 'E15g @1.50+-.01 FWHM=.16<.18+-.02 AMP=1.'

         chsimu(8)= 'I13a @1.27+-.01 FWHM=.03<.04+-.005 AMP=1.'
         chsimu(9)= 'I13b @1.27+-.01 FWHM=.06<.07+-.005 AMP=1.'
         chsimu(10)= 'I13c @1.27+-.01 FWHM=.13<.14+-.005 AMP=1.'
         chsimu(11)= 'I13d @1.27+-.01 FWHM=.20<.21+-.005 AMP=1.'
         chsimu(12)= 'I13e @1.29+-.01 FWHM=.09<.11+-.02 AMP=1.'
         chsimu(13)= 'I13f @1.25+-.01 FWHM=.09<.11+-.02 AMP=1.'

         chsimu(14)= 'E11a @1.09+-.01 FWHM=.05<.06+-.005 AMP=1.'
         chsimu(15)= 'E11b @1.09+-.01 FWHM=.11<.12+-.005 AMP=1.'
         chsimu(16)= 'E11c @1.09+-.01 FWHM=.17<.18+-.005 AMP=1.'

         chsimu(17)= 'I09a @.89+-.01 FWHM=.05<.06+-.005 AMP=1.'
         chsimu(18)= 'I09b @.89+-.01 FWHM=.12<.13+-.005 AMP=1.'
         chsimu(19)= 'I09c @.89+-.01 FWHM=.19<.20+-.005 AMP=1.'

         chsimu(20)= 'E23a @2.33+-.03 FWHM=.04<.05+-.005 AMP=.5 
     1                @1.8 FWHM=.59 AMP=.5'
         chsimu(21)= 'E23b @2.33+-.03 FWHM=.11<.12+-.005 AMP=.5 
     1                @1.8 FWHM=.59 AMP=.5'
         chsimu(22)= 'E23c @2.33+-.03 FWHM=.18<.19+-.005 AMP=.5 
     1                @1.8 FWHM=.69 AMP=.5'
         chsimu(23)= 'E23d @2.33+-.03 FWHM=.25<.26+-.005 AMP=.5 
     1                @1.8 FWHM=.79 AMP=.5'
         chsimu(24)= 'E23e @2.33+-.03 FWHM=.32<.33+-.005 AMP=.5 
     1                @1.8 FWHM=.79 AMP=.5'

         chsimu(25)= 'I21a @2.15+-.03 FWHM=.04<.05+-.005 AMP=.556
     1                @1.8 FWHM=.59 AMP=.444'
         chsimu(26)= 'I21b @2.15+-.03 FWHM=.11<.12+-.005 AMP=.556
     1                @1.8 FWHM=.59 AMP=.444'
         chsimu(27)= 'I21c @2.15+-.03 FWHM=.18<.19+-.005 AMP=.556
     1                @1.8 FWHM=.69 AMP=.444'
         chsimu(28)= 'I21d @2.15+-.03 FWHM=.25<.26+-.005 AMP=.556
     1                @1.8 FWHM=.79 AMP=.444'
         chsimu(29)= 'I21e @2.15+-.03 FWHM=.32<.33+-.005 AMP=.556
     1                @1.8 FWHM=.79 AMP=.444'

         chsimu(30)= 'cr39a @3.93+-.02 FWHM=.025<.035+-.005  AMP=1.'
         chsimu(31)= 'cr39b @3.93+-.02 FWHM=.055<.065+-.005  AMP=1.'

         chsimu(32)= 'cr28a @2.90+-.004 FWHM=.03<.04+-.004 AMP=1.'
         chsimu(33)= 'cr28b @2.86+-.004 FWHM=.03<.04+-.004 AMP=1.'
         chsimu(34)= 'cr28c @2.82+-.004 FWHM=.03<.04+-.004 AMP=1.'
         chsimu(35)= 'cr28d @2.78+-.004 FWHM=.03<.04+-.004 AMP=1.'
         chsimu(36)= 'cr28e @2.74+-.004 FWHM=.03<.04+-.004 AMP=1.'

         chsimu(37)= 'cho1 @3.24+-.003 FWHM=.025<.03+-.003  AMP=1.'
         chsimu(38)= 'cho2 @3.21+-.003 FWHM=.025<.03+-.003  AMP=1.'
         chsimu(39)= 'cho3 @3.18+-.003 FWHM=.025<.03+-.003  AMP=1.'
         chsimu(40)= 'cho4 @3.15+-.003 FWHM=.025<.03+-.003  AMP=1.'

         chsimu(41)= 'Cr1 @3.03+-.02 FWHM=.025<.035+-.005  AMP=1.'
         chsimu(42)= 'Cr2 @3.03+-.02 FWHM=.055<.065+-.005  AMP=1.'

         chsimu(43)= 'tau1 @3.58+-.006 FWHM=.04<.06+-.01 AMP=1.'
         chsimu(44)= 'tau2 @3.52+-.006 FWHM=.04<.06+-.01 AMP=1.'
         chsimu(45)= 'tau3 @3.46+-.006 FWHM=.04<.06+-.01 AMP=1.'
         chsimu(46)= 'tau4 @3.40+-.006 FWHM=.04<.06+-.01 AMP=1.'
         chsimu(47)= 'tau5 @3.34+-.006 FWHM=.04<.06+-.01 AMP=1.'

         chsimu(48)= 'E55a @5.52+-.01 FWHM=.06<.07+-.005 AMP=1.'
         chsimu(49)= 'E55b @5.52+-.01 FWHM=.13<.14+-.005 AMP=1.'
         chsimu(50)= 'E55c @5.52+-.01 FWHM=.20<.21+-.005 AMP=1.'
         chsimu(51)= 'E55d @5.52+-.01 FWHM=.27<.28+-.005 AMP=1.'

         chsimu(52)= 'I53a @5.30+-.01 FWHM=.06<.07+-.005 AMP=1.'
         chsimu(53)= 'I53b @5.30+-.01 FWHM=.13<.14+-.005 AMP=1.'
         chsimu(54)= 'I53c @5.30+-.01 FWHM=.20<.21+-.005 AMP=1.'
         chsimu(55)= 'I53d @5.30+-.01 FWHM=.27<.28+-.005 AMP=1.'

         if (roomt) then
            chsimu(56)= 'W1 @4.83+-.02 FWHM=.03<.04+-.005 AMP=1.'
            chsimu(57)= 'W2 @4.83+-.02 FWHM=.06<.07+-.005 AMP=1.'
            chsimu(58)= 'W3 @4.83+-.02 FWHM=.14<.15+-.005 AMP=1.'
            chsimu(59)= 'W4 @4.83+-.02 FWHM=.22<.23+-.005 AMP=1.'
            chsimu(60)= 'W5 @4.83+-.02 FWHM=.30<.31+-.005 AMP=1.'
            chsimu(61)= 'W6 @4.83+-.02 FWHM=.38<.39+-.005 AMP=1.'
            chsimu(62)= 'W7 @4.83+-.02 FWHM=.46<.47+-.005 AMP=1.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(56)= 'W1 @4.65+-.02 FWHM=.03<.04+-.005 AMP=1.'
            chsimu(57)= 'W2 @4.65+-.02 FWHM=.06<.07+-.005 AMP=1.'
            chsimu(58)= 'W3 @4.65+-.02 FWHM=.14<.15+-.005 AMP=1.'
            chsimu(59)= 'W4 @4.65+-.02 FWHM=.22<.23+-.005 AMP=1.'
            chsimu(60)= 'W5 @4.65+-.02 FWHM=.30<.31+-.005 AMP=1.'
            chsimu(61)= 'W6 @4.65+-.02 FWHM=.38<.39+-.005 AMP=1.'
            chsimu(62)= 'W7 @4.65+-.02 FWHM=.46<.47+-.005 AMP=1.'
         end if

         dkntmn(2)= 2.4
         endpha= .true.
         fwhmmx= .21
         fwhmst= .2
         gauss_rt2= .false.
         idgppm= 2
         incsmx= 1
         isdbol= 5
         mrepha(1)= 1
         namrel= 'Cr'
         ncombi= 15
         ngrsh= 11
         nnot1= 16
         nnot2= 4
         nobasi= .true.
         nratio= 0
         nrefpk(2)= 1
         nsidmn= 1
         nsidmx= 1
         nsimul= 62
         onlyco= .true.
         ppmref(1,2)= 3.03
         if (ppmst .lt. -9998.) ppmst = 3.8
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         prmnmx(1,1)= .5
         prmnmx(2,1)= .95
         ralinc= 1.4
         rrt2mq= .25
         sdgrsh(1)= .002
         sdgrsh(2)= .002
         sdgrsh(3)= .002
         sdgrsh(4)= .002
         sdgrsh(5)= .002
         sdgrsh(6)= .004
         sdgrsh(7)= .004
         sdgrsh(8)= .002
         sdgrsh(9)= .002
         sdgrsh(10)= .002
         sdgrsh(11)= .007
         useany= .true.
         vitro= .true.
         wconc= .5
      end if

c     44444444444444444444444444444444444444444444444444444444444444444444444
c     muscle-4 or muscle-5.  Estimate E55, I53 & Water.
c     NPOWER(JMETAB) = SIFWSD
C     FRACT_POWER_SD(JPOWER,JMETAB) = SIFWMN(JGAU,JMETAB) (when SIPPM=999)
c     muscle-5 uses ASYMLP=T
c     44444444444444444444444444444444444444444444444444444444444444444444444
      if (sptype(:8) .eq. 'muscle-4'   .or.
     1    sptype(:8) .eq. 'muscle-5') then
         is_sptype= .true.
         if (sptype(:8) .eq. 'muscle-5') asymlp= .true.
         atth2o= 1.

         chcom2(1)= 'IMCL13'
         chcom2(2)= 'EMCL15'
         chcom2(3)= 'Cr'
         chcom2(4)= 'I09'
         chcom2(5)= 'E11'
         chcom2(6)= 'Water'
         chcom2(7)= 'I53'
         chcom2(8)= 'E55'
         chcom2(9)= 'I53+E55'
         chcom2(10)= 'cho'
         chcom2(11)= 'I21'
         chcom2(12)= 'E23'
         chcom2(13)= 'tau'
         chcom2(14)= 'cr39'
         chcom2(15)= 'cr28'

         chcomb(1)= 'xIMC13+xIM13a'
         chcomb(2)= 'xEMC15+xEM15a'
         chcomb(3)= 'Cr'
         chcomb(4)= 'xI09'
         chcomb(5)= 'xE11'
         chcomb(6)= 'Water+Water2'
         chcomb(7)= 'yI53'
         chcomb(8)= 'yE55'
         chcomb(9)= 'yI53+yE55'
         chcomb(10)= 'cho4+cho3+cho2+cho1'
         chcomb(11)= 'I21'
         chcomb(12)= 'E23'
         chcomb(13)= 'tau5+tau4+tau3+tau2+tau1'
         chcomb(14)= 'cr39'
         chcomb(15)= 'cr28e+cr28d+cr28c+cr28b+cr28a'

         chgrsh(1)= 'x'
         chgrsh(2)= 'y'

         chsimu(1)= 'xIMC13 @1.27+-.025 FWHM=.03<9.+-2. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'
         chsimu(2)= 'xIM13a @1.28+-.025 FWHM=.03<9.+-2. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'
         chsimu(3)= 'xEMC15 @1.52+-.025 FWHM=.03<9.+-2. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'
         chsimu(4)= 'xEM15a @1.51+-.025 FWHM=.03<9.+-2. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'

         if (.not.asymlp) then
            chcomb(1)= 'xIMC13'
            chcomb(2)= 'xEMC15'
            chomit(1)= 'xIM13a'
            chomit(2)= 'xEM15a'
            chsimu(1)= 'xIMC13 @1.27+-.025 FWHM=.03<9.+-3. AMP=1.
     1                  @999. FWHM=.5 AMP=-1.
     2                  @999. FWHM=.05 AMP=-1.
     3                  @999. FWHM=.1 AMP=-1.'
            chsimu(3)= 'xEMC15 @1.52+-.025 FWHM=.03<9.+-3. AMP=1.
     1                  @999. FWHM=.5 AMP=-1.
     2                  @999. FWHM=.05 AMP=-1.
     3                  @999. FWHM=.1 AMP=-1.'
            nomit= 2
         end if

         chsimu(5)= 'Cr @3.03+-.02 FWHM=.02<9.+-2. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'
         chsimu(6)= 'xI09 @.89+-.025 FWHM=.05<9.+-1. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'
         chsimu(7)= 'xE11 @1.09+-.025 FWHM=.05<9.+-1. AMP=1.
     1               @999. FWHM=.5 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.
     3               @999. FWHM=.1 AMP=-1.'

         chsimu(8)= 'cho1 @3.24+-.003 FWHM=.03<9.+-1. AMP=1.
     1               @999. FWHM=.05 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.'
         chsimu(9)= 'cho2 @3.21+-.003 FWHM=.03<9.+-1. AMP=1.
     1               @999. FWHM=.05 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.'
         chsimu(10)= 'cho3 @3.18+-.003 FWHM=.03<9.+-1. AMP=1.
     1               @999. FWHM=.05 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.'
         chsimu(11)= 'cho4 @3.15+-.003 FWHM=.03<9.+-1. AMP=1.
     1               @999. FWHM=.05 AMP=-1.
     2               @999. FWHM=.05 AMP=-1.'

         chsimu(12)= 'I21 @2.15+-.03 FWHM=.04<9.+-1. AMP=.556
     1                @1.6 FWHM=.50 AMP=.444
     2                @999. FWHM=1. AMP=-1.
     3                @999. FWHM=.05 AMP=-1.
     4                @999. FWHM=.1 AMP=-1.'
         chsimu(13)= 'E23 @2.33+-.03 FWHM=.04<9.+-1. AMP=.556
     1                @1.8 FWHM=.50 AMP=.444
     2                @999. FWHM=1. AMP=-1.
     3                @999. FWHM=.05 AMP=-1.
     4                @999. FWHM=.1 AMP=-1.'

         chsimu(14)= 'tau1 @3.58+-.006 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(15)= 'tau2 @3.52+-.006 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(16)= 'tau3 @3.46+-.006 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(17)= 'tau4 @3.40+-.006 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(18)= 'tau5 @3.34+-.006 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'

         chsimu(19)= 'cr39 @3.93+-.02 FWHM=.025<9.+-1. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.
     3                @999. FWHM=.1 AMP=-1.'

         chsimu(20)= 'cr28a @2.90+-.004 FWHM=.04<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(21)= 'cr28b @2.86+-.004 FWHM=.04<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(22)= 'cr28c @2.82+-.004 FWHM=.04<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(23)= 'cr28d @2.78+-.004 FWHM=.04<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'
         chsimu(24)= 'cr28e @2.74+-.004 FWHM=.04<9.+-1. AMP=1.
     1                @999. FWHM=.05 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.'

         chsimu(25)= 'yE55 @5.52+-.01 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.
     3                @999. FWHM=.1 AMP=-1.'
         chsimu(26)= 'yI53 @5.30+-.01 FWHM=.06<9.+-1. AMP=1.
     1                @999. FWHM=.5 AMP=-1.
     2                @999. FWHM=.05 AMP=-1.
     3                @999. FWHM=.1 AMP=-1.'
         if (roomt) then
            chsimu(27)= 'Water @4.83+-.02 FWHM=.03<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.
     3                   @999. FWHM=.1 AMP=-1.'
            chsimu(28)= 'Water2 @4.82+-.02 FWHM=.03<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.
     3                   @999. FWHM=.1 AMP=-1.'
            ppmcen= 4.83
            ppmref(1,1) = 4.83
         else
            chsimu(27)= 'Water @4.65+-.02 FWHM=.03<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.
     3                   @999. FWHM=.1 AMP=-1.'
            chsimu(28)= 'Water2 @4.66+-.02 FWHM=.03<9.+-2. AMP=1.
     1                   @999. FWHM=.5 AMP=-1.
     2                   @999. FWHM=.1 AMP=-1.
     3                   @999. FWHM=.1 AMP=-1.'
         end if
         chuse1(1)= 'xIMC13'
         chuse1(2)= 'xEMC15'
         chuse1(3)= 'Cr'
         chuse1(4)= 'I21'
         chuse1(5)= 'E23'
         chuse1(6)= 'cr39'
         chuse1(7)= 'Water'
         chuse1(8)= 'yE55'
         chuse1(9)= 'yI53'
         dkntmn(2)= 2.4
         endpha= .true.
         fixshf= .true.
         fwhmmx= .21
         fwhmst= .2
         gauss_rt2= .false.
         idgppm= 2
         imethd= 3
         incsmx= 1
         isdbol= 5
         mpower= 2
         if (.not.asymlp) mpower= 3
         mrepha(1)= 1
         namrel= 'Cr'
         ncombi= 15
         ngrsh= 2
         nobasi= .true.
         nratio= 0
         nrefpk(2)= 4
         nsidmn= 1
         nsidmx= 1
         nsimul= 28
         onlyco= .true.
         ppmref(1,2)= 3.2
         ppmref(2,2)= 3.03
         ppmref(3,2)= 1.52
         ppmref(4,2)= 1.27
         if (ppmst .lt. -9998.) ppmst = 3.8
         if (ppmend .gt. 9998.) ppmend = -2.
         if (ppmst .lt. 4.) then
            nuse1 = 5
         else if (ppmst .lt. 6.) then
            nuse1 = 6
         else
            nuse1 = 9
         end if
         if (ppmst - ppmend .ge. 6.4) then
            nsubtk= 5
            xstep= 1.
         end if
         prmnmx(1,1)= .05
         prmnmx(2,1)= .9
         ralinc= 1.4
         rpowmq= .1
         sdgrsh(1)= .004
         sdgrsh(2)= .004
         useany= .true.
         vitro= .true.
         wconc= .5
      end if
