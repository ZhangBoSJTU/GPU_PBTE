!  GPU_PBTE, a solver for the Phonon Boltzmann Transport Equation
!  Copyright (C) 2019-2021 Zhang Bo <zhangbo.sjtu@qq.com>
!  Copyright (C) 2019-2021 Xiaokun Gu <xiaokun.gu@sjtu.edu.cn>

!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with GPU_PBTE.  If not, see <http://www.gnu.org/licenses/>.

!  By using the software, you agree to fully comply with the terms 
!  and conditions of the HPC SDK Software License Agreement. see 
!  <https://docs.nvidia.com/hpc-sdk/eula/index.html>.

!  This file is part of GPU_PBTE.
!  Calculate the joint density of states

subroutine Cal_DOS_3PH
   use mesh
   use crystal
   use symm
   use smearingFunctions
   implicit none

   integer::i,is0,iq0,is1,iq1,nfreq
   real(8)::freq_hi,freq_lo,ff,ff0,ff1,sigma1,sigma2,V0(3),V1(3),dV(3)
   real(8),allocatable::freqplot(:),PDOS1(:),PDOS2(:)
   freq_lo=0d0
   freq_hi=1.d0
   nfreq=2500
   allocate(freqplot(0:nfreq),PDOS1(0:nfreq),PDOS2(0:nfreq))   
   PDOS1(:)=0d0
   PDOS2(:)=0d0
   do i=0,nfreq
      ff=freq_lo+(freq_hi-freq_lo)*dble(i)/dble(nfreq)
      freqplot(i)=ff
      do iq0=1,NQtot
      do is0=1,nBand
         ff0=freq(iq0,is0)
         V0=(/grpvelx(iq0,is0),grpvely(iq0,is0),grpvelz(iq0,is0)/)
         do iq1=1,NQtot
         do is1=1,nBand
            ff1=freq(iq1,is1)
            V1=(/grpvelx(iq1,is1),grpvely(iq1,is1),grpvelz(iq1,is1)/)
            dV=V1-V0
            sigma1 = sig3(dV)
            dV=V1+V0
            sigma2 = sig3(dV)         
            if (abs((ff-ff0+ff1)/sigma1) .le. 2.0) then 
               PDOS1(i)=PDOS1(i)+1/dble(NQtot*NQtot)*dG(ff-ff0+ff1,sigma1)
            endif
            if (abs((ff-ff0-ff1)/sigma2) .le. 2.0) then
               PDOS2(i)=PDOS2(i)+1/dble(NQtot*NQtot)*dG(ff-ff0-ff1,sigma2)
            endif
         enddo
         enddo
      enddo
      enddo
   enddo
   open(10,file='dos3.dat')
   do i=0,nfreq
      write(10,*) freqplot(i),PDOS1(i),PDOS2(i)
   enddo
   return
end subroutine Cal_DOS_3PH