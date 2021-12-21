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
!  Calculate the density of states.

subroutine Cal_DOS
   use Mesh
   use Crystal
   use symm
   use SmearingFunctions
   implicit none
   
   integer::i,is,iq,nfreq
   real(8)::freq_hi,freq_lo,ff,ff0,sigma,dV(3)
   real(8),allocatable::freqplot(:),PDOS(:)
   
   freq_lo=0d0
   freq_hi=0.5d0
   nfreq=2500
   allocate(freqplot(0:nfreq),PDOS(0:nfreq))
   PDOS(i)=0d0
   do i=0,nfreq
      ff=freq_lo+(freq_hi-freq_lo)*dble(i)/dble(nfreq)
      freqplot(i)=ff
      do iq=1,NQtot
         do is=1,nBand
            ff0=freq(iq,is)
            dV=(/grpvelx(iq,is),grpvely(iq,is),grpvelz(iq,is)/)
            sigma = sig3(dV)
            if (abs((-ff+ff0)/sigma) .le. 2.0) then 
               PDOS(i)=PDOS(i)+1/dble(NQtot)*dG(-ff+ff0,sigma)
            endif
         enddo
      enddo
   enddo
   open(10,file='dos.dat')
   do i=0,nfreq
      write(10,*) freqplot(i),PDOS(i)
   enddo
   return
end subroutine Cal_DOS
  
   
   
   