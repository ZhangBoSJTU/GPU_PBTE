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

subroutine Cal_DOS_4PH(iq0_used)
   use mesh
   use crystal
   use smearingFunctions
   use symm
   implicit none

   integer::iq0_used,iq0,iq1,iq2,iq3,is0,is1,is2,is3,nB01,nB11,nB21,nB31
   integer::G(3),q12List(NQtot,2)
   real(8)::f11,f22,f33,f00,ff0(nBand),ff3(nBand)
   real(8)::dgg,sigma,sigmamin
   real(8),allocatable,dimension(:)::GG
   real(8)::V1(3),V2(3),dV(3)
   character(50)::s1
   
   sigmamin=1d-8
   allocate(GG(nBand))
   GG=0d0
   iq0=List_used(iq0_used)
   ff0=freq(iq0,:)
   nB01=1

   if (iq0 .eq. iGamma) nB01=4
   do iq3=1,NQtot
      q12List=0
      ff3=freq(iq3,:)
      do iq1=1,NQtot
         G=ind1(iq0,:)+ind1(iq1,:)-ind1(iq3,:)
         G=(/modulo(G(1),nQx),modulo(G(2),nQy),modulo(G(3),nQz)/)
         iq2=ind3(G(1),G(2),G(3))
         q12List(iq1,1)=iq1
         q12List(iq1,2)=iq2
      enddo

      do is0=nB01,nBand
         f00=ff0(is0) 
         nB31=1
         if (iq3 .eq. iGamma) nB31=4
         do is3=nB31,nBand
            f33=ff3(is3)
            if ( .not. (iq0 .eq. iq3 .and. abs(f00-f33) .le. 1e-5)) then
               do iq1=1,NQtot
                  iq2=q12List(iq1,2)
                  nB11=1
                  nB21=1
                  if ( iq1 .eq. iGamma ) nB11=4
                  if ( iq2 .eq. iGamma)  nB21=4
                  do is1=nB11,nBand
                     f11=freq(iq1,is1)
                     V1=(/grpvelx(iq1,is1),grpvely(iq1,is1),grpvelz(iq1,is1)/)
                     do is2=nB21,nBand
                        f22=freq(iq2,is2)
                        V2=(/grpvelx(iq2,is2),grpvely(iq2,is2),grpvelz(iq2,is2)/)
                        dV=V1+V2
                        sigma = sig4(dV)
                        if (abs((f00+f11+f22-f33)/sigma) .le. 2.0 .and. sigma .ge. sigmamin) then
                           dgg=dG(f00+f11+f22-f33,sigma)
                           GG(is0)=GG(is0)+dgg
                        endif
                        if (abs((f00-f11-f22-f33)/sigma) .le. 2.0 .and. sigma .ge. sigmamin) then
                           dgg=dG(f00-f11-f22-f33,sigma)
                           GG(is0)=GG(is0)+dgg
                        endif
                        dV=V1-V2
                        sigma = sig4(dV)
                        if (abs((f00+f11-f22-f33)/sigma) .le. 2.0 .and. sigma .ge. sigmamin) then
                           dgg=dG(f00+f11-f22-f33,sigma)
                           GG(is0)=GG(is0)+dgg
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
   call system('mkdir dos4')
   write(s1,'(i7)') iq0
   open(10,file='dos4/'//trim(adjustl(s1))//'')
   do is0=1,nBand
      write(10,*) iq0,freq(iq0,is0),GG(is0)
   enddo
   close(10)
   deallocate(GG)
end subroutine Cal_DOS_4PH