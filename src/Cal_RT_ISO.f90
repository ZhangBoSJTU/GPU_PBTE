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
!  Calculate isotope scattering

subroutine Cal_RT_ISO(iq0_used)
use mesh
use crystal
use symm
use tmpSpace
use quantum
use smearingFunctions
implicit none

integer::iq0_used,iq0,iq1,j1,j2,iT
integer::isDone,G(3)
real(8)::f11,f00,Qiso,Qisot
real(8)::ff0(nBand)
real(8)::GG,sigma
real(8)::nn0,nn1
real(8)::qq0(3),qq1(3),dV(3)
real(8),dimension(:,:),allocatable::A
integer::is0,is1

iq0=List_used(iq0_used)
call IS_DONE(iq0_used,isDone)
if (isDone .eq. 1) return

if (.not. RTA) then
   allocate(A(nBand*nQtot,nTemp))
endif

qq0=Qpoint(iq0,:)
ff0=freq(iq0,:)
do is0=1,nBand
   if (.not. RTA) A=0
   f00=ff0(is0) 
   do iT=1,nTemp
      temp=TempList(iT)
      nn0=BE(f00,temp)
      do iq1=1,nQtot
         do is1=1,nBand
            f11=freq(iq1,is1) 
            nn1=BE(f11,temp)
            dV=(/grpvelx(iq1,is1),grpvely(iq1,is1),grpvelz(iq1,is1)/)
            sigma = sig3(dV)
            if (abs((-f11+f00)/sigma) .le. 2.0) then
               Qiso=0d0
               do j1=1,nBasis
                  Qisot=0d0
                  do j2=1,3
                     Qisot=Qisot+conjg(ev(iq0,(j1-1)*3+j2,is0))*ev(iq1,(j1-1)*3+j2,is1)
                  enddo
                  Qiso=Qiso+(abs(Qisot))**2d0*g2(j1)
               enddo
               GG=Qiso*PI/2d0/dble(nQtot)*f11*f00*(nn0*nn1+(nn1+nn0)/2d0)*dG(-f11+f00,sigma)
               if (.not. RTA) then
                  A((iq1-1)*nBand+is1,iT)=A((iq1-1)*nBand+is1,iT)-GG
               endif
               QQQ(is0,iT)=QQQ(is0,iT)+GG
               QQQiso(is0,iT)=QQQiso(is0,iT)+GG
            endif
         enddo
      enddo
      if (.not. RTA)  AAA(:,is0,iT)=AAA(:,is0,iT)+A(:,iT)
   enddo
enddo

if (.not. rta) then
   deallocate(A)
endif

end subroutine Cal_RT_ISO