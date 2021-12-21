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
!  Generate a grid for dividing the Brillouin zone.

subroutine Gen_QPoint
use mesh
use Crystal
use symm

implicit none
integer::i,j,k
real(8)::qq(3)
dx=1d0/nQx
dy=1d0/nQy
dz=1d0/nQz
allocate(ind1(nQx*nQy*nQz,3),Qpoint(nQx*nQy*nQz,3))
allocate(ind3(0:nQx-1,0:nQy-1,0:nQz-1))

NQtot=0
do i=0,nQx-1
   do j=0,nQy-1
      do k=0,nQz-1
         qq=i*dx*b1+j*dy*b2+k*dz*b3
         NQtot=NQtot+1
         if (i.eq.0 .and. j .eq.0 .and. k.eq.0) iGamma=NQtot
         ind1(NQtot,:)=(/modulo(i,nQx),modulo(j,nQy),modulo(k,nQz)/)
         ind3(modulo(i,nQx),modulo(j,nQy),modulo(k,nQz))=NQtot
         Qpoint(NQtot,:)=qq
      enddo
   enddo
enddo
call Find_Symm
allocate(Nequi_ibz(NQtot),List_ibz(NQtot),IdEqui_ibz(NQtot), &
   &  ALLEquiList_ibz(Nsymm_rot_ibz,NQtot),TypeofSymmetry_ibz(Nsymm_rot_ibz,NQtot))
call Find_EQP(N0_ibz,Nequi_ibz,List_ibz,IdEqui_ibz,AllEquiList_ibz, & 
   &  TypeofSymmetry_ibz,Nsymm_rot_ibz,qrotations_ibz,NQtot)
if (.not. RTA) then
   allocate(Nequi_8fbz(NQtot),List_8fbz(NQtot),IdEqui_8fbz(NQtot), &
      &  ALLEquiList_8fbz(Nsymm_rot_8fbz,NQtot),TypeofSymmetry_8fbz(Nsymm_rot_8fbz,NQtot))
   call Find_EQP(N0_8fbz,Nequi_8fbz,List_8fbz,IdEqui_8fbz,AllEquiList_8fbz, & 
      &   TypeofSymmetry_8fbz,Nsymm_rot_8fbz,qrotations_8fbz,NQtot)
endif

write(*,*) NQtot
allocate(freq(NQtot,nBand),ev(NQtot,nBand,nBand))
allocate(grpvelx(NQtot,nBand),grpvely(NQtot,nBand),grpvelz(NQtot,nBand),Q(NQtot,nBand))

return
end subroutine Gen_QPoint