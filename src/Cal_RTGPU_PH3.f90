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
!  Calculate three-phonon scattering processes.

subroutine Cal_RTGPU_PH3(iq0_used)
   use mesh
   use Crystal
   use SmearingFunctions
   use symm
   use TmpSpace,only: AAA,QQQ,QQQ3ph
   use Quantum
   use IFC
   use cudafor
   use Accrelate_V3
   implicit none

   integer::iq0_used,iq0,iq1,iq2,is0,is1,is2,i3,nB01,nB11,nB21,isscattcount
   integer::isDone,G(3),q12List(NQtot,2)
   integer::isscatt(NQtot*nBand*nBand,3)
   real(8)::f00
   real(8),dimension(:,:),allocatable::A
   complex*16::ee0(nBand)
   logical::isScattList(NQtot,nBand,nBand)
   type(dim3)::dimGrid, dimBlock
   integer::tPB=32
   integer::istat,nDevices,ierrSync,ierrAsync,ierr
   integer,dimension(:,:),device::isscattd(NQtot*nBand*nBand,3)
   integer,dimension(:,:),device::q12Listd(NQtot,2)
   integer,dimension(:,:),device::CFCIndextaild(nTriple,3)
   real(8),dimension(:),device::TempListd(nTemp)
   real(8),dimension(:),device::c1d(3),c2d(3),c3d(3)
   real(8),dimension(:,:),device::freqd(NQtot,nBand)
   real(8),dimension(:,:),device::Qpointd(nQx*nQy*nQz,3)
   real(8),dimension(:,:),device::Ad(nBand*NQtot,nTemp)
   real(8),dimension(:,:,:,:),device::CFCmd(nTriple,3,3,3)
   real(8),dimension(:,:,:),device::CFCIndexRd(nTriple,2,3)
   real(8),dimension(:,:,:,:),device::sigmad(NQtot,nBand,nBand,2)
   real(8),dimension(:,:),device::QQQ3phd(nBand,nTemp)
   real(8),dimension(:,:),device::grpvelxd(NQtot,nBand),grpvelyd(NQtot,nBand),grpvelzd(NQtot,nBand) 
   complex*16,dimension(:),device::ee0d(nBand)
   complex*16,dimension(:,:,:),device::evd(NQtot,nBand,nBand)
   logical,dimension(:,:,:),device::isScattList1d(NQtot,nBand,nBand),isScattList2d(NQtot,nBand,nBand),isScattListd(NQtot,nBand,nBand)

   integer:: int1, int2

   c1d=c1
   c2d=c2
   c3d=c3
   TempListd=TempList
   Qpointd=Qpoint(1:nQx*nQy*nQz,1:3)
   freqd=freq(1:NQtot,1:nBand)
   evd=ev(1:NQtot,1:nBand,1:nBand)
   grpvelxd=grpvelx(1:NQtot,1:nBand)
   grpvelyd=grpvely(1:NQtot,1:nBand)
   grpvelzd=grpvelz(1:NQtot,1:nBand)
   CFCIndextaild(1:nTriple,1:3) = CFCIndextail(1:nTriple,1:3)
   CFCmd(1:nTriple,1:3,1:3,1:3) = CFCm(1:nTriple,1:3,1:3,1:3)
   CFCIndexRd(1:nTriple,1:2,1:3) = CFCIndexR(1:nTriple,1:2,1:3)
   
   iq0=List_used(iq0_used)   
   call IS_DONE(iq0_used,isDone)
   if (isDone .eq. 1) return
   if (.not. RTA) then
      allocate(A(nBand*NQtot,nTemp))
   endif
   do iq1=1,NQtot
      G=-ind1(iq0,:)-ind1(iq1,:)
      G=(/modulo(G(1),nQx),modulo(G(2),nQy),modulo(G(3),nQz)/)
      iq2=ind3(G(1),G(2),G(3))
      q12List(iq1,1)=iq1
      q12List(iq1,2)=iq2
   enddo
   q12Listd=q12List
   nB01=1
   if (iq0 .eq. iGamma) nB01=4
   do is0=nB01,nBand
      if (.not. RTA) then 
         A=0d0
         Ad(1:(nBand*NQtot),1:nTemp)=A(1:(nBand*NQtot),1:nTemp)
      endif
      f00=freq(iq0,is0) 
      ee0=ev(iq0,:,is0)
      ee0d(1:nBand)=ee0(1:nBand)
      dimGrid = dim3(ceiling(real(NQtot)/tPB),nBand,nBand)
      dimBlock = dim3( tPB, 1, 1 )
      call system_clock(int1)
      call Cal_IsScat_V3<<<dimGrid,dimBlock>>>(NQtot,nBand,nQx,nQy,nQz,tpSmearing,scoeff,f00,&
         q12Listd,freqd,grpvelxd,grpvelyd,grpvelzd,c1d,c2d,c3d,isScattList1d,isScattList2d,isScattListd,sigmad)
      isscattcount=0
      isScattList(1:NQtot,1:nBand,1:nBand)=isScattListd(1:NQtot,1:nBand,1:nBand)
      do i3=1,NQtot
         iq1=q12List(i3,1)
         iq2=q12List(i3,2)
         nB11=1
         nB21=1
         if ( iq1 .eq. iGamma ) nB11=4
         if ( iq2 .eq. iGamma ) nB21=4
         do is1=nB11,nBand
         do is2=nB21,nBand
            if ( isScattList(i3,is1,is2)) then
               isscattcount=isscattcount+1
               isscatt(isscattcount,:)=(/i3,is1,is2/)
            endif
         enddo
         enddo
      enddo
      isscattd(1:NQtot*nBand*nBand,1:3)=isscatt(1:NQtot*nBand*nBand,1:3)
      dimGrid = dim3(ceiling(real(isscattcount)/tPB),1,1)
      dimBlock = dim3(tPB, 1, 1 )
      call Cal_V3<<<dimGrid,dimBlock>>>(NQtot,nTriple,nBand,nQx,nQy,nQz,isscattcount,nTemp,tpSmearing,RTA,f00,is0,&
         ee0d,isscattd,q12Listd,Qpointd,freqd,isScattList1d,isScattList2d,TempListd,sigmad,evd,CFCIndexRd,CFCmd,QQQ3phd,Ad,CFCIndextaild)
      QQQ3ph(1:nBand,1:nTemp)=QQQ3phd(1:nBand,1:nTemp)
      if (.not. RTA) then
         A(1:(nBand*NQtot),1:nTemp)=Ad(1:(nBand*NQtot),1:nTemp)
         AAA(:,is0,:)=AAA(:,is0,:)+A(:,:)
      endif
   enddo
   QQQ(1:nBand,1:nTemp)=QQQ(1:nBand,1:nTemp)+QQQ3ph(1:nBand,1:nTemp)
   if (.not. rta) then
      deallocate(A)
   endif
end subroutine Cal_RTGPU_PH3