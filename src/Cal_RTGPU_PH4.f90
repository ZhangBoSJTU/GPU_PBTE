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
!  Calculate four-phonon scattering processes.

subroutine Cal_RTGPU_PH4(iq0_used)
   use mesh
   use Crystal
   use SmearingFunctions
   use symm
   use TmpSpace,only: AAA,QQQ,QQQ4ph
   use Quantum
   use IFC
   use cudafor
   use Accrelate_V4
   implicit none
   
   integer::iq0_used,iq0,iq1,iq2,iq3,is0,is1,is2,is3,nB01,nB11,nB21,nB31,isscattcount
   integer::isDone,G(3),q12List(NQtot,2)
   integer::isscatt(NQtot*nBand*nBand,3)
   real(8)::f33,f00
   real(8)::qq3(3)
   real(8),dimension(:,:,:),allocatable::AA
   complex(8)::ee0(nBand),ee3(nBand)
   logical::isScattList(NQtot,nBand,nBand)
   type(dim3)::dimGrid, dimBlock
   integer::tPB=32
   integer::istat,nDevices,ierrSync,ierrAsync,ierr
   integer,device::isscattd(NQtot*nBand*nBand,3)
   integer,dimension(:,:),device::q12Listd(NQtot,2)
   integer,dimension(:,:),device::QFCIndextaild(nQuart,4)
   real(8),dimension(:),device::qq3d(3)
   real(8),dimension(:),device::TempListd(nTemp)
   real(8),dimension(:),device::c1d(3),c2d(3),c3d(3)
   real(8),dimension(:,:),device::freqd(NQtot,nBand)
   real(8),dimension(:,:),device::Qpointd(nQx*nQy*nQz,3)
   real(8),dimension(:,:,:),device::QFCIndexRd(nQuart,3,3)
   real(8),dimension(:,:,:),device::AAd(nBand*NQtot,nBand,nTemp)   
   real(8),dimension(:,:,:,:),device::sigmad(NQtot,nBand,nBand,2)
   real(8),dimension(:,:),device::QQQ4phd(nBand,nTemp)
   real(8),dimension(:,:),device::grpvelxd(NQtot,nBand),grpvelyd(NQtot,nBand),grpvelzd(NQtot,nBand)
   complex(8),dimension(:,:),device::QFC_complex_sumd(9,nQuart)   
   complex*16,dimension(:,:,:),device::evd(NQtot,nBand,nBand)
   logical,dimension(:,:,:),device::isScattList1d(NQtot,nBand,nBand),isScattList2d(NQtot,nBand,nBand),isScattList3d(NQtot,nBand,nBand),isScattListd(NQtot,nBand,nBand)
   
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
   QFCIndextaild(1:nQuart,1:4)=QFCIndextail(1:nQuart,1:4)
   QFCIndexRd(1:nQuart,1:3,1:3)=QFCIndexR(1:nQuart,1:3,1:3)
   iq0=List_used(iq0_used)
   call IS_DONE(iq0_used,isDone)
   if (isDone .eq. 1) return
   if (.not. RTA) then
      allocate(AA(nBand*NQtot,nBand,nTemp))
      AA=0d0
      AAd(1:(nBand*NQtot),1:nBand,1:nTemp)=AA(1:(nBand*NQtot),1:nBand,1:nTemp)
   endif
   nB01=1
   if (iq0 .eq. iGamma) nB01=4
   do iq3=1,NQtot
      write(*,*) iq0_used,iq3,QQQ4ph(1,1)
      q12List=0
      qq3=-Qpoint(iq3,:)
      qq3d=qq3
      do iq1=1,NQtot
         q12List(iq1,1)=iq1
         G=ind1(iq0,:)+ind1(iq1,:)-ind1(iq3,:)
         q12List(iq1,2)=ind3(modulo(G(1),nQx),modulo(G(2),nQy),modulo(G(3),nQz))
      enddo
      q12Listd=q12List(1:NQtot,1:2)
      do is0=nB01,nBand
         f00=freq(iq0,is0)
         ee0=ev(iq0,:,is0)
         nB31=1
         if (iq3 .eq. iGamma) nB31=4
         do is3=nB31,nBand
            f33=freq(iq3,is3)
            if ( .not. (iq0 .eq. iq3 .and. abs(f00-f33) .le. 1e-5)) then
               ee3=conjg(ev(iq3,:,is3))
               call Cal_V4_pre(ee0,ee3)
               QFC_complex_sumd(1:9,1:nQuart)=QFC_complex_sum(1:9,1:nQuart)
               
               dimGrid = dim3(ceiling(real(NQtot)/tPB),nBand,nBand)
               dimBlock = dim3(tPB, 1, 1 )
               call system_clock(int1)
               call Cal_IsScat_V4<<<dimGrid,dimBlock>>>(NQtot,nBand,nQx,nQy,nQz,tpSmearing,f00,f33,scoeff4,&
                  q12Listd,freqd,grpvelxd,grpvelyd,grpvelzd,c1d,c2d,c3d,isScattList1d,isScattList2d,isScattList3d,isScattListd,sigmad)
               isscattcount=0
               isScattList(1:NQtot,1:nBand,1:nBand)=isScattListd(1:NQtot,1:nBand,1:nBand)
               do iq1=1,NQtot
                  iq2=q12List(iq1,2)
                  nB11=1
                  nB21=1
                  if ( iq1 .eq. iGamma ) nB11=4
                  if ( iq2 .eq. iGamma ) nB21=4
                  do is1=nB11,nBand
                  do is2=nB21,nBand 
                     if (isScattList(iq1,is1,is2)) then
                        isscattcount=isscattcount+1
                        isscatt(isscattcount,:)=(/iq1,is1,is2/)
                     endif
                  enddo
                  enddo
               enddo
               isscattd=isscatt(1:NQtot*nBand*nBand,1:3)
               dimGrid = dim3(ceiling(real(isscattcount)/tPB),1, 1 )
               dimBlock = dim3(tPB, 1, 1 )
               call Cal_V4<<<dimGrid,dimBlock>>>(NQtot,nQuart,nBand,nQx,nQy,nQz,isscattcount,tpSmearing,is0,iq3,is3,f00,f33,nTemp,RTA,&
                  isscattd,evd,QFCIndexRd,QFC_complex_sumd,QFCIndextaild,q12Listd,qq3d,Qpointd,freqd,isScattList1d,isScattList2d,isScattList3d,AAd,QQQ4phd,TempListd,sigmad)
               QQQ4ph(1:nBand,1:nTemp)=QQQ4phd(1:nBand,1:nTemp)
            endif
         enddo
      enddo
   enddo
   QQQ(1:nBand,1:nTemp) = QQQ(1:nBand,1:nTemp)+QQQ4ph(1:nBand,1:nTemp)
   if (.not. RTA) then
      AA(1:nBand*NQtot,1:nBand,1:nTemp)=AAd(1:nBand*NQtot,1:nBand,1:nTemp)
      AAA=AAA+AA
      deallocate(AA)
   endif
end subroutine Cal_RTGPU_PH4