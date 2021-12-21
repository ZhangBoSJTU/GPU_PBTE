!  PBTE, a solver for the Phonon Boltzmann Transport Equation
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
!  along with GPU_PBTE.If not, see <http://www.gnu.org/licenses/>.

!  By using the software, you agree to fully comply with the terms 
!  and conditions of the HPC SDK Software License Agreement. see 
!  <https://docs.nvidia.com/hpc-sdk/eula/index.html>.

!  This file is part of GPU_PBTE.
!  Suroutine to read interatomic force constant HFC.txt CFC.txt and QFC.txt.

module IFC
   integer:: nHarm
   integer,dimension(:,:),allocatable::HFCIndex
   real(8),dimension(:),allocatable::HFC
   integer::nTriple
   integer,dimension(:,:),allocatable::CFCIndex
   real(8),dimension(:,:,:,:),allocatable::CFC
   real(8),dimension(:,:,:,:),allocatable::CFCm
   real(8),dimension(:,:,:),allocatable::CFCIndexR
   integer,dimension(:,:),allocatable::CFCIndextail
   integer::nQuart
   integer,dimension(:,:),allocatable::QFCIndex,QFCI
   real(8),dimension(:,:),allocatable::QFCR
   real(8),dimension(:,:,:,:,:),allocatable::QFC
   real(8),dimension(:,:,:),allocatable::QFCIndexR
   real(8),dimension(:,:,:,:,:),allocatable::QFCm
   complex(8),dimension(:,:,:,:,:),allocatable::QFC_complex
   complex(8),dimension(:,:),allocatable::QFC_complex_sum
   integer,dimension(:,:),allocatable::QFCIndextail
   contains
   
   subroutine Read_HFC
      implicit none
      integer::i
      open(10,file='../input/HFC.dat',status='unknown')
      read(10,*) nHarm
      allocate(HFCIndex(nHarm,7),HFC(nHarm))
      do i=1,nHarm
         read(10,*) HFCIndex(i,1:7),HFC(i)
      enddo
      close(10)
      return
   end subroutine Read_HFC
   
   subroutine Read_CFC
      use Crystal,only:a1,a2,a3,mass
      implicit none
      integer::i,j,k,l
      open(10,file='../input/CFC.dat',status='unknown')
      read(10,*) nTriple
      allocate(CFC(nTriple,3,3,3),CFCIndex(nTriple,9))
      allocate(CFCm(nTriple,3,3,3),CFCIndexR(nTriple,2,3),CFCIndextail(nTriple,3))
      do i=1,nTriple
         read(10,*) CFCIndex(i,1:9)
         read(10,*) CFC(i,1,1,1),CFC(i,1,1,2),CFC(i,1,1,3), &
            & CFC(i,1,2,1),CFC(i,1,2,2),CFC(i,1,2,3), &
            & CFC(i,1,3,1),CFC(i,1,3,2),CFC(i,1,3,3), &
            & CFC(i,2,1,1),CFC(i,2,1,2),CFC(i,2,1,3), &
            & CFC(i,2,2,1),CFC(i,2,2,2),CFC(i,2,2,3), &
            & CFC(i,2,3,1),CFC(i,2,3,2),CFC(i,2,3,3), &
            & CFC(i,3,1,1),CFC(i,3,1,2),CFC(i,3,1,3), &
            & CFC(i,3,2,1),CFC(i,3,2,2),CFC(i,3,2,3), &
            & CFC(i,3,3,1),CFC(i,3,3,2),CFC(i,3,3,3)
      enddo
      do i=1,nTriple
         CFCIndexR(i,1,:)=CFCIndex(i,1)*a1+CFCIndex(i,2)*a2+CFCIndex(i,3)*a3
         CFCIndexR(i,2,:)=CFCIndex(i,4)*a1+CFCIndex(i,5)*a2+CFCIndex(i,6)*a3
         CFCIndextail(i,1)=(CFCIndex(i,7)-1)*3
         CFCIndextail(i,2)=(CFCIndex(i,8)-1)*3
         CFCIndextail(i,3)=(CFCIndex(i,9)-1)*3  
         do j=1,3
         do k=1,3
         do l=1,3
            CFCm(i,j,k,l)=CFC(i,j,k,l)/sqrt(mass(CFCIndex(i,7)) &
               &  *mass(CFCIndex(i,8))*mass(CFCIndex(i,9)))
         enddo
         enddo
         enddo
      enddo      
      close(10)
      return
   end subroutine Read_CFC
   
   subroutine Read_QFC
      use crystal,only:a1,a2,a3,mass
      implicit none
      integer::i,j,k,l,m
      integer::i2,j1,j2,j3,j4
      real(8)::coeff3
      open(10,file='../input/QFC.dat',status='unknown')
      read(10,*) nQuart
      allocate(QFC(nQuart,3,3,3,3),QFCIndex(nQuart,13))
      allocate(QFCIndexR(nQuart,3,3),QFCm(nQuart,3,3,3,3),QFC_complex(nQuart,3,3,3,3),&
         &  QFC_complex_sum(9,nQuart),QFCIndextail(nQuart,4))
      do i=1,nQuart
         read(10,*) QFCIndex(i,1:13)
         do j=1,3
         do k=1,3
         do l=1,3
         do m=1,3
            read(10,*) QFC(i,j,k,l,m)
         enddo
         enddo
         enddo
         enddo
      enddo
      do i2=1,nQuart
         QFCIndextail(i2,1)=(QFCIndex(i2,10)-1)*3
         QFCIndextail(i2,2)=(QFCIndex(i2,11)-1)*3
         QFCIndextail(i2,3)=(QFCIndex(i2,12)-1)*3
         QFCIndextail(i2,4)=(QFCIndex(i2,13)-1)*3
         QFCIndexR(i2,1,:)=QFCIndex(i2,1)*a1+QFCIndex(i2,2)*a2+QFCIndex(i2,3)*a3
         QFCIndexR(i2,2,:)=QFCIndex(i2,4)*a1+QFCIndex(i2,5)*a2+QFCIndex(i2,6)*a3
         QFCIndexR(i2,3,:)=QFCIndex(i2,7)*a1+QFCIndex(i2,8)*a2+QFCIndex(i2,9)*a3
         coeff3=1d0/sqrt(mass(QFCIndex(i2,10))*mass(QFCIndex(i2,11))&
            &  *mass(QFCIndex(i2,12))*mass(QFCIndex(i2,13)))
         do j1=1,3  
         do j2=1,3
         do j3=1,3
         do j4=1,3
            QFCm(i2,j1,j2,j3,j4)=QFC(i2,j1,j2,j3,j4)*coeff3
         enddo
         enddo
         enddo
         enddo
      enddo
      close(10)
      return
   end subroutine Read_QFC
      
end module IFC