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
!  Calculate the phonon dispersion.

module Disp
   implicit none
   contains
   subroutine Cal_Disp
      use mesh,only:Qpoint,ev,freq,NQtot
      use Crystal,only:nBand
      use IFC
      implicit none
      integer::i
      real(8)::qq(3),ff(nBand),ww(nBand)
      complex*16::eigenv(nBand,nBand),DM(nBand,nBand)
      do i=1,NQtot
         qq=Qpoint(i,:)
         call Cal_Disp_point(qq,ff,eigenv,ww,DM)
         freq(i,:)=ff
         ev(i,:,:)=eigenv
      enddo
      return
   end subroutine Cal_Disp
   
   subroutine Cal_Disp_point(qq,ff,eigenv,ww,DM)
      use Crystal,only:nBand
      use IFC
      implicit none
      integer::j
      real(8)::qq(3),ff(nBand),ww(nBand)
      complex*16::eigenv(nBand,nBand),DM(nBand,nBand)
      call Cal_DM(qq,DM,nBand)
      call Cdiagh(nBand,DM,ww,eigenv)
      do j=1,nBand
         ff(j)=sqrt(abs(ww(j)))*sign(1d0,ww(j))
      enddo
      return
   end subroutine Cal_Disp_point
   
   subroutine Cal_DM(qq,DM,nDim)
      use Crystal
      use IFC
      implicit none
      integer::nDim
      real(8),dimension(3)::qq
      complex*16::DM(nDim,nDim)
      integer::i,nx,ny,nz,nB1,nB2,ix1,ix2
      real(8)::f
      DM=0
      do i=1,nHarm
         nx=HFCIndex(i,1);
         ny=HFCIndex(i,2);
         nz=HFCIndex(i,3);
         nB1=HFCIndex(i,4);
         nB2=HFCIndex(i,5);
         ix1=HFCIndex(i,6);
         ix2=HFCIndex(i,7);
         f=HFC(i);
         DM((nB1-1)*3+ix1,(nB2-1)*3+ix2)=DM((nB1-1)*3+ix1,(nB2-1)*3+ix2) &
            &  +1d0/sqrt(mass(nB1)*mass(nB2))*f*exp(2d0*PI*um*dot_product(nx*a1+ny*a2+nz*a3,qq))
      enddo
      DM=(DM+transpose(conjg(DM)))/2d0
      return
   end subroutine Cal_DM

   subroutine Cal_DDM(qq,DDM,nDim)
      use Crystal,only: mass,PI,um,a1,a2,a3
      use IFC
      implicit none
      integer::nDim
      real(8)::qq(3),rr(3)
      complex*16::DDM(nDim,nDim,3)
      integer::i,j,nx,ny,nz,nB1,nB2,ix1,ix2
      real(8)::f
      DDM=0
      do i=1,nHarm
         nx=HFCIndex(i,1);
         ny=HFCIndex(i,2);
         nz=HFCIndex(i,3);
         nB1=HFCIndex(i,4);
         nB2=HFCIndex(i,5);
         ix1=HFCIndex(i,6);
         ix2=HFCIndex(i,7);
         f=HFC(i);
         rr=nx*a1+ny*a2+nz*a3
         do j=1,3
            DDM((nB1-1)*3+ix1,(nB2-1)*3+ix2,j)=DDM((nB1-1)*3+ix1,(nB2-1)*3+ix2,j) &
               & +um*rr(j)/sqrt(mass(nB1)*mass(nB2))*f*exp(2d0*PI*um*dot_product(rr,qq))
         enddo
      enddo
      return
   end subroutine Cal_DDM
end module Disp