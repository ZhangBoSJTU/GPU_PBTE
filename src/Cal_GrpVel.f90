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
!  Calculate the group velocity of phonons.

subroutine Cal_GrpVel
   use Crystal,only:nBand,a0
   use mesh,only:grpvelx,grpvely,grpvelz,freq,ev,NQtot,Qpoint
   use symm,only:nsymm_rot_ibz,qrotations_ibz,crotations_ibz
   use Disp

   implicit none
   integer::i,j,is,kk
   real(8)::qq(3)
   real(8)::newvelocity(nBand,3),velocity(nBand,3)
   complex*16,dimension(:,:,:),allocatable::DDM
   integer::ID_equi(nsymm_rot_ibz,NQtot)

   allocate(DDM(nBand,nBand,3))
   grpvelx=0
   grpvely=0
   grpvelz=0
   do i=1,NQtot
      qq=Qpoint(i,:)
      call Cal_DDM(qq,DDM,nBand)
      do j=1,nBand
         grpvelx(i,j)=real(dot_product(ev(i,:,j),matmul(DDM(:,:,1),ev(i,:,j))))/2d0/freq(i,j)
         grpvely(i,j)=real(dot_product(ev(i,:,j),matmul(DDM(:,:,2),ev(i,:,j))))/2d0/freq(i,j)
         grpvelz(i,j)=real(dot_product(ev(i,:,j),matmul(DDM(:,:,3),ev(i,:,j))))/2d0/freq(i,j)
      enddo
   enddo
   grpvelx=grpvelx*a0
   grpvely=grpvely*a0
   grpvelz=grpvelz*a0
   call Equ_List(ID_equi,nsymm_rot_ibz,qrotations_ibz)
   do i=1,NQtot
      velocity=0
      do is=1,nBand
         velocity(is,1:3)=(/grpvelx(i,is),grpvely(i,is),grpvelz(i,is)/)
      enddo
      newvelocity=0
      kk=0
      do j=1,nsymm_rot_ibz
         if (ID_equi(j,i).eq.i) then   
            newvelocity = newvelocity + transpose(&
                & matmul(crotations_ibz(:,:,j),transpose(velocity)))
            kk=kk+1
         endif
      enddo
      if (kk .ge. 1) then
         velocity=newvelocity/kk
      endif
      do is=1,nBand
         grpvelx(i,is)=velocity(is,1)
         grpvely(i,is)=velocity(is,2)
         grpvelz(i,is)=velocity(is,3)
      enddo
   enddo
   return
end subroutine Cal_GrpVel