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
!  GPU_PBTE data.

module Crystal
    integer:: nBasis,nBand,NIt,nTemp
    real(8):: a0,vol,temp,LBoundary0
    real(8),dimension(:),allocatable::TempList
    logical::RTA,isKx,isKy,isKz,isCubic,isQuart,isIso
    real(8),dimension(3),parameter::e1=(/1.d0,0.d0,0.d0/)
    real(8),dimension(3),parameter::e2=(/0.d0,1.d0,0.d0/)
    real(8),dimension(3),parameter::e3=(/0.d0,0.d0,1.d0/)
    real(8),dimension(3)::a1,a2,a3,b1,b2,b3,c1,c2,c3
    real(8),allocatable::mass(:),g2(:),pos(:,:)
    integer,allocatable:: types(:)
    real(8),parameter::hb=0.033123406258297d0,kb=6.333623318249880d-006,PI=3.141592653589793d0
    complex*16:: um=(0d0,1d0)
    real(8):: wmin=1d-5,wmin0=1d-9
end module crystal

module mesh
    real(8):: dx,dy,dz
    integer:: NQtot,nQx,nQy,nQz,iGamma,nIntp
    real(8),dimension(:,:),allocatable::Qpoint
    real(8),dimension(:,:),allocatable:: freq
    complex*16,dimension(:,:,:),allocatable::ev
    real(8),dimension(:,:),allocatable::grpvelx,grpvely,grpvelz
    integer,dimension(:,:,:),allocatable::ind3
    integer,dimension(:,:),allocatable::ind1
    real(8),dimension(:,:),allocatable::Q
end module mesh

module Task
    integer::JobTask,nJob,nNode
    integer,dimension(:),allocatable::jobList
end module Task
 
module TmpSpace
    real(8),allocatable::AAA(:,:,:),QQQ(:,:)
    real(8),allocatable::QQQ3ph(:,:),QQQ4ph(:,:),QQQiso(:,:)
end module TmpSpace