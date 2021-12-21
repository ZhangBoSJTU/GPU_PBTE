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
!  Symmetry Operation.

module Symm
   real(8),allocatable :: qrotations_ibz(:,:,:),crotations_ibz(:,:,:),rotations_ibz(:,:,:)
   integer::nsymm_rot_ibz,N0_ibz
   integer(kind=4),allocatable :: Nequi_ibz(:),List_ibz(:),IdEqui_ibz(:)
   integer(kind=4),allocatable :: ALLEquiList_ibz(:,:),TypeofSymmetry_ibz(:,:)
   real(8),allocatable ::qrotations_8fbz(:,:,:),crotations_8fbz(:,:,:),rotations_8fbz(:,:,:)
   integer::nsymm_rot_8fbz,N0_8fbz
   integer(kind=4),allocatable :: Nequi_8fbz(:),List_8fbz(:),IdEqui_8fbz(:)
   integer(kind=4),allocatable :: ALLEquiList_8fbz(:,:),TypeofSymmetry_8fbz(:,:)
   real(8),allocatable::qrotations_used(:,:,:),crotations_used(:,:,:),rotations_used(:,:,:)
   integer::nsymm_rot_used,N0_used
   integer(kind=4),allocatable :: Nequi_used(:),List_used(:),IdEqui_used(:)
   integer(kind=4),allocatable :: ALLEquiList_used(:,:),TypeofSymmetry_used(:,:)
   contains

   subroutine Choose_Symm
      use Crystal,only:RTA
      use mesh,only: NQtot
      implicit none

      if (RTA) then
         N0_used=N0_ibz
         nsymm_rot_used=nsymm_rot_ibz 
         allocate(Nequi_used(NQtot),List_used(NQtot),IdEqui_used(NQtot), &
            & ALLEquiList_used(Nsymm_rot_used,NQtot),TypeofSymmetry_used(Nsymm_rot_used,NQtot))
         Nequi_used=Nequi_ibz
         List_used=List_ibz
         IdEqui_used=IdEqui_ibz
         ALLEquiList_used=ALLEquiList_ibz
         TypeofSymmetry_used=TypeofSymmetry_ibz
         deallocate(Nequi_ibz,List_ibz,IdEqui_ibz,ALLEquiList_ibz,TypeofSymmetry_ibz)
      else
         N0_used=N0_8fbz
         nsymm_rot_used=nsymm_rot_8fbz 
         allocate(Nequi_used(NQtot),List_used(NQtot),IdEqui_used(NQtot), &
            & ALLEquiList_used(Nsymm_rot_used,NQtot),TypeofSymmetry_used(Nsymm_rot_used,NQtot), &
            & crotations_used(3,3,Nsymm_rot_used),qrotations_used(3,3,Nsymm_rot_used))
         Nequi_used=Nequi_8fbz
         List_used=List_8fbz
         IdEqui_used=IdEqui_8fbz
         ALLEquiList_used=ALLEquiList_8fbz
         TypeofSymmetry_used=TypeofSymmetry_8fbz
         crotations_used=crotations_8fbz
         qrotations_used=qrotations_8fbz
         deallocate(Nequi_ibz,List_ibz,IdEqui_ibz,ALLEquiList_ibz,TypeofSymmetry_ibz)
         deallocate(Nequi_8fbz,List_8fbz,IdEqui_8fbz,ALLEquiList_8fbz,TypeofSymmetry_8fbz)
      endif
   endsubroutine Choose_Symm
end module Symm