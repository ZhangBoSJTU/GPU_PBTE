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

! Thin, specialized wrapper around spglib, a library by Atsushi Togo.

module symmetry
  use iso_c_binding
  implicit none
  real(kind=C_DOUBLE),parameter :: symprec=1d-5
  contains
  function get_num_operations(lattice,natoms,types,positions)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: natoms
    integer(kind=4),dimension(natoms),intent(in) :: types
    real(kind=8),dimension(3,natoms),intent(in) :: positions
    integer(kind=4) :: get_num_operations
    real(kind=C_DOUBLE),dimension(3,3) :: clattice
    integer(kind=C_INT) :: cnatoms
    integer(kind=C_INT),dimension(natoms) :: ctypes
    real(kind=C_DOUBLE),dimension(3,natoms) :: cpositions
    integer(kind=C_INT) :: num

    clattice=transpose(lattice)
    cnatoms=natoms
    ctypes=types
    cpositions=positions
    call spg_get_multiplicity(num,clattice,&
         cpositions,ctypes,cnatoms,symprec)
    get_num_operations=num
  end function get_num_operations

  subroutine get_operations(lattice,natoms,types,positions,nops,&
       rotations,translations,international)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: natoms
    integer(kind=4),dimension(natoms),intent(in) :: types
    real(kind=8),dimension(3,natoms),intent(in) :: positions
    integer(kind=4),intent(inout) :: nops
    integer(kind=4),dimension(3,3,nops),intent(out) :: rotations
    real(kind=8),dimension(3,nops),intent(out) :: translations
    character(len=10),intent(out) :: international

    integer(kind=C_INT) ::  newnops,i
    real(kind=C_DOUBLE),dimension(3,3) :: clattice
    integer(kind=C_INT) :: cnatoms
    integer(kind=C_INT),dimension(natoms) :: ctypes
    real(kind=C_DOUBLE),dimension(3,natoms) :: cpositions
    integer(kind=C_INT) :: cnops
    integer(kind=C_INT),dimension(3,3,nops) :: crotations
    real(kind=C_DOUBLE),dimension(3,nops) :: ctranslations
    character(len=11,kind=C_CHAR) :: intertmp

    clattice=transpose(lattice)
    cnatoms=natoms
    ctypes=types
    cpositions=positions
    cnops=nops

    call spg_get_symmetry(newnops,crotations,ctranslations,cnops,&
         clattice,cpositions,ctypes,cnatoms,symprec)
    call spg_get_international(i,intertmp,clattice,&
         cpositions,ctypes,cnatoms,symprec)
    international=intertmp(1:10)
    nops=newnops
    do i=1,nops
       rotations(:,:,i)=transpose(crotations(:,:,i))
    end do
    translations=ctranslations
  end subroutine get_operations

  subroutine get_cartesian_operations(lattice,nops,rotations,translations,&
       crotations,ctranslations)
    implicit none
    real(kind=8),dimension(3,3),intent(in) :: lattice
    integer(kind=4),intent(in) :: nops
    integer(kind=4),dimension(3,3,nops),intent(in) :: rotations
    real(kind=8),dimension(3,nops),intent(in) :: translations
    real(kind=8),dimension(3,3,nops),intent(out) :: crotations
    real(kind=8),dimension(3,nops),intent(out) :: ctranslations

    integer(kind=4) :: i,info
    integer(kind=4),dimension(3) :: P
    real(kind=8),dimension(3,3) :: tmp1,tmp2

    ctranslations=matmul(lattice,translations)
    do i=1,nops
       tmp1=transpose(lattice)
       tmp2=transpose(rotations(:,:,i))
       call dgesv(3,3,tmp1,3,P,tmp2,3,info)
       crotations(:,:,i)=transpose(matmul(tmp2,transpose(lattice)))
    end do
  end subroutine get_cartesian_operations
end module symmetry