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
!  Vector functions.

module VectorFunctions
  implicit none
contains    
  function cross_product(a, b)
    real(8), dimension(3) :: cross_product
    real(8), dimension(3), intent(in) :: a, b
 
    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product
   function s3_product(a, b, c)
    real(8) :: s3_product
    real(8), dimension(3), intent(in) :: a, b, c
 
    s3_product = dot_product(a, cross_product(b, c))
  end function s3_product
  function v3_product(a, b, c)
    real(8), dimension(3) :: v3_product
    real(8), dimension(3), intent(in) :: a, b, c
 
    v3_product = cross_product(a, cross_product(b, c))
  end function v3_product
end module VectorFunctions