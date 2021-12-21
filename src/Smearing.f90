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
!  Smearing functions.
module SmearingFunctions
   implicit none
   integer::tpSmearing
   real(8)::scoeff,scoeff4
   contains
   
   real(8) function sig3(v)
      use Crystal, only: c1,c2,c3
      use mesh, only: nQx, nQy, nQz
      implicit none
      real(8)::v(3)
      sig3=0.
      if (tpSmearing .eq. 0) then  
         sig3=scoeff
      elseif (tpSmearing .eq. 1) then 
         sig3=sig3+(dot_product(c1,v)/nQx)**2
         sig3=sig3+(dot_product(c2,v)/nQy)**2
         sig3=sig3+(dot_product(c3,v)/nQz)**2
         sig3=sqrt(sig3/6d0)*scoeff
      endif
      return
   end function sig3

   real(8) function sig4(v)
      use Crystal, only: c1,c2,c3
      use mesh, only: nQx,nQy,nQz
      implicit none
      real(8)::v(3)
      sig4=0.
      if (tpSmearing .eq. 0) then  
         sig4=scoeff4
      elseif (tpSmearing .eq. 1) then 
         sig4=sig4+(dot_product(c1,v)/nQx)**2
         sig4=sig4+(dot_product(c2,v)/nQy)**2
         sig4=sig4+(dot_product(c3,v)/nQz)**2
         sig4=sqrt(sig4/6d0)*scoeff4
      endif
      return
   end function sig4

   real(8) function dG(x,ee)
      use Crystal, only: PI
      real(8) x,ee
      dG=1d0/sqrt(PI)/ee*exp(-(x/ee)**2d0)
      return
   end function dG
end module SmearingFunctions