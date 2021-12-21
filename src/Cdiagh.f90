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
!  Math function.

subroutine Cdiagh(n,h,e,v)
   implicit none
   integer::n
   real(8)::e(n)  
   complex*16::h(n,n)
   complex*16::v(n,n)
   call cdiagh_lapack(v,e)
   return  
   contains
   subroutine cdiagh_lapack(v,e)
      implicit none
      real(8)::e(n) 
      complex*16::v(n,n)
      integer::lwork,nb,info
      real(8),allocatable::rwork(:)
      complex*16,allocatable::work(:)
      integer,external::ILAENV
      nb=ILAENV(1,'ZHETRD','U',n,-1,-1,-1)
      if (nb .lt. 1 .or. nb .ge. n)then
         lwork=2*n
      else
         lwork=(nb+1)*n
      endif
      v=h
      allocate(work(lwork))
      allocate(rwork(3*n-2))
      call ZHEEV('V','U',n,v,n,e,work,lwork,rwork,info)
      deallocate(work)
      deallocate(rwork)
      return
   end subroutine cdiagh_lapack
end subroutine Cdiagh