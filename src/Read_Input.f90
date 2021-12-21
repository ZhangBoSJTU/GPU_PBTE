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
!  Suroutine to read input.txt.


subroutine Read_Input
   use Crystal
   use mesh
   use VectorFunctions
   use SmearingFunctions,only:scoeff,scoeff4,tpSmearing
   implicit none
   integer::i
   open(10,file='../input/input.txt')
   read(10,*) RTA,NIt
   if (RTA .eqv. .false.) write(*,*) 'Iterative calculation'
   if (RTA .eqv. .true.) write(*,*) 'RTA calculation'
   read(10,*) isCubic    
   read(10,*) isQuart      
   read(10,*) isIso        
   read(10,*) a0           
   read(10,*) a1
   read(10,*) a2 
   read(10,*) a3
   b1=cross_product(a2,a3)/s3_product(a1,a2,a3)
   b2=cross_product(a3,a1)/s3_product(a1,a2,a3)
   b3=cross_product(a1,a2)/s3_product(a1,a2,a3)
   c1=b1*2*PI/a0
   c2=b2*2*PI/a0
   c3=b3*2*PI/a0
   vol=abs(s3_product(a1,a2,a3)*a0**3.0)
   read(10,*) nBasis     
   allocate(mass(nBasis),g2(nBasis),pos(nBasis,3),types(nBasis))
   read(10,*) mass(:)      
   read(10,*) types(:)     
   read(10,*) g2(:)        
 
   do i=1,nBasis
      read(10,*) pos(i,:)  
   enddo
   read(10,*) nQx,nQy,nQz  
   read(10,*) nTemp        
   allocate(TempList(nTemp))
   read(10,*) TempList(:)
   read(10,*) isKx,isKy,isKz
   read(10,*) LBoundary0    
   read(10,*) tpSmearing,scoeff,scoeff4
   nBand=3*nBasis
   return
end subroutine Read_Input