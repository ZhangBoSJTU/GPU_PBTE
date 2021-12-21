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

subroutine IS_DONE(iq0_used,isdone)
   use mesh
   use Crystal
   use symm
   implicit none
   integer::iq0_used,isDone
   character(30)::s1,s2,ss1,ss2
   write(s1,'(i7)') (iq0_used-1)/100+1
   write(ss1,*) '../result/'//trim(adjustl(s1))//'/'
   write(s2,'(i7)') iq0_used
   ss2 = trim(adjustl(ss1))//trim(adjustl(s2))//'.z'
   open(10,file=ss2)
   read(10,*) isDone 
   close(10)
   return
end subroutine IS_DONE

subroutine Cal_RT_Init(iq0_used)
   use mesh, only: NQtot
   use Crystal, only: nBand,RTA,nTemp
   use symm
   use TmpSpace
   implicit none
   integer::iq0_used,iq0
   integer::isDone
   iq0=List_used(iq0_used)
   write(*,*) NQtot,N0_used,iq0_used,iq0
   call IS_DONE(iq0_used,isDone)
   if (isDone .ne. 1) then
      allocate(QQQ(nBand,nTemp),QQQ3ph(nBand,nTemp),QQQ4ph(nBand,nTemp),QQQiso(nBand,nTemp))
      if (.not. RTA) then
         allocate(AAA(nBand*NQtot,nBand,nTemp))
         AAA=0d0
      endif
      QQQ=0
      QQQ3ph=0
      QQQ4ph=0
      QQQiso=0
   endif
   return                  
end subroutine Cal_RT_Init

subroutine Write_RT(iq0_used)
   use mesh
   use Crystal
   use symm
   use TmpSpace
   use Quantum
   implicit none
   integer::iq0_used,iq0
   integer::isDone,is0,iT
   real(8)::nn0
   character(30)::s1,s2,s3,ss1,ss2
   iq0=List_used(iq0_used)
   call IS_DONE(iq0_used,isDone)
   if (isDone .eq. 1) return
   write(s1,'(i7)') (iq0_used-1)/100+1
   write(ss1,*) '../result/'//trim(adjustl(s1))//'/'
   write(s2,'(i7)') iq0_used
   if (.not. RTA) then
      do iT=1,nTemp
         write(s3,'(i7)') iT
         if (iq0 .eq. iGamma ) then
            AAA(:,1:3,iT) = 0d0
         endif
         write(ss2,*) trim(adjustl(ss1))//'A_'//trim(adjustl(s2))//'_'//trim(adjustl(s3))//'.dat'
         open(1000,file=ss2)
         write(1000,*) AAA(:,:,iT)
         close(1000)
      enddo
      deallocate(AAA)
   endif
   do iT=1,nTemp
      write(s3,'(i7)') iT
      write(ss2,*) trim(adjustl(ss1))//'Q_'//trim(adjustl(s2))//'_'//trim(adjustl(s3))// '.dat'
      open(10,file=ss2)
      do is0=1,nBand  
         nn0=BE(freq(iq0,is0),TempList(iT))
         write(10,*) QQQ(is0,iT),QQQ(is0,iT)/nn0/(nn0+1)
      enddo
      write(10,*)
      do is0=1,nBand 
         nn0=BE(freq(iq0,is0),TempList(iT))
         write(10,*) QQQ3ph(is0,iT),QQQ3ph(is0,iT)/nn0/(nn0+1)
      enddo
      write(10,*)
      do is0=1,nBand 
         nn0=BE(freq(iq0,is0),TempList(iT))
         write(10,*) QQQ4ph(is0,iT),QQQ4ph(is0,iT)/nn0/(nn0+1)
      enddo
      write(10,*)
      do is0=1,nBand 
         nn0=BE(freq(iq0,is0),TempList(iT))
         write(10,*) QQQiso(is0,iT),QQQiso(is0,iT)/nn0/(nn0+1)
      enddo
      close(10)
   enddo
   deallocate(QQQ,QQQ3ph,QQQ4ph,QQQiso)
   ss2 = trim(adjustl(ss1))//trim(adjustl(s2))//'.z'
   open(10,file=ss2)
   write(10,*) 1
   close(10)
end subroutine Write_RT