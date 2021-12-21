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
!  Write & Read functions.
module WriteRead
   implicit none
   contains

   subroutine Write_Data
      use Crystal,only:nBand
      use mesh,only:NQtot,ind1,Qpoint,grpvelx,grpvely,grpvelz,ev,freq
      implicit none
      integer::i,j,k
      open(10,file='../result/mesh.dat')
      write(10,*) NQtot
      do i=1,NQtot
         write(10,*) i,ind1(i,:)
         write(10,'(1x,d20.13,1x,d20.13,1x,d20.13)') Qpoint(i,:)
         do j=1,nBand
            if (j .lt. nBand) then
               write(10,'(1x,d20.13$)') freq(i,j)
            else
               write(10,'(1x,d20.13)') freq(i,j)
            endif
         enddo
         do k=1,nBand
            do j=1,nBand
               if (j .lt. nBand) then
                  write(10,'(1x,d15.8,1x,d15.8$)') dble(ev(i,j,k)),aimag(ev(i,j,k))
               else
                  write(10,'(1x,d15.8,1x,d15.8)') dble(ev(i,j,k)),aimag(ev(i,j,k))
               endif
            enddo
         enddo
      enddo
      close(10)
      open(40,file='../result/grpvelx.dat')
      open(50,file='../result/grpvely.dat')
      open(60,file='../result/grpvelz.dat')
      do i=1,NQtot
         do j=1,nBand-1
            write(40,'(1x,d20.13$)') grpvelx(i,j)
            write(50,'(1x,d20.13$)') grpvely(i,j)
            write(60,'(1x,d20.13$)') grpvelz(i,j)
         enddo
         write(40,'(1x,d20.13)') grpvelx(i,nBand)
         write(50,'(1x,d20.13)') grpvely(i,nBand)
         write(60,'(1x,d20.13)') grpvelz(i,nBand)
      enddo
      close(40)
      close(50)
      close(60)
      return
   end subroutine Write_Data

   subroutine Write_Identifier
      use crystal, only: RTA
      use symm,only:N0_ibz,N0_8fbz
      implicit none
      integer::i,j,N0
      character(50)::ss1,ss2,s1,s2
      if (RTA) then
         N0=N0_ibz
      else
         N0=N0_8fbz
      endif
      do i=1,N0
         if (mod(i,100).eq.1) then
            j=(i-1)/100+1
            write(s2,'(i7)') j
            write(ss1,*) '../result/'//trim(adjustl(s2))//'/'
            call system('mkdir '//trim(adjustl(ss1)))
         endif
         write(s1,'(i7)') i
         ss2 = trim(adjustl(ss1))//trim(adjustl(s1))//'.z'
         open(10,file=ss2)
         write(10,*) 0
         close(10)
      enddo
   end subroutine Write_Identifier

   subroutine Read_Data
      use Crystal,only:nBand
      use mesh,only:freq,ev,grpvelx,grpvely,grpvelz,NQtot 
      implicit none
      integer::i,j,k,l1,l2
      integer::ss
      real(8),dimension(:),allocatable::cmp
      allocate(cmp(nBand*2))
      open(10,file='../result/mesh.dat')
      read(10,*) NQtot
      do i=1,NQtot
         read(10,*) ss,ss,ss,ss
         read(10,*)
         read(10,*) freq(i,:)
         do j=1,nBand
            read(10,*) cmp(:)
            do l1=1,nBand/3
            do l2=1,3
               k=(l1-1)*3+l2
               ev(i,k,j)=cmplx(cmp(2*k-1),cmp(2*k),8)
            enddo
            enddo
         enddo
      enddo
      close(10)
      open(20,file='../result/grpvelx.dat')
      do i=1,NQtot
         read(20,*) grpvelx(i,:)
      enddo
      close(20)
      open(30,file='../result/grpvely.dat')
      do i=1,NQtot
         read(30,*) grpvely(i,:)
      enddo
      close(30)
      open(40,file='../result/grpvelz.dat')
      do i=1,NQtot
         read(40,*) grpvelz(i,:)
      enddo
      close(40)
      return
   end subroutine Read_Data
end module WriteRead