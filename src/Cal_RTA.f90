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
!  RTA calculation

subroutine Cal_RTA(iT)
use Crystal,only:vol,kb,hb,temp,nBand,wmin,LBoundary0,TempList
use mesh,only:freq,grpvelx,grpvely,grpvelz,NQtot,ind1
use symm,only:Idequi_used
use Quantum
implicit none
integer::iq0_used,iq0,is0,iT
character(50)::ss,s1,s2,s3
real(8)::n,Kxx,Kyy,Kzz,Kxx1,Kyy1,Kzz1,cc,sr,QQQ,LBoundary
real(8),dimension(:),allocatable::ff


temp=TempList(iT)
allocate(ff(nBand))
write(s3,'(i7)') iT
LBoundary=LBoundary0*18.89726132885643d0
open(20,file='../result/rt_'//trim(adjustl(s3))//'.dat')
Kxx=0
Kyy=0
Kzz=0
Kxx1=0
Kyy1=0
Kzz1=0
do iq0=1,NQtot
      ff=freq(iq0,:)      
      iq0_used=Idequi_used(iq0)
      write(ss,'(i7)') iq0_used
      write(s2,'(i7)') (iq0_used-1)/100+1
      s1='../result/'// trim(adjustl(s2))//'/Q_'//trim(adjustl(ss))//'_'//trim(adjustl(s3))//'.dat'
      open(10,file=s1,status='unknown') 
      do is0=1,nBand
         read(10,*) QQQ,sr
         n=BE(freq(iq0,is0),temp)
         if (freq(iq0,is0) .ge. wmin) then
            Kxx=Kxx+grpvelx(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/sr
            Kyy=Kyy+grpvely(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/sr
            Kzz=Kzz+grpvelz(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/sr
            Kxx1=Kxx1+grpvelx(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/(sr+abs(grpvelx(iq0,is0))*2d0/LBoundary)
            Kyy1=Kyy1+grpvely(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/(sr+abs(grpvely(iq0,is0))*2d0/LBoundary)
            Kzz1=Kzz1+grpvelz(iq0,is0)**2d0*ff(is0)*ff(is0)*n*(n+1)/(sr+abs(grpvelz(iq0,is0))*2d0/LBoundary)
         endif
         if (sr .lt. 1e-30) then
            cc=1d30
         else
            cc=1d0/sr
         endif
         write(20,'(I3,1x,I3,1x,I3,1x,I6,1x,I4,1x,e20.13,1x,e20.13)') ind1(iq0,1),ind1(iq0,2),ind1(iq0,3),iq0,is0,ff(is0),cc
      enddo
      close(10)
enddo

close(20)

Kxx=Kxx/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
Kyy=Kyy/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
Kzz=Kzz/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
Kxx1=Kxx1/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
Kyy1=Kyy1/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
Kzz1=Kzz1/vol/NQtot/kB/temp/temp*hb*hb*2.179872e-18/1.460528847358342e-15/5.2917720859e-11

write(*,*) Kxx,Kxx1
write(*,*) Kyy,Kyy1
write(*,*) Kzz,Kzz1
return
end subroutine Cal_RTA
