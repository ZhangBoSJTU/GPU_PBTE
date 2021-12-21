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
!  Calculate the thermal conductivity

subroutine Cal_TC(iK,iTemp)
use, intrinsic :: ISO_FORTRAN_ENV
use Crystal
use mesh
use Task
use symm
use Quantum
implicit none
integer::i,j,k,kk,is0,iq0,iIt,iq0t,iq0_used,iK,iTemp
real(8)::nn0,TCx,lamda,LBoundary
real(8)::sr,sign1,grpvel(NQtot,nBand)
real(8)::QQ(NQtot,nBand),FF0(NQtot,nBand),FF(NQtot,nBand)
character(len=40)::ss,s2,filename,ss1,ss2,ss3,s3

integer,allocatable::ID_equi(:,:),isCount(:)
real(8)::rou0,alf,omg0,rou1,omg1,beta
real(8),dimension(:),allocatable::AoutMat,SqAoutMat,bMat,FFFSMA,FFF0,FFF,bMatTmp
real(8),dimension(:),allocatable::r0Mat,r1Mat,r0hat,v0Mat,v1Mat,p0Mat,p1Mat,sMat,tMat
real(8),dimension(:,:),allocatable::AAA,AinMat,AAATmp

allocate(AoutMat(nBand*NQtot),SqAoutMat(nBand*NQtot),bMat(nBand*NQtot),FFFSMA(nBand*NQtot),&
   &  FFF0(nBand*NQtot),AinMat((nBand*NQtot),nBand),FFF(nBand*NQtot))
allocate(r0Mat(nBand*NQtot),r1Mat(nBand*NQtot),r0hat(nBand*NQtot),v0Mat(nBand*NQtot),&
   &  p0Mat(nBand*NQtot),v1Mat(nBand*NQtot),p1Mat(nBand*NQtot))
allocate(sMat(nBand*NQtot),tMat(nBand*NQtot),bmattmp(nBand*NQtot))
allocate(ID_equi(nsymm_rot_used,NQtot),isCount(NQtot))
allocate(AAA(nBand*N0_used,nBand*NQtot),AAATmp(nBand,(nBand*NQtot)))

if (iK .eq.1) grpvel=grpvelx
if (iK .eq.2) grpvel=grpvely
if (iK .eq.3) grpvel=grpvelz
write(s3,'(I3)') iK

QQ=0
FF0=0
FF=0
temp=TempList(iTemp)
write(ss3,'(I7)') iTemp
write(ss1,'(I7)') floor(temp)
if (LBoundary0 .lt. 1e9) then
   write(ss2,'(I10)') floor(LBoundary0)
else
   ss2='inf'
endif
write(filename,*) '../TC_',trim(adjustl(ss1)),'_',trim(adjustl(ss2)),'.dat'
open(20,file=filename)

LBoundary=LBoundary0*18.89726132885643d0
do iq0=1,NQtot
    iq0_used=Idequi_used(iq0)
    write(ss,'(I7)') iq0_used
    write(s2,'(I7)') (iq0_used-1)/100+1
    write(filename,*) '../result/'//trim(adjustl(s2))//'/Q_'//trim(adjustl(ss))//'_'//trim(adjustl(ss3))//'.dat'
    open(100,file=filename)
    do is0=1,nBand
        read(100,*) Q(iq0,is0),sr
    enddo
    close(100)
enddo

do iq0=1,NQtot
   if (sqrt(Qpoint(iq0,1)**2d0+Qpoint(iq0,2)**2d0+Qpoint(iq0,3)**2d0).gt.1d-5) then
       do is0=1,nBand
         nn0=BE(freq(iq0,is0),temp)
         FF0(iq0,is0)=grpvel(iq0,is0)*hb*freq(iq0,is0)*nn0*(nn0+1d0)/temp**2d0/kb/Q(iq0,is0)
         TCx=TCx+hb*freq(iq0,is0)*grpvel(iq0,is0)*nn0*(nn0+1d0)*FF0(iq0,is0)
       enddo
   endif
enddo
write(20,*) TCx/NQtot/vol*2.179872e-18/1.460528847358342e-15/5.2917720859e-11
write(*,*) TCx/NQtot/vol*2.179872e-18/1.460528847358342e-15/5.2917720859e-11

do iq0=1,NQtot
   do is0=1,nBand
      if (freq(iq0,is0).gt.1e-5) then
         nn0=BE(freq(iq0,is0),temp)
         Q(iq0,is0)=Q(iq0,is0)+nn0*(nn0+1d0)/(LBoundary/2d0)*abs(grpvel(iq0,is0))
      endif
   enddo
enddo

AoutMat=0
do iq0=1,NQtot
   do is0=1,nBand
      if (iq0 .ne.iGamma) then
         AoutMat((iq0-1)*nBand+is0)=Q(iq0,is0)
      else
         AoutMat((iq0-1)*nBand+is0)=1
      endif
   enddo
enddo

bmattmp=0
do iq0_used=1,N0_used
   iq0=List_used(iq0_used)
   do is0=1,nBand
      if (iq0 .ne.iGamma) then
         nn0=BE(freq(iq0,is0),temp)
         bmattmp((iq0_used-1)*nBand+is0)=-grpvel(iq0,is0)*hb*freq(iq0,is0)*nn0*(nn0+1)
      else
         bmattmp((iq0_used-1)*nBand+is0)=0d0
      endif
   enddo
enddo


do iq0_used=1,N0_used
   do i=1,Nequi_used(iq0_used)
      iq0=ALLEquiList_used(i,iq0_used)
      sign1=crotations_used(iK,iK,typeofsymmetry_used(i,iq0_used))
      do is0=1,nBand
         bmat((iq0-1)*nBand+is0)=bmattmp((iq0_used-1)*nBand+is0)*sign1
      enddo
   enddo
enddo
deallocate(bmattmp)

do i=1,NQtot*nBand
   FFFSMA(i)=bmat(i)/AoutMat(i) 
enddo

TCx=0
lamda=1d0/NQtot/vol*2.179872e-18/1.460528847358342e-15/5.2917720859e-11/kb/temp**2d0
do iq0=1,NQtot
   do is0=1,nBand
      if (freq(iq0,is0) .gt.1d-5) then
         TCx=TCx+lamda*bmat((iq0-1)*nBand+is0)*FFFSMA((iq0-1)*nBand+is0)
      endif
   enddo
enddo
write(20,*) 'SMRTA = ',TCx, ' W/mK'
write(*,*) 'SMRTA = ',TCx, ' W/mK'

SqAoutMat=sqrt(AoutMat)
FFF0=FFFSMA*SqAoutMat
bmat=bmat/SqAoutMat

do iq0_used=1,N0_used
   iq0=List_used(iq0_used)
   AinMat=0
   if (mod(i,1000) .eq.0)  write(*,*) iq0_used
   write(ss,'(I7)') iq0_used
   write(s2,'(I7)') (iq0_used-1)/100+1
   write(filename,*) '../result/'//trim(adjustl(s2))//'/A_'//trim(adjustl(ss))//'_'//trim(adjustl(ss3))//'.dat'
   open(10,file=filename)
   read(10,*) AinMat
   do k=1,nBand
      AinMat(:,k)=AinMat(:,k)/SqAoutMat
      AinMat(:,k)=AinMat(:,k)/SqAoutMat((iq0-1)*nBand+k)
      AinMat((iq0-1)*nBand+k,k)=AinMat((iq0-1)*nBand+k,k)+1
   enddo
   close(10)
   AAA((iq0_used-1)*nBand+1:iq0_used*nBand,:)=transpose(AinMat)
enddo

deallocate(AinMat,SqAoutMat)

call Equ_List(ID_equi,nsymm_rot_used,qrotations_used)

isCount=0
do iq0_used=1,N0_used
   iq0=List_used(iq0_used)
   do j=1,nsymm_rot_used
      iq0t=ID_equi(j,iq0)
      if (isCount(iq0t) .eq.1) cycle
      isCount(iq0t)=1
      do k=1,NQtot
         kk=ID_equi(j,k)
         AAAtmp(:,(kk-1)*nBand+1:kk*nBand)=AAA((iq0_used-1)*nBand+1:iq0_used*nBand,(k-1)*nBand+1:k*nBand)
      enddo
      do is0=1,nBand
         i=(iq0t-1)*nBand+is0
         r0Mat(i)=bMat(i)-dot_product(AAAtmp(is0,:),FFF0)
      enddo
   enddo
enddo
r1Mat=r0Mat
r0Hat=r0Mat
rou0=1
alf=1
omg0=1
v0Mat=0
v1Mat=v0Mat
p0Mat=0
p1Mat=p0Mat


do iIt=1,NIt
   rou1=dot_product(r0Hat,r0Mat)
   beta=(rou1/rou0)*(alf/omg0)
   p1Mat=r0Mat+beta*(p0Mat-omg0*v0Mat)
   isCount=0
   do iq0_used=1,N0_used
      iq0=List_used(iq0_used)
      do j=1,nsymm_rot_used
         iq0t=ID_equi(j,iq0)
         if (isCount(iq0t) .eq.1) cycle
         isCount(iq0t)=1
         do k=1,NQtot
            kk=ID_equi(j,k)
            AAAtmp(:,(kk-1)*nBand+1:kk*nBand)=AAA((iq0_used-1)*nBand+1:iq0_used*nBand,(k-1)*nBand+1:k*nBand)
         enddo
         do is0=1,nBand
            i=(iq0t-1)*nBand+is0
            v1Mat(i)=dot_product(AAAtmp(is0,:),p1Mat)
         enddo
      enddo
   enddo


   alf=rou1/dot_product(r0Hat,v1Mat)
   sMat=r0Mat-alf*v1Mat
   isCount=0
   do iq0_used=1,N0_used
      iq0=List_used(iq0_used)
      do j=1,nsymm_rot_used
         iq0t=ID_equi(j,iq0)
         if (isCount(iq0t) .eq.1) cycle
         isCount(iq0t)=1
         do k=1,NQtot
            kk=ID_equi(j,k)
            AAAtmp(:,(kk-1)*nBand+1:kk*nBand)=AAA((iq0_used-1)*nBand+1:iq0_used*nBand,(k-1)*nBand+1:k*nBand)
         enddo
         do is0=1,nBand
            i=(iq0t-1)*nBand+is0
            tMat(i)=dot_product(AAAtmp(is0,:),sMat)
         enddo
      enddo
   enddo
   omg1=dot_product(tmat,smat)/dot_product(tmat,tmat)
   FFF=FFF0+alf*p1Mat+omg1*sMat
   r1Mat=sMat-omg1*tMat

   rou0=rou1
   p0Mat=p1Mat
   v0Mat=v1Mat
   omg0=omg1
   FFF0=FFF
   r0Mat=r1Mat

   TCx=0
   do iq0=1,NQtot
   do is0=1,nBand
      if (freq(iq0,is0) .gt.1d-5) then
         TCx=TCx+lamda*bmat((iq0-1)*nBand+is0)*FFF0((iq0-1)*nBand+is0)
      endif
   enddo
   enddo
   write(20,*) 'step ',iIt,': ',TCx
   write(*,*) 'step ',iIt,': ',TCx

   write(filename,*) '../result/FF_',trim(adjustl(ss1)),'_',trim(adjustl(s3)),'_',trim(adjustl(ss2)),'_',trim(adjustl(ss3)),'.dat'
   open(30,file=filename)
   do iq0=1,NQtot
   do is0=1,nBand
      write(30,*) FFF0((iq0-1)*nBand+is0)
   enddo
   enddo
   close(30)
enddo
deallocate(ID_equi,tMat,sMat,r0Mat,r1Mat,r0hat,v0mat,p0Mat,v1Mat,p1Mat)
deallocate(bMat,FFFSMA,FFF0,FFF,AAA,AAAtmp,isCount)
return
end subroutine Cal_TC