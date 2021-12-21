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

subroutine Find_EQP(N0,Nequi,List,IdEqui,AllEquiList,TypeofSymmetry,Nsymm_rot,qrotations,NQtot)
implicit none

integer,intent(in)::Nsymm_rot,NQtot
real(8),intent(in)::qrotations(3,3,Nsymm_rot)
integer,intent(out)::N0,Nequi(NQtot),List(NQtot),IdEqui(NQtot)
integer,intent(out)::ALLEquiList(Nsymm_rot,NQtot),TypeofSymmetry(Nsymm_rot,NQtot)

integer(kind=4) :: ID_equi(Nsymm_rot,NQtot),ii,jj,ll,iaux,Ntot
integer(kind=4) :: NAllList,AllList(NQtot),EquiList(Nsymm_rot)

call Equ_List(Id_equi,nsymm_rot,qrotations)

IdEqui=0
AllEquiList=0
TypeofSymmetry=0
Ntot=0
N0=0
Nequi=0
NAllList=0
List=0
AllList=0
do ii=1,NQtot
   iaux=1
   if (ii.ne.1) then
      do jj=1,NAllList
         if (ii.eq.AllList(jj)) then
            iaux=0
            exit
         endif
      end do
   end if
   if (iaux.eq.1) then
      N0=N0+1
      List(N0)=ii
      Nequi(N0)=0
      EquiList(:)=ID_equi(:,ii)
      if(any(EquiList.eq.-1)) then
         write(*,*) "Error: and incompatible symmetry was not removed"
      end if
      do ll=1,Nsymm_rot
         iaux=1
         if (ll.ne.1) then
            do jj=1,Nequi(N0)
               If(EquiList(ll).eq.EquiList(jj)) then
                  iaux=0
                  exit
               end if
            end do
         end if
         if (iaux.eq.1) then
            Nequi(N0)=Nequi(N0)+1
            EquiList(Nequi(N0))=EquiList(ll)
            NAllList=NAllList+1
            AllList(NAllList)=EquiList(Nequi(N0))
            ALLEquiList(Nequi(N0),N0)=EquiList(Nequi(N0))
            TypeofSymmetry(Nequi(N0),N0)=ll
          end if
      end do
       Ntot=Ntot+Nequi(N0)
   end if
end do

write(*,*) "Info: Ntot =",Ntot
write(*,*) "Info: N0 =",N0
if(ntot.ne.NQtot) then
   write(*,*) "Error: ntot!=nptk"
end if

do ii=1,N0
   do jj=1,Nequi(ii)
      IdEqui(ALLEquiList(jj,ii))= ii
   enddo
enddo

end subroutine Find_EQP

subroutine Find_Symm
use mesh
use Crystal
use symm
use symmetry
implicit none
integer::i
integer::nsymm,nsymm_rot,P(3),info,ll,ii,jj,kk
logical,allocatable :: valid(:)
real(8)::lattvec(3,3),positions(3,nBasis)
integer,allocatable :: rotations(:,:,:),rtmp(:,:,:),ID_equi(:,:)
real(8),allocatable :: qrotations(:,:,:),crotations(:,:,:),crtmp(:,:,:),qrtmp(:,:,:)
real(8),allocatable :: translations(:,:),ctranslations(:,:)
character(len=10) :: international
real(8)::tmp1(3,3),tmp2(3,3),tmp3(3,3)


   lattvec(:,1)=a1
   lattvec(:,2)=a2
   lattvec(:,3)=a3
   do i=1,nBasis
      positions(:,i)=pos(i,:)
   enddo
   nsymm=get_num_operations(lattvec,nBasis,types,positions)
   nsymm_rot=2*nsymm
   allocate(translations(3,nsymm),rotations(3,3,nsymm_rot),&
         ctranslations(3,nsymm),crotations(3,3,nsymm_rot),&
         qrotations(3,3,nsymm_rot))
   call get_operations(lattvec,nBasis,types,positions,nsymm,&
        rotations,translations,international)
   write(*,*) "Info: symmetry group ",trim(international)," detected"
   write(*,*) "Info: ",nsymm," symmetry operations"
   call get_cartesian_operations(lattvec,nsymm,rotations,translations,&
      crotations,ctranslations)
 
   ! Transform the rotation matrices to the reciprocal-space basis.
   do i=1,nsymm
      tmp1=matmul(transpose(lattvec),lattvec)
      tmp2=transpose(rotations(:,:,i))
      tmp3=tmp1
      call dgesv(3,3,tmp1,3,P,tmp2,3,info)
      qrotations(:,:,i)=transpose(matmul(tmp2,tmp3))
   end do
   rotations(:,:,nsymm+1:2*nsymm)=-rotations(:,:,1:nsymm)
   qrotations(:,:,nsymm+1:2*nsymm)=-qrotations(:,:,1:nsymm)
   crotations(:,:,nsymm+1:2*nsymm)=-crotations(:,:,1:nsymm)
   allocate(ID_Equi(nsymm_rot,NQtot),valid(nsymm_rot))
   valid=.TRUE.
   ll=0
   do ii=2,nsymm_rot
      do jj=1,ii-1
         if(.not.valid(jj))cycle
         if(all(rotations(:,:,ii).eq.rotations(:,:,jj))) then
            valid(ii)=.FALSE.
            ll=ll+1
            exit
         end if
      end do
   end do
   if(ll.ne.0)  write(*,*) "Info:",ll,"duplicated rotations will be discarded"
   call Equ_List(ID_equi,nsymm_rot,qrotations)
   jj=0
   do ii=1,nsymm_rot
      if(valid(ii).and.any(ID_equi(ii,:).eq.-1)) then
         valid(ii)=.FALSE.
         jj=jj+1
      end if
   end do
   if(jj.ne.0)write(*,*) "Info:",jj,"rotations are incompatible with the q-point grid and will be discarded"
   if(ll+jj.ne.0) then
      call move_alloc(rotations,rtmp)
      call move_alloc(crotations,crtmp)
      call move_alloc(qrotations,qrtmp)
      allocate(rotations(3,3,nsymm_rot-ll-jj),crotations(3,3,nsymm_rot-ll-jj),&
          qrotations(3,3,nsymm_rot-ll-jj))
      kk=0
      do ii=1,nsymm_rot
         if(valid(ii)) then
            kk=kk+1
            rotations(:,:,kk)=rtmp(:,:,ii)
            crotations(:,:,kk)=crtmp(:,:,ii)
            qrotations(:,:,kk)=qrtmp(:,:,ii)
         end if
      end do
      nsymm_rot=nsymm_rot-ll-jj
      deallocate(rtmp,crtmp,qrtmp)
   end if

   nsymm_rot_ibz=nsymm_rot
   allocate(qrotations_ibz(3,3,nsymm_rot_ibz),crotations_ibz(3,3,nsymm_rot_ibz), &
       & rotations_ibz(3,3,nsymm_rot_ibz))
   qrotations_ibz=qrotations
   crotations_ibz=crotations
   rotations_ibz=rotations

   deallocate(ID_Equi,valid)
   write(*,*) "Info:",nsymm_rot_ibz,"rotations will be used for ibz"

   allocate(valid(nsymm_rot_ibz))
   valid=.TRUE.
   nsymm_rot_8fbz=0
   if ( .not. RTA) then
      do ii=1,nsymm_rot_ibz
         if (abs(abs(crotations(1,1,ii))-1d0) .gt.1d-2) valid(ii)=.FALSE.
         if (abs(abs(crotations(2,2,ii))-1d0) .gt.1d-2) valid(ii)=.FALSE.
         if (abs(abs(crotations(3,3,ii))-1d0) .gt.1d-2) valid(ii)=.FALSE.
         if (valid(ii) .eqv. .true.) nsymm_rot_8fbz=nsymm_rot_8fbz+1
      enddo
      if (nsymm_rot_8fbz .ne. 8) write(*,*) 'warning for symmetry for iterative calculation'
      allocate(qrotations_8fbz(3,3,nsymm_rot_8fbz),crotations_8fbz(3,3,nsymm_rot_8fbz), &
         & rotations_8fbz(3,3,nsymm_rot_8fbz))
      jj=0
      do ii=1,nsymm_rot_ibz
         if (valid(ii)) then
            jj=jj+1
            qrotations_8fbz(:,:,jj)=qrotations_ibz(:,:,ii)
            crotations_8fbz(:,:,jj)=crotations_ibz(:,:,ii)
            rotations_8fbz(:,:,jj)=rotations_ibz(:,:,ii)
         endif
      enddo
      deallocate(valid)
   endif
   deallocate(qrotations,crotations,rotations)
   deallocate(translations,ctranslations)
return
end subroutine Find_Symm

subroutine Equ_List(ID_equi,nsymm_rot,qrotations)
use mesh, only:NQtot,ind1,ind3,nQx,nQy,nQz
implicit none

integer,intent(in)::nsymm_rot
real(8),intent(in)::qrotations(3,3,nsymm_rot)
integer,intent(out)::ID_equi(nsymm_rot,NQtot)

integer::i,ii,ivec(3)
real(8)::vec(3)
   do i=1,NQtot
      do ii=1,nsymm_rot
         vec=(/nQx,nQy,nQz/)*matmul(qrotations(:,:,ii),dble(ind1(i,:))/(/nQx,nQy,nQz/))
         ivec=nint(vec)
         if(sqrt(dot_product(vec-dble(ivec),vec-dble(ivec))) .gt. 1e-2) then
            ID_equi(ii,i)=-1
         else
            ID_equi(ii,i)=ind3(modulo(ivec(1),nQx),modulo(ivec(2),nQy),modulo(ivec(3),nQz))
         end if
      enddo
   enddo
return
end subroutine Equ_List