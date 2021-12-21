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
!  Accelerate the calculation of three-phonon scattering process.

module Accrelate_V3
   
   contains
   attributes(global)subroutine Cal_IsScat_V3(NQtot,nBand,nQx,nQy,nQz,tpSmearing,scoeff,f00,&
        q12List,freq,grpvelx,grpvely,grpvelz,c1,c2,c3,isScattList1,isScattList2,isScattList,sigma)
   implicit none
   integer,value::NQtot,nBand,nQx,nQy,nQz,tpSmearing
   real(8),value::scoeff,f00
   integer,dimension(:),intent(in)::q12List(NQtot,2)
   real(8),dimension(:,:),intent(in):: freq(NQtot,nBand)
   real(8),dimension(:,:),intent(in)::grpvelx(NQtot,nBand),grpvely(NQtot,nBand),grpvelz(NQtot,nBand)
   real(8),dimension(:),intent(in)::c1(3),c2(3),c3(3)
   integer::iq1,iq2,i3,is1,is2
   real(8)::f11,f22,sigma1,sigma2
   real(8)::V1(3),V2(3),dV(3)
   real(8),parameter,constant::sigmamin=1d-8
   real(8),dimension(:,:,:,:,:),intent(out)::sigma(NQtot,nBand,nBand,2)
   logical,dimension(:,:,:,:),intent(out)::isScattList1(NQtot,nBand,nBand),&
      &  isScattList2(NQtot,nBand,nBand),isScattList(NQtot,nBand,nBand)
   
   i3=blockDim%x*(blockIdx%x-1)+threadIdx%x
   is1=blockIdx%y
   is2=blockIdx%z
   if ( i3 .le. NQtot ) then
      iq1=q12List(i3,1)
      iq2=q12List(i3,2) 
      f11=freq(iq1,is1)
      f22=freq(iq2,is2)
      V1=(/grpvelx(iq1,is1),grpvely(iq1,is1),grpvelz(iq1,is1)/)
      V2=(/grpvelx(iq2,is2),grpvely(iq2,is2),grpvelz(iq2,is2)/)
      dV=V1+V2
      if (tpSmearing .eq. 0) then
         sigma1=scoeff
         elseif (tpSmearing .eq. 1) then
             sigma1=0. 
             sigma1=sigma1+(dot_product(c1,dV)/nQx)**2+(dot_product(c2,dV)/nQy)**2+(dot_product(c3,dV)/nQz)**2
             sigma1=sqrt(sigma1/6d0)*scoeff
         endif
      sigma(i3,is1,is2,1)=sigma1
      dV=V1-V2
      if (tpSmearing .eq. 0) then 
         sigma2=scoeff
         elseif (tpSmearing .eq. 1) then 
             sigma2=0.
             sigma2=sigma2+(dot_product(c1,dV)/nQx)**2+(dot_product(c2,dV)/nQy)**2+(dot_product(c3,dV)/nQz)**2
             sigma2=sqrt(sigma2/6d0)*scoeff
         endif
      sigma(i3,is1,is2,2)=sigma2
      isScattList1(i3,is1,is2)=abs(( f11-f22+f00)/sigma1) .LE. 2.0 .AND. sigma1 .GE. sigmamin
      isScattList2(i3,is1,is2)=abs((-f11-f22+f00)/sigma2) .LE. 2.0 .AND. sigma1 .GE. sigmamin
      isScattList(i3,is1,is2)=isScattList1(i3,is1,is2) .OR. isScattList2(i3,is1,is2)
   endif
   return
   end subroutine Cal_IsScat_V3
   
   attributes(global) subroutine Cal_V3(NQtot,nTriple,nBand,nQx,nQy,nQz,isscattcount,nTemp,tpSmearing,RTA,f00,is0,&
      ee0,isscatt,q12List,Qpoint,freq,isScattList1,isScattList2,TempList,sigma,ev,CFCIndexR,CFCm,QQQ3ph,A,CFCIndextail)
   implicit none
   integer,value::NQtot,nTriple,nBand,nQx,nQy,nQz,isscattcount,nTemp,tpSmearing,is0
   logical,value::RTA
   real(8),value::f00
   integer,dimension(:,:),intent(in)::q12List(NQtot,2)
   integer,dimension(:,:),intent(in)::CFCIndextail(nTriple,3)
   integer,dimension(:,:),intent(in)::isscatt(NQtot*nBand*nBand,3)
   real(8),dimension(:),intent(in)::TempList(nTemp)
   real(8),dimension(:,:),intent(in):: freq(NQtot,nBand)
   real(8),dimension(:,:),intent(in)::Qpoint(nQx*nQy*nQz,3)
   real(8),dimension(:,:,:,:),intent(in)::CFCm(nTriple,3,3,3)
   real(8),dimension(:,:,:),intent(in)::CFCIndexR(nTriple,2,3)
   real(8),dimension(:,:,:,:),intent(in)::sigma(NQtot,nBand,nBand,2)
   complex*16,dimension(:),intent(in)::ee0(nBand)   
   complex*16,dimension(:,:,:),intent(in)::ev(NQtot,nBand,nBand)
   logical,dimension(:,:,:),intent(in)::isScattList1(NQtot,nBand,nBand),isScattList2(NQtot,nBand,nBand)
   integer::i1,i3,iq1,iq2,is1,is2,i2,iT,istat,j1
   real(8)::qq1(3),qq2(3)
   real(8)::f11,f22,temp,ato,GG1,GG2,GG3,dgg,nn0,nn1,nn2,sigma1,sigma2
   complex*16::coeff1,coeff2,cval
   real(8),parameter,constant::sigmamin=1d-8
   real(8),parameter,constant::hb=0.033123406258297d0,kb=6.333623318249880d-006,PI=3.141592653589793d0
   complex(8),parameter,constant:: um=(0d0,1d0)
   real(8),dimension(:,:),intent(out)::QQQ3ph(nBand,nTemp)
   real(8),dimension(:,:),intent(out)::A(nBand*NQtot,nTemp)
   i1=blockDim%x*(blockIdx%x-1)+threadIdx%x
   if( i1 .le. isscattcount ) then
      i3 =isscatt(i1,1)
      is1=isscatt(i1,2)
      is2=isscatt(i1,3)
      iq1=q12List(i3,1)
      iq2=q12List(i3,2)
      f11=freq(iq1,is1)
      f22=freq(iq2,is2)
      qq1=Qpoint(iq1,:)
      qq2=Qpoint(iq2,:)
      cval=0d0
      do i2=1,nTriple
         coeff1=exp(2d0*PI*um*( dot_product(qq1,CFCIndexR(i2,1,:))+dot_product(qq2,CFCIndexR(i2,2,:))))
      do j1=1,3
         coeff2=ee0(CFCIndextail(i2,1)+j1)    
         cval=cval+coeff1*coeff2*ev(iq1,CFCIndextail(i2,2)+1,is1)*&
         &  (CFCm(i2,j1,1,1)*ev(iq2,CFCIndextail(i2,3)+1,is2)+CFCm(i2,j1,1,2)*ev(iq2,CFCIndextail(i2,3)+2,is2)+&
         &  CFCm(i2,j1,1,3)*ev(iq2,CFCIndextail(i2,3)+3,is2))
         cval=cval+coeff1*coeff2*ev(iq1,CFCIndextail(i2,2)+2,is1)*&
         &  (CFCm(i2,j1,2,1)*ev(iq2,CFCIndextail(i2,3)+1,is2)+CFCm(i2,j1,2,2)*ev(iq2,CFCIndextail(i2,3)+2,is2)+&
         &  CFCm(i2,j1,2,3)*ev(iq2,CFCIndextail(i2,3)+3,is2))
         cval=cval+coeff1*coeff2*ev(iq1,CFCIndextail(i2,2)+3,is1)*&
         &  (CFCm(i2,j1,3,1)*ev(iq2,CFCIndextail(i2,3)+1,is2)+CFCm(i2,j1,3,2)*ev(iq2,CFCIndextail(i2,3)+2,is2)+&
         &  CFCm(i2,j1,3,3)*ev(iq2,CFCIndextail(i2,3)+3,is2))
      enddo
      enddo
      sigma1=sigma(i3,is1,is2,1)
      sigma2=sigma(i3,is1,is2,2)
      do iT=1,nTemp
         temp=TempList(iT)
         nn0=1d0/(exp(hb*f00/kb/temp)-1d0)
         nn1=1d0/(exp(hb*f11/kb/temp)-1d0)
         nn2=1d0/(exp(hb*f22/kb/temp)-1d0)
         IF (sigma1 .GE. sigmamin) THEN
            dgg=1d0/sqrt(PI)/sigma1*exp(-((f11-f22+f00)/sigma1)**2d0)                  
            GG1=(PI*hb/4d0/dble(NQtot)/(f11*f22*f00))*(abs(cval))**2d0*(nn0*nn1*(nn2+1d0))*dgg
            GG2=(PI*hb/4d0/dble(NQtot)/(f11*f22*f00))*(abs(cval))**2d0*((nn0+1d0)*(nn1+1d0)*nn2)*dgg
         ENDIF
         IF (sigma2 .GE. sigmamin) THEN
            dgg=1d0/sqrt(PI)/sigma2*exp(-((-f11-f22+f00)/sigma2)**2d0)
            GG3=(PI*hb/4d0/dble(NQtot)/(f11*f22*f00))*(abs(cval))**2d0*(nn0*(nn1+1d0)*(nn2+1d0))*dgg
         ENDIF
         if (.not. RTA) then
            istat=atomicadd(A((iq1-1)*nBand+is1,iT),GG1*dble(isScattList1(i3,is1,is2))+GG3*dble(isScattList2(i3,is1,is2)))
            istat=atomicadd(A((iq2-1)*nBand+is2,iT),GG2*dble(isScattList1(i3,is1,is2)))
         endif
         ato=GG1*dble(isScattList1(i3,is1,is2))+GG3*0.5d0*dble(isScattList2(i3,is1,is2))
         istat=atomicadd(QQQ3ph(is0,iT),ato)   
      enddo
   endif
return
end subroutine Cal_V3
end module Accrelate_V3