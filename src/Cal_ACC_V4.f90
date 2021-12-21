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
!  Accelerate the calculation of four-phonon scattering process.

module Accrelate_V4
    contains
    subroutine Cal_V4_pre(ee0,ee3)
        use IFC,only:nQuart,QFC_complex,QFCm,QFCIndextail,QFC_complex_sum
        use Crystal,only:nBand
        implicit none
        complex(8)::ee0(nBand),ee3(nBand)
        integer::i,j1,j2,j3,j4
        do i=1,nQuart
            do j1=1,3 
            do j2=1,3
            do j3=1,3
            do j4=1,3
                QFC_complex(i,j1,j2,j3,j4)=QFCm(i,j1,j2,j3,j4)& 
                    &   *ee0(QFCIndextail(i,1)+j1)*ee3(QFCIndextail(i,4)+j4)
            enddo
            enddo
            enddo
            enddo
            QFC_complex_sum(1,i)=sum(QFC_complex(i,:,1,1,1:3))
            QFC_complex_sum(2,i)=sum(QFC_complex(i,:,1,2,1:3))
            QFC_complex_sum(3,i)=sum(QFC_complex(i,:,1,3,1:3))
            QFC_complex_sum(4,i)=sum(QFC_complex(i,:,2,1,1:3))
            QFC_complex_sum(5,i)=sum(QFC_complex(i,:,2,2,1:3))
            QFC_complex_sum(6,i)=sum(QFC_complex(i,:,2,3,1:3))
            QFC_complex_sum(7,i)=sum(QFC_complex(i,:,3,1,1:3))
            QFC_complex_sum(8,i)=sum(QFC_complex(i,:,3,2,1:3))
            QFC_complex_sum(9,i)=sum(QFC_complex(i,:,3,3,1:3))
        enddo
        return
    end subroutine Cal_V4_pre

    attributes(global) subroutine Cal_IsScat_V4(NQtot,nBand,nQx,nQy,nQz,tpSmearing,f00,f33,scoeff4,&
        q12List,freq,grpvelx,grpvely,grpvelz,c1,c2,c3,isScattList1,isScattList2,isScattList3,isScattList,sigma)
    implicit none
    integer,value::NQtot,nBand,nQx,nQy,nQz,tpSmearing
    real(8),value::f00,f33,scoeff4
    integer,dimension(:),intent(in)::q12List(NQtot,2)
    real(8),dimension(:),intent(in)::c1(3),c2(3),c3(3)
    real(8),dimension(:,:),intent(in):: freq(NQtot,nBand)
    real(8),dimension(:,:),intent(in)::grpvelx(NQtot,nBand),grpvely(NQtot,nBand),grpvelz(NQtot,nBand)
    integer::iq1,iq2,is1,is2
    real(8)::V1(3),V2(3),dV(3)
    real(8)::f11,f22,sigma1,sigma2
    real(8),parameter,constant::sigmamin=1d-8
    real(8),dimension(:,:,:,:),intent(out)::sigma(NQtot,nBand,nBand,2)
    logical,dimension(:,:,:),intent(out)::isScattList1(NQtot,nBand,nBand),&
        &   isScattList2(NQtot,nBand,nBand),isScattList3(NQtot,nBand,nBand),isScattList(NQtot,nBand,nBand)
    iq1=blockDim%x*(blockIdx%x-1)+threadIdx%x
    is1=blockIdx%y
    is2=blockIdx%z
    if ( iq1 .le. NQtot ) then
        iq2=q12List(iq1,2)
        f11=freq(iq1,is1)
        V1=(/grpvelx(iq1,is1),grpvely(iq1,is1),grpvelz(iq1,is1)/)
        f22=freq(iq2,is2)
        V2=(/grpvelx(iq2,is2),grpvely(iq2,is2),grpvelz(iq2,is2)/)
        dV=V1+V2 
        if (tpSmearing .eq. 0) then
         sigma1=scoeff4
         elseif (tpSmearing .eq. 1) then
             sigma1=0. 
             sigma1=sigma1+(dot_product(c1,dV)/nQx)**2+(dot_product(c2,dV)/nQy)**2+(dot_product(c3,dV)/nQz)**2
             sigma1=sqrt(sigma1/6d0)*scoeff4
         endif
        sigma(iq1,is1,is2,1)=sigma1
        dV=V1-V2
        if (tpSmearing .eq. 0) then 
         sigma2=scoeff4
         elseif (tpSmearing .eq. 1) then 
             sigma2=0.
             sigma2=sigma2+(dot_product(c1,dV)/nQx)**2+(dot_product(c2,dV)/nQy)**2+(dot_product(c3,dV)/nQz)**2
             sigma2=sqrt(sigma2/6d0)*scoeff4
         endif
        sigma(iq1,is1,is2,2)=sigma2
        isScattList1(iq1,is1,is2) = abs((f00+f11+f22-f33)/sigma1) .le. 2.0 .and. sigma1 .ge. sigmamin
        isScattList2(iq1,is1,is2) = abs((f00-f11-f22-f33)/sigma1) .le. 2.0 .and. sigma1 .ge. sigmamin
        isScattList3(iq1,is1,is2) = abs((f00+f11-f22-f33)/sigma2) .le. 2.0 .and. sigma2 .ge. sigmamin
        isScattList(iq1,is1,is2)=isScattList1(iq1,is1,is2) .or. isScattList2(iq1,is1,is2) .or. isScattList3(iq1,is1,is2) 
    endif                 
    return
    end subroutine Cal_IsScat_V4

    attributes(global) subroutine Cal_V4(NQtot,nQuart,nBand,nQx,nQy,nQz,isscattcount,tpSmearing,is0,iq3,is3,f00,f33,nTemp,RTA,&
    isscatt,ev,QFCIndexR,QFC_complex_sum,QFCIndextail,q12List,qq3,Qpoint,freq,isScattList1,isScattList2,isScattList3,AA,QQQ4ph,temp0,sigma)
    implicit none
    integer,value::NQtot,nBand,nQuart,nQx,nQy,nQz,isscattcount,nTemp,tpSmearing,is0,iq3,is3
    real(8),value::f00,f33
    logical,value::RTA
    integer,dimension(:,:),intent(in)::q12List(NQtot,2)
    integer,dimension(:,:),intent(in)::QFCIndextail(nQuart,4)
    integer,dimension(:,:),intent(in)::isscatt(NQtot*nBand*nBand,3)
    real(8),dimension(:),intent(in)::qq3(3)    
    real(8),dimension(:),intent(in)::temp0(nTemp)
    real(8),dimension(:,:),intent(in):: freq(NQtot,nBand)
    real(8),dimension(:,:),intent(in)::Qpoint(nQx*nQy*nQz,3)
    real(8),dimension(:,:,:),intent(in)::QFCIndexR(nQuart,3,3)
    real(8),dimension(:,:,:,:),intent(in)::sigma(NQtot,nBand,nBand,2)
    complex(8),dimension(:,:),intent(in)::QFC_complex_sum(9,nQuart)
    complex*16,dimension(:,:,:),intent(in)::ev(NQtot,nBand,nBand)
    logical,dimension(:,:,:),intent(in)::isScattList1(NQtot,nBand,nBand),&
        &   isScattList2(NQtot,nBand,nBand),isScattList3(NQtot,nBand,nBand)
    integer::i1,iq1,iq2,is1,is2,i2,iT,istat
    real(8)::qq1(3),qq2(3)
    real(8)::f11,f22,temp,ato,GG1,GG2,GG3,dgg,nn0,nn1,nn2,nn3,sigma1,sigma2
    complex(8)::coeff1,cval
    real(8),parameter,constant::hb=0.033123406258297d0,kb=6.333623318249880d-006,PI=3.141592653589793d0
    complex(8),parameter,constant:: um=(0d0,1d0)
    real(8),dimension(:,:),intent(out)::QQQ4ph(nBand,nTemp)
    real(8),dimension(:,:,:),intent(out)::AA(nBand*NQtot,nBand,nTemp)
    
    i1=blockDim%x*(blockIdx%x-1)+threadIdx%x
    if( i1 .le. isscattcount ) then
        iq1=isscatt(i1,1)
        iq2=q12List(iq1,2)
        is1=isscatt(i1,2)
        is2=isscatt(i1,3)
        f11=freq(iq1,is1)
        f22=freq(iq2,is2)
        qq1=Qpoint(iq1,:)
        qq2=-Qpoint(iq2,:) 
        cval=0d0
        do i2=1,nQuart
            coeff1=exp(2d0*PI*um*(dot_product(qq1,QFCIndexR(i2,1,:))+dot_product(qq2,QFCIndexR(i2,2,:))+dot_product(qq3,QFCIndexR(i2,3,:))))
            cval=cval+coeff1*(ev(iq1,QFCIndextail(i2,2)+1,is1)*(conjg(ev(iq2,QFCIndextail(i2,3)+1,is2))*QFC_complex_sum(1,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+2,is2))*QFC_complex_sum(2,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+3,is2))*QFC_complex_sum(3,i2))+ &
            ev(iq1,QFCIndextail(i2,2)+2,is1)*(conjg(ev(iq2,QFCIndextail(i2,3)+1,is2))*QFC_complex_sum(4,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+2,is2))*QFC_complex_sum(5,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+3,is2))*QFC_complex_sum(6,i2))+ &
            ev(iq1,QFCIndextail(i2,2)+3,is1)*(conjg(ev(iq2,QFCIndextail(i2,3)+1,is2))*QFC_complex_sum(7,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+2,is2))*QFC_complex_sum(8,i2)+conjg(ev(iq2,QFCIndextail(i2,3)+3,is2))*QFC_complex_sum(9,i2))  )
        enddo

        sigma1=sigma(iq1,is1,is2,1)
        sigma2=sigma(iq1,is1,is2,2)
        do iT=1,nTemp
            temp=temp0(iT)
            nn0=1d0/(exp(hb*f00/kb/temp)-1d0)
            nn3=1d0/(exp(hb*f33/kb/temp)-1d0)
            nn1=1d0/(exp(hb*f11/kb/temp)-1d0)
            nn2=1d0/(exp(hb*f22/kb/temp)-1d0)
            IF (sigma1 .NE. 0d0) THEN
            dgg=1d0/sqrt(PI)/sigma1*exp(-((f00+f11+f22-f33)/sigma1)**2d0)
            GG1=(PI*hb*hb/8d0/dble(NQtot*NQtot)/(f00*f11*f22*f33))*(abs(cval))**2d0*(nn0*nn1*nn2*(nn3+1d0))*dgg
            dgg=1d0/sqrt(PI)/sigma1*exp(-((f00-f11-f22-f33)/sigma1)**2d0)
            GG2=(PI*hb*hb/8d0/dble(NQtot*NQtot)/(f00*f11*f22*f33))*(abs(cval))**2d0*(nn0*(nn1+1d0)*(nn2+1d0)*(nn3+1d0))*dgg
            ENDIF
            IF (sigma2 .NE. 0d0) THEN
            dgg=1d0/sqrt(PI)/sigma2*exp(-((f00+f11-f22-f33)/sigma2)**2d0)
            GG3=(PI*hb*hb/8d0/dble(NQtot*NQtot)/(f00*f11*f22*f33))*(abs(cval))**2d0*(nn0*nn1*(nn2+1d0)*(nn3+1d0))*dgg
            ENDIF
            ato=GG1/2d0*dble(isScattList1(iq1,is1,is2))+GG2/6d0*dble(isScattList2(iq1,is1,is2))+GG3/2d0*dble(isScattList3(iq1,is1,is2))
            if (.not. RTA) then
            istat=atomicadd(AA((iq1-1)*nBand+is1,is0,iT),ato)
            istat=atomicadd(AA((iq2-1)*nBand+is2,is0,iT),-ato)
            istat=atomicadd(AA((iq3-1)*nBand+is3,is0,iT),-ato)
            endif
            istat=atomicadd(QQQ4ph(is0,iT),ato)
        enddo
    endif
    return
    end subroutine Cal_V4
end module Accrelate_V4