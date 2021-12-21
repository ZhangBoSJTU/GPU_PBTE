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
!  Main program.

program GPU_PBTE
   use crystal,only:isKx,isKy,isKz,isCubic,isQuart,isISO,RTA,nTemp
   use Task
   use IFC
   use Disp
   use WriteRead
   use Symm
   use Job_Task
   implicit none
   integer::i,iT

   open(10,file='../input/jobtype.txt')
   read(10,*) JobTask
   close(10) 
   call Read_Input
   call Gen_QPoint
   call Read_HFC
   if (isCubic .eqv. .true.) call Read_CFC
   if (isQuart .eqv. .true.) call Read_QFC
   if (JobTask .eq. 1) then
      call Cal_Disp
      call Cal_GrpVel
      call Write_Data
      call Write_Identifier
   elseif (JobTask .eq. 2) then
      call Choose_Symm
      call Job_Task_RT
      call Read_Data
      do i=1,nJob
         call Cal_RT_Init(jobList(i))
         if (isCubic .eqv. .true. ) call Cal_RTGPU_PH3(jobList(i))
         if (isQuart .eqv. .true. ) call Cal_RTGPU_PH4(jobList(i))
         if (isIso .eqv. .true. ) call Cal_RT_iso(jobList(i))
         call Write_RT(jobList(i))
      enddo
   elseif (JobTask .eq. 3) then
      call Choose_Symm
      call Read_Data
      do iT=1,nTemp
         call Cal_RTA(iT)
         if (.not. RTA) then
            if (isKx) call Cal_TC(1,iT)
            if (isKy) call Cal_TC(2,iT)
            if (isKz) call Cal_TC(3,iT)
         endif
      enddo
   elseif (JobTask .eq. 4) then
      call Read_Data
      call Cal_DOS
   elseif (JobTask .eq. 5) then
      call Read_Data
      call Cal_DOS_3PH
   elseif (JobTask .eq. 6) then
      call Choose_symm
      call Read_Data
      call Cal_DOS_4PH
   endif

end program GPU_PBTE