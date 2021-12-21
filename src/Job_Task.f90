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
!  Jod tasks to parallel the Calculation.

module Job_Task
   implicit none
   contains
   subroutine Job_Task_RT
      use task
      use Crystal,only: RTA
      use symm,only:N0_8fbz,N0_ibz
      implicit none
      integer::N,nJobInNode,jobnumber,njob1,njob2,i,j,jobCount
      integer,dimension(:,:),allocatable::jobIndex
      integer,dimension(:),allocatable::QList
      if (RTA) then
         N=N0_ibz
      else
         N=N0_8fbz
      endif
      open(10,file='jobid.txt')
      read(10,*) nNode
      read(10,*) jobNumber
      close(10)
      allocate(QList(N))
      do i=1,N
         QList(i)=i
      enddo
      nJobInNode=N/nNode+1
      allocate(jobIndex(nJobInNode,nNode),jobList(nJobInNode))
      jobIndex=0
      jobList=0
      do i=1,N
         njob1=(i-1)/nNode+1
         njob2=i-(njob1-1)*nNode
         jobIndex(njob1,njob2)=1
      enddo
      jobCount=0
      do j=1,nNode
         do i=1,nJobInNode
            if (jobIndex(i,j) .eq. 1) then
               jobCount=jobCount+1
               jobIndex(i,j)=QList(jobCount)
            endif
         enddo
      enddo
      jobCount=0
      do i=1,nJobInNode
         if (jobIndex(i,jobNumber) .ne. 0) then
            jobList(i)=jobIndex(i,jobNumber)
            jobCount=jobCount+1
         endif
      enddo
      deallocate(jobList)
      allocate(jobList(jobCount))
      do i=1,jobCount
        jobList(i)=jobIndex(i,jobNumber)
      enddo
      nJob=jobCount
   end subroutine Job_Task_RT
end module Job_Task