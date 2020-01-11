C     shown to work on 1/10/2020
C     this version locked from any other changes 
      PROGRAM ParGauss
      include 'mpif.h'
      integer ProcNum, ProcRank, ierr, root, i, j, i2, j2 , j3 , i6, j6 
      integer Size ! size of the square matrix A
      integer RowNum ! Number of  the matrix rows 
      integer RestRows !Number of rows, that haven't been distributed 
      integer count, count_rate, count_max !clock stuff
      integer PivotPos, Iter, RowIndex, temp1

      integer, dimension (:), allocatable :: pParallelPivotPos
      integer, dimension (:), allocatable :: pProcPivotIter
      integer, dimension (:), allocatable :: pProcInd
      integer, dimension (:), allocatable :: pProcNum


      double precision , dimension (:), allocatable :: pMatrix
      double precision , dimension (:), allocatable :: pVector
      double precision , dimension (:), allocatable :: pResult
      double precision , dimension (:), allocatable :: pProcRows
      double precision , dimension (:), allocatable :: pProcVector
      double precision , dimension (:), allocatable :: pProcResult
      double precision , dimension (:), allocatable ::  pPivotRow
C     Buffer for storing the vector, that is a result of multiplication 
C     of the linear system matrix by the vector of unknowns */
      double precision , dimension (:), allocatable :: pRightPartVector
      double precision , dimension (:), allocatable :: ResultVector

      integer equal       
      double precision Accuracy ! Comparison accuracy

      integer IterProcRank, IterPivotPos      
      double precision IterResult, val

      integer, dimension (:), allocatable :: pSendNum     ! Number of the elements sent to the process
      integer, dimension (:), allocatable :: pSendInd     ! Index of the first data element sent 
                      ! to the process
C      struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;
      double precision in (2), out(2), MaxValue, multiplier


      data root/0/       
C     size of the matrix 
      Size = 6
      call MPI_INIT( ierr ) 
      call MPI_COMM_RANK( MPI_COMM_WORLD, ProcRank, ierr ) 
C     Number of the available processes // Rank of the current process
      call MPI_COMM_SIZE( MPI_COMM_WORLD, ProcNum, ierr ) 
C      write(*,*) ProcRank
      if (ProcRank .eq. 0) then
         write(*,*) 'Parallel Gauss algorithm for solving linear 
     & systems\n'
      end if 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Memory allocation and data initialization 
C      call ProcessInitialization (pMatrix, pVector, 
C     & pResult, pProcRows, pProcVector, pProcResult, Size, RowNum,
C     & pParallelPivotPos, pProcPivotIter, pProcInd, pProcNum, ProcNum)

C     Subroutine for memory allocation and data initialization
C     subroutine ProcessInitialization
C      Loop variable
C      write(*,*) 'in allocation '
      RestRows = Size
      do 70 i=0,(ProcRank-1)
         RestRows = RestRows-RestRows/(ProcNum-i)
   70 continue
      RowNum = RestRows/(ProcNum-ProcRank)
C      write(*,*) RowNum, ProcRank
      allocate (pProcRows(RowNum*Size))
      allocate (pProcVector(RowNum))
      allocate(pProcResult (RowNum))
      allocate(pParallelPivotPos(Size))      
      allocate(pProcPivotIter(RowNum))
  
      allocate (pProcInd(ProcNum))
      allocate(pProcNum(ProcNum))

      do 71 i=1, (RowNum)
         pProcPivotIter(i) = -1
   71 continue 

      if (ProcRank .eq. 0)  then
         allocate (pMatrix(Size*Size))
         allocate (pVector(Size))
         allocate (pResult(Size))
         allocate (ResultVector(Size))
C     DummyDataInitialization (pMatrix, pVector, Size);
C         call RandomDataInitialization(pMatrix, pVector, Size)
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Function for random definition of matrix and vector elements 
C     Loop variables
         CALL SYSTEM_CLOCK(count, count_rate, count_max)
         CALL srand(count)
         do 80 i=1,Size
            pVector(i) = rand(0)*10.0
            do 90 j= 1, Size
               if (j .le. i) then
                val = (rand(0)*10.0)
                i2=(i-1)*Size + j 
                pMatrix (i2) = val
               else
                  pMatrix((i-1)*Size+j) = 0.0
               end if 
   90       continue
   80    continue 
      do 151 i6=1,Size
          write(*,*) pMatrix( (Size*(i6-1)+1): (Size*(i6-1)+Size) ) 
          write(*,*) ' '
  151 continue 
      write(*,*) pVector
      end if 
C end of initalization 
      write(*,*) 'process initializated', ProcRank
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      call DataDistribution(pMatrix, pProcRows, pVector, pProcVector, 
C     & Size, RowNum, pProcInd, ProcNum)

      RestRows=Size ! Number of rows, that have not been 
C     distributed yet
C      Alloc memory for temporary objects
      allocate (pSendInd (ProcNum))
      allocate (pSendNum(ProcNum))

C      Define the disposition of the matrix rows for the current process
      RowNum = (Size/ProcNum)
      pSendNum(1) = RowNum*Size
      pSendInd(1) = 0
      do 100 i=1, (ProcNum-1)
         RestRows = RestRows - RowNum
         RowNum = RestRows/(ProcNum-i)
         pSendNum(i+1) = RowNum*Size
         pSendInd(i+1) = pSendInd(i)+pSendNum(i)
  100 continue 

C      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
C     write (*,*) pSendNum, 'Send Rank ' ,ProcRank
C      write (*,*) pSendInd, 'Ind Rank ' ,ProcRank
C      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

C      Scatter the rows
      call MPI_Scatterv(pMatrix,pSendNum,pSendInd,MPI_DOUBLE_PRECISION, 
     & pProcRows,pSendNum(ProcRank+1), MPI_DOUBLE_PRECISION, 
     & 0, MPI_COMM_WORLD, ierr)
C      Define the disposition of the matrix rows for current process
      RestRows = Size
      pProcInd(1) = 0
      pProcNum(1) = Size/ProcNum
      do 110 i=1, (ProcNum-1)
         RestRows = RestRows- pProcNum(i)
         pProcNum(i+1) = RestRows/(ProcNum-i)
         pProcInd(i+1) = pProcInd(i)+pProcNum(i)
  110 continue 

      CALL MPI_Scatterv(pVector,pProcNum,pProcInd, MPI_DOUBLE_PRECISION, 
     & pProcVector,pProcNum(ProcRank+1),MPI_DOUBLE_PRECISION,
     & 0,MPI_COMM_WORLD, ierr)

C      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
C      write (*,*) pProcRows, 'Rank ' ,ProcRank
C      Write(*,*) '  '
C      write (*,*) pProcVector, ' Vector Rank ' ,ProcRank
C      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

C Free the memory
      deallocate(pSendNum)
      deallocate(pSendInd)
      write(*,*) 'data distributed', ProcRank
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     The execution of the parallel Gauss algorithm 

C     call ParallelResultCalculation(pProcRows, pProcVector, 
C     & pProcResult, Size, RowNum, pProcPivotIter, pParallelPivotPos)

CCCCCCCC      
C     call ParallelGaussianElimination (pProcRows, pProcVector, Size, 
C     & RowNum, pProcPivotIter, pParallelPivotPos, pProcInd)
C     pPivotRow is used for storing the pivot row and the corresponding 
      allocate(pPivotRow(Size+1))


C     The iterations of the Gaussian elimination stage 
      do 10 i=1, (Size)
         MaxValue = 0.0
C      Calculating the local pivot row 
         do 20 j=1, RowNum
            if ((pProcPivotIter(j) .eq. -1) .AND. (MaxValue .lt. 
     & abs( pProcRows((j-1)*Size+i) ) )) then 
               MaxValue = abs(pProcRows((j-1)*Size+i)) 
               PivotPos = j;
            end if
   20    continue   


         in(1) = MaxValue
         in(2) = ProcRank
C         write(*,*) i, 'int',in(1), in(2)
C        Finding the pivot process (process with the maximum value of MaxValue) 
         call MPI_Allreduce(in, out, 1 , MPI_2DOUBLE_PRECISION, 
     & MPI_MAXLOC, MPI_COMM_WORLD, ierr)
C      write(*,*) i, 'out',out(1), out(2)
C        Broadcasting the pivot row 
         if ( ProcRank .eq. out(2) ) then
C           iteration number 
            pProcPivotIter(PivotPos)= i 
            pParallelPivotPos(i)= pProcInd(ProcRank+1) + PivotPos
C            write(*,*) 'pPPP',pParallelPivotPos(i),pProcInd(ProcRank+1)
         end if 
         temp1 = out(2)
         call MPI_Bcast(pParallelPivotPos(i), 1, MPI_INTEGER, temp1, 
     & MPI_COMM_WORLD, ierr)
C         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
C         write(*,*) pParallelPivotPos(i),'rank',ProcRank,'pP',PivotPos
C         call MPI_BARRIER(MPI_COMM_WORLD, ierr)

         if ( ProcRank .eq. temp1 ) then
C        Fill the pivot row 
            do 30 j=1, Size
               pPivotRow(j) = pProcRows((PivotPos-1)*Size + j)
   30       continue 
         end if 

         pPivotRow(Size+1) = pProcVector(PivotPos)
         call MPI_Bcast(pPivotRow, Size+1,MPI_DOUBLE_PRECISION, temp1, 
     & MPI_COMM_WORLD, ierr)


CCCCCC
C         call ParallelEliminateColumns(pProcRows, pProcVector, 
C     & pPivotRow, Size, RowNum, i)
C     Fuction for the column elimination
         do 160 i2=1, RowNum
            if (pProcPivotIter(i2) .eq. -1) then
               multiplier = pProcRows((i2-1)*Size+i)/pPivotRow(i)
C               write(*,*) 'multi',multiplier, i2, 'Rank', ProcRank
               do 170 j2= i, Size
                  pProcRows((i2-1)*Size+j2)=pProcRows((i2-1)*Size+j2)- 
     & pPivotRow(j2)*multiplier
  170          continue
               pProcVector(i2) = pProcVector(i2)-
     & pPivotRow(Size+1)*multiplier
            end if            
  160    continue
C      write(*,*)' Vector Rank ' ,ProcRank, i, '--------'
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
CCCCCCCCCC

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         write (*,*) 'Rank ' ,ProcRank, 'pPR', pProcRows, 
     & 'pPV', pProcVector, i
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   10 continue
C      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c      write (*,*) 'Rank ' ,ProcRank, 'pPR', pProcRows, 
c     & 'pPV', pProcVector
c      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

C     end of the subroutine 
      write(*,*) 'Gaussian elimination done'
      


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     call ParallelBackSubstitution (pProcRows, pProcVector, 
C     & pProcResult, Size, RowNum, pParallelPivotPos, pProcPivotIter)
C      write(*,*) 'I amhere',pProcInd,'Rank',ProcRank, pParallelPivotPos
C     Iterations of the back substitution stage 
      do 40 i= Size, 1, -1 
C     Calculating the rank of the process, which holds the pivot row 
CCCCCCCCCCCCCC
C        call FindBackPivotRow(pParallelPivotPos(i), Size, 
C     & IterProcRank, IterPivotPos) 
         RowIndex = pParallelPivotPos(i)
C    Subroutine for finding the pivot row of the back substitution
C      subroutine FindBackPivotRow(RowIndex, Size, IterProcRank,
C    & IterPivotPos, pProcInd)

         do 60 i2=1, (ProcNum-1)
C            write (*,*) RowIndex, pProcInd(i2)
            if (pProcInd(i2) .lt. RowIndex) then
               if (RowIndex .le. pProcInd(i2+1)) then
                  IterProcRank = i2
C                  write(*,*) 'in', RowIndex, i2
               end if  
            end if
   60    continue
            
         if (RowIndex .gt. pProcInd(ProcNum)) then
            IterProcRank = ProcNum
         end if

C         write(*,*) 'I am here', IterProcRank , ProcRank
         IterPivotPos = RowIndex - pProcInd(IterProcRank)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
C         write(*,*) 'IPR ', IterProcRank, 'IPP', IterPivotPos, i
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
CCCCCCCCCC
C     Calculating the unknown 
         if (ProcRank .eq. (IterProcRank-1)) then
            IterResult = pProcVector(IterPivotPos)/pProcRows
     & ((IterPivotPos-1)*Size+i)
            pProcResult(IterPivotPos) = IterResult
C            write(*,*) 'I am here', ProcRank, IterResult
C            write(*,*) 'I am here', pProcVector(IterPivotPos),pProcRows
C     & ((IterPivotPos-1)*Size+i)
         end if 
C     Broadcasting the value of the current unknown 
         call MPI_Bcast(IterResult, 1, MPI_DOUBLE_PRECISION, 
     & (IterProcRank-1), MPI_COMM_WORLD, ierr)
C         write(*,*) ' IterResult',IterResult,RowIndex         
C     Updating the values of the vector b 
         do 50 j=1, RowNum
C            write(*,*) ' pPP j',pProcPivotIter(j) ,j  
            if ( pProcPivotIter(j) .lt. i ) then
               val = pProcRows((j-1)*Size + i) * IterResult
               pProcVector(j)=pProcVector(j) - val
C               write(*,*) ' Val',val,j  
            end if 
   50    continue      
   40 continue
C     end of the subroutine 

      write(*,*) 'BackSubstitution done'
C      write(*,*) 'Result calculated', pProcResult, ProcRank
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     call ResultCollection(pProcResult, pResult)
C     Gather the whole result vector on every processor
      call MPI_Gatherv(pProcResult, pProcNum(ProcRank+1),  
     & MPI_DOUBLE_PRECISION,pResult, pProcNum, pProcInd,
     &  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      write(*,*) 'Result gathered'
      write(*,*) pResult
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     call TestResult(pMatrix, pVector, pResult, Size, 
C    & pParallelPivotPos)

C     Function for testing the result

      equal = 0        
      Accuracy = 1.e-6 ! Comparison accuracy

      if (ProcRank .eq. 0) then
         allocate (pRightPartVector(Size))
         do 120 i=1, Size
            pRightPartVector(i) = 0
            do 130 j=1, Size
               pRightPartVector(i) = pRightPartVector(i) + 
     & pMatrix((i-1)*Size+j)*pResult(pParallelPivotPos(j))
               ResultVector(j) = pResult(pParallelPivotPos(j))
  130       continue
  120    continue 
      end if 

      do 140 i=1,Size
C         write(*,*) abs(pRightPartVector(i)-pVector(i))
         if ( abs(pRightPartVector(i)-pVector(i)) .gt. Accuracy) then
            equal = 1
         end if 
  140 continue

      if (equal .eq. 1) then
         write(*,*) 'The result of the parallel Gauss algorithm is 
     & NOT correct. Check your code'
      else
         write(*,*) 'The result of the parallel Gauss algorithm 
     & is correct'
         write(*,*) ResultVector
         deallocate(pRightPartVector)
      end if 

C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

CC    Computational process termination
C     call ProcessTermination(pMatrix, pVector, pResult, pProcRows,
C    & pProcVector, pProcResult, pProcPivotIter, 
C    & pParallelPivotPos, pProcInd)
      if (ProcRank .eq. 0)  then
         deallocate(pMatrix)
         deallocate(pVector)
         deallocate(pResult)
      end if 

      deallocate(pProcRows)
      deallocate(pProcVector)
      deallocate(pProcResult)

      deallocate(pParallelPivotPos)
      deallocate(pProcPivotIter)

      deallocate(pProcInd)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call MPI_Finalize(ierr)

      end program 



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC CUSTOM SUBROUTINES for Parallelization CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 


C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



         


C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    





C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

