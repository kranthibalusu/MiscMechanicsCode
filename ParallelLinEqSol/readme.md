# Parallel solution for a system of linear equations, dense matrices

Contains 
1) ParGaussV2.for
2) inpScript.sh

This code can solve a system of linear equations, parallel algorithm. 
Gaussian elimination procedure. Applicable to all square matrices, including dense matrices. 
Distributed memory type parallel implementation. 
Used MPI on intel MPI/2017 compiler. 
Implemented in Fortran77. No other Math libraries (such as ScaLapack) needed. 

Time saving for large number of equations. 

Inspired by algorithm by Gergel V.P.  
C++ code in http://www.hpcc.unn.ru/mskurs/RUS/PPT/CODES/ParallelGauss.html#main


Feel free to contact me at kranthibalusu@gmail.com for future code enhancment. 
