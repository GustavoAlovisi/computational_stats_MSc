
############this script tests whether 1000 symmetrical NxN normal(0,1) matrixes are positive-definite or not
##########elapsed time for the 1000 matrixes tests are computed
library(Matrix)
library(tictoc)


##symmetrical NxN data generation function
gen_rnorm_sym_matrix <- function(N){
  normal_mat <- matrix(rnorm(N*N,mean=0,sd=1), N, N) 
  normal_mat <- as.matrix(Matrix::forceSymmetric(normal_mat))
  return(normal_mat)
}

###prop: a symmetrical NxN matrix is positive definite i.f.f all eigenvalues of matrix > 0.
##this is the most straightforward way to test for positive-definiteness.
##however, computing NxN eigenvalues scales very badly as matrix size grows.
#eigen(matrix)$values


###prop: a symmetrical NxN matrix is pos.def. i.f.f det(each Leading Principal Minor) > 0. 
##this alternative method uses determinant, is more efficient and scales not so bad as eigenvalues.

##leading principal minor method for a matrix 
LPMinors_det <- function(mat){
  for(i in 1:nrow(mat)){
    if(det(as.matrix(mat[1:i,1:i])) <= 0) return(FALSE) #if at least 1 det <= 0, matrix is not pos.def.
  }
  return(TRUE) #if every det > 0, matrix is pos.def. 
}

##########test: lets use our built function to test for a known pos.def. matrix: 
posdefMat <- matrix(0,2,2)
posdefMat[1,1] <- posdefMat[2,2] <- 5

LPMinors_det(posdefMat)
#TRUE

#########now lets test simulated symmetrical normal matrixes: 
posdefprop <- numeric(4) #initiating an array for proportion of pos.def. matrixes over 1000 matrixes generated 
#for each case

#######testing for 10x10 matrixes
tic('10x10') #initiating tic time counter
matrix_test <- integer(1000) #an array to save whether the j matrix is posdef (1) or not (0). 
for(j in 1:1000){
  matrixx <- gen_rnorm_sym_matrix(10)
  matrix_test[j] <- LPMinors_det(matrixx)
}
posdefprop[1] <- sum(matrix_test)/1000 #proportion of posdef (1) over 1k matrixes. 
toc(log = T) #time counter 

#######testing for 100x100 matrixes
tic('100x100')
matrix_test <- integer(1000)
for(j in 1:1000){
  matrixx <- gen_rnorm_sym_matrix(100)
  matrix_test[j] <- LPMinors_det(matrixx)
}
posdefprop[2] <- sum(matrix_test)/1000
toc(log = T)

#######testing for 1000x1000 matrixes
tic('1000x1000')
matrix_test <- integer(1000)
for(j in 1:1000){
  matrixx <- gen_rnorm_sym_matrix(1000)
  matrix_test[j] <- LPMinors_det(matrixx)
}
posdefprop[3] <- sum(matrix_test)/1000
toc(log = T)

#######testing for 10000x10000 matrixes
tic('10000x10000')
matrix_test <- integer(1000)
for(j in 1:1000){
  matrixx <- gen_rnorm_sym_matrix(10000)
  matrix_test[j] <- LPMinors_det(matrixx)
}
posdefprop[4] <- sum(matrix_test)/1000
toc(log = T)


#######proportion of pos.def. matrixes in each case. There were no pos.def. matrixes in simulated data. 
print(posdefprop)
#[1] 0 0 0 0

tic.log(format = T) #elapsed time log for each test 

#######below, elapsed time for a run is presented. As we can see, difference in time for 
#######1000x1000 and 10000x10000 test is almost 12000% bigger. 
#[[4]]
#[1] "10x10: 0.14 sec elapsed"

#[[5]]
#[1] "100x100: 0.61 sec elapsed"

#[[6]]
#[1] "1000x1000: 47.25 sec elapsed"

#[[7]]
#[1] "10000x10000: 5672.14 sec elapsed"

#####we ran optimizations using a PC with 16gb ram, RX5700XT graphic card and amd ryzen 5 3600