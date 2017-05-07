!Subroutine comes here 

!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
     use prec
     IMPLICIT NONE
     !Declarations
     INTEGER, INTENT(IN) :: n
     INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
     REAL(dp), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
     REAL(dp), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
     
     LOGICAL :: FLAG = .TRUE.
     INTEGER :: i, j, k, l
     REAL(dp) :: m
     REAL(dp), DIMENSION(n,2*n) :: augmatrix !augmented matrix
     
     !Augment input matrix with an identity matrix
     DO i = 1, n
        DO j = 1, 2*n
            IF (j <= n ) THEN
                  augmatrix(i,j) = matrix(i,j)
                ELSE IF ((i+n) == j) THEN
                  augmatrix(i,j) = 1d0
                Else
                 augmatrix(i,j) = 0d0
             ENDIF
         END DO
     END DO
     
     !Reduce augmented matrix to upper traingular form
     DO k =1, n-1
        IF (augmatrix(k,k) == 0d0) THEN
             FLAG = .FALSE.
             DO i = k+1, n
                IF (augmatrix(i,k) /= 0d0) THEN
                    DO j = 1,2*n
                       augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                    END DO
                    FLAG = .TRUE.
                    EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                    PRINT*, "Matrix is non - invertible"
                    inverse = 0d0
                    errorflag = -1
                    return
                ENDIF
             END DO
        ENDIF
        DO j = k+1, n
              m = augmatrix(j,k)/augmatrix(k,k)
              DO i = k, 2*n
                     augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
              END DO
        END DO
     END DO
     
     !Test for invertibility
     DO i = 1, n
         IF (augmatrix(i,i) == 0d0) THEN
                PRINT*, "Matrix is non - invertible"
                inverse = 0d0
                errorflag = -1
                return
        ENDIF
     END DO
     
     !Make diagonal elements as 1
     DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)
               augmatrix(i,j) = (augmatrix(i,j) / m)
        END DO
     END DO
     
     !Reduced right side half of augmented matrix to identity matrix
     DO k = n-1, 1, -1
        DO i =1, k
           m = augmatrix(i,k+1)
                DO j = k, (2*n)
                   augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                END DO
        END DO
     END DO
     
     !store answer
     DO i =1, n
        DO j = 1, n
          inverse(i,j) = augmatrix(i,j+n)
        END DO
     END DO
     errorflag = 0
END SUBROUTINE FINDinv
