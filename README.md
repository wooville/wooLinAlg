# wooLinAlg
A linear algebra library for C++ created while following the excellent series on the YouTube channel QuantitativeBytes, with additions and changes (ex inverse from adjoint and related functions):

www.youtube.com/c/QuantitativeBytes

https://www.youtube.com/playlist?list=PL3WoIG-PLjSv9vFx2dg0BqzDZH_6qzF8-

Functions are implemented in header files to allow for generic template data type T.

# Functions
wooPCA.h
Implementation of Principal Component Analysis (PCA).

https://youtu.be/ifxUSa5r_Ls

wooQR.h
Function to perform QR decomposition on the given matrix, returning an orthogonal matrix, Q, and an upper-triangular matrix, R. Uses the method of Householder reflections to perform the decomposition.

https://youtu.be/MR54VHqhROw

wooEIG.h
Functions for computing the eigenvectors and eigenvalues for a given matrix. Contains an implementation of the power iteration method for computing the dominant eigenvector, the inverse-power-iteration method and an implementation of the QR algorithm to estimate eigenvalue / eigenvector pairs for a given symmetric matrix.

https://youtu.be/hnLyWa2_hd8

https://youtu.be/tYqOrvUOMFc

wooLinSolve.h
Function for solving systems of linear equations. Uses an implementation of Gaussian elimination and back-substitution.

https://youtu.be/GKkUU4T6o08

https://youtu.be/NJIv0xH-S0I

https://youtu.be/KsrlAnEmRNE

wooLSQ.h
Function for computing the linear least squares solution to an over-determined system of linear equations.

https://youtu.be/4UVPXs3vIHk

https://youtu.be/fG1JXf7WSQw

wooMatrix.h
Class for handling matrices. Implements a number of useful functions:

wooMatrix - inverse()
Compute the inverse of the matrix using the Gauss-Jordan elimination method.

https://youtu.be/wOlG_fnd3v8

https://youtu.be/AEuNHdgn-R8

https://youtu.be/JWM8Y8b1ZVQ

wooMatrix - inverseFromAdj()
Compute the inverse of the matrix using the adjoint and determinant.

wooMatrix - rowEchelon()
Convert the matrix to row echelon form.

wooMatrix - transpose()
Transpose the matrix.

wooMatrix - determinant()
Compute the determinant of the matrix.

https://youtu.be/YVk0nYrwBb0

wooVector.h
Class for handling vectors. Implements a number of useful functions:

https://youtu.be/YfWX-EsvX7c

https://youtu.be/c5AB5T7LBCI

wooVector - normalized()
Returns a normalized copy of the vector.

wooVector - normalize()
Normalizes the vector 'in-place'.

wooVector - dot()
Computes the vector dot product.

wooVector - cross()
Computes the vector cross product.
