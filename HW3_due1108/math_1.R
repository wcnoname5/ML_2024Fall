# Create the matrix with 10 rows and 3 columns
# features in columns
matx <- matrix(c(1,2,3,
                 4,8,5,
                 3,12,9,
                 1,8,5,5,14,2,7,4,1,9,8,9,3,8,1,11,5,6,10,11,7),
               byrow = TRUE,
               nrow = 10) 
# Perform PCA using prcomp()
# Note: data were centered
pca_prcomp <- prcomp(matx, center = TRUE, scale. = FALSE)

# 1.  Principal Axes (in columns)
pca_prcomp$rotation |> xtable::xtable(digits = 3)

# 2.  principal components
t(pca_prcomp$x) |> xtable::xtable(digits = 3)
# 3. reconstruction error
matx.reconstruct <-
  pca_prcomp$x[, 1:2] %*% t(pca_prcomp$rotation[,1:2])
matx.colmean <- matrix(rep(colMeans(matx), 10), nrow = 10, byrow = T)
matx.reconstruct <- matx.reconstruct + matx.colmean
## Error
sum((matx - matx.reconstruct)^2)


# 3.

W <- matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
              1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
              0, 1, 1, 0, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
              0, 0, 0, 0, 0, 0, 1, 0, 1, 0),
            byrow = TRUE, nrow = 10)
# D <- diag(rowSums(W))
D <- diag(c(3, 3, 2, 2, 2, 1, 2, 3, 2, 2))
# D inverse
# D_inv <- diag(1/rowSums(W))
D_inv.sqrt <- diag(sqrt(1/(rowSums(W))))
D.sqrt <- diag(sqrt((rowSums(W))))
 
L <- D - W

eigen_decom  <- eigen(solve(D) %*% L)  # Solving D^(-1) * L
# Sort the eignevalues indices ascending
sorted_indices <- order(eigen_decom$values)

n_col <- ncol(eigen_decom$vectors)
# 3 smallest eigenvectors
psi <- eigen_decom$vectors[,sorted_indices[1:3]]

## (4)

psi2 <- eigen_decom$vectors[,sorted_indices[2:4]]

# z_df <- as.data.frame(psi2)

sum(diag(t(psi2) %*% L %*% psi2))

t(psi2) %*% D %*% psi2

# Orthonormalize Psi with respect to D
# Using the Cholesky decomposition of D for D^-1/2
D_chol <- chol(D)
D_inv_half <- solve(D_chol)

# Transform Psi to satisfy Psi^T D Psi = I_3
Psi <- D_inv_half %*% psi2
Psi <- qr.Q(qr(Psi))  # QR decomposition to ensure orthonormality

# Check conditions
trace_value <- sum(diag(t(Psi) %*% L %*% Psi))
identity_check <- t(Psi) %*% D %*% Psi  # Should be close to I_3




