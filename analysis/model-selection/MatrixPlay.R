### Create H, C, L

N <- 2 # Number of grid points
k <- 4 # Number of predictor variables

n <- N^k # Size of matrices

iA <- cbind(1:n, 1+mod(0:(n-1), N))

mtx <- data.matrix(
        expand.grid(
          h=seq(from=1, to=2, length.out=N), # h = h_apical.x
          w=seq(from=3, to=4, length.out=N), # w = herb_avg.x
          c=seq(from=5, to=6, length.out=N),
          l=seq(from=7, to=8, length.out=N)
          )
        )

# Initialize with zeros
H <- array(0, dim=c(n, N))
C <- array(0, dim=c(n, N))
L <- array(0, dim=c(n, N))

H[iA] <- mtx[,"w"]
C[iA] <- mtx[,"c"]
L[iA] <- mtx[,"l"]


### Create Survival, Flowering, Pods

# Commented out because function won't really work (more of a prototype)

# imtx <- data.matrix( expand.grid(h=1:N, w=1:N, c=1:N, l=1:N) )  # Matrix of indices
# K[imtx] <- apply(mtx, 1, function(x) predict(h_apical = x["h"], herb_avg = x["w"], card = x["c"], lma = x["l"]))
# c(K) will give you ONE row of the matrix that you want.  Replicate this row N times
