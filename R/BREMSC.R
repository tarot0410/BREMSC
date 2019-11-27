## usethis namespace: start
#' @useDynLib BREMSC, .registration = TRUE
## usethis namespace: end
NULL
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

# Initialize alpha
#' @import dirmult
#' @import stats
EM_initial_alpha <- function(data, clusters_initial, method_alpha_intial)
{
  num_cluster <- length(unique(clusters_initial))
  cluster <- matrix(,num_cluster,1)
  for (i in 1:num_cluster){
    cluster[i,] <-  length(which(clusters_initial==i))/ncol(data)
  }
  cluster <- as.numeric(cluster)
  sort <- rank(cluster,ties.method="random")
  label.new <- clusters_initial
  for(j in 1:num_cluster){
    label.new[which(clusters_initial==j)] <- sort[j]
  }

  location <- list()
  for (m in 1:num_cluster){
    location[[m]] <- which(label.new==m)
  }
  newdata <- list()
  for (m in 1:num_cluster){
    newdata[[m]] <- data[,location[[m]]]
  }

  p <- matrix(,nrow(data),num_cluster)
  for (m in 1:num_cluster){
    sum <- sum(as.vector(newdata[[m]]))
    for ( i in 1:nrow(data)){
      p[i,m] <- sum(newdata[[m]][i,])/sum
    }
  }

  new.data <- list()
  for (m in 1:num_cluster){
    new.data[[m]] <- newdata[[m]][rowSums(newdata[[m]]) != 0, colSums(newdata[[m]]) != 0]
  }
  C <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    C[m,1] <- ncol(new.data[[m]])
  }
  G <- matrix(,num_cluster,1)
  for (m in 1:num_cluster){
    G[m,1] <- nrow(new.data[[m]])
  }

  ppp <- list()
  for (m in 1:num_cluster){
    ppp[[m]] <- matrix(,G[m,1],C[m,1])
    for (j in 1:C[m,1]){
      tmp = sum(new.data[[m]][,j])
      ppp[[m]][,j]=new.data[[m]][,j]/tmp
    }
  }

  pp <- list()
  for (m in 1:num_cluster){
    pp[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      pp[[m]][i,1] <- mean(ppp[[m]][i,])
    }
  }

  v <- list()
  for (m in 1:num_cluster){
    v[[m]] <- matrix(,G[m,1],1)
    for ( i in 1:G[m,1]){
      v[[m]][i,1] <- var(ppp[[m]][i,])
    }
    v[[m]][which(v[[m]]==0),1] <- mean(v[[m]],na.rm=T)
  }

  s <- matrix(,num_cluster,1)
  for(m in 1:num_cluster){
    sum <- 0
    for (i in 1:(G[m,1]-1)){
      tmp <- log( ( pp[[m]][i,1]*(1-pp[[m]][i,1])/ v[[m]][i,1] ) -1 )
      sum <- sum+tmp
    }
    s[m,1] <- exp( (1/(G[m,1]-1))*sum )
  }

  new.alpha <- list()
  for(m in 1:num_cluster){
    new.alpha[[m]] <- s[m,1]*p[,m]
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }

  new_alpha_R <- matrix(,num_cluster,nrow(data))
  for (m in 1:num_cluster){
    new_alpha_R[m,] <- new.alpha[[m]]
  }

  new.alpha <- list()
  for (m in 1:num_cluster){
    mom <- weirMoM(t(new.data[[m]]), se=FALSE)
    if (mom <= 0) {mom <- 0.005}
    initscalar <- (1 - mom)/mom
    new.alpha[[m]] <- p[,m]*initscalar
    new.alpha[[m]][which(new.alpha[[m]]==0)] <- 0.000001
  }

  new_alpha_W <- matrix(,num_cluster,nrow(data))
  for (m in 1:num_cluster){
    new_alpha_W[m,] <- new.alpha[[m]]
  }

  if(method_alpha_intial == "Ronning"){
    return(new_alpha_R)
  } else{
    return(new_alpha_W)
  }
}

# Update the matrix of log tau*pdf, a key interium product
# Return a matrix of nrow = number of cells (n), ncol = number of clusters (K)
# tauVec is a vector of length equal to number of cluster, indicating cluster proportions
# data1/data2, alphamtx1/alphamtx2 are source specific
# bVec is a vector (of length J) of random effect, each element for a cell
#' @import Rfast
#' @import stats
logDMultDirMtx_R = function(dataMtx, alphaVec, bVec, J, P) {
  alphaMtx <- matrix(rep(alphaVec, J), nrow = P)
  alphaBmtx <- alphaMtx * t(matrix(rep(bVec, P), nrow = J))

  part1 <- colSums(Lgamma(dataMtx +  alphaBmtx) - Lgamma(alphaBmtx))
  part2 <- Lgamma(colSums(alphaBmtx)) - Lgamma(colSums(dataMtx) + colSums(alphaBmtx))

  return(part1 + part2)
}

#' @import stats
updateLogTauF_R = function(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec, J, K) {
  P1 <- nrow(data1)
  P2 <- nrow(data2)
  taufMtxLog = matrix(nrow = J, ncol = K)

  for (k in 1:K) {
    alphaVec1 <- alphaMtx1[k,]
    alphaVec2 <- alphaMtx2[k,]

    logD1 <- logDMultDirMtx_R(data1, alphaVec1, bVec, J, P1)
    logD2 <- logDMultDirMtx_R(data2, alphaVec2, bVec, J, P2)

    taufMtxLog[,k] <- logD1 + logD2 + log(tauVec[k])
  }

  return(taufMtxLog)
}

#' @import stats
updateLogTauFCompare_R = function(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec, bVecNew, J, K) {
  P1 <- nrow(data1)
  P2 <- nrow(data2)
  taufMtxLog <- matrix(nrow = J, ncol = K)
  taufMtxLog2 <- taufMtxLog

  for (k in 1:K) {
    alphaVec1 <- alphaMtx1[k,]
    alphaVec2 <- alphaMtx2[k,]

    logD1 <- logDMultDirMtx_R(data1, alphaVec1, bVec, J, P1)
    logD2 <- logDMultDirMtx_R(data2, alphaVec2, bVec, J, P2)
    taufMtxLog[,k] <- logD1 + logD2 + log(tauVec[k])

    logD1New <- logDMultDirMtx_R(data1, alphaVec1, bVecNew, J, P1)
    logD2New <- logDMultDirMtx_R(data2, alphaVec2, bVecNew, J, P2)
    taufMtxLog2[,k] <- logD1New + logD2New + log(tauVec[k])
  }

  return(list(taufOld = taufMtxLog, taufNew = taufMtxLog2))
}

# Update the matrix of cluster-specific posteriror probabilities, a key product
# Return a matrix of nrow = number of cells (J), ncol = number of clusters (K)
#' @import stats
updateZmtx_R = function(taufMtxLog, J, K) {
  zMtx <- matrix(nrow = J, ncol = K)
  for (k in 1:K) {
    taufVec <- taufMtxLog[,k]
    zMtx[,k] <- 1/rowSums(exp(taufMtxLog-taufVec))
  }
  return(zMtx)
}

#' @import stats
updateZmtxCompare_R = function(taufMtxLogOld, taufMtxLogNew, J, K) {
  zMtx <- matrix(nrow = J, ncol = K)
  zMtx2 <- zMtx
  for (k in 1:K) {
    taufVec <- taufMtxLogOld[,k]
    zMtx[,k] <- 1/rowSums(exp(taufMtxLogOld-taufVec))

    taufVec2 <- taufMtxLogNew[,k]
    zMtx2[,k] <- 1/rowSums(exp(taufMtxLogNew-taufVec2))
  }
  return(list(zMtxOld = zMtx, zMtxNew = zMtx2))
}

# Compute complete log-likelihood
# sigmaB is the estiamted sd of random effect (a positive number)
#' @import stats
getLogL_R = function(zMtx, taufMtxLog, bVec, sigmaB) {
  logL1 <- sum(zMtx * taufMtxLog)
  logL2 <- sum(-log(bVec)-(log(bVec))^2/(2*sigmaB^2))
  logL3 <- nrow(zMtx)*(-log(sigmaB))
  return(logL1+logL2+logL3)
}

# Updated B vector for random effect (should be a vector not a matrix, name is misleading)
# z is a vector (of length eqault to number of cells) of cell-specific labels (drawn from posterior prob)
# sd_b reflects the jump distance in random walk
# sigmaB is the previous estimated sd of random effect
#' @import stats
#' @import doParallel
#' @import foreach
updateBMtxOne_R = function(data1, data2, z, tauVec, alphaMtx1, alphaMtx2, sd_b, bVec, sigmaB, J, K, nCores, distList) {
  exportFuns <- c("updateLogTauFCompare_R", "updateZmtxCompare_R", "getLogL_R", "logDMultDirMtx_R")

  out = foreach(i = 1:nCores, .combine = 'c', .export = exportFuns, .packages = "Rfast") %dopar% {
    cellVec <- distList[[i]]
    for (j in cellVec) {
      bVecNew <- bVec
      b_old <- bVec[j]
      b_new <- rnorm(1, b_old, sd=sd_b)
      if (b_new<0)   {b_new<- -b_new}
      bVecNew[j] <- b_new

      taufDual <- updateLogTauFCompare_R(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec, bVecNew, J, K)
      zMtxDual <- updateZmtxCompare_R(taufDual$taufOld, taufDual$taufNew, J, K)

      likeOld <- getLogL_R(zMtxDual$zMtxOld, taufDual$taufOld, bVec, sigmaB)
      likeNew <- getLogL_R(zMtxDual$zMtxNew, taufDual$taufNew, bVecNew, sigmaB)

      acc_p <- min(0, likeNew-likeOld)
      u <- runif(1,0,1)
      if(log(u) < acc_p){bVec[j] <- b_new }
    }
    bVec[cellVec]
  }

  return(out)
}

# Try to avoid zero counts for a cluster
# z is the cell label, zMtx is the matrix of cluster-specific posterior prob for each cell, K is number of cluster
# Always force a cluster to have at least 1 cell
#' @import stats
avoidZero = function(z, zMtx, K) {
  indexVec <- rep(NA, K)

  if (K == 1) {
    return(z)
  } else {
    indexVec[1] <- which.max(zMtx[,1])
    z[indexVec[1]] <- 1
    for (k in 2:K) {
      max0 <- max(zMtx[-indexVec[1:(k-1)],k])
      indexVec[k] <- setdiff(which(max0 == zMtx[,k]), indexVec[1:(k-1)])[1]
      z[indexVec[k]] <- k
    }
    return(z)
  }
}

# Final function to be called to return a list of objects we need (parallel version)
#' @import stats
#' @import parallel
#' @import doParallel
startMCMC_para = function(data1, data2, sd_alpha, sd_b, sigmaB, K, nMCMC, nCores) {
  cl <- makeCluster(nCores)
  registerDoParallel(cl)

  J <- ncol(data1)
  P1 <- nrow(data1)
  P2 <- nrow(data2)

  # Evenly distribute workers
  cellIndex <- 1:J
  coreIndex <- as.numeric(cut(cellIndex, nCores))
  distList <- list()
  for (l in 1:nCores) {
    distList[[l]] <- cellIndex[coreIndex == l]
  }

  ### Initial values set up
  clusters_initial0 <- kmeans(t(rbind(as.matrix(log2(data1+1)), as.matrix(log2(data2+1)))),K)$cluster

  clusters_initial1 <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  if(length(as.numeric(table(clusters_initial1))) < K | min(as.numeric(table(clusters_initial1))) == 1) {
    clusters_initial1 <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  }
  alpha_initial1 <-  EM_initial_alpha(data1,clusters_initial1,"Ronning")

  clusters_initial2 <- kmeans(t(as.matrix(log2(data2+1))),K)$cluster
  if(length(as.numeric(table(clusters_initial2))) < K | min(as.numeric(table(clusters_initial2))) == 1) {
    clusters_initial2 <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  }
  alpha_initial2 <-  EM_initial_alpha(data2,clusters_initial2,"Ronning")

  #b_initial = exp(rnorm(J, 0, sigmaB))
  b_initial <- rep(1, J)
  # tau vector
  tauVec0 <- as.numeric(table(clusters_initial0) / J)
  # zMtx
  taufMtxLog <- updateLogTauF_R(data1, data2, tauVec0, alpha_initial1, alpha_initial2, b_initial, J, K)
  zMtx <- updateZmtx_R(taufMtxLog, J, K)

  z <- rep(NA, J)
  for (i in 1:J) {
    z[i] <- sample(1:K, 1, prob = zMtx[i,])
  }
  z <- avoidZero(z, zMtx, K)
  # Initial log likelihood
  logL0 <- getLogL_R(zMtx, taufMtxLog, b_initial, sigmaB)
  # Updated njVec
  njVec <- updateNjVec(zMtx)
  # Updated tauVec
  tauVec <- updateTauVec(njVec, J)

  ### MCMC ###
  alphaMtx1 <- alpha_initial1
  alphaMtx2 <- alpha_initial2
  bVec <- b_initial

  alphaMtxMCMC1 <- matrix(nrow = nMCMC, ncol = K*P1)
  alphaMtxMCMC2 <- matrix(nrow = nMCMC, ncol = K*P2)
  bMtxMCMC <- matrix(nrow = nMCMC, ncol = J)
  sigmaBMCMC <- rep(NA, nMCMC)
  logLMtx <- rep(NA, nMCMC)
  tauOut <- matrix(nrow = nMCMC, ncol = K)
  zMtxOut <- matrix(nrow = nMCMC, ncol = J*K)
  for (i in 1:nMCMC) {
    if (i %% 50 == 1) {
      cat(paste0("MCMC: ", i, " starts..."), "\n")
    }

    # Update alpha1
    alphaMtx_temp1 <- as.numeric(updateAlphaMtxOne(data1, z, alphaMtx1, sd_alpha, bVec, K, P1))
    alphaMtxMCMC1[i,] <- alphaMtx_temp1
    alphaMtx1 <- matrix(alphaMtx_temp1, nrow = K, ncol = P1)
    # Update alpha2
    alphaMtx_temp2 <- as.numeric(updateAlphaMtxOne(data2, z, alphaMtx2, sd_alpha, bVec, K, P2))
    alphaMtxMCMC2[i,] <- alphaMtx_temp2
    alphaMtx2 <- matrix(alphaMtx_temp2, nrow = K, ncol = P2)
    # tau vector
    tauVec <- as.numeric(table(z) / J)
    tauOut[i,] <- tauVec
    # zMtx
    taufMtxLog <- updateLogTauF_R(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec, J, K)
    zMtx <- updateZmtx_R(taufMtxLog, J, K)
    zMtxOut[i,] <- as.numeric(zMtx)
    for (j in 1:J) {
      z[j] <- sample(1:K, 1, prob = zMtx[j,])
    }
    # Add a procedure to make sure no cluster has 0 cell
    z <- avoidZero(z, zMtx, K)
    # Update b
    #start.time = proc.time()
    bVec <- updateBMtxOne_R(data1, data2, z, tauVec, alphaMtx1, alphaMtx2, sd_b, bVec, sigmaB, J, K, nCores, distList)
    bMtxMCMC[i,] <- bVec
    #t0 = proc.time() - start.time
    #cat("Update bVec: ", as.numeric(t0[3]), "\n")

    # update sigma_b
    sigmaB<-sd(log(bVec))
    sigmaBMCMC[i]<-sigmaB
    # log-likelihood
    logL1 <- getLogL_R(zMtx, taufMtxLog, bVec, sigmaB)
    logLMtx[i] <- logL1
  }
  return(list(tau = tauOut, # proportion of each cluster
              alphaOut1 = alphaMtxMCMC1,
              alphaOut2 = alphaMtxMCMC2,
              zOut = zMtxOut, # updated p for cell i belongs to cluster j
              bOut = bMtxMCMC,
              sigB = sigmaBMCMC,
              logLik = logLMtx)) # final updated log-likelihood
}

# Final function to be called to return a list of objects we need (single core version)
#' @import stats
startMCMC_nonPara = function(data1, data2, sd_alpha, sd_b, sigmaB, K, nMCMC) {
  n <- ncol(data1)
  P1 <- nrow(data1)
  P2 <- nrow(data2)

  ### Initial values set up
  clusters_initial0 <- kmeans(t(rbind(as.matrix(log2(data1+1)), as.matrix(log2(data2+1)))),K)$cluster

  clusters_initial1 <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  alpha_initial1 <-  EM_initial_alpha(data1,clusters_initial1,"Ronning")

  clusters_initial2 <- kmeans(t(as.matrix(log2(data2+1))),K)$cluster
  alpha_initial2 <-  EM_initial_alpha(data2,clusters_initial2,"Ronning")

  #b_initial = exp(rnorm(n, 0, sigmaB))
  b_initial <- rep(1, n)
  # tau vector
  tauVec0 <- as.numeric(table(clusters_initial0) / n)
  # zMtx
  taufMtxLog <- updateLogTauF(data1, data2, tauVec0, alpha_initial1, alpha_initial2, b_initial, n, K)
  zMtx <- updateZmtx(taufMtxLog, n, K)

  z <- rep(NA, n)
  for (i in 1:n) {
    z[i] <- sample(1:K, 1, prob = zMtx[i,])
  }

  z <- avoidZero(z, zMtx, K)

  # Initial log likelihood
  logL0 <- getLogL(zMtx, taufMtxLog, b_initial, sigmaB)
  # Updated njVec
  njVec <- updateNjVec(zMtx)
  # Updated tauVec
  tauVec <- updateTauVec(njVec, n)

  ### MCMC ###
  alphaMtx1 <- alpha_initial1
  alphaMtx2 <- alpha_initial2
  bVec <- b_initial

  alphaMtxMCMC1 <- matrix(nrow = nMCMC, ncol = K*P1)
  alphaMtxMCMC2 <- matrix(nrow = nMCMC, ncol = K*P2)
  bMtxMCMC <- matrix(nrow = nMCMC, ncol = n)
  sigmaBMCMC <- rep(NA, nMCMC)
  logLMtx <- rep(NA, nMCMC)
  tauOut <- matrix(nrow = nMCMC, ncol = K)
  zMtxOut <- matrix(nrow = nMCMC, ncol = n*K)

  for (i in 1:nMCMC) {
    if (i %% 50 == 1) {
      cat(paste0("MCMC: ", i, " starts..."), "\n")
    }
    # Update alpha1
    alphaMtx_temp1 <- as.numeric(updateAlphaMtxOne(data1, z, alphaMtx1, sd_alpha, bVec, K, P1))
    alphaMtxMCMC1[i,] <- alphaMtx_temp1
    alphaMtx1 <- matrix(alphaMtx_temp1, nrow = K, ncol = P1)

    # Update alpha2
    alphaMtx_temp2 <- as.numeric(updateAlphaMtxOne(data2, z, alphaMtx2, sd_alpha, bVec, K, P2))
    alphaMtxMCMC2[i,] <- alphaMtx_temp2
    alphaMtx2 <- matrix(alphaMtx_temp2, nrow = K, ncol = P2)

    # tau vector
    tauVec <- as.numeric(table(z) / n)
    tauOut[i,] <- tauVec

    # zMtx
    taufMtxLog <- updateLogTauF(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec, n, K)
    zMtx <- updateZmtx(taufMtxLog, n, K)
    zMtxOut[i,] <- as.numeric(zMtx)
    for (j in 1:n) {
      z[j] <- sample(1:K, 1, prob = zMtx[j,])
    }

    # Add a procedure to make sure no cluster has 0 cell
    z <- avoidZero(z, zMtx, K)

    # Update b
    bVec <- updateBMtxOne(data1, data2, z, tauVec, alphaMtx1, alphaMtx2, sd_b, bVec, sigmaB, n, K)
    bMtxMCMC[i,] <- bVec

    # update sigma_b
    sigmaB<-sd(log(bVec))
    sigmaBMCMC[i]<-sigmaB

    # log-likelihood
    logL1 <- getLogL(zMtx, taufMtxLog, bVec, sigmaB)
    logLMtx[i] <- logL1

  }

  return(list(tau = tauOut, # proportion of each cluster
              alphaOut1 = alphaMtxMCMC1,
              alphaOut2 = alphaMtxMCMC2,
              zOut = zMtxOut, # updated p for cell i belongs to cluster j
              bOut = bMtxMCMC,
              sigB = sigmaBMCMC,
              logLik = logLMtx)) # final updated log-likelihood

}

#' BREMSC function
#' @name BREMSC
#' @aliases BREMSC
#' @param dataProtein a D*C matrix with D proteins and C cells
#' @param dataRNA a G*C matrix with G genes and C cells
#' @param K number of desired clusters
#' @param nCores number of CPUs to use
#' @param nMCMC number of MCMCs
#' @param useGene number of genes to keep
#' @param diagnostic whether store a dignostic plot
#' @param sd_alpha random walk depth for alpha
#' @param sd_b random walk depth for b
#' @param sigmaB initial estimate of SD of random effects
#' @return BREMSC returns a list object containing:
#' @return \itemize{
#'   \item clusterID: estimated label for each cell
#'   \item posteriorProb: a C*K matrix with probability that each cell belongs to each cluster
#'   \item sdRF: estimated SD of random effects
#'   \item logLik: estimated log likelihood
#'   \item vecB: estimated random effects for each cell
#'   \item alphaMtxProtein: a K*D matrix of alpha estimates for protein source
#'   \item alphaMtxRNA: a K*G matrix of alpha estimates for RNA source
#' }
#' @author Xinjun Wang <xiw119@pitt.edu>, Zhe Sun <zhs31@pitt.edu>, Wei Chen <wei.chen@chp.edu>.
#' @references Xinjun Wang, Zhe Sun, Yanfu Zhang, Heng Huang, Kong Chen, Ying Ding, Wei Chen. BREM-SC: A Bayesian Random Effects Mixture Model for Joint Clustering Single Cell Multi-omics Data. Submitted 2019.
#' @examples
#' # Load the example data data_DIMMSC
#' data("dataADT")
#' data("dataRNA")
#' # Test run of BREMSC: use small number of MCMC to save time
#' testRun = BREMSC(dataADT, dataRNA, K=4, nCores=5, nMCMC=20)
#' # End
#' @import stats
#' @export
BREMSC = function(dataProtein, dataRNA, K, nCores = 1, nMCMC = 500, useGene = 100, diagnostic = T,
                  sd_alpha = 1, sd_b = 0.5, sigmaB = 0.8) {
  # Format input data
  cat(paste0("Start loading data..."), "\n")
  data1 <- data.matrix(dataProtein)
  data2 <- data.matrix(dataRNA)
  n = ncol(data1) # compute number of cells

  if (ncol(data1) != ncol(data2)) {
    stop("Dimension of two data sources don't match. Check if each source of data is feature by cell,
         and the cells should match in two sources.")
  }

  # Keep only top genes (G = useGene)
  cat(paste0("Start selecting top genes..."), "\n")
  sd<-apply(data2,1,sd)
  or<-order(sd)
  or<-or[dim(data2)[1]:1]
  or<-or[1:useGene]
  list<-sort(or)
  data2<-data2[list,]
  # Decide to use single core coding or multi-core coding
  cat(paste0("Start running BREMSC..."), "\n")
  if (nCores > 1) {
    BREM_rslt <- startMCMC_para(data1, data2, sd_alpha, sd_b, sigmaB, K, nMCMC, nCores)
  } else {
    BREM_rslt <- startMCMC_nonPara(data1, data2, sd_alpha, sd_b, sigmaB, K, nMCMC)
  }
  # Save dignostic plot
  if (diagnostic == T) {
    cat(paste0("Start collecting results..."), "\n")
    jpeg('Diagnostic_LogLik.jpg')
    plot(BREM_rslt$logLik, type = "l", xlab = "MCMC", ylab = "logLik")
    dev.off()
  }
  # Collect posterior probability of cluster for cell (N*K)
  postPMtx <- matrix(BREM_rslt$zOut[nMCMC,], ncol = K)
  # Collect final clustering information
  label_out <- apply(postPMtx,1,which.max)
  # Collect final estimate of sigB
  estSigmaB <- BREM_rslt$sigB[nMCMC]
  # Collect final estimate of log likelihood
  logLikLast <- BREM_rslt$logLik[nMCMC]
  # Collect final estimate of log likelihood
  bVec <- BREM_rslt$bOut[nMCMC,]
  # Collect final estimate of log likelihood
  alphaMtx1 <- matrix(BREM_rslt$alphaOut1[nMCMC,], nrow = K)
  # Collect final estimate of log likelihood
  alphaMtx2 <- matrix(BREM_rslt$alphaOut2[nMCMC,], nrow = K)

  cat(paste0("All done!"), "\n")
  return(list(clusterID = label_out,
              posteriorProb = postPMtx,
              sdRF = estSigmaB,
              logLik = logLikLast,
              vecB = bVec,
              alphaMtxProtein = alphaMtx1,
              alphaMtxRNA = alphaMtx2))
}

