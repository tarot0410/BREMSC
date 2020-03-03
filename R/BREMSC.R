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

# Initialize alpha
#' @import Rfast
#' @import stats
# Compute log density of multinomial dirichlet distribution for cell j
logCondAlphaGibbs = function(yK_vec, sumYK_vec, alpha, alphaK_vec, bK_vec) {
  alphaB_vec = alpha*bK_vec
  part1 = Lgamma(yK_vec + alphaB_vec) - Lgamma(alphaB_vec)
  sumAlpha = sum(alphaK_vec)
  sumAlphaB_vec = sumAlpha*bK_vec
  part2 = Lgamma(sumAlphaB_vec) - Lgamma(sumYK_vec + sumAlphaB_vec)

  log_comb = part1+part2
  return(sum(log_comb))
}

### Updated alphaMtx ###
# yK_mtx is input y matrix (gene by cell, all gene with cells belongs to cluster K)
# bK_vec is input random effects of cells belong to cluster K
# alphaK_vec is input alpha vector corresponding to cluster K
# returns the conditional log likelihood for all alpha in cluster K
#' @import stats
updateAlphaMtxOneGibbs = function(data0, z, alphaMtx, sd_alpha, bVec, K, P) {
  sumY_vec = colSums(data0)
  for(k in 1:K){
    type = which(z==k)
    yK_mtx = data0[, type]
    sumYK_vec = sumY_vec[type]
    bK_vec = bVec[type]

    alphaK_poolOld<-alphaMtx[k,]
    alphaK_vecNew = rnorm(P, alphaK_poolOld, sd_alpha)
    alphaK_vecNew[alphaK_vecNew>500] = 1000 - alphaK_vecNew[alphaK_vecNew>500]
    alphaK_vecNew[alphaK_vecNew<0] = -alphaK_vecNew[alphaK_vecNew<0]

    for(p in 1:P){
      yK_vec = yK_mtx[p,]
      alphaOld = alphaK_poolOld[p]
      alphaNew = alphaK_vecNew[p]
      alphaK_poolNew <- alphaK_poolOld
      alphaK_poolNew[p]<- alphaNew

      sumLikOld = logCondAlphaGibbs(yK_vec, sumYK_vec, alphaOld, alphaK_poolOld, bK_vec)
      sumLikNew = logCondAlphaGibbs(yK_vec, sumYK_vec, alphaNew, alphaK_poolNew, bK_vec)

      acc_p <- min(0, sumLikNew - sumLikOld)

      u <- runif(1,0,1)
      if(log(u) < acc_p) {
        alphaK_poolOld[p] <- alphaNew
      }
    }
    alphaMtx[k,] = alphaK_poolOld
  }
  return(alphaMtx)
}

#' @import Rfast
#' @import stats
logCondBCore = function(RNA_K_mtx, ADT_K_mtx, alphaRNA_K_vec, alphaADT_K_vec, bK_vec, P_RNA, P_ADT) {
  tempB_RNA_mtx = matrix(rep(bK_vec, P_RNA), nrow=P_RNA, byrow=T)
  tempB_ADT_mtx = matrix(rep(bK_vec, P_ADT), nrow=P_ADT, byrow=T)

  alphaB_RNA_mtx = tempB_RNA_mtx*alphaRNA_K_vec
  alphaB_ADT_mtx = tempB_ADT_mtx*alphaADT_K_vec

  part1RNA_mtx = Lgamma(RNA_K_mtx + alphaB_RNA_mtx) - Lgamma(alphaB_RNA_mtx)
  part1ADT_mtx = Lgamma(ADT_K_mtx + alphaB_ADT_mtx) - Lgamma(alphaB_ADT_mtx)

  sumAlphaB_RNA_vec = bK_vec*sum(alphaRNA_K_vec)
  part2RNA_vec = Lgamma(sumAlphaB_RNA_vec) - Lgamma(sumAlphaB_RNA_vec + colSums(RNA_K_mtx))
  sumAlphaB_ADT_vec = bK_vec*sum(alphaADT_K_vec)
  part2ADT_vec = Lgamma(sumAlphaB_ADT_vec) - Lgamma(sumAlphaB_ADT_vec + colSums(ADT_K_mtx))

  sum_vec = colSums(part1RNA_mtx) + part2RNA_vec + colSums(part1ADT_mtx) + part2ADT_vec
  return(sum_vec)
}

# zMtx is N*K
#' @import stats
updateZ_mtx = function(data1, data2, alpha1, alpha2, bVec, tau, K, P1, P2, N) {
  log_tempMtx = matrix(nrow = N, ncol = K)
  outMtx = log_tempMtx
  for (k in 1:K) {
    log_tempMtx[,k] = logCondBCore(data1, data2, alpha1[k, ], alpha2[k, ], bVec, P1, P2)
  }
  for (k in 1:K) {
    tempMtx2 = log_tempMtx - log_tempMtx[,k]
    tempMtx3 = exp(tempMtx2) %*% matrix(tau, ncol=1)
    outMtx[,k] = tau[k] / tempMtx3
  }
  return(outMtx)
}

# get loglik
#' @import stats
getLogL = function(dataRNA, dataADT, z, alphaRNA, alphaADT, b_vec, sigmaB, K, P_RNA, P_ADT, N) {
  coreOldVec = rep(NA, N)
  for (k in 1:K) {
    type <- which(z==k)
    num_cell<-length(type)

    RNA_K_mtx = dataRNA[, type]
    ADT_K_mtx = dataADT[, type]

    alphaRNA_K_vec = alphaRNA[k,]
    alphaADT_K_vec = alphaADT[k,]

    bK_vec = b_vec[type]
    coreOldVec[type] = logCondBCore(RNA_K_mtx, ADT_K_mtx, alphaRNA_K_vec, alphaADT_K_vec, bK_vec, P_RNA, P_ADT)
  }

  logLikOld_vec = coreOldVec - log(b_vec) - (log(b_vec))^2 / (2*sigmaB^2)
  return(sum(logLikOld_vec))
}

### Updated bMtx ###
#' @import stats
updateBMtxOne = function(dataRNA, dataADT, z, alphaRNA, alphaADT, b_vec, sd_b, sigmaB, K, P_RNA, P_ADT, N) {
  b_vecNew = rnorm(N, b_vec, sd_b)
  b_vecNew[b_vecNew<0] = -b_vecNew[b_vecNew<0]

  coreOldVec = rep(NA, N)
  coreNewVec = coreOldVec
  for (k in 1:K) {
    type <- which(z==k)
    num_cell<-length(type)

    RNA_K_mtx = dataRNA[, type]
    ADT_K_mtx = dataADT[, type]

    alphaRNA_K_vec = alphaRNA[k,]
    alphaADT_K_vec = alphaADT[k,]

    bK_vec = b_vec[type]
    bK_vecNew = b_vecNew[type]
    coreOldVec[type] = logCondBCore(RNA_K_mtx, ADT_K_mtx, alphaRNA_K_vec, alphaADT_K_vec, bK_vec, P_RNA, P_ADT)
    coreNewVec[type] = logCondBCore(RNA_K_mtx, ADT_K_mtx, alphaRNA_K_vec, alphaADT_K_vec, bK_vecNew, P_RNA, P_ADT)
  }

  logLikOld_vec = coreOldVec - log(b_vec) - (log(b_vec))^2 / (2*sigmaB^2)
  logLikNew_vec = coreNewVec - log(b_vecNew) - (log(b_vecNew))^2 / (2*sigmaB^2)

  acc_p_vec <- pmin(0, logLikNew_vec-logLikOld_vec)
  u_vec <- runif(N,0,1)

  updateInd = which(log(u_vec) < acc_p_vec)
  b_vec[updateInd] = b_vecNew[updateInd]

  return(b_vec)
}

# Try to avoid zero counts for a cluster
# z is the cell label, zMtx is the matrix of cluster-specific posterior prob for each cell, N*K
# Always force a cluster to have at least 2 cells
#' @import stats
avoidZero = function(z, K) {
  N = length(z)
  k_null = setdiff(1:K, unique(z))
  z_vec = as.numeric(names(table(z)))
  k_only1 = z_vec[as.numeric(table(z)) == 1]

  maxK = z_vec[which.max(as.numeric(table(z)))]
  randPool = (1:N)[z==maxK]
  for (k in k_null) {
    rand2 = sample(randPool,2)
    z[rand2] = k
    randPool = setdiff(randPool, rand2)
  }
  for (k in k_only1) {
    rand1 = sample(randPool,1)
    z[rand1] = k
    randPool = setdiff(randPool, rand1)
  }
  return(z)
}

# Final function to be called to return a list of objects we need (single core version)
#' @import stats
startMCMC_singleGibbs = function(data1, data2, sd_alpha, sd_b, sigmaB, K, nMCMC) {
  n = ncol(data1)
  P1 = nrow(data1)
  P2 = nrow(data2)

  ### Initial values set up
  clusters_initial0 <- kmeans(t(rbind(as.matrix(log2(data1+1)), as.matrix(log2(data2+1)))),K)$cluster

  clusters_initial1 <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  alpha_initial1 <-  EM_initial_alpha(data1,clusters_initial1,"Ronning")

  clusters_initial2 <- kmeans(t(as.matrix(log2(data2+1))),K)$cluster
  alpha_initial2 <-  EM_initial_alpha(data2,clusters_initial2,"Ronning")

  #b_initial = exp(rnorm(n, 0, sigmaB))
  b_initial = rep(1, n)

  # tau vector
  tauVec0 = as.numeric(table(clusters_initial0) / n)
  # zMtx
  zMtx = updateZ_mtx(data1, data2, alpha_initial1, alpha_initial2, b_initial, tauVec0, K, P1, P2, n)

  z = rep(NA, n)
  for (i in 1:n) {
    z[i] = sample(1:K, 1, prob = zMtx[i,])
  }

  # Add a procedure to make sure all clusters have >= 2 cells
  if (!(all(as.numeric(table(z))>=2) & length(unique(z)) == K)) {
    z = avoidZero(z, K)
  }

  # Initial log likelihood
  logL0 = getLogL(data1, data2, z, alpha_initial1, alpha_initial2, b_initial, sigmaB, K, P1, P2, n)

  ### MCMC ###
  alphaMtx1 = alpha_initial1
  alphaMtx2 = alpha_initial2
  bVec = b_initial

  sigmaBMCMC = rep(NA, nMCMC)
  logLMtx = rep(NA, nMCMC)

  for (i in 1:nMCMC) {
    if (i %% 100 == 1) {
      cat(paste0("MCMC: ", i, " starts..."), "\n")
    }
    # Update alpha1
    alphaMtx1 = updateAlphaMtxOneGibbs(data1, z, alphaMtx1, sd_alpha, bVec, K, P1)
    # Update alpha2
    alphaMtx2 = updateAlphaMtxOneGibbs(data2, z, alphaMtx2, sd_alpha, bVec, K, P2)
    # tau vector
    tauVec = as.numeric(table(z) / n)
    #tauOut[i,] = tauVec
    # zMtx
    zMtx = updateZ_mtx(data1, data2, alphaMtx1, alphaMtx2, bVec, tauVec, K, P1, P2, n)

    for (j in 1:n) {
      z[j] = sample(1:K, 1, prob = zMtx[j,])
    }

    # Add a procedure to make sure all clusters have >= 2 cells
    if (!(all(as.numeric(table(z))>=2) & length(unique(z)) == K)) {
      z = avoidZero(z, K)
    }

    # Update b
    bVec = updateBMtxOne(data1, data2, z, alphaMtx1, alphaMtx2, bVec, sd_b, sigmaB, K, P1, P2, n)

    # update sigma_b
    sigmaB<-sd(log(bVec))
    sigmaBMCMC[i]<-sigmaB

    # log-likelihood
    logL1 = getLogL(data1, data2, z, alphaMtx1, alphaMtx2, bVec, sigmaB, K, P1, P2, n)
    logLMtx[i] = logL1
  }

  return(list(tau = tauVec, # proportion of each cluster
              alphaOut1 = alphaMtx1,
              alphaOut2 = alphaMtx2,
              zOut = zMtx, # updated p for cell i belongs to cluster j
              bOut = bVec,
              sigB = sigmaBMCMC,
              logLik = logLMtx)) # final updated log-likelihood

}

#' @import stats
#' @import foreach
#' @import doParallel
#' @import parallel
multiChainGibbs = function(data1, data2, sd_alpha_range, sd_b_range, sigmaB, K, nMCMC, nChain) {
  sd_alpha_vec = seq(sd_alpha_range[1], sd_alpha_range[2], (sd_alpha_range[2] - sd_alpha_range[1]) / (nChain-1))
  sd_b_vec = seq(sd_b_range[1], sd_b_range[2], (sd_b_range[2] - sd_b_range[1]) / (nChain-1))

  cl <- makeCluster(nChain)
  registerDoParallel(cl)

  exportFuns <- c("EM_initial_alpha", "logCondAlphaGibbs",
                  "updateAlphaMtxOneGibbs", "logCondBCore", "updateZ_mtx", "getLogL", "updateBMtxOne",
                  "avoidZero", "startMCMC_singleGibbs")

  out = foreach(i = 1:nChain, .export = exportFuns, .packages = c("Rfast", "dirmult")) %dopar% {
    startMCMC_singleGibbs(data1, data2, sd_alpha_vec[i], sd_b_vec[i], sigmaB, K, nMCMC)
  }

  stopCluster(cl)

  logLikVec = rep(NA, nChain)
  for (j in 1:nChain) {
    logLikVec[j] = out[[j]]$logLik[nMCMC]
  }

  highest = which.max(logLikVec)

  return(list(tau = out[[highest]]$tau, # proportion of each cluster
              alphaOut1 = out[[highest]]$alphaOut1,
              alphaOut2 = out[[highest]]$alphaOut2,
              zOut = out[[highest]]$zOut, # updated p for cell i belongs to cluster j
              bOut = out[[highest]]$bOut,
              sigB = out[[highest]]$sigB,
              logLik = out[[highest]]$logLik)) # final updated log-likelihood
}


#' A novel joint clustering method for scRNA-seq and CITE-seq data
#'
#' BREMSC is developed to joint cluster scRNA-seq and CITE-seq data with the introduction of random effects to incorporate the two data sources
#' @name BREMSC
#' @param dataProtein a D*C matrix with D proteins and C cells
#' @param dataRNA a G*C matrix with G genes and C cells
#' @param K number of desired clusters
#' @param nChains number of MCMC chains to run in parallel, each using a CPU
#' @param nMCMC number of MCMCs
#' @param sd_alpha_range range of random walk depth for alpha
#' @param sd_b_range range of random walk depth for b
#' @param sigmaB initial estimate of SD of random effects
#' @return BREMSC returns a list object containing:
#' @return \itemize{
#'   \item clusterID: estimated label for each cell
#'   \item posteriorProb: a C*K matrix with probability that each cell belongs to each cluster
#'   \item sdRF_final: estimated SD of random effects
#'   \item logLik_final: estimated log likelihood
#'   \item vecB: estimated random effects for each cell
#'   \item alphaMtxProtein: a K*D matrix of alpha estimates for protein source
#'   \item alphaMtxRNA: a K*G matrix of alpha estimates for RNA source
#'   \item vecLogLik: estimated log likelihood for each MCMC
#'   \item vecSDRF: estimated SD of random effects for each MCMC
#' }
#' @author Xinjun Wang <xiw119@pitt.edu>, Zhe Sun <zhs31@pitt.edu>, Wei Chen <wei.chen@chp.edu>.
#' @references Xinjun Wang, Zhe Sun, Yanfu Zhang, Zhongli Xu, Heng Huang, Richard H Duerr, Kong Chen, Ying Ding, Wei Chen. BREM-SC: A Bayesian Random Effects Mixture Model for Joint Clustering Single Cell Multi-omics Data. Submitted 2019.
#' @examples
#' # Load the example data data_DIMMSC
#' data("dataADT")
#' data("dataRNA")
#' # Test run of BREMSC: use small number of MCMC to save time
#' testRun = BREMSC(dataADT, dataRNA, K=4)
#' # End
#' @import stats
#' @export
BREMSC = function(dataProtein, dataRNA, K, nChains = 3, nMCMC = 1000,
                  sd_alpha_range = c(0.5, 1.5), sd_b_range = c(0.2, 1), sigmaB = 0.8) {
  # Format input data
  cat(paste0("Start loading data..."), "\n")
  data1 <- data.matrix(dataProtein)
  data2 <- data.matrix(dataRNA)

  # Add pseudo gene names to RNA data for seurat to select top genes
  if (is.null(rownames(data2))) {
    rownames(data2) = paste0("v", 1:nrow(data2))
  }

  n = ncol(data1) # compute number of cells

  if (ncol(data1) != ncol(data2)) {
    stop("Dimension of two data sources don't match. Check if each source of data is feature by cell,
         and the cells should match in two sources.")
  }

  # Keep only top genes (G = useGene)
  if (nrow(dataRNA) > 1000) {
    cat(paste0("Consider to select a smaller number of RNA features if running too slow!"), "\n")
  }

  # Decide to use single core coding or multi-core coding
  cat(paste0("Start running BREMSC..."), "\n")
  if (nChains > 1) {
    BREM_rslt <- multiChainGibbs(data1, data2, sd_alpha_range, sd_b_range, sigmaB, K, nMCMC, nChains)
  } else {
    BREM_rslt <- startMCMC_singleGibbs(data1, data2, sum(sd_alpha_range)/2, sum(sd_b_range)/2, sigmaB, K, nMCMC)
  }

  # Collect posterior probability of cluster for cell (N*K)
  postPMtx <- BREM_rslt$zOut
  # Collect final clustering information
  label_out <- apply(postPMtx,1,which.max)
  # Collect final estimate of sigB
  estSigmaB <- BREM_rslt$sigB[nMCMC]
  # Collect final estimate of log likelihood
  logLikLast <- BREM_rslt$logLik[nMCMC]
  # Collect final estimate of log likelihood
  bVec <- BREM_rslt$bOut
  # Collect final estimate of log likelihood
  alphaMtx1 <- BREM_rslt$alphaOut1
  # Collect final estimate of log likelihood
  alphaMtx2 <- BREM_rslt$alphaOut2

  cat(paste0("All done!"), "\n")
  return(list(clusterID = label_out,
              posteriorProb = postPMtx,
              sdRF_final = estSigmaB,
              logLik_final = logLikLast,
              alphaMtxProtein = alphaMtx1,
              alphaMtxRNA = alphaMtx2,
              vecLogLik = BREM_rslt$logLik,
              vecSDRF = BREM_rslt$sigB,
              vecB = bVec))
}


# EM algorithm for DIMMSC_Joint
#' @import stats
#' @import Rfast
EM_multinomial = function(dataProtein, dataRNA, K, alpha_adt, alpha, maxiter, tol, lik.tol)
{
  pie <- rdirichlet(1,rep(1,K))
  J <- dim(dataRNA)[2]
  G <- dim(dataRNA)[1]
  n = nrow(dataProtein)
  p = ncol(dataProtein)


  differ=1
  iter=0
  loglik=1
  dif=100

  while ( (differ > tol | dif > lik.tol) &  iter < maxiter ) {
    ## E-step: compute omega:
    ## ADT
    num1_adt <- matrix(0,J,K)
    num2_adt <- matrix(,J,K)
    for (j in 1:J)
      for (k in 1:K){
        num1_adt[j,k] <- sum(Lgamma(dataProtein[,j]+alpha_adt[k,]) - Lgamma(alpha_adt[k,]))
        num2_adt[j,k] <- Lgamma(sum(alpha_adt[k,]))-Lgamma(sum(alpha_adt[k,])+sum(dataProtein[,j]))
      }

    delta_tmp_adt <- matrix(,J,K)
    for (j in 1:J){
      for (k in 1:K){
        delta_tmp_adt[j,k] <- num1_adt[j,k]+num2_adt[j,k]
      }
    }

    ## RNA
    num1 <- matrix(0,J,K)
    num2 <- matrix(,J,K)
    for (j in 1:J)
      for (k in 1:K){
        num1[j,k] <- sum(Lgamma(dataRNA[,j]+alpha[k,]) - Lgamma(alpha[k,]))
        num2[j,k] <- Lgamma(sum(alpha[k,]))-Lgamma(sum(alpha[k,])+sum(dataRNA[,j]))
      }

    delta_tmp <- matrix(,J,K)
    for (j in 1:J){
      for (k in 1:K){
        delta_tmp[j,k] <- num1[j,k]+num2[j,k]+log(pie[k])
      }
    }

    ########
    new_delta_tmp<-matrix(,J,K)
    for (j in 1:J){
      for (k in 1:K){
        new_delta_tmp[j,k]<-delta_tmp[j,k] + delta_tmp_adt[j,k]
      }
    }
    delta <- matrix(,J,K)
    for (j in 1:J){
      for (k in 1:K) {
        M <- c(1:K)[-k]
        sum <- 0
        for (m in 1:length(M)){
          sum_tmp <- exp(new_delta_tmp[j,M[m]]-new_delta_tmp[j,k])
          sum <- sum_tmp+sum
        }
        delta[j,k] <- 1/(1+sum)
      }
    }

    ## M-step: update pie and alpha
    # cat("start M-step\n")
    pie.new <- matrix(,1,K)
    for (k in 1:K){
      pie.new[1,k] <- sum(delta[,k])/J
    }

    alpha_new <- matrix(,K,G)
    for (k in 1:K){
      den <- 0
      for (j in 1:J){
        den_tmp <- delta[j,k]*(sum(dataRNA[,j])/(sum(dataRNA[,j])-1+sum(alpha[k,])))
        den <- den+den_tmp
      }
      for ( i in 1:G){
        num <- 0
        num <- sum(delta[,k]*((dataRNA[i,]+0.000001)/(dataRNA[i,]+0.000001-1+alpha[k,i])))
        alpha_new[k,i] <- alpha[k,i]*num/den
      }
    }
    ####
    alpha_adt_new <- matrix(,K,n)
    for (k in 1:K){
      den <- 0
      for (j in 1:J){
        den_tmp <- delta[j,k]*(sum(dataProtein[,j])/(sum(dataProtein[,j])-1+sum(alpha_adt[k,])))
        den <- den+den_tmp
      }
      for ( i in 1:n){
        num <- 0
        num <- sum(delta[,k]*((dataProtein[i,]+0.000001)/(dataProtein[i,]+0.000001-1+alpha_adt[k,i])))
        alpha_adt_new[k,i] <- alpha_adt[k,i]*num/den
      }
    }
    ####
    mem <- matrix(,J,1)
    for (i in 1:J){
      mem[i,1] <- which(delta[i,]==max(delta[i,]))
    }

    sort <- rank(pie)
    res <- mem
    for(k in 1:K){
      res[which(mem==k)] <- sort[k]
    }
    mem <- res

    num <- num1+num2+num1_adt+num2_adt
    lik <- matrix(,J,1)
    for(j in 1:J){
      lik[j,1] <- num[j,mem[j,]]
    }
    new.loglik <- sum(lik)

    ## calculate diff to check convergence
    dif <- abs((new.loglik-loglik)/loglik*100)

    sumd <- 0
    for (k in 1:K){
      diff <- (pie.new[k]-pie[k])^2
      sumd=diff+sumd
    }

    differ=sqrt(abs(sumd))
    #differ=1;
    pie=pie.new;
    alpha=alpha_new;
    alpha_adt=alpha_adt_new;
    loglik <- new.loglik

    for (k in 1:k){
      alpha[k,which(alpha[k,]<=0)] <- 0.000001
    }
    for (k in 1:k){
      alpha_adt[k,which(alpha_adt[k,]<=0)] <- 0.000001
    }

    iter=iter+1;
    #cat("Iter", iter, ", differ=", differ, "\n")
  }

  mem <- matrix(,J,1)
  for (i in 1:J){
    if (length(which(delta[i,]==max(delta[i,])) )>1)
    {mem[i,1] <- which(delta[i,]==max(delta[i,]))[1]  }
    if (length(which(delta[i,]==max(delta[i,])) )==1)
    {mem[i,1] <- which(delta[i,]==max(delta[i,]))   }
  }

  return(list(pie=pie,delta=delta,alpha=alpha, alpha_adt=alpha_adt, mem=mem,loglik=loglik))
}


#' A faster version of algorithm for joint clustering scRNA-seq and CITE-seq data
#'
#' jointDIMMSC is developed as a direct extension of DIMMSC without the introduction of random effects to incorporate the two data sources
#' @rdname BREMSC
#' @param dataProtein a D*C matrix with D proteins and C cells
#' @param dataRNA a G*C matrix with G genes and C cells
#' @param K number of desired clusters
#' @param maxiter maximum number of iterations, with default value 100
#' @param tol a convergence tolerance for the difference of vector pie between iterations, with default value 1e-4
#' @param lik.tol a convergence tolerance for the difference of log-likelihoods between iterations, with default value 1e-2
#' @return jointDIMMSC returns a list object containing:
#' @return \itemize{
#'   \item clusterID: estimated label for each cell
#'   \item posteriorProb: a C*K matrix with probability that each cell belongs to each cluster
#'   \item logLik: estimated log likelihood
#'   \item alphaMtxProtein: a K*D matrix of alpha estimates for protein source
#'   \item alphaMtxRNA: a K*G matrix of alpha estimates for RNA source
#' }
#' # Load the example data data_DIMMSC
#' data("dataADT")
#' data("dataRNA")
#' # Test run of BREMSC: use small number of MCMC to save time
#' testRun = jointDIMMSC(dataADT, dataRNA, K=4)
#' # End
#' @import stats
#' @export
jointDIMMSC = function(dataProtein, dataRNA, K, maxiter = 100, tol = 1e-4, lik.tol = 1e-2) {
  # Format input data
  cat(paste0("Start loading data..."), "\n")
  data1 <- data.matrix(dataProtein)
  data2 <- data.matrix(dataRNA)
  n = ncol(data1) # compute number of cells

  if (ncol(data1) != ncol(data2)) {
    stop("Dimension of two data sources don't match. Check if each source of data is feature by cell,
         and the cells should match in two sources.")
  }

  # Initialize jointDIMMSC
  cat(paste0("Start initializing jointDIMMSC..."), "\n")
  # adt_data
  clusters_initial_adt <- kmeans(t(as.matrix(log2(data1+1))),K)$cluster
  alpha_adt <- EM_initial_alpha(data1,clusters_initial_adt,"Ronning")
  # rna_data
  clusters_initial <- kmeans(t(as.matrix(log2(data2+1))),K)$cluster
  alpha <- EM_initial_alpha(data2,clusters_initial,"Ronning")

  # Start running jointDIMMSC
  cat(paste0("Start running jointDIMMSC..."), "\n")
  jointDIMMSC_rslt <- EM_multinomial(data1, data2, K, alpha_adt, alpha, maxiter, tol, lik.tol)

  cat(paste0("All done!"), "\n")
  return(list(clusterID = jointDIMMSC_rslt$mem,
              posteriorProb = jointDIMMSC_rslt$delta,
              logLik = jointDIMMSC_rslt$loglik,
              alphaMtxProtein = jointDIMMSC_rslt$alpha_adt,
              alphaMtxRNA = jointDIMMSC_rslt$alpha))
}




