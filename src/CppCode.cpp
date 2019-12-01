#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


// Estimating clusters for droplet-based single cell data using EM algorithm
//
// @name BREMSC
// @aliases BREMSC
// @param dataProtein a D*C matrix with D proteins and C cells
// @param dataRNA a G*C matrix with G genes and C cells
// @param K number of desired clusters
// @param nCores number of CPUs to use
// @param nMCMC number of MCMCs
// @param useGene number of genes to keep
// @param diagnostic whether store a dignostic plot
// @param sd_alpha random walk depth for alpha
// @param sd_b random walk depth for b
// @param sigmaB initial estimate of SD of random effects
// @return BREMSC returns a list object containing:
// @return \itemize{
//   \item clusterID: estimated label for each cell
//   \item posteriorProb: a C*K matrix with probability that each cell belongs to each cluster
//   \item sdRF: estimated SD of random effects
//   \item logLik: estimated log likelihood
//   \item vecB: estimated random effects for each cell
//   \item alphaMtxProtein: a K*D matrix of alpha estimates for protein source
//   \item alphaMtxRNA: a K*G matrix of alpha estimates for RNA source
// Author: Xinjun Wang and Zhe Sun

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
arma::rowvec Arma_colSums(arma::mat& x) {
  return arma::sum(x, 0);
}

arma::colvec Arma_rowSums(arma::mat& x) {
  return arma::sum(x, 1);
}

//So...alright no comment is allowed between export line and the function

double logDMultDirFull(arma::vec& yVec, arma::vec& alphaVec, double bj) {
   double log_part2 = arma::sum(lgamma(yVec+alphaVec*bj) - lgamma(alphaVec*bj));
   double log_part3 = lgamma(arma::sum(alphaVec*bj)) - lgamma(arma::sum(yVec) + arma::sum(alphaVec*bj));
   double out = log_part2 + log_part3;
   return out;
}

arma::mat genRepMtx(arma::colvec& x, int P, int J) {
  arma::mat out = arma::zeros<arma::mat>(P,J);
  out.each_col() = x;
  return out;
}

arma::mat genRepMtxTran(arma::colvec& x, int P, int J) {
  arma::mat out = arma::zeros<arma::mat>(P,J);
  arma::rowvec xTran = x.t();
  out.each_row() = xTran;
  return out;
}

arma::rowvec logDMultDirMtx(arma::mat& dataMtx, arma::colvec& alphaVec, arma::colvec& bVec, int J, int P) {
  arma::mat alphaMtx = genRepMtx(alphaVec, P, J);
  arma::mat bMtx = genRepMtxTran(bVec, P, J);
  arma::mat alphaBmtx = alphaMtx % bMtx;

  arma::mat mtx1 = lgamma(dataMtx + alphaBmtx) - lgamma(alphaBmtx);
  arma::rowvec vec1 = Arma_colSums(mtx1);
  arma::rowvec vec2 = lgamma(Arma_colSums(alphaBmtx)) - lgamma(Arma_colSums(dataMtx) + Arma_colSums(alphaBmtx));

  return vec1+vec2;
}

// [[Rcpp::export]]
arma::mat updateLogTauF(arma::mat& data1, arma::mat& data2, arma::vec& tauVec, arma::mat& alphaMtx1, arma::mat& alphaMtx2,  arma::vec& bVec, int J, int K) {
  int P1 = data1.n_rows;
  int P2 = data2.n_rows;

  arma::mat taufMtxLog;
  taufMtxLog.zeros(J, K);

  for (int k=0; k<K; k++) {
    arma::colvec alphaVec1 = alphaMtx1.row(k).t();
    arma::colvec alphaVec2 = alphaMtx2.row(k).t();
    arma::colvec bVecCol = bVec;

    arma::rowvec logD1 = logDMultDirMtx(data1, alphaVec1, bVecCol, J, P1);
    arma::rowvec logD2 = logDMultDirMtx(data2, alphaVec2, bVecCol, J, P2);

    arma::rowvec sumVec = logD1 + logD2 + log(tauVec(k));
    taufMtxLog.col(k) = sumVec.t();
  }

  return taufMtxLog;
}

arma::cube updateLogTauFCompare(arma::mat& data1, arma::mat& data2, arma::vec& tauVec, arma::mat& alphaMtx1, arma::mat& alphaMtx2, arma::vec& bVec, arma::vec& bVecNew, int J, int K) {

   arma::cube taufMtxLogCube=arma::zeros<arma::cube>(J,K,2);
	int P1, P2;
	arma::colvec alphaVec1, alphaVec2, bVecCol, bVecNewCol;
   arma::rowvec logD1, logD2, sumVec, logD1New, logD2New, sumVecNew;

   P1 = data1.n_rows;
   P2 = data2.n_rows;

   for (int k=0; k<K; k++) {
     alphaVec1 = alphaMtx1.row(k).t();
     alphaVec2 = alphaMtx2.row(k).t();
     bVecCol = bVec;
     bVecNewCol = bVecNew;

     logD1 = logDMultDirMtx(data1, alphaVec1, bVecCol, J, P1);
     logD2 = logDMultDirMtx(data2, alphaVec2, bVecCol, J, P2);
     sumVec = logD1 + logD2 + log(tauVec(k));
     taufMtxLogCube.slice(0).col(k) = sumVec.t();

     logD1New = logDMultDirMtx(data1, alphaVec1, bVecNewCol, J, P1);
     logD2New = logDMultDirMtx(data2, alphaVec2, bVecNewCol, J, P2);
     sumVecNew = logD1New + logD2New + log(tauVec(k));
     taufMtxLogCube.slice(1).col(k) = sumVecNew.t();
   }
   return taufMtxLogCube;
}

// [[Rcpp::export]]
arma::mat updateZmtx(arma::mat& taufMtxLog, int J, int K) {
  arma::mat zMtx;
  zMtx.zeros(J, K);

  for (int k=0; k<K; k++) {
    arma::colvec taufVec = taufMtxLog.col(k);
    arma::mat subtractMtx = genRepMtx(taufVec, J, K);
    arma::mat tempMtx1 = taufMtxLog-subtractMtx;
    arma::mat tempMtx2 = exp(tempMtx1);
    arma::colvec invVec = Arma_rowSums(tempMtx2);
    zMtx.col(k) = 1/invVec;
  }

  return zMtx;
}

arma::cube updateZmtxCompare(arma::cube& taufMtxLogCube, int J, int K) {
  arma::cube zMtxCube = arma::zeros<arma::cube>(J,K,2);
  arma::colvec taufVec, taufVecNew, invVec, invVecNew;
  arma::mat subtractMtx, tempMtx1, tempMtx2, subtractMtxNew, tempMtx1New, tempMtx2New;

  for (int k=0; k<K; k++) {
    taufVec = taufMtxLogCube.slice(0).col(k);
    subtractMtx = genRepMtx(taufVec, J, K);
    tempMtx1 = taufMtxLogCube.slice(0)-subtractMtx;
    tempMtx2 = exp(tempMtx1);
    invVec = Arma_rowSums(tempMtx2);
    zMtxCube.slice(0).col(k) = 1/invVec;

    taufVecNew = taufMtxLogCube.slice(1).col(k);
    subtractMtxNew = genRepMtx(taufVecNew, J, K);
    tempMtx1New = taufMtxLogCube.slice(1)-subtractMtxNew;
    tempMtx2New = exp(tempMtx1New);
    invVecNew = Arma_rowSums(tempMtx2New);
    zMtxCube.slice(1).col(k) = 1/invVecNew;
  }
  return zMtxCube;
}

// [[Rcpp::export]]
arma::rowvec updateNjVec(arma::mat& zMtx) {
  return(Arma_colSums(zMtx));
}

// [[Rcpp::export]]
arma::rowvec updateTauVec(arma::rowvec njVec, int K) {
  return njVec/K;
}

// [[Rcpp::export]]
double getLogL(arma::mat& zMtx, arma::mat& taufMtxLog, arma::vec& bVec, double sigmaB) {
  double logL1 = arma::accu(zMtx % taufMtxLog);
  double logL2 = arma::sum( -log(bVec) - pow(log(bVec),2)/(2*pow(sigmaB,2)) );
  double logL3 = zMtx.n_rows*(-log(sigmaB));
  double logLtot = logL1+logL2+logL3;
  return logLtot;
}


// [[Rcpp::export]]
arma::mat updateAlphaMtxOne(arma::mat& data0, arma::vec& z, arma::mat& alphaMtx, double sd_alpha, arma::vec& bVec, int K, int P) {
  for(int k=0; k<K; k++){
    int cl0 = k+1;
    arma::uvec typeInd = find(z >= cl0 && z <= cl0); // Find indices
    int numCell = typeInd.size();

    for(int p=0; p<P; p++){
      arma::vec alphaAll = alphaMtx.row(k).t();
      arma::vec likeNew = arma::zeros<arma::vec>(numCell);
      arma::vec likeOld = arma::zeros<arma::vec>(numCell);

      double m = alphaMtx(k,p);
      double alphaNew = Rcpp::as<double>(Rcpp::rnorm(1, m, sd_alpha));

      if (alphaNew > 500) {
         alphaNew = 1000 - alphaNew;
      }
      if (alphaNew < 0) {
         alphaNew = -alphaNew;
      }

      arma::vec alphaAllNew = alphaAll;
      alphaAllNew(p) = alphaNew;

      for (int j=0; j<numCell; j++) {
        int jNew = typeInd(j);
        double bj = bVec(jNew);

        arma::vec yVec = data0.col(jNew);

        likeOld(j) = logDMultDirFull(yVec, alphaAll, bj);
        likeNew(j) = logDMultDirFull(yVec, alphaAllNew, bj);
      }

      double likDiff = arma::sum(likeNew)-arma::sum(likeOld);
      double accP = std::min(0.0, likDiff);
      double u = Rcpp::as<double>(Rcpp::runif(1));

      if (log(u) < accP) {
         alphaMtx(k,p) = alphaNew;
      }
    }
  }
  return alphaMtx;
}

// [[Rcpp::export]]
arma::vec updateBMtxOne(arma::mat& data1, arma::mat& data2, arma::vec& z, arma::vec& tauVec, arma::mat& alphaMtx1, arma::mat& alphaMtx2, double sdB, arma::vec& bVec, double sigmaB, int J, int K) {
  for (int j=0; j<J; j++) {
    arma::vec bVecNew = bVec;
    double bOld = bVec(j);
    double bNew = Rcpp::as<double>(Rcpp::rnorm(1, bOld, sdB));
    if (bNew<0) {
      bNew = -bNew;
    }
    bVecNew(j) = bNew;

    arma::cube taufMtxLogCube = updateLogTauFCompare(data1, data2, tauVec, alphaMtx1, alphaMtx2, bVec,bVecNew, J, K);
    arma::cube zMtxCube=updateZmtxCompare(taufMtxLogCube,J,K);
    double likeOld = getLogL(zMtxCube.slice(0), taufMtxLogCube.slice(0), bVec, sigmaB);
    double likeNew = getLogL(zMtxCube.slice(1), taufMtxLogCube.slice(1), bVecNew, sigmaB);

    double accP = std::min(0.0, likeNew-likeOld);
    double u = Rcpp::as<double>(Rcpp::runif(1));
    if (log(u) < accP) {
      bVec(j) = bNew;
    }
  }
  return bVec;
}



