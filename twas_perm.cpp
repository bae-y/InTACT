#define ARMA_NO_DEBUG
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec baseSVD(arma::mat &X)
{
    arma::mat U, V;
    arma::vec S;
    arma::svd(U, S, V, X);
    return S;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat cSum(arma::vec U, arma::mat CovS)
{
    
    //const int n = X1.n_rows;
    const int n = U.size();
    
    double pval;
    arma::mat testStat;
    
    testStat = arma::sum(pow(U, 1));
    
    arma::vec a;
    a.ones(n);
    arma::mat denom = sqrt(a.t() * CovS * a);
    
    double testStat2 = testStat(0) / denom(0, 0);
    pval = R::pchisq(pow(testStat2, 2), 1.0, 0, FALSE);
    
    arma::vec res(2);
    res.fill(0);
    res(0) = testStat2;
    res(1) = pval;
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
double cSSU(arma::vec U, arma::mat CovS)
{
    
    //const int n = X1.n_rows;
    const int n = U.size();
    
    double pval;
    arma::vec testStat;
    
    if (n == 1) {
        arma::mat tmpres = cSum(U,CovS);
        pval = tmpres(1);
    } else {
        testStat = arma::sum(pow(U, 2));
        arma::vec cr = baseSVD(CovS);
        
        arma::vec alpha1 = arma::sum(pow(cr, 3)) / arma::sum(pow(cr, 2));
        
        arma::vec beta1 = arma::sum(cr) - pow(arma::sum(pow(cr, 2)), 2) / arma::sum(pow(cr, 3));
        arma::vec d1 = pow(arma::sum(pow(cr, 2)), 3) / pow(arma::sum(pow(cr, 3)), 2);
        
        double d12 = d1(0);
        testStat = (testStat - beta1) / alpha1;
        double test2 = testStat(0);
        pval = R::pchisq(test2, d12, 0, FALSE);
    }
    
    return pval;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0W2(arma::mat &CvSqrt, arma::mat &wgt, Rcpp::List &indx_mat, Rcpp::List &LD, arma::mat &MultMat, double MultDf, int nperm)
{
    
    //const int n = X1.n_rows;
    const int n_weight = wgt.n_cols;
    const int k = CvSqrt.n_rows;
    
    // containers
    arma::mat T0s(nperm, n_weight - 1);
    T0s.fill(0);
    
    arma::mat T1s(nperm, 3);
    T1s.fill(0);
    
    for (int i = 0; i < nperm; i++)
    {
        arma::mat U00 = arma::randn(k, 1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int j = 0; j < (n_weight - 1); j++)
        {
            arma::vec wgtused = wgt.col(j);
            arma::vec U0tmp = wgtused % U0;
            
            arma::vec rowidx_in = indx_mat[j];
            arma::uvec rowidx = arma::conv_to<arma::uvec>::from(rowidx_in) - 1;
            
            /*arma::mat s1 = CvSqrt.submat(rowidx, rowidx);
             arma::vec wgtused2 = wgtused.elem(rowidx);
             arma::mat LDused = arma::diagmat(wgtused2) * s1 *  arma::diagmat(wgtused2).t();*/
            arma::vec Utmp = U0tmp.elem(rowidx);
            arma::mat LDused = LD[j];
            arma::vec Res = cSum(Utmp, LDused);
            T0s(i, j) = Res(0);
        }
        
        int j = n_weight - 1;
        arma::vec wgtused = wgt.col(j);
        arma::vec U0tmp = wgtused % U0;
        
        arma::vec rowidx_in = indx_mat[j];
        arma::uvec rowidx = arma::conv_to<arma::uvec>::from(rowidx_in) - 1;
        
        arma::vec Utmp = U0tmp.elem(rowidx);
        arma::mat LDused = LD[j];
        arma::vec Res = cSum(Utmp, LDused);
        
        T1s(i, 0) = Res(1);
        T1s(i, 1) = cSSU(Utmp, LDused);
        
        arma::mat MultStat = T0s.row(i) * MultMat * T0s.row(i).t();
        
        double tmpStat = MultStat(0, 0);
        
        T1s(i, 2) = R::pchisq(tmpStat, MultDf, 0, FALSE);
    }
    
    Rcpp::List res;
    res["T0s"] = T0s;
    res["T1s"] = T1s;
    
    return (res);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0Wsingle(arma::mat &CvSqrt, arma::mat &wgt, Rcpp::List &indx_mat, Rcpp::List &LD, int nperm)
{
    
    //const int n = X1.n_rows;
    const int n_weight = wgt.n_cols;
    const int k = CvSqrt.n_rows;
    
    // containers
    arma::mat T0s(nperm, n_weight + 1);
    T0s.fill(0);
    
    arma::mat T1s(nperm, 3);
    T1s.fill(0);
    
    for (int i = 0; i < nperm; i++)
    {
        arma::mat U00 = arma::randn(k, 1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int j = 0; j < n_weight; j++)
        {
            arma::vec wgtused = wgt.col(j);
            arma::vec U0tmp = wgtused % U0;
            
            arma::vec rowidx_in = indx_mat[j];
            arma::uvec rowidx = arma::conv_to<arma::uvec>::from(rowidx_in) - 1;
            
            /*arma::mat s1 = CvSqrt.submat(rowidx, rowidx);
             arma::vec wgtused2 = wgtused.elem(rowidx);
             arma::mat LDused = arma::diagmat(wgtused2) * s1 *  arma::diagmat(wgtused2).t();*/
            arma::vec Utmp = U0tmp.elem(rowidx);
            arma::mat LDused = LD[j];
            arma::vec Res = cSum(Utmp, LDused);
            if(j == 0 ) {
                T0s(i, j) = Res(1);
            } else {
                T0s(i, j) = Res(1);
                T0s(i, j+1) = cSSU(Utmp, LDused);
            }
        }
    }
    
    Rcpp::List res;
    res["T0s"] = T0s;
    
    return (res);
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec calcMultiVar(arma::mat &LD, arma::mat &wgt, Rcpp::List &indx_mat)
{
    
    //const int n = X1.n_rows;
    const int n_weight = wgt.n_cols;
    const int k = LD.n_rows;
    
    // containers
    arma::vec var(n_weight);
    var.fill(0);
    
    for (int i = 0; i < n_weight; i++)
    {
        
        arma::vec idx_in = indx_mat[i];
        int n_nonzero = idx_in.n_elem;
        
        double tmp_res = 0;
        for (int k1 = 0; k1 <n_nonzero; k1++) {
            for (int k2 =0; k2 <n_nonzero; k2++) {
                
                int tmp1 = idx_in[k1];
                int tmp2 = idx_in[k2];
                tmp_res = tmp_res + wgt(tmp1,i) * wgt(tmp2,i) * LD(tmp1, tmp2);
            }
        }
        var(i) = tmp_res;
    }
    
    return (var);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat calcMultiCov(arma::mat &LD, arma::mat &wgt, Rcpp::List &indx_mat,arma::vec &var)
{
    //const int n = X1.n_rows;
    const int n_weight = wgt.n_cols;
    const int k = LD.n_rows;
    
    // containers
    arma::mat T0s(n_weight,n_weight);
    T0s.fill(1.0);
    
    for (int i = 0; i < n_weight; i++)
    {
        for(int j = (i+1); j< n_weight; j++) {
            arma::vec idxi = indx_mat[i];
            int n_nonzeroi = idxi.n_elem;
            
            arma::vec idxj = indx_mat[j];
            int n_nonzeroj = idxj.n_elem;
            
            double tmp_res = 0;
            for (int k1 = 0; k1 < n_nonzeroi; k1++) {
                for (int k2 =0; k2 < n_nonzeroj; k2++) {
                    
                    int tmp1 = idxi[k1];
                    int tmp2 = idxj[k2];
                    tmp_res = tmp_res + wgt(tmp1,i) * wgt(tmp2,j) * LD(tmp1, tmp2);
                }
            }
            tmp_res = tmp_res / sqrt(var(i) * var(j));
            T0s(i,j) = tmp_res;
            T0s(j,i) = tmp_res;
        }
    }
    
    return (T0s);
}

