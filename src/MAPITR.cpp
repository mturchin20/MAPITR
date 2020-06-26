// load Rcpp
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return X.t()*X/p;
}

////////////////////////////////////////////////////////////////////////////

//Below is a function for MAPITR looking for interaction effects for pathways

////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List MAPITRBase(arma::mat Y,arma::mat X,List regions,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int nsnp = X.n_rows;
    const int p = regions.size();
    //const int q = Z.n_rows;
    
//    Rcout << "nana1" << endl;

    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
  

    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
	//Extract phenotype
	vec y = Y.col(i);

        //Pre-compute the Linear GSM
        uvec j = regions[i];
        
        //Compute K covariance matrices
        mat K = (GSM*nsnp-GetLinearKernel(X.rows(j))*j.n_elem)/(nsnp-j.n_elem-1);
        mat G = GetLinearKernel(X.rows(j))%K;
        
        //Transform K and G using projection M
        mat b = zeros(n);
	b.col(0) = ones<vec>(n); 
        mat btb_inv = inv(b.t()*b);
        mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(3); //Create k-vector q to save
        mat S = zeros(3,3); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Kc*yc);
        q(1) = as_scalar(yc.t()*Gc*yc);
        q(2) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Kc%Kc));
        S(0,1) = as_scalar(accu(Kc%Gc));
        S(0,2) = as_scalar(accu(Kc%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(Gc%Gc));
        S(1,2) = as_scalar(accu(Gc%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(2,0) = S(0,2);
        S(2,1) = S(1,2);
        S(2,2) = as_scalar(accu((eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Record nu^2, and tau^2 under the null hypothesis
        vec q_sub = zeros(2);
        mat S_sub = zeros(2,2);
        
        q_sub(0)=q(0);
        q_sub(1)=q(2);
        
        S_sub(0,0)=S(0,0);
        S_sub(0,1)=S(0,2);
        
        S_sub(1,0)=S(2,0);
        S_sub(1,1)=S(2,2);
        
        //Compute P and P^{1/2} matrix
        vec delta_null = inv(S_sub)*q_sub;
        
        vec eigval;
        mat eigvec;
        
        eig_sym(eigval,eigvec,delta_null(0)*Kc+delta_null(1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()));
        
        //Find the eigenvalues of the projection matrix
        vec evals;
        
        eig_sym(evals, (eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0))))*(Sinv(1,0)*Kc+Sinv(1,1)*Gc+Sinv(1,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0)))));
        Lambda.col(i) = evals;
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(1);
        
        //Compute the PVE
        pve(i) = delta(1)/accu(delta);
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}
