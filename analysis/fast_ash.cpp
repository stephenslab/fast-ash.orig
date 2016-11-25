#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix my_mmult(const NumericMatrix & m , const NumericVector & v , int a, int b){
  
  if( ! (m.ncol() == v.size()) ) stop("Non-conformable arrays") ;
  
  NumericMatrix out(b-a+1 , m.ncol()) ;
  
  for (int j = 0; j < m.ncol(); j++) {
    for (int i = 0; i < b-a+1; i++) {
      out(i,j) = m(i+a,j) * v[j];
    }
  }
  
  return out ;
}


// directly sum the responsibilities between a and b
// [[Rcpp::export]]
NumericVector wsum_direct(NumericVector pi,const NumericMatrix & matrix_lik,int a,int b){
    int k=matrix_lik.ncol();
    NumericVector wsum(k);
    NumericVector classprob_rowsum(b-a+1);
    NumericMatrix classprob = my_mmult(matrix_lik, pi, a, b);
    
    for (int i=0;i<k;i++){
      classprob_rowsum=classprob_rowsum+classprob.column(i);
    }
    
    
    for (int i=0;i<k;i++){
      classprob.column(i) = classprob.column(i)/classprob_rowsum;
      wsum[i] = sum(classprob.column(i));
    }
    
    return(wsum);
}
  

// compute responsibility (posterior class prob) for row i
// [[Rcpp::export]]
NumericVector resp(const NumericVector & pi,const NumericMatrix & lik, int i){
  return(pi*lik.row(i)/sum(pi*lik.row(i)));
}

// [[Rcpp::export]]
double maxdiff(const NumericVector & a,const NumericVector & b){
  return(max(abs(b-a)));
}

// computes sum of responsibilities from ath to bth row (inclusive)
// performs binary search strategy, and assumes responsibilities between
// a and b are all similar if resp(a) is within tol of resp(b)
// [[Rcpp::export]]
NumericVector wsum(const NumericVector & pi,const NumericMatrix & matrix_lik, int a, int b, NumericVector wa, NumericVector wb,
                   double tol, double ntol=1){
  
  if((b-a) < ntol){return(wsum_direct(pi, matrix_lik, a, b));}  
 
  if(wa.length()<2){ 
    wa = resp(pi,matrix_lik,a);
  } 
  if(wb.length()<2){
    wb = resp(pi, matrix_lik,b); 
  }
  if(maxdiff(wa,wb)<tol){ // if within tolerance, just average the two and multiply by number of rows
    //Rcout << a << "," << b << "\n";
    return (b-a+1) * 0.5*(wa+wb);
  } else { // split interval [a,b] in half and sum each half
    int c = trunc(0.5*(a+b));
    return(wsum(pi, matrix_lik, a, c, wa, 0, tol, ntol) + wsum(pi, matrix_lik, c+1, b, 0, wb, tol, ntol));
  }
}


