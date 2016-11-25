#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix my_mmult( NumericMatrix m , NumericVector v , int a, int b){
  
  if( ! (m.ncol() == v.size()) ) stop("Non-conformable arrays") ;
  
  NumericMatrix out(b-a+1 , m.ncol()) ;
  
  for (int j = 0; j < m.ncol(); j++) {
    for (int i = 0; i < b-a+1; i++) {
      out(i,j) = m(i+a,j) * v[j];
    }
  }
  
  return out ;
}

// [[Rcpp::export]]
NumericMatrix mmult( NumericMatrix m , NumericVector v , bool byrow = true ){
  if(byrow){
    if( ! (m.ncol() == v.size()) ) stop("Non-conformable arrays") ;
  }
  if( ! byrow ){
    if( ! (m.nrow() == v.size()) ) stop("Non-conformable arrays") ;
  }
  
  NumericMatrix out(m) ;
  
  if( byrow ){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if( ! byrow ){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
}

// wsum_current uses essentially the code from the current fixedpointfunction
// for summing up the responsibilities
// [[Rcpp::export]]
NumericVector wsum_current(NumericVector pi,const NumericMatrix & matrix_lik,int a,int b){
  int n=matrix_lik.nrow(), k=matrix_lik.ncol();
  NumericVector wsum(k);
  NumericMatrix m(n,k);
  NumericMatrix classprob(m);
  NumericVector m_rowsum(n);
  //IntegerVector subset(prior);
  for (int i=0;i<k;i++){
    m.column(i)=pi[i]*matrix_lik.column(i);
    m_rowsum=m_rowsum+m.column(i);
  }
  for (int i=0;i<k;i++){
    classprob.column(i)=classprob.column(i)/m_rowsum;
  }
  
  for (int i=0;i<k;i++){
    wsum[i]=sum(classprob.column(i));
  }
  return(wsum);
}
  
    
// an attempt to optimize/streamline the current procedure
// [[Rcpp::export]]
NumericVector wsum_current_opt(NumericVector pi,const NumericMatrix & matrix_lik,int a,int b){
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
  
// this is the code from the current fixptfn (does more than just compute wsum)
// [[Rcpp::export]]
List fixptfn(NumericVector pi_est,const NumericMatrix & matrix_lik, NumericVector prior){
  int n=matrix_lik.nrow(), k=matrix_lik.ncol();
  NumericVector pi_new(k);
  for (int i=0;i<k;i++){
    pi_est[i]=std::max(0.0,pi_est[i]);
  }
  pi_est=pi_est/sum(pi_est); //normalize pi
  
  double loglik,lpriordens=0.0;
  NumericMatrix m(n,k);
  NumericMatrix classprob(m);
  NumericVector m_rowsum(n);
  //IntegerVector subset(prior);
  for (int i=0;i<k;i++){
    m.column(i)=pi_est[i]*matrix_lik.column(i);
    m_rowsum=m_rowsum+m.column(i);
  }
  for (int i=0;i<k;i++){
    classprob.column(i)=classprob.column(i)/m_rowsum;
  }
  //calculating objective value--probability
  loglik=sum(log(m_rowsum));
  
  for (int i=0;i<k;i++){
    if(prior[i]!=1.0){
      lpriordens +=(prior[i]-1.0)*log(pi_est[i]);
    }
    
  }
  //generating new pi
  for (int i=0;i<k;i++){//set any estimates that are less than zero, which can happen with prior<1, to 0
    pi_new[i]=std::max(0.0,sum(classprob.column(i))+prior[i]-1.0);
  }
  
  pi_new=pi_new/sum(pi_new); //normalize pi
  
  return(List::create(Named("fixedpointvector")=pi_new,
                      Named("objfn")=-loglik-lpriordens));
}


// [[Rcpp::export]]
NumericVector resp(NumericVector pi,NumericVector lik){
  return(pi*lik/sum(pi*lik));
}

// [[Rcpp::export]]
NumericVector resp2(const NumericVector & pi,const NumericMatrix & lik, int i){
  return(pi*lik.row(i)/sum(pi*lik.row(i)));
}

// [[Rcpp::export]]
NumericVector wsum_naive(NumericVector pi,NumericMatrix & matrix_lik,int a,int b){
  NumericVector sum = resp(pi, matrix_lik.row(a));

  for(int i = a+1; i<=b; i++){
    sum += resp(pi, matrix_lik.row(i));
  }
  return(sum);
}

// [[Rcpp::export]]
NumericVector wsum_naive_transposed(NumericVector pi,const NumericMatrix & matrix_lik,int a,int b){
  NumericVector sum = resp(pi, matrix_lik.column(a));
  
  for(int i = a+1; i<=b; i++){
    sum += resp(pi, matrix_lik.column(i));
  }
  return(sum);
}

// [[Rcpp::export]]
double maxdiff(const NumericVector & a,const NumericVector & b){
  return(max(abs(b-a)));
}


// [[Rcpp::export]]
NumericVector wsum(const NumericVector & pi,const NumericMatrix & matrix_lik, int a, int b, NumericVector wa, NumericVector wb,
                   double tol){
  
  if(wa.length()<2){
    wa = resp2(pi,matrix_lik,a);
  }
  
//if((a-b) < ntol){return(wsum_current_opt(pi, matrix_lik, a, b);}  
  if(a==b){return(wa);}  
  
  if(wb.length()<2){
    wb = resp2(pi, matrix_lik,b); //resp(pi,matrix_lik.row(b));
  }
  if(maxdiff(wa,wb)<tol){
    //Rcout << a << "," << b << "\n";
    return (b-a+1) * 0.5*(wa+wb);
  } else {
    int c = trunc(0.5*(a+b));
    return(wsum(pi, matrix_lik, a, c, wa, 0, tol) + wsum(pi, matrix_lik, c+1, b, 0, wb, tol));
  }
}

// [[Rcpp::export]]
List fixptfn2(NumericVector pi_est,NumericMatrix matrix_lik, NumericVector prior){
  int n=matrix_lik.nrow(), k=matrix_lik.ncol();
  NumericVector pi_new(k);
  for (int i=0;i<k;i++){
    pi_est[i]=std::max(0.0,pi_est[i]);
  }
  pi_est=pi_est/sum(pi_est); //normalize pi
  
  double loglik,lpriordens=0.0;
  NumericMatrix m(n,k);
  NumericMatrix classprob(m);
  NumericVector m_rowsum(n);
  //IntegerVector subset(prior);
  for (int i=0;i<k;i++){
    m.column(i)=pi_est[i]*matrix_lik.column(i);
    m_rowsum=m_rowsum+m.column(i);
  }
  for (int i=0;i<k;i++){
    classprob.column(i)=classprob.column(i)/m_rowsum;
  }
  //calculating objective value--probability
  loglik=sum(log(m_rowsum));
  
  for (int i=0;i<k;i++){
    if(prior[i]!=1.0){
      lpriordens +=(prior[i]-1.0)*log(pi_est[i]);
    }
    
  }
  //generating new pi
  for (int i=0;i<k;i++){//set any estimates that are less than zero, which can happen with prior<1, to 0
    pi_new[i]=std::max(0.0,sum(classprob.column(i))+prior[i]-1.0);
  }
  
  pi_new=pi_new/sum(pi_new); //normalize pi
  
  return(List::create(Named("fixedpointvector")=pi_new,
                      Named("objfn")=-loglik-lpriordens));
}

