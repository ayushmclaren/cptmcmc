#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int cpt_uni_modf(NumericVector chain)
{
  int n = chain.length();
  NumericVector z(n);
  z[0] = 0;
  int localmax_idx = 0;

  for(int i=0; i< n-1; i++)
  {
    z[i+1] = z[i] + (chain[i] - Rcpp::mean(chain[Rcpp::Range(i+1, n)]))/std::sqrt(n);

    if(std::fabs(z[i])>std::fabs(z[i+1]))
    {
      localmax_idx = i+1;
      break;
    }
  }
//Rcpp::Rcout << "\nModified LocalMax method:\nThe estimated break pt. is at:" << localmax_idx << "\n";
  return(localmax_idx);
}

// [[Rcpp::export]]
NumericVector cpt_multi_modf(NumericMatrix chain)
{
  int n = chain.cols();
  NumericVector breaks(n);
  for(int i=0; i<n; i++)
  {
    breaks[i] = cpt_uni_modf(chain(_,i));
  }
  return(breaks);
}
