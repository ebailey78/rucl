#include <Rcpp.h> 
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List ple(NumericVector x, LogicalVector d) {

  std::vector<double> v;

  int n = x.size();
  int vstart = 2;

  for(int i = 0; i < n; i++) 
    if(d[i] == true)
      v.insert(v.end(), x[i]);

  std::sort(v.begin(), v.end());
  std::reverse(v.begin(), v.end());

  double m = *std::min_element(x.begin(), x.end());

  if(v[v.size()] != m) {
    vstart = 3;
    v.insert(v.end(), m);
  }
    
  double it;  
  it = std::distance(v.begin(), std::unique(v.begin(), v.end()));

  v.resize(it);
  int vn = v.size();
  
  std::vector<double> B(vn);
  std::vector<double> C(vn);
  std::vector<double> trec(vn);
  std::vector<double> brec(vn);
  double pple = 1;
  double mean = 0;
  double var = 0;
  
  for(int i = 0; i < vn; i++) {

    for(int j = 0; j < n; j++) {
    
      if(x[j] <= v[i]) {
        B[i]++;
      }
      if(x[j] == v[i]) {
        if(d[j] == true) {
          C[i]++;
        }
      }
    }

    pple = (B[i] - C[i]) / B[i] * pple;
    brec[i] = pple * (v[i] - v[i+1]);

    mean = mean + (1 - pple) * (v[i] - v[i+1]);

  }

  double prev = 0;

  for(int i = vn-vstart; i >= 0; i--) {
    
    prev = brec[i] + prev;
    var = var + ((prev * prev) / (B[i] * (B[i] - C[i]))) * C[i];
    
  }

  mean = mean + v[vn-1];
  var = var * ((double)vn / ((double)vn-1));
  double se = sqrt(var);
  
  return Rcpp::List::create(Rcpp::Named("mean") = mean,
                            Rcpp::Named("var") = var,
                            Rcpp::Named("se") = se);
  
}

// [[Rcpp::export]]
NumericVector BSmean(NumericVector x, int N) {
  
  NumericVector means(N);
  int n = x.size();

  double bss;
  
  for(int i = 0; i < N; i++) {

    bss = 0;

    for(int j = 0; j < x.size(); j++) {

      bss = bss + x[rand() % n];

    }

    means[i] = bss / n;

  }
  
  return means;
  
}

// [[Rcpp::export]]
NumericVector BSple(NumericVector x, LogicalVector d, int N) {
  
  NumericVector means(N);
  int n = x.size();
  int r;
  NumericVector m;
  
  NumericVector bsx(n);
  LogicalVector bsd(n);
  
  for(int i = 0; i < N; i++) {

    for(int j = 0; j < x.size(); j++) {

      r = rand() % n;
      bsx[j] = x[r];
      bsd[j] = d[r];

    }

    m = ple(bsx, bsd)[1];
  
    means[i] = m[1];

  }
  
  return means;
  
}