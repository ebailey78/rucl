#include <Rcpp.h> 
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ple(NumericVector x, LogicalVector d) {

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
  
  return Rcpp::NumericVector::create(Rcpp::Named("mean") = mean,
                                     Rcpp::Named("var") = var,
                                     Rcpp::Named("se") = se);
  
}

double sumcpp(NumericVector x) {
 
  double s = 0;
  int n = x.size();

  for(int i = 0; i < n; i++) { 
    s = s + x[i];
  }

  return s;
  
}

double meancpp(NumericVector x) {
  
  return sumcpp(x)/double(x.size());
  
}

// [[Rcpp::export]]
double tic(NumericVector x, LogicalVector d, double m) {
  
  int n = x.size();
  NumericVector p = ple(x, d);
  double ti = sqrt((double)n) * ((p["mean"] - m) / (p["se"] * sqrt(n)));

  return ti;
  
}

// [[Rcpp::export]]
double tiu(NumericVector x, double m) {
  
  int n = x.size();
  double mi = meancpp(x);
  double sd = 0;
  for(int i = 0; i < n; i++) sd = sd + pow(x[i] - mi, 2);
  double si = sqrt(sd)/((double)n - 1);
  
  double ti = sqrt((double)n) * ((mi - m) / si);
  
  return ti;
  
}

// [[Rcpp::export]]
NumericVector BSmean(NumericVector x, int N) {

  NumericVector means(N);
  int n = x.size();

  double bss;
  
  for(int i = 0; i < N; i++) {

    bss = 0;

    for(int j = 0; j < n; j++) {

      bss = bss + x[rand() % n];

    }

    means[i] = bss / n;

  }
  
  return means;
  
}

// [[Rcpp::export]]
NumericVector BSple(NumericVector x, LogicalVector d, int N) {
  
  int n = x.size();
  int r;
  int ds = 0;
  NumericVector m;
  NumericVector z;
    
  NumericVector bsx(n);
  LogicalVector bsd(n);
  NumericVector means(N);
  
  for(int i = 0; i < N; i++) {

    for(int j = 0; j < n; j++) {

      r = rand() % n;
      bsx[j] = x[r];
      bsd[j] = d[r];
      if(d[r] == true) {
        ds++;
      }

    }

    if(ds >= 3) {

      m = ple(bsx, bsd);
    
      z = m["mean"];
    
      means[i] = m["mean"];
      
    } else {
      i--;
    }

  }
  
  return means;
  
}