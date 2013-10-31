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

double meancpp(NumericVector x) {
  
  double s = 0;
  int n = x.size();

  for(int i = 0; i < n; i++) { 
    s = s + x[i];
  }

  return s/(double)n;
  
}

double sd(NumericVector x) {
  
  int n = x.size();
  double m = meancpp(x);
  double s = 0;
  
  for(int i = 0; i < n; i++) s = s + pow(x[i] - m, 2);
  
  return sqrt(s / ((double)n - 1));
  
}

// [[Rcpp::export]]
double u3(NumericVector x) {
  
  int n = x.size();
  double m = meancpp(x);
  double u = 0;
  for(int i = 0; i < n; i++) u = u + pow(x[i] - m, 3);
    
  return u / (((double)n - 1) * ((double)n - 2));
  
}

// [[Rcpp::export]]
double skew(NumericVector x) {
  
  double n = x.size();
  double u = u3(x);
  double s = sd(x);
  
  return n * u / pow(s, 3);
  
}

double tvalue(NumericVector x, LogicalVector d = NA_LOGICAL, double mean = NA_REAL) {
  
  int n = x.size();
  double mi;
  double si;
  
  if(d[0] == NA_LOGICAL) {
    mi = meancpp(x);
    double sd = 0;
    for(int i = 0; i < n; i++) sd = sd + pow(x[i] - mi, 2);
    si = sqrt(sd);
  } else {
    NumericVector p = ple(x, d);
    mi = p[0];
    si = p[2];
  }
  
  double ti = sqrt((double)n) * ((mi - mean) / (si * sqrt((double)n)));
  
  return ti;
  
  
}

double HallW(NumericVector x, double mean) {
  
  int n = x.size();
  int mi = meancpp(x);
  double sd = 0;
  for(int i = 0; i < n; i++) sd = sd + pow(x[i] - mi, 2);
  double si = sqrt(sd) / (double(n) - 1);
  double wi = (mi - mean) / si;
  double k3 = skew(x);
  
  return wi + k3 * pow(wi,2) / 3 + pow(k3, 2) * pow(wi, 3) / 27 + k3 / (6 * (double)n);  
  
}

List bootstrapSample(NumericVector x, LogicalVector d = NA_LOGICAL) {
  
  int n = x.size();
  int r;
  bool cen = true;
  int ds = 0;
  NumericVector bsx(n);
  LogicalVector bsd(n);
  if(d[0] == NA_LOGICAL) cen = false;

  do {

    ds = 0;

    for(int j = 0; j < n; j++) { // Inner Loop for Bootstrap Samples
  
        r = rand() % n; 
        bsx[j] = x[r];
        if(cen) {
          bsd[j] = d[r];
          if(d[r]) ds++;
        }
  
    }
  
  } while (cen && ds <= 3);
  
  if(cen) {
    return List::create(Named("x") = bsx,
                        Named("d") = bsd);
  } else {
    return List::create(Named("x") = bsx,
                        Named("d") = NA_LOGICAL);
  }
  
}

// [[Rcpp::export]]
NumericVector bootstrap(NumericVector x, LogicalVector d = NA_LOGICAL, int N = 2000) {
  
  bool cen = (d[0] != NA_LOGICAL);
  double mean;
  NumericVector means(N);
  
  List bs;
  
  for(int i = 0; i < N; i++) {
   
    bs = bootstrapSample(x, d);
    if(cen) {
      mean = ple(bs["x"], bs["d"])[0];
    } else {
      mean = meancpp(bs["x"]);
    }
   
    means[i] = mean;
   
  }
  
  return means;
  
}

NumericVector boottvalue(NumericVector x, LogicalVector d = NA_LOGICAL, int N = 2000) {
  
  bool cen = (d[0] != NA_LOGICAL);
  double t;
  NumericVector tvalues(N);
  double mean;
  if(cen) {
    mean = ple(x,d)[0];
  } else {
    mean = meancpp(x);
  }
  
  List bs;
  
  for(int i = 0; i < N; i++) {
   
    bs = bootstrapSample(x, d);
    t = tvalue(bs["x"], bs["d"], mean);
   
    tvalues[i] = t;
   
  }
  
  return tvalues;
  
}

NumericVector bootHallw(NumericVector x, int N = 2000) {
  
  double w;
  NumericVector wvalues(N);
  double mean = meancpp(x);

  List bs;
  
  for(int i = 0; i < N; i++) {
   
    bs = bootstrapSample(x);
    w = HallW(bs["x"], mean);
   
    wvalues[i] = w;
   
  }
  
  return wvalues;
  
}