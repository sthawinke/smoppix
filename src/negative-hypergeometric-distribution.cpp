// Taken from extraDistr
#include <Rcpp.h>
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

inline bool is_large_int(double x) {
  if (x > std::numeric_limits<int>::max())
    return true;
  return false;
}

inline double to_dbl(int x) {
  return static_cast<double>(x);
}

inline int to_pos_int(double x) {
  if (x < 0.0 || ISNAN(x))
    Rcpp::stop("value cannot be coerced to integer");
  if (is_large_int(x))
    Rcpp::stop("value out of integer range");
  return static_cast<int>(x);
}

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
//inline int to_pos_int(double x);
//inline double to_dbl(int x);
//inline bool is_large_int(double x);

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;

bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

std::vector<double> nhyper_table(
    double n, double m, double r,
    bool cumulative = false
) {
  
  if (n < 0.0 || m < 0.0 || r < 0.0 || r > m)
    Rcpp::stop("inadmissible values");
  
  double j, N, start_eps;
  int ni = to_pos_int(n);
  N = m+n;
  
  std::vector<double> t(ni), h(ni), p(ni+1);
  start_eps = 1e-200;
  h[0] = start_eps * r*n/(N-r);
  t[0] = start_eps + h[0];
  
  for (int i = 1; i <= ni-1; i++) {
    j = to_dbl(i) + r;
    h[i] = h[i-1] * j*(n+r-j)/(N-j)/(j+1.0-r);
    t[i] = t[i-1] + h[i];
  }
  
  p[0] = start_eps / t[ni-1];
  
  if (cumulative) {
    for (int i = 1; i < ni; i++)
      p[i] = t[i-1] / t[ni-1];
    p[ni] = 1.0;
  } else {
    for (int i = 1; i <= ni; i++)
      p[i] = h[i-1] / t[ni-1];
  }
  
  return p;
}

// [[Rcpp::export]]
NumericVector cpp_pnhyper(
    const NumericVector& x,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  if (std::min({x.length(), n.length(),
               m.length(), r.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    n.length(),
    m.length(),
    r.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 100 == 0)
      Rcpp::checkUserInterrupt();
    
#ifdef IEEE_754
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(n, i)) ||
        ISNAN(GETV(m, i)) || ISNAN(GETV(r, i))) {
      p[i] = GETV(x, i) + GETV(n, i) + GETV(m, i) + GETV(r, i);
      continue;
    }
#endif
    
    if (GETV(r, i) > GETV(m, i) || GETV(n, i) < 0.0 ||
        GETV(m, i) < 0.0 || GETV(r, i) < 0.0 ||
        !isInteger(GETV(n, i), false) ||
        !isInteger(GETV(m, i), false) ||
        !isInteger(GETV(r, i), false)) {
        throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < GETV(r, i)) {
      p[i] = 0.0;
    } else if (GETV(x, i) >= (GETV(n, i) + GETV(r, i))) {
      p[i] = 1.0;
    } else if (is_large_int(GETV(x, i))) {
      p[i] = NA_REAL;
      Rcpp::warning("NAs introduced by coercion to integer range");
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(
        static_cast<int>(i % n.length()),
        static_cast<int>(i % m.length()),
        static_cast<int>(i % r.length())
      )];
      
      if (!tmp.size()) {
        tmp = nhyper_table(GETV(n, i), GETV(m, i), GETV(r, i), true);
      }
      p[i] = tmp[to_pos_int( GETV(x, i) - GETV(r, i) )];
      
    }
  } 
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}

