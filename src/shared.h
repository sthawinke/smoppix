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

inline int to_pos_int(double x);
inline double to_dbl(int x);
#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
inline bool is_large_int(double x);