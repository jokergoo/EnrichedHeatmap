#include <Rcpp.h>
using namespace Rcpp;


float closeness(NumericVector x1, NumericVector x2) {
	float y = 0;
	int y_n = 0;
	int n = x1.size();

	for(int i = 0; i < n; i ++) {
		for(int j = 0; j < n; j ++) {
			if( (std::abs(x1[i] - 0) < 1e-6) | (std::abs(x2[j] - 0) < 1e-6) ) {
				continue;
			}
			y += std::abs(i - j);
			y_n += 1;
		}
	}

	if(y_n == 0) {
		return 0;
	} else {
		return y/y_n;
	}
}

// [[Rcpp::export]]
NumericMatrix dist_by_closeness(NumericMatrix mat) {
	int n = mat.nrow();

	NumericMatrix dist(n, n);

	for(int i = 0; i < n - 1; i ++) {
		for(int j = i+1; j < n; j ++) {
			dist(i, j) = closeness(mat(i, _), mat(j, _));
			dist(j, i) = dist(i, j);
		}
	}
	for(int i = 0; i < n; i ++) {
		dist(i, i) = 0;
	}
	return dist;
}

