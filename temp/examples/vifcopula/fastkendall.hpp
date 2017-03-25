#ifndef _GOODMAN_FAST_KENDALL_H_
#define _GOODMAN_FAST_KENDALL_H_

#include <algorithm>
#include <vector>
#include <stdint.h>

#include "goodman_kruskal_gamma.hpp"

/*
Package: fastkendall
Title: Fast Calculation of Kendall's tau and the Goodman-Kruskal gamma
Version: 0.99
Date: 2012-12-13
Author: Karsten Looschen
Maintainer: Karsten Looschen <karsten.looschen@gmail.com>
Description: Fast Calculation of Kendall's tau and the Goodman-Kruskal gamma in O(n log(n))
License: GPL-3
URL: https://github.com/looschen/fastkendall

*/

template<typename compare_t = std::less<double>, typename T = std::vector<double> >
struct compare_other{
  // Functor for use with std::sort.
  // Don't compare directly but use i, j as indices for another vector/array
  compare_other(const T& other_): other(&other_), compare(compare_t()) {}

  bool operator()(size_t i, size_t j) const{
    return compare((*other)[i], (*other)[j]);
  }

private:
  const T* other;
  compare_t compare;
};

template <typename T_u, typename T_v>
double fast_kendall(const T_u& u, const T_v& v, size_t N){
    // prepare for and do concordance_count() from goodman_kruskal_gamma.hpp,
    // returns numeric vector
    // c(concordant, discordant, extraX, extraY, spare)

    double* Xptr = u;
    double* Yptr = v;

    // for(size_t i = 0; i != N; ++i)
    //   std::cout << Xptr[i] << " ";
    // std::cout << "\n";
    //   for(size_t i = 0; i != N; ++i)
    // 	std::cout << Yptr[i] << " ";
    //   std::cout << "\n";

    // preparation
    // sort X and Y
    std::vector<size_t> order_X(N); // use uint32_t if memory becomes a problem
    // get order with compare_other
    for(size_t i = 0; i != order_X.size(); ++i)
      order_X[i] = i;
    sort(order_X.begin(), order_X.end(), compare_other<std::less<double>, double*>(Xptr));

    // reorder X and Y
    std::vector<double> X_ord(N);
    std::vector<double> Y_ord(N);
    for(size_t i = 0; i != N; ++i)
      X_ord[i] = Xptr[order_X[i]];
    for(size_t i = 0; i != N; ++i)
      Y_ord[i] = Yptr[order_X[i]];


    // sort Y for equal values of X
    secondary_sort(X_ord.begin(), X_ord.end(), Y_ord.begin(), Y_ord.end());


    uint64_t tmp[5];

    // calculate result
    concordance_count(X_ord.begin(), X_ord.end(), Y_ord.begin(), Y_ord.end(),
		      tmp[0], tmp[1], tmp[2], tmp[3], tmp[4]);


    double kendal = kendall_tau(tmp[0], tmp[1], tmp[2], tmp[3]);

    return kendal;
  }

#endif
