/*
 * R random generator for scrm, using R's C API (no Rcpp).
 * Based on scrm's r_random_generator.h (GPL-3).
 */

#ifndef scrmr_r_random_generator
#define scrmr_r_random_generator

#include <R.h>
#include <Rmath.h>
#include <memory>
#include "scrm/random/random_generator.h"

class RRandomGenerator : public RandomGenerator {
 public:
  RRandomGenerator(std::shared_ptr<FastFunc> ff):RandomGenerator(ff) {
    this->initializeUnitExponential();
  };
  virtual ~RRandomGenerator() {};

  void initialize() {};

  double sample() {
    return unif_rand();
  }

  double sampleUnitExponential() {
    return exp_rand();
  }
};

#endif
