/** @file pairwise.h
 *  @brief Shared definitions for the morpho-pairwise package
 */

#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <morpho.h>

#define PAIRWISE_CLASSNAME                   "Pairwise"

#define PAIRWISE_POTENTIAL_PROPERTY          "potential"
#define PAIRWISE_CUTOFF_PROPERTY             "cutoff"
#define PAIRWISE_PERIODIC_PROPERTY           "box"
#define PAIRWISE_SIGMA_PROPERTY              "sigma"
#define PAIRWISE_CENTER_PROPERTY             "center"

#define PAIRWISE_VALUE_METHOD                "value"
#define PAIRWISE_DERIVATIVE_METHOD           "derivative"

#define PAIRWISE_PRP                         "PrwsPrp"
#define PAIRWISE_PRP_MSG                     "Pairwise properties."

/* Shared interned symbols used across the package */
extern value pairwise_potentialproperty;
extern value pairwise_cutoffproperty;
extern value pairwise_periodicproperty;
extern value pairwise_sigmaproperty;
extern value pairwise_centerproperty;
extern value pairwise_valuemethod;
extern value pairwise_derivativemethod;

void pairwise_initialize(void);

#endif
