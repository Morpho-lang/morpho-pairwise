#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>

#define PAIRWISE_CLASSNAME                   "Pairwise"

#define PAIRWISE_POTENTIAL_PROPERTY          "potential"
#define PAIRWISE_CUTOFF_PROPERTY             "cutoff"

#define COULOMB_CLASSNAME                    "CoulombPotential"

#define PAIRWISE_VALUE_METHOD                "value"
#define PAIRWISE_DERIVATIVE_METHOD           "derivative"

#define PAIRWISE_PRP                         "PrwsPrp"
#define PAIRWISE_PRP_MSG                     "Pairwise properties."

void pairwise_initialize(void);