#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>

#define PAIRWISE_CLASSNAME                   "Pairwise"

#define PAIRWISE_POTENTIAL_PROPERTY          "potential"
#define PAIRWISE_CUTOFF_PROPERTY             "cutoff"
#define PAIRWISE_PERIODIC_PROPERTY           "box"
#define PAIRWISE_SIGMA_PROPERTY              "sigma"
#define PAIRWISE_CENTER_PROPERTY             "center"
#define LJ_SIGMA_PROPERTY                    "sigma"

#define GRAVITY_CLASSNAME                    "GravityPotential"
#define COULOMB_CLASSNAME                    "CoulombPotential"
#define HERTZIAN_CLASSNAME                   "HertzianPotential"
#define LENNARDJONES_CLASSNAME               "LJPotential"

#define SPHEROCYLINDER_CLASSNAME             "SpherocylinderOverlap"

#define PAIRWISE_VALUE_METHOD                "value"
#define PAIRWISE_DERIVATIVE_METHOD           "derivative"

#define PAIRWISE_PRP                         "PrwsPrp"
#define PAIRWISE_PRP_MSG                     "Pairwise properties."

#define SPHEROCYLINDER_FLD                   "SphrCylFld"
#define SPHEROCYLINDER_FLD_MSG               "Spherocylinder requires a Field to define orientation."

#define SPHEROCYLINDER_DIM                   "SphrCylDim"
#define SPHEROCYLINDER_DIM_MSG               "Dimension of vector field must match dimension of space."

void functional_vecsub_periodic(unsigned int n, double *a, double *b, double box, double *out);

void pairwise_initialize(void);
