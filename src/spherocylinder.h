/** @file spherocylinder.h
 *  @brief Spherocylinder overlap functional
 */

#ifndef SPHEROCYLINDER_H
#define SPHEROCYLINDER_H

#include <morpho.h>

#define SPHEROCYLINDER_CLASSNAME             "SpherocylinderOverlap"

#define SPHEROCYLINDER_FLD                   "SphrCylFld"
#define SPHEROCYLINDER_FLD_MSG               "Spherocylinder requires a Field to define orientation."

#define SPHEROCYLINDER_DIM                   "SphrCylDim"
#define SPHEROCYLINDER_DIM_MSG               "Dimension of vector field must match dimension of space."

void spherocylinder_initialize(value objclass);

#endif
