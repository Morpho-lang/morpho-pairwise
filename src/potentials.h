/** @file potentials.h
 *  @brief Common pairwise potential classes
 */

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <morpho.h>

#define GRAVITY_CLASSNAME                    "GravityPotential"
#define COULOMB_CLASSNAME                    "CoulombPotential"
#define HERTZIAN_CLASSNAME                   "HertzianPotential"
#define LENNARDJONES_CLASSNAME               "LJPotential"

void potentials_initialize(value objclass);

#endif
