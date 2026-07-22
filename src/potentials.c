/** @file potentials.c
 *  @brief Common pairwise potential classes
 */

#define MORPHO_INCLUDE_LINALG

#include <stdio.h>
#include <math.h>
#include <morpho.h>
#include <builtin.h>
#include <classes.h>

#include "pairwise.h"
#include "potentials.h"

/* ----------------------------------------------
 * Gravity potential
 * ---------------------------------------------- */

value Gravity_value(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(-1/r);
        }
    }
    return out;
}

value Gravity_deriv(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(2/(r*r));
        }
    }
    return out;
}

MORPHO_BEGINCLASS(Gravity)
MORPHO_METHOD(PAIRWISE_VALUE_METHOD, Gravity_value, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_DERIVATIVE_METHOD, Gravity_deriv, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

/* ----------------------------------------------
 * Coulomb potential
 * ---------------------------------------------- */

value Coulomb_value(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(1/r);
        }
    }
    return out;
}

value Coulomb_deriv(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(-1/(r*r));
        }
    }
    return out;
}

MORPHO_BEGINCLASS(Coulomb)
MORPHO_METHOD(PAIRWISE_VALUE_METHOD, Coulomb_value, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_DERIVATIVE_METHOD, Coulomb_deriv, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

/* ----------------------------------------------
 * Hertzian potential
 * ---------------------------------------------- */

value Hertzian_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));

    if (nargs>0 && MORPHO_ISNUMBER(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, pairwise_sigmaproperty, MORPHO_GETARG(args, 0));
    } else {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    return MORPHO_NIL;
}

bool hertzian_getsigma(value obj, double *sigma) {
    value val;
    return (objectinstance_getproperty(MORPHO_GETINSTANCE(obj), pairwise_sigmaproperty, &val) &&
            morpho_valuetofloat(val, sigma));
}

value Hertzian_value(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double sigma, r;
        if (hertzian_getsigma(MORPHO_SELF(args), &sigma) &&
            morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            if (r<sigma) {
                double u = 1-r/sigma;
                out = MORPHO_FLOAT(pow(u, 2.5));
            } else {
                out = MORPHO_FLOAT(0.0);
            }
        }
    }
    return out;
}

value Hertzian_deriv(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double sigma, r;
        if (hertzian_getsigma(MORPHO_SELF(args), &sigma) &&
            morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)){
            if (r<sigma) {
                double u = 1-r/sigma;
                out = MORPHO_FLOAT(-2.5*pow(u, 1.5)/sigma);
            } else {
                out = MORPHO_FLOAT(0.0);
            }
        }
    }
    return out;
}

MORPHO_BEGINCLASS(Hertzian)
MORPHO_METHOD(MORPHO_INITIALIZER_METHOD, Hertzian_init, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_VALUE_METHOD, Hertzian_value, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_DERIVATIVE_METHOD, Hertzian_deriv, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

/* ----------------------------------------------
 * Lennard Jones
 * ---------------------------------------------- */

static value lj_sigmaproperty;

value LennardJones_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));

    if (nargs>0 && MORPHO_ISNUMBER(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, lj_sigmaproperty, MORPHO_GETARG(args, 0));
    } else {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    return MORPHO_NIL;
}

bool lennardjones_getsigma(value obj, double *sigma) {
    value val;
    return (objectinstance_getproperty(MORPHO_GETINSTANCE(obj), lj_sigmaproperty, &val) &&
            morpho_valuetofloat(val, sigma));
}

value LennardJones_value(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double sigma, r, val;
        if (lennardjones_getsigma(MORPHO_SELF(args), &sigma) &&
            morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            val = 4*(pow((sigma/r), 12) - pow((sigma/r), 6));
            out = MORPHO_FLOAT(val) ;
        }
    }
    return out;
}

value LennardJones_deriv(vm *v, int nargs, value *args) {
    value out = MORPHO_NIL;
    if (nargs==1) {
        double sigma, r, val;
        if (lennardjones_getsigma(MORPHO_SELF(args), &sigma) &&
            morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)){
            val = -24 * (2 * pow((sigma/r),12) - pow((sigma/r),6)) / r;
            out = MORPHO_FLOAT(val);
        }
    }
    return out;
}

MORPHO_BEGINCLASS(LennardJones)
MORPHO_METHOD(MORPHO_INITIALIZER_METHOD, LennardJones_init, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_VALUE_METHOD, LennardJones_value, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_DERIVATIVE_METHOD, LennardJones_deriv, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

void potentials_initialize(value objclass) {
    lj_sigmaproperty=builtin_internsymbolascstring(LJ_SIGMA_PROPERTY);

    builtin_addclass(COULOMB_CLASSNAME, MORPHO_GETCLASSDEFINITION(Coulomb), objclass);
    builtin_addclass(GRAVITY_CLASSNAME, MORPHO_GETCLASSDEFINITION(Gravity), objclass);
    builtin_addclass(HERTZIAN_CLASSNAME, MORPHO_GETCLASSDEFINITION(Hertzian), objclass);
    builtin_addclass(LENNARDJONES_CLASSNAME, MORPHO_GETCLASSDEFINITION(LennardJones), objclass);
}
