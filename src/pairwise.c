#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>
#include <morpho/value.h>
#include <morpho/veneer.h>
#include <morpho/functional.h>

#include "pairwise.h"

/* ----------------------------------------------
 * Some common pairwise potentials
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
 * Pairwise class
 * ---------------------------------------------- */

value pairwise_potentialproperty;
value pairwise_cutoffproperty;

value pairwise_valuemethod; 
value pairwise_derivativemethod; 

typedef struct {
    value potential; 
    value valuemethod; 
    value derivmethod; 
} pairwiseref;

/** Prepares the reference structure from the object's properties */
bool pairwise_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, pairwiseref *ref) {
    bool success=false;
    
    if (objectinstance_getproperty(self, pairwise_potentialproperty, &ref->potential) && 
        MORPHO_ISOBJECT(ref->potential) && 
        morpho_lookupmethod(ref->potential, pairwise_valuemethod, &ref->valuemethod) &&
        morpho_lookupmethod(ref->potential, pairwise_derivativemethod, &ref->derivmethod)) {

        success=true; 
    }

    return success;
}

/** Calculate pairwise interaction */
bool pairwise_integrand(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, double *out) {
    pairwiseref *eref = (pairwiseref *) ref; 
    double *x0, *x1, s[mesh->dim], sum = 0.0;

    matrix_getcolumn(mesh->vert, id, &x0);

    for (int j=0; j<id; j++) {
        // Compute separation
        matrix_getcolumn(mesh->vert, j, &x1);
        functional_vecsub(mesh->dim, x0, x1, s);
        double r = functional_vecnorm(mesh->dim, s);

        // Call potential function 
        value rval = MORPHO_FLOAT(r), ret;
        if (!morpho_invoke(v, eref->potential, eref->valuemethod, 1, &rval, &ret)) return false; 

        // Add to sum 
        double val; 
        if (morpho_valuetofloat(ret, &val)) sum+=val; 
    }

    *out=sum; 

    return true;
}

/** Calculate scaled gradient */
bool pairwise_gradient(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, objectmatrix *frc) {
    pairwiseref *eref = (pairwiseref *) ref; 
    double *x0, *x1, s[mesh->dim];

    matrix_getcolumn(mesh->vert, id, &x0);

    for (int j=0; j<id; j++) {
        // Compute separation
        matrix_getcolumn(mesh->vert, j, &x1);
        functional_vecsub(mesh->dim, x0, x1, s);
        double r = functional_vecnorm(mesh->dim, s);

        // Call potential derivative function 
        value rval = MORPHO_FLOAT(r), ret;
        if (!morpho_invoke(v, eref->potential, eref->derivmethod, 1, &rval, &ret)) return false; 

        // Add to sum 
        double val; 
        if (morpho_valuetofloat(ret, &val)) {
            matrix_addtocolumn(frc, id, val/r, s);
            matrix_addtocolumn(frc, j, -val/r, s);
        }
    }

    return true;
}

/** Initialize a Pairwise object */
value Pairwise_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));

    if (nargs>0 && MORPHO_ISOBJECT(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, pairwise_potentialproperty, MORPHO_GETARG(args, 0));
    } else {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    return MORPHO_NIL;
}

FUNCTIONAL_METHOD(Pairwise, integrand, MESH_GRADE_VERTEX, pairwiseref, pairwise_prepareref, functional_mapintegrand, pairwise_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

FUNCTIONAL_METHOD(Pairwise, total, MESH_GRADE_VERTEX, pairwiseref, pairwise_prepareref, functional_sumintegrand, pairwise_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

value Pairwise_gradient(vm *v, int nargs, value *args) {
    functional_mapinfo info;
    pairwiseref ref;
    value out=MORPHO_NIL;

    if (functional_validateargs(v, nargs, args, &info)) {
        if (pairwise_prepareref(MORPHO_GETINSTANCE(MORPHO_SELF(args)), info.mesh, MESH_GRADE_LINE, info.sel, &ref)) {
            info.g=MESH_GRADE_VERTEX;
            info.integrand=pairwise_integrand;
            info.grad=pairwise_gradient;
            info.ref=&ref;
            functional_mapgradient(v, &info, &out);
        } else morpho_runtimeerror(v, PAIRWISE_PRP);
    }
    if (!MORPHO_ISNIL(out)) morpho_bindobjects(v, 1, &out);
    return out;
}


MORPHO_BEGINCLASS(Pairwise)
MORPHO_METHOD(MORPHO_INITIALIZER_METHOD, Pairwise_init, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_INTEGRAND_METHOD, Pairwise_integrand, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_TOTAL_METHOD, Pairwise_total, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_GRADIENT_METHOD, Pairwise_gradient, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

void pairwise_initialize(void) { 
    pairwise_potentialproperty=builtin_internsymbolascstring(PAIRWISE_POTENTIAL_PROPERTY);
    pairwise_cutoffproperty=builtin_internsymbolascstring(PAIRWISE_CUTOFF_PROPERTY);

    pairwise_valuemethod=builtin_internsymbolascstring(PAIRWISE_VALUE_METHOD);
    pairwise_derivativemethod=builtin_internsymbolascstring(PAIRWISE_DERIVATIVE_METHOD);

    objectstring objclassname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objclassname));
    builtin_addclass(PAIRWISE_CLASSNAME, MORPHO_GETCLASSDEFINITION(Pairwise), objclass);

    builtin_addclass(COULOMB_CLASSNAME, MORPHO_GETCLASSDEFINITION(Coulomb), objclass);

    morpho_defineerror(PAIRWISE_PRP, ERROR_HALT, PAIRWISE_PRP_MSG);
}
