#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>
#include <morpho/value.h>
#include <morpho/veneer.h>
#include <morpho/functional.h>

#include <math.h>
#include "pairwise.h"

/** Calculate the difference of two vectors assuming a periodicity length `box` in each dimension */
void functional_vecsub_periodic(unsigned int n, double *a, double *b, double box, double *out) {
    for (unsigned int i=0; i<n; i++) {
        out[i]=a[i]-b[i];
        if (out[i]>0.5*box) {
            out[i] = fmod(out[i] + box/2, box) - box/2;
        }
        else if (out[i]<-0.5*box) {
            out[i] = fmod(out[i] - box/2, box) + box/2;
        }
    }
}

/* ----------------------------------------------
 * Some common pairwise potentials
 * ---------------------------------------------- */

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

value pairwise_sigmaproperty;

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

value lj_sigmaproperty;

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

/* ----------------------------------------------
 * Pairwise class
 * ---------------------------------------------- */

value pairwise_potentialproperty;
value pairwise_cutoffproperty;
value pairwise_periodicproperty;

value pairwise_valuemethod; 
value pairwise_derivativemethod; 

typedef struct {
    value potential; 
    value valuemethod; 
    value derivmethod; 
    bool cutoff; 
    double cutoffdist; 
    bool periodic; 
    double box; 
} pairwiseref;

/** Prepares the reference structure from the object's properties */
bool pairwise_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, pairwiseref *ref) {
    bool success=false;
    value cutoff; 
    value box;

    ref->cutoff=(objectinstance_getproperty(self, pairwise_cutoffproperty, &cutoff) && 
                 morpho_valuetofloat(cutoff, &ref->cutoffdist));
    ref->periodic=(objectinstance_getproperty(self, pairwise_periodicproperty, &box) && 
                 morpho_valuetofloat(box, &ref->box));

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
        if (eref->periodic) {
            functional_vecsub_periodic(mesh->dim, x0, x1, eref->box, s);
        }
        double r = functional_vecnorm(mesh->dim, s);

        if (eref->cutoff && r > eref->cutoffdist) continue; 

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
        if (eref->periodic) {
            functional_vecsub_periodic(mesh->dim, x0, x1, eref->box, s);
        }
        double r = functional_vecnorm(mesh->dim, s);

        if (eref->cutoff && r > eref->cutoffdist) continue; 

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
    int nfixed=nargs; 
    value cutoff = MORPHO_NIL;
    value box = MORPHO_NIL;

    if (builtin_options(v, nargs, args, &nfixed, 2, pairwise_cutoffproperty, &cutoff, pairwise_periodicproperty, &box) && 
        nfixed>0 && MORPHO_ISOBJECT(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, pairwise_potentialproperty, MORPHO_GETARG(args, 0));
        objectinstance_setproperty(self, pairwise_cutoffproperty, cutoff);
        objectinstance_setproperty(self, pairwise_periodicproperty, box);
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

/* ----------------------------------------------
 * SpherocylinderOverlap 
 * ---------------------------------------------- */

value pairwise_centerproperty;

typedef struct {
    value tangent; 
    objectfield *tangentfield; 
    value potential; 
    value valuemethod; 
    value derivmethod; 
    bool usesigma; 
    double sigma; 
} spherocylinderref;

/* 1D subproblem for Spherocylinder distance calculations. 
   Minimizes the polynomial p x^2 + q x + r for xlower <= x <= xupper. 
   Outputs the value of the polynomial in sep and the value of x in xout
   Returns true on success */
bool spherocylinder_distance1d(double p, double q, double r, double xlower, double xupper, double *sep, double *xout) {
    //printf("p: %g q: %g r: %g\n", p, q, r);
    double x, s; 
    
    x = -0.5*q/p;  // Solve unconstrained 1D problem 
    s = -0.25*q*q/p + r; 

    //printf("1d unconstrained: x: %g s: %g\n", x, s);

    if (x<xlower) { // If we are outside the bounds, evaluate the solution at the appropriate bound
        x = xlower; 
        s = xlower*(p*xlower + q) + r;
    } else if (x>xupper) {
        x = xupper; 
        s = xupper*(p*xupper + q) + r;
    }

    *xout = x; 
    *sep = s; 

    //printf("1d constrained: x: %g s: %g\n", x, s);

    return true;
}

/** Spherocylinder distance function. Each spherocylinder is set by a point and a vector
 * @param[in] dim - dimension of space
 * @param[in] x0  - } Location of spherocylinders
 * @param[in] x1  - }
 * @param[in] t0  - ] Tangent along spherocylinders
 * @param[in] t1  - ] 
 * @param[in] center  - true if x0 is the center of the spherocylinder, or false otherwise 
 * @param[out] dist - Shortest distance squared
 * @returns true on success, false otherwise */
bool spherocylinder_distance(unsigned int dim, double *x0, double *x1, double *t0, double *t1, bool center, double *dist) {
    double deltax[dim]; // Separation between two points
    functional_vecsub(dim, x1, x0, deltax); // x1 - x0 

    double dxdx = functional_vecdot(dim, deltax, deltax);

    double t00 = functional_vecdot(dim, t0, t0), // Dot products between tangent vectors
           t01 = functional_vecdot(dim, t0, t1),
           t11 = functional_vecdot(dim, t1, t1); 

    double dxt0 = functional_vecdot(dim, deltax, t0), // Dot products between separation and tangent vectors 
           dxt1 = functional_vecdot(dim, deltax, t1);

    double u,v,s; // Points of closest contact and closest contact
    double ulower = (center ? -1 : 0), uupper = 1; // Upper and lower bounds for u 
    double vlower = (center ? -1 : 0), vupper = 1; // Upper and lower bounds for u

    double discrim = t01*t01 - t00*t11; // Alignment discriminant

    if (fabs(discrim)>MORPHO_EPS) { 
        // Solve unconstrained quadratic problem
        u = (t01*dxt1 - t11*dxt0)/discrim;
        v = (t00*dxt1 - t01*dxt0)/discrim;
        s = u*(t00*u - t01*v) + v*(t11*v - t01*u) - 2*(dxt0*u - dxt1*v) + dxdx; 
        // printf("nondegen s: %g u: %g v: %g\n", s, u, v);

    } else { // If t0 and t1 are aligned, we must solve a degenerate problem
        double unum = (dxt0 + dxt1), udemon = t00 + t11 + 2*t01;
        u = unum/udemon; 
        v = -u; 
        s = dxdx - u*unum; 
        // printf("degen s: %g u: %g v: %g\n", s, u, v);
    }

    if (u<ulower || u>uupper || v<vlower || v>vupper) { // Solve 1D subproblems
        double uu[4] = { ulower, uupper, 0, 0 }, vv[4] = { 0, 0, vlower, vupper}, ss[4]; // Solutions to 1D subproblem
        if (spherocylinder_distance1d(t11, 2*(dxt1 - t01*ulower), (t00*ulower - 2*dxt0)*ulower + dxdx, ulower, uupper, ss, vv+0) && // Min on v at u=ulower
            spherocylinder_distance1d(t11, 2*(dxt1 - t01*uupper), (t00*uupper - 2*dxt0)*uupper + dxdx, ulower, uupper, ss+1, vv+1) && // Min on v at u=uupper
            spherocylinder_distance1d(t00, -2*(dxt0 + t01*vlower), (t11*vlower + 2*dxt1)*vlower + dxdx, vlower, vupper, ss+2, uu+2) && // Min on u at v=vlower
            spherocylinder_distance1d(t00, -2*(dxt0 + t01*vupper), (t11*vupper + 2*dxt1)*vupper + dxdx, vlower, vupper, ss+3, uu+3) // Min on u at v=vupper
        ) {
            u = uu[0]; v = vv[0]; s = ss[0]; // Find smallest separation
            for (int i=1; i<4; i++) if (ss[i]<s) { s = ss[i]; u = uu[i]; v = vv[i]; }
            //for (int i=0; i<4; i++) printf("1D subproblem: %g %g %g\n", uu[i], vv[i], ss[i]);
        } else return false; 
    }

    *dist = s;

    return true;
}

/** Prepares the reference structure from the object's properties */
bool spherocylinder_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, spherocylinderref *ref) {
    bool success=false;
    value sigma; 
    value box;

    ref->usesigma=(objectinstance_getproperty(self, pairwise_sigmaproperty, &sigma) && 
                 morpho_valuetofloat(sigma, &ref->sigma));

    ref->potential = MORPHO_NIL; 
    if (objectinstance_getproperty(self, pairwise_potentialproperty, &ref->potential) && 
        MORPHO_ISOBJECT(ref->potential) && 
        morpho_lookupmethod(ref->potential, pairwise_valuemethod, &ref->valuemethod) &&
        morpho_lookupmethod(ref->potential, pairwise_derivativemethod, &ref->derivmethod)) {
    }

    if (objectinstance_getproperty(self, functional_fieldproperty, &ref->tangent) && 
        MORPHO_ISFIELD(ref->tangent)) {
        ref->tangentfield = MORPHO_GETFIELD(ref->tangent);
        success=true; 
    } 

    return success;
}

/** Calculate pairwise interaction */
bool spherocylinder_integrand(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, double *out) {
    spherocylinderref *eref = (spherocylinderref *) ref;  
    double *x0, *x1, *t0, *t1, s[mesh->dim], sum = 0.0;
    unsigned int nel; 

    matrix_getcolumn(mesh->vert, id, &x0);
    field_getelementaslist(eref->tangentfield, MESH_GRADE_VERTEX, id, 0, &nel, &t0); 
    if (nel!=mesh->dim) return false; 

    for (int j=0; j<id; j++) {
        double r; 

        matrix_getcolumn(mesh->vert, j, &x1);
        field_getelementaslist(eref->tangentfield, MESH_GRADE_VERTEX, j, 0, &nel, &t1);

        if (!spherocylinder_distance(nel, x0, x1, t0, t1, true, &r)) return false; 

        r = sqrt(r);
        
        if (eref->usesigma && r > eref->sigma) continue; 

        // Call potential function 
        value rval = MORPHO_FLOAT(r), ret;
        if (!MORPHO_ISNIL(eref->potential)) {
            if (!morpho_invoke(v, eref->potential, eref->valuemethod, 1, &rval, &ret)) return false; 
            if (!morpho_valuetofloat(ret, &r)) return false; 
        }
        sum+=r;
    }

    *out=sum; 

    return true;
}

value SpherocylinderOverlap_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));
    value field = MORPHO_NIL; 
    value potential = MORPHO_NIL; 
    value sigma = MORPHO_NIL; 

    for (int i=0; i<nargs; i++) {
        value arg = MORPHO_GETARG(args, i);
        if (MORPHO_ISFIELD(arg)) field = arg; 
        else if (morpho_isnumber(arg)) sigma = arg;
        else if (MORPHO_ISOBJECT(arg)) potential = arg; 
        else morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    if (MORPHO_ISNIL(field)) {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    } else {
        objectinstance_setproperty(self, functional_fieldproperty, field);
        objectinstance_setproperty(self, pairwise_sigmaproperty, sigma);
        objectinstance_setproperty(self, pairwise_potentialproperty, potential);
    }

    return MORPHO_NIL;
}

FUNCTIONAL_METHOD(SpherocylinderOverlap, integrand, MESH_GRADE_VERTEX, spherocylinderref, spherocylinder_prepareref, functional_mapintegrand, spherocylinder_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

FUNCTIONAL_METHOD(SpherocylinderOverlap, total, MESH_GRADE_VERTEX, spherocylinderref, spherocylinder_prepareref, functional_sumintegrand, spherocylinder_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

MORPHO_BEGINCLASS(SpherocylinderOverlap)
MORPHO_METHOD(MORPHO_INITIALIZER_METHOD, SpherocylinderOverlap_init, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_INTEGRAND_METHOD, SpherocylinderOverlap_integrand, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_TOTAL_METHOD, SpherocylinderOverlap_total, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

void pairwise_initialize(void) { 
    pairwise_potentialproperty=builtin_internsymbolascstring(PAIRWISE_POTENTIAL_PROPERTY);
    pairwise_cutoffproperty=builtin_internsymbolascstring(PAIRWISE_CUTOFF_PROPERTY);
    pairwise_periodicproperty=builtin_internsymbolascstring(PAIRWISE_PERIODIC_PROPERTY);
    pairwise_sigmaproperty=builtin_internsymbolascstring(PAIRWISE_SIGMA_PROPERTY);
    pairwise_centerproperty=builtin_internsymbolascstring(PAIRWISE_CENTER_PROPERTY);
    lj_sigmaproperty=builtin_internsymbolascstring(LJ_SIGMA_PROPERTY);

    pairwise_valuemethod=builtin_internsymbolascstring(PAIRWISE_VALUE_METHOD);
    pairwise_derivativemethod=builtin_internsymbolascstring(PAIRWISE_DERIVATIVE_METHOD);

    objectstring objclassname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objclassname));
    builtin_addclass(PAIRWISE_CLASSNAME, MORPHO_GETCLASSDEFINITION(Pairwise), objclass);

    builtin_addclass(COULOMB_CLASSNAME, MORPHO_GETCLASSDEFINITION(Coulomb), objclass);
    builtin_addclass(HERTZIAN_CLASSNAME, MORPHO_GETCLASSDEFINITION(Hertzian), objclass);
    builtin_addclass(LENNARDJONES_CLASSNAME, MORPHO_GETCLASSDEFINITION(LennardJones), objclass);
    builtin_addclass(SPHEROCYLINDER_CLASSNAME, MORPHO_GETCLASSDEFINITION(SpherocylinderOverlap), objclass);

    morpho_defineerror(PAIRWISE_PRP, ERROR_HALT, PAIRWISE_PRP_MSG);
    morpho_defineerror(SPHEROCYLINDER_FLD, ERROR_HALT, SPHEROCYLINDER_FLD_MSG);
    morpho_defineerror(SPHEROCYLINDER_DIM, ERROR_HALT, SPHEROCYLINDER_DIM_MSG);
}
