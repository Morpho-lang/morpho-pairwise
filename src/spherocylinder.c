/** @file spherocylinder.c
 *  @brief Spherocylinder overlap functional
 */

#define MORPHO_INCLUDE_LINALG
#define MORPHO_INCLUDE_SPARSE
#define MORPHO_INCLUDE_GEOMETRY

#include <math.h>
#include <morpho.h>
#include <builtin.h>
#include <classes.h>
#include <geometry.h>
#include <functional.h>

#include "pairwise.h"
#include "spherocylinder.h"

typedef struct {
    value tangent;
    objectfield *field;
    value potential;
    value valuemethod;
    value derivmethod;
    bool usesigma;
    double sigma;
    bool center;
} spherocylinderref;

/* 1D subproblem for Spherocylinder distance calculations.
   Minimizes the polynomial p x^2 + q x + r for xlower <= x <= xupper.
   Outputs the value of the polynomial in sep and the value of x in xout
   Returns true on success */
bool spherocylinder_distance1d(double p, double q, double r, double xlower, double xupper, double *sep, double *xout) {
    double x, s;

    x = -0.5*q/p;  // Solve unconstrained 1D problem
    s = -0.25*q*q/p + r;

    if (x<xlower) { // If we are outside the bounds, evaluate the solution at the appropriate bound
        x = xlower;
        s = xlower*(p*xlower + q) + r;
    } else if (x>xupper) {
        x = xupper;
        s = xupper*(p*xupper + q) + r;
    }

    *xout = x;
    *sep = s;

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
 * @param[out] uout  - } Parameter values at which shortest distance occurs
 * @param[out] vout  - }
 * @returns true on success, false otherwise */
bool spherocylinder_distance(unsigned int dim, double *x0, double *x1, double *t0, double *t1, bool center, double *dist, double *uout, double *vout) {
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
    double vlower = (center ? -1 : 0), vupper = 1; // Upper and lower bounds for v

    double discrim = t01*t01 - t00*t11; // Alignment discriminant

    if (fabs(discrim)>MORPHO_EPS) {
        // Solve unconstrained quadratic problem
        u = (t01*dxt1 - t11*dxt0)/discrim;
        v = (t00*dxt1 - t01*dxt0)/discrim;
        s = u*(t00*u - t01*v) + v*(t11*v - t01*u) - 2*(dxt0*u - dxt1*v) + dxdx;

    } else { // If t0 and t1 are aligned, we must solve a degenerate problem
        double unum = (dxt0 + dxt1), udemon = t00 + t11 + 2*t01;
        u = unum/udemon;
        v = -u;
        s = dxdx - u*unum;
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
        } else return false;
    }

    *dist = fabs(s);
    if (uout) *uout = u;
    if (vout) *vout = v;

    return true;
}

/** Prepares the reference structure from the object's properties */
bool spherocylinder_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, spherocylinderref *ref) {
    bool success=false;
    value sigma, center, pot;

    ref->usesigma=(objectinstance_getproperty(self, pairwise_sigmaproperty, &sigma) &&
                 morpho_valuetofloat(sigma, &ref->sigma));

    ref->center = true;
    if (objectinstance_getproperty(self, pairwise_centerproperty, &center) &&
        MORPHO_ISBOOL(center)) {
            ref->center=MORPHO_GETBOOLVALUE(center);
    }

    ref->potential = MORPHO_NIL;
    ref->valuemethod = MORPHO_NIL;
    ref->derivmethod = MORPHO_NIL;
    if (objectinstance_getproperty(self, pairwise_potentialproperty, &pot) &&
        MORPHO_ISOBJECT(pot) &&
        morpho_lookupmethod(pot, pairwise_valuemethod, &ref->valuemethod) &&
        morpho_lookupmethod(pot, pairwise_derivativemethod, &ref->derivmethod)) {
        ref->potential = pot;
    }

    if (objectinstance_getproperty(self, functional_fieldproperty, &ref->tangent) &&
        MORPHO_ISFIELD(ref->tangent)) {
        ref->field = MORPHO_GETFIELD(ref->tangent);
        success=true;
    }

    return success;
}

/** Clones the spherocylinder reference with a given substitute field */
void *spherocylinder_cloneref(void *ref, objectfield *field, objectfield *sub) {
    spherocylinderref *nref = (spherocylinderref *) ref;
    spherocylinderref *clone = MORPHO_MALLOC(sizeof(spherocylinderref));

    if (clone) {
        *clone = *nref;
        if (clone->field==field) clone->field=sub;
    }

    return clone;
}

/** Calculate pairwise interaction */
bool spherocylinder_integrand(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, double *out) {
    spherocylinderref *eref = (spherocylinderref *) ref;
    double *x0, *x1, *t0, *t1, sum = 0.0;
    unsigned int nel;

    if (matrix_getcolumnptr(mesh->vert, id, &x0)!=LINALGERR_OK) return false;
    field_getelementaslist(eref->field, MESH_GRADE_VERTEX, id, 0, &nel, &t0);
    if (nel!=mesh->dim) {
        morpho_runtimeerror(v, SPHEROCYLINDER_DIM);
        return false;
    }

    for (int j=0; j<id; j++) {
        double r;

        if (matrix_getcolumnptr(mesh->vert, j, &x1)!=LINALGERR_OK) return false;
        field_getelementaslist(eref->field, MESH_GRADE_VERTEX, j, 0, &nel, &t1);

        if (!spherocylinder_distance(nel, x0, x1, t0, t1, eref->center, &r, NULL, NULL)) return false;

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

/** Calculate gradient */
bool spherocylinder_gradient(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, objectmatrix *frc) {
    spherocylinderref *eref = (spherocylinderref *) ref;
    double *x0, *x1, *t0, *t1, s[mesh->dim];
    unsigned int nel;

    if (matrix_getcolumnptr(mesh->vert, id, &x0)!=LINALGERR_OK) return false;
    field_getelementaslist(eref->field, MESH_GRADE_VERTEX, id, 0, &nel, &t0);
    if (nel!=mesh->dim) {
        morpho_runtimeerror(v, SPHEROCYLINDER_DIM);
        return false;
    }

    for (int j=0; j<id; j++) {
        double rsq, r, uu, vv, dv=1.0;

        if (matrix_getcolumnptr(mesh->vert, j, &x1)!=LINALGERR_OK) return false;
        field_getelementaslist(eref->field, MESH_GRADE_VERTEX, j, 0, &nel, &t1);

        if (!spherocylinder_distance(nel, x0, x1, t0, t1, eref->center, &rsq, &uu, &vv)) return false;

        r = sqrt(rsq);

        if (eref->usesigma && r > eref->sigma) continue;
        if (fabs(r)<MORPHO_EPS) continue;

        // Call potential derivative function
        value rval = MORPHO_FLOAT(r), ret;
        if (!MORPHO_ISNIL(eref->potential)) {
            if (!morpho_invoke(v, eref->potential, eref->derivmethod, 1, &rval, &ret)) return false;
            if (!morpho_valuetofloat(ret, &dv)) return false;
        }

        // Grad_x0 s^2 = 2*((x0-x1) + u*t0 - v*t1 ) / (2 r)
        functional_vecsub(mesh->dim, x0, x1, s);
        functional_vecaddscale(mesh->dim, s, uu, t0, s);
        functional_vecaddscale(mesh->dim, s, -vv, t1, s);
        if (matrix_addtocolumnptr(frc, id, dv/r, s)!=LINALGERR_OK) return false;

        // Grad_x1 s^2 = -Grad_x0 s^2
        if (matrix_addtocolumnptr(frc, j, -dv/r, s)!=LINALGERR_OK) return false;
    }

    return true;
}

/** Calculate field gradient */
bool spherocylinder_fieldgradient(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, objectfield *frc) {
    spherocylinderref *eref = (spherocylinderref *) ref;
    double *x0, *x1, *t0, *t1, *ft0, *ft1, s[mesh->dim];
    unsigned int nel, fnel;

    if (matrix_getcolumnptr(mesh->vert, id, &x0)!=LINALGERR_OK) return false;
    field_getelementaslist(eref->field, MESH_GRADE_VERTEX, id, 0, &nel, &t0);
    field_getelementaslist(frc, MESH_GRADE_VERTEX, id, 0, &fnel, &ft0);
    if (nel!=mesh->dim || fnel!=mesh->dim) {
        morpho_runtimeerror(v, SPHEROCYLINDER_DIM);
        return false;
    }

    for (int j=0; j<id; j++) {
        double rsq, r, uu, vv, dv=1.0;

        if (matrix_getcolumnptr(mesh->vert, j, &x1)!=LINALGERR_OK) return false;
        field_getelementaslist(eref->field, MESH_GRADE_VERTEX, j, 0, &nel, &t1);

        if (!spherocylinder_distance(nel, x0, x1, t0, t1, eref->center, &rsq, &uu, &vv)) return false;

        r = sqrt(rsq);

        if (eref->usesigma && r > eref->sigma) continue;
        if (fabs(r)<MORPHO_EPS) continue;

        // Call potential derivative function
        value rval = MORPHO_FLOAT(r), ret;
        if (!MORPHO_ISNIL(eref->potential)) {
            if (!morpho_invoke(v, eref->potential, eref->derivmethod, 1, &rval, &ret)) return false;
            if (!morpho_valuetofloat(ret, &dv)) return false;
        }

        // Grad_t0 s^2 = 2*u*((x0-x1) + u*t0 - v*t1 ) / (2 r)
        functional_vecsub(mesh->dim, x0, x1, s);
        functional_vecaddscale(mesh->dim, s, uu, t0, s);
        functional_vecaddscale(mesh->dim, s, -vv, t1, s);

        functional_vecaddscale(mesh->dim, ft0, dv*uu/r, s, ft0);

        // Grad_t1 s^2 = - 2*v*((x0-x1) + u*t0 - v*t1 ) / (2 r)
        field_getelementaslist(frc, MESH_GRADE_VERTEX, j, 0, &fnel, &ft1);
        functional_vecaddscale(mesh->dim, ft1, -dv*vv/r, s, ft1);
    }

    return true;
}

value SpherocylinderOverlap_init(vm *v, int nargs, value *args) {
    int nfixed;
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));
    value field = MORPHO_NIL;
    value potential = MORPHO_NIL;
    value sigma = MORPHO_NIL;
    value center = MORPHO_TRUE;

    if (builtin_options(v, nargs, args, &nfixed, 1, pairwise_centerproperty, &center) &&
        MORPHO_ISBOOL(center)) {
        objectinstance_setproperty(self, pairwise_centerproperty, center);
    } else {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    for (int i=0; i<nfixed; i++) {
        value arg = MORPHO_GETARG(args, i);
        if (MORPHO_ISFIELD(arg)) field = arg;
        else if (morpho_isnumber(arg)) sigma = arg;
        else if (MORPHO_ISOBJECT(arg)) potential = arg;
        else morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    if (MORPHO_ISNIL(field)) {
        morpho_runtimeerror(v, SPHEROCYLINDER_FLD);
    } else {
        objectinstance_setproperty(self, functional_fieldproperty, field);
        objectinstance_setproperty(self, pairwise_sigmaproperty, sigma);
        objectinstance_setproperty(self, pairwise_potentialproperty, potential);
    }

    return MORPHO_NIL;
}

static bool spherocylinder_mapfieldgradient(vm *v, functional_mapinfo *info, value *out) {
    info->fieldgrad = spherocylinder_fieldgradient;
    return functional_mapfieldgradient(v, info, out);
}

FUNCTIONAL_MD_REF_BIND(SpherocylinderOverlap, spherocylinderref, spherocylinder_prepareref, spherocylinder_integrand, SPHEROCYLINDER_FLD)
FUNCTIONAL_MD_REF_INTEGRAND(SpherocylinderOverlap, spherocylinderref, MESH_GRADE_VERTEX)
FUNCTIONAL_MD_REF_TOTAL(SpherocylinderOverlap, spherocylinderref, MESH_GRADE_VERTEX)
FUNCTIONAL_MD_REF_GRADIENT(SpherocylinderOverlap, spherocylinderref, MESH_GRADE_VERTEX, spherocylinder_gradient, SYMMETRY_NONE)
FUNCTIONAL_MD_REF_FIELDGRADIENT_MAP(SpherocylinderOverlap, spherocylinderref, MESH_GRADE_VERTEX, spherocylinder_mapfieldgradient, spherocylinder_cloneref, NULL)

MORPHO_BEGINCLASS(SpherocylinderOverlap)
MORPHO_METHOD_SIGNATURE(MORPHO_INITIALIZER_METHOD, "(...)", SpherocylinderOverlap_init, MORPHO_FN_MUTATES|MORPHO_FN_OPTARGS),

FUNCTIONAL_MD_INTEGRAND_METHODS(SpherocylinderOverlap),
FUNCTIONAL_MD_TOTAL_METHODS(SpherocylinderOverlap),
FUNCTIONAL_MD_GRADIENT_METHODS(SpherocylinderOverlap),
FUNCTIONAL_MD_FIELDGRADIENT_METHODS(SpherocylinderOverlap)
MORPHO_ENDCLASS

void spherocylinder_initialize(value objclass) {
    builtin_addclass(SPHEROCYLINDER_CLASSNAME, MORPHO_GETCLASSDEFINITION(SpherocylinderOverlap), objclass);

    morpho_defineerror(SPHEROCYLINDER_FLD, ERROR_HALT, SPHEROCYLINDER_FLD_MSG);
    morpho_defineerror(SPHEROCYLINDER_DIM, ERROR_HALT, SPHEROCYLINDER_DIM_MSG);
}
