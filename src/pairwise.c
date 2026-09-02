/** @file pairwise.c
 *  @brief Pairwise mesh functional
 */

#define MORPHO_INCLUDE_LINALG
#define MORPHO_INCLUDE_SPARSE
#define MORPHO_INCLUDE_GEOMETRY

#include <stdio.h>
#include <math.h>
#include <morpho.h>
#include <builtin.h>
#include <classes.h>
#include <geometry.h>
#include <functional.h>

#include "pairwise.h"
#include "potentials.h"
#include "spherocylinder.h"

value pairwise_potentialproperty;
value pairwise_cutoffproperty;
value pairwise_periodicproperty;
value pairwise_sigmaproperty;
value pairwise_centerproperty;
value pairwise_valuemethod;
value pairwise_derivativemethod;

/** Calculate the difference of two vectors assuming a periodicity length `box` in each dimension */
static void functional_vecsub_periodic(unsigned int n, double *a, double *b, double box, double *out) {
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

typedef struct {
    value potential;
    value valuemethod;
    value derivmethod;
    bool cutoff;
    double cutoffdist;
    bool periodic;
    double box;
    grade g;
    objectsparse *conn; // Connectivity matrix for grade g
} pairwiseref;

/** Prepares the reference structure from the object's properties */
bool pairwise_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, pairwiseref *ref) {
    bool success=false;
    value cutoff;
    value box;
    value gradeprop;

    ref->cutoff=(objectinstance_getproperty(self, pairwise_cutoffproperty, &cutoff) &&
                 morpho_valuetofloat(cutoff, &ref->cutoffdist));
    ref->periodic=(objectinstance_getproperty(self, pairwise_periodicproperty, &box) &&
                 morpho_valuetofloat(box, &ref->box));

    ref->g=MESH_GRADE_VERTEX;
    ref->conn=NULL;
    if (objectinstance_getproperty(self, functional_gradeproperty, &gradeprop)) {
        morpho_valuetoint(gradeprop, &ref->g);
    }
    if (ref->g > MESH_GRADE_VERTEX) {
        ref->conn = mesh_getconnectivityelement(mesh, 0, ref->g);
    }

    if (objectinstance_getproperty(self, pairwise_potentialproperty, &ref->potential) &&
        MORPHO_ISOBJECT(ref->potential) &&
        morpho_lookupmethod(ref->potential, pairwise_valuemethod, &ref->valuemethod) &&
        morpho_lookupmethod(ref->potential, pairwise_derivativemethod, &ref->derivmethod)) {

        success=true;
    }

    return success;
}

/** Compute the average vertex position
 * @param[in] mesh - the mesh object
 * @param[in] nv   - number of vertices
 * @param[in] vid  - vertex ids
 * @param[out] xmean - the mean vertex position
 * @returns true on success
 */
bool pairwise_averagevertexposition(objectmesh *mesh, int nv, int *vid, double *xmean) {
    double *x0;
    for (int i=0; i<mesh->dim; i++) xmean[i]=0.0;
    for (int i=0; i<nv; i++) {
        if (matrix_getcolumnptr(mesh->vert, vid[i], &x0)!=LINALGERR_OK) return false;
        functional_vecadd(mesh->dim, x0, xmean, xmean);
    }
    functional_vecscale(mesh->dim, 1.0/nv, xmean, xmean);
    return true;
}

/** Calculate pairwise interaction */
bool pairwise_integrand(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, double *out) {
    pairwiseref *eref = (pairwiseref *) ref;
    double *x0, *x1, s[mesh->dim], sum = 0.0;
    double x0mean[mesh->dim], x1mean[mesh->dim];

    // Extract x0
    if (nv==1) {
        if (matrix_getcolumnptr(mesh->vert, id, &x0)!=LINALGERR_OK) return false;
    } else { // Compute average position from vertices
        if (!pairwise_averagevertexposition(mesh, nv, vid, x0mean)) return false;
        x0 = x0mean;
        if (!eref->conn) UNREACHABLE("Connectivity matrix not available in Pairwise_integrand");
    }

    for (int j=0; j<id; j++) {
        // Extract x1
        if (nv==1) {
            if (matrix_getcolumnptr(mesh->vert, j, &x1)!=LINALGERR_OK) return false;
        } else {
            int nvj, *vidj;
            if (!sparseccs_getrowindices(&eref->conn->ccs, j, &nvj, &vidj)) return false;
            if (!pairwise_averagevertexposition(mesh, nvj, vidj, x1mean)) return false;
            x1 = x1mean;
        }

        // Compute separation
        functional_vecsub(mesh->dim, x0, x1, s);
        if (eref->periodic) {
            functional_vecsub_periodic(mesh->dim, x0, x1, eref->box, s);
        }
        double r = functional_vecnorm(mesh->dim, s);

        if (eref->cutoff && r > eref->cutoffdist) continue;

        // Call potential function
        value rval = MORPHO_FLOAT(r), ret;
        if (!morpho_invoke(v, eref->potential, eref->valuemethod, 1, &rval, &ret)) return false;

        double val;
        if (morpho_valuetofloat(ret, &val)) {
            sum+=val;
        } else return false;
    }

    *out=sum;

    return true;
}

/** Calculate scaled gradient */
bool pairwise_gradient(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, objectmatrix *frc) {
    pairwiseref *eref = (pairwiseref *) ref;
    double *x0, *x1, s[mesh->dim];
    double x0mean[mesh->dim], x1mean[mesh->dim];
    int nvj, *vidj;

    // Extract x0
    if (nv==1) {
        if (matrix_getcolumnptr(mesh->vert, id, &x0)!=LINALGERR_OK) return false;
    } else { // Compute average position from vertices
        if (!pairwise_averagevertexposition(mesh, nv, vid, x0mean)) return false;
        x0 = x0mean;
        if (!eref->conn) UNREACHABLE("Connectivity matrix not available in Pairwise_integrand");
    }

    for (int j=0; j<id; j++) {
        // Extract x1
        if (nv==1) {
            if (matrix_getcolumnptr(mesh->vert, j, &x1)!=LINALGERR_OK) return false;
        } else {
            if (!sparseccs_getrowindices(&eref->conn->ccs, j, &nvj, &vidj)) return false;
            if (!pairwise_averagevertexposition(mesh, nvj, vidj, x1mean)) return false;
            x1 = x1mean;
        }

        // Compute separation
        functional_vecsub(mesh->dim, x0, x1, s);
        if (eref->periodic) {
            functional_vecsub_periodic(mesh->dim, x0, x1, eref->box, s);
        }
        double r = functional_vecnorm(mesh->dim, s);

        if (eref->cutoff && r > eref->cutoffdist) continue;
        if (fabs(r)<MORPHO_EPS) continue;

        // Call potential derivative function
        value rval = MORPHO_FLOAT(r), ret;
        if (!morpho_invoke(v, eref->potential, eref->derivmethod, 1, &rval, &ret)) return false;

        // Add to sum
        double val;
        if (morpho_valuetofloat(ret, &val)) {
            if (nv==1) {
                if (matrix_addtocolumnptr(frc, id, val/r, s)!=LINALGERR_OK) return false;
                if (matrix_addtocolumnptr(frc, j, -val/r, s)!=LINALGERR_OK) return false;
            } else {
                double nnv = (double) nv;
                for (int i=0; i<nv; i++) {
                    if (matrix_addtocolumnptr(frc, vid[i], val/r/nnv, s)!=LINALGERR_OK) return false;
                }
                for (int i=0; i<nvj; i++) {
                    if (matrix_addtocolumnptr(frc, vidj[i], -val/r/nnv, s)!=LINALGERR_OK) return false;
                }
            }
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
    value grade = MORPHO_INTEGER(0);

    if (builtin_options(v, nargs, args, &nfixed, 3, functional_gradeproperty, &grade, pairwise_cutoffproperty, &cutoff, pairwise_periodicproperty, &box) &&
        nfixed>0 && MORPHO_ISOBJECT(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, pairwise_potentialproperty, MORPHO_GETARG(args, 0));
        objectinstance_setproperty(self, functional_gradeproperty, grade);
        objectinstance_setproperty(self, pairwise_cutoffproperty, cutoff);
        objectinstance_setproperty(self, pairwise_periodicproperty, box);
    } else {
        morpho_runtimeerror(v, PAIRWISE_PRP);
    }

    return MORPHO_NIL;
}

FUNCTIONAL_MD_REF_BIND(Pairwise, pairwiseref, pairwise_prepareref, pairwise_integrand, PAIRWISE_PRP)
FUNCTIONAL_MD_REF_INTEGRAND(Pairwise, pairwiseref, ref.g)
FUNCTIONAL_MD_REF_TOTAL(Pairwise, pairwiseref, ref.g)
FUNCTIONAL_MD_REF_GRADIENT(Pairwise, pairwiseref, ref.g, pairwise_gradient, SYMMETRY_NONE)

MORPHO_BEGINCLASS(Pairwise)
MORPHO_METHOD_SIGNATURE(MORPHO_INITIALIZER_METHOD, "(...)", Pairwise_init, MORPHO_FN_MUTATES|MORPHO_FN_OPTARGS),

FUNCTIONAL_MD_INTEGRAND_METHODS(Pairwise),
FUNCTIONAL_MD_TOTAL_METHODS(Pairwise),
FUNCTIONAL_MD_GRADIENT_METHODS(Pairwise)
MORPHO_ENDCLASS

void pairwise_initialize(void) {
    pairwise_potentialproperty=builtin_internsymbolascstring(PAIRWISE_POTENTIAL_PROPERTY);
    pairwise_cutoffproperty=builtin_internsymbolascstring(PAIRWISE_CUTOFF_PROPERTY);
    pairwise_periodicproperty=builtin_internsymbolascstring(PAIRWISE_PERIODIC_PROPERTY);
    pairwise_sigmaproperty=builtin_internsymbolascstring(PAIRWISE_SIGMA_PROPERTY);
    pairwise_centerproperty=builtin_internsymbolascstring(PAIRWISE_CENTER_PROPERTY);

    pairwise_valuemethod=builtin_internsymbolascstring(PAIRWISE_VALUE_METHOD);
    pairwise_derivativemethod=builtin_internsymbolascstring(PAIRWISE_DERIVATIVE_METHOD);

    objectstring objclassname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objclassname));
    builtin_addclass(PAIRWISE_CLASSNAME, MORPHO_GETCLASSDEFINITION(Pairwise), objclass);

    potentials_initialize(objclass);
    spherocylinder_initialize(objclass);

    morpho_defineerror(PAIRWISE_PRP, ERROR_HALT, PAIRWISE_PRP_MSG);
}
