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
 * Hudson Gravity Attempt
 * ---------------------------------------------- */

value Gravity_value(vm *v, int nargs, value *args) { 
    value out = MORPHO_NIL; 
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(-1/(r*r));
        } 
    }
    return out; 
}

value Gravity_deriv(vm *v, int nargs, value *args) { 
    value out = MORPHO_NIL; 
    if (nargs==1) {
        double r;
        if (morpho_valuetofloat(MORPHO_GETARG(args, 0), &r)) {
            out = MORPHO_FLOAT(2/(r*r*r));
        } 
    }
    return out; 
}

MORPHO_BEGINCLASS(Gravity)
MORPHO_METHOD(PAIRWISE_VALUE_METHOD, Gravity_value, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(PAIRWISE_DERIVATIVE_METHOD, Gravity_deriv, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

/*
    END HUDSON'S GRAVITY ATTEMPT
*/

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

value lj_sigmaproperty;

value LennardJones_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));

    if (nargs>0 && MORPHO_ISNUMBER(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, lj_sigmaproperty, MORPHO_GETARG(args, 0));
    } else {
        printf("Error here!\n");
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
        if (!matrix_getcolumn(mesh->vert, vid[i], &x0)) return false; 
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
    double w0=1.0, w1=1.0;

    // Extract x0 
    if (nv==1) {
        matrix_getcolumn(mesh->vert, id, &x0);
    } else { // Compute average position from vertices
        pairwise_averagevertexposition(mesh, nv, vid, x0mean);
        x0 = x0mean; 
        //for (int i=0; i<mesh->dim; i++) printf("%g ", x0mean[i]);
        //printf("\n");
        if (!eref->conn) UNREACHABLE("Connectivity matrix not available in Pairwise_integrand");
        functional_elementsize(v, mesh, eref->g, id, nv, vid, &w0);
    }

    for (int j=0; j<id; j++) {
        // Extract x1
        if (nv==1) {
            matrix_getcolumn(mesh->vert, j, &x1);
        } else {
            int nvj, *vidj;
            //printf("conn: %p\n", (void *) eref->conn);
            if (!sparseccs_getrowindices(&eref->conn->ccs, j, &nvj, &vidj)) return false; 
            //printf("nvj: %i\n", nvj);
            pairwise_averagevertexposition(mesh, nvj, vidj, x1mean);
            x1 = x1mean; 
            //printf("x1: %i\n", j);
            //for (int i=0; i<mesh->dim; i++) printf("%g ", x1mean[i]);
            //printf("\n");
            functional_elementsize(v, mesh, eref->g, j, nvj, vidj, &w1);
            //printf("Element sizes [%g %g]\n", w0, w1);
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

        //printf("Value returned by potential function ");
        //morpho_printvalue(ret);
        //printf("\n");

        w0=1.0; w1=1.0;  // Disable including the area 

        double val;
        if (morpho_valuetofloat(ret, &val)) {
            sum+=w0*w1*val; 
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
    double w0=1.0, w1=1.0;
    int nvj, *vidj;

    //printf("GRADIENT\n");

    // Extract x0 
    if (nv==1) {
        matrix_getcolumn(mesh->vert, id, &x0);
    } else { // Compute average position from vertices
        pairwise_averagevertexposition(mesh, nv, vid, x0mean);
        x0 = x0mean; 
        //for (int i=0; i<mesh->dim; i++) printf("%g ", x0mean[i]);
        //printf("\n");
        if (!eref->conn) UNREACHABLE("Connectivity matrix not available in Pairwise_integrand");
        functional_elementsize(v, mesh, eref->g, id, nv, vid, &w0);
    }

    for (int j=0; j<id; j++) {
        // Extract x1
        if (nv==1) {
            matrix_getcolumn(mesh->vert, j, &x1);
        } else {
            
            //printf("conn: %p\n", (void *) eref->conn);
            if (!sparseccs_getrowindices(&eref->conn->ccs, j, &nvj, &vidj)) return false; 
            //printf("nvj: %i\n", nvj);
            pairwise_averagevertexposition(mesh, nvj, vidj, x1mean);
            x1 = x1mean; 
            //printf("x1: %i\n", j);
            //for (int i=0; i<mesh->dim; i++) printf("%g ", x1mean[i]);
            //printf("\n");
            functional_elementsize(v, mesh, eref->g, j, nvj, vidj, &w1);
            //printf("Element sizes [%g %g]\n", w0, w1);
        }

        // Compute separation
        functional_vecsub(mesh->dim, x0, x1, s);
        if (eref->periodic) {
            functional_vecsub_periodic(mesh->dim, x0, x1, eref->box, s);
        }
        double r = functional_vecnorm(mesh->dim, s);

        if (eref->cutoff && r > eref->cutoffdist) continue; 

        // Call potential derivative function 
        value rval = MORPHO_FLOAT(r), ret;
        if (!morpho_invoke(v, eref->potential, eref->derivmethod, 1, &rval, &ret)) return false; 

        // printf("Value returned by potential derivative function ");
        // morpho_printvalue(ret);
        // printf("\n");

        w0=1.0; w1=1.0;  // Disable including the area 

        // Add to sum 
        double val; 
        if (morpho_valuetofloat(ret, &val)) {
            if (nv==1) {
                matrix_addtocolumn(frc, id, w0*w1*val/r, s);
                matrix_addtocolumn(frc, j, -w0*w1*val/r, s);
            } else {
                double nnv = (double) nv; 
                for (int i=0; i<nv; i++) matrix_addtocolumn(frc, vid[i], w0*w1*val/r/nnv, s);
                for (int i=0; i<nvj; i++) matrix_addtocolumn(frc, vidj[i], -w0*w1*val/r/nnv, s);
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

FUNCTIONAL_METHOD(Pairwise, integrand, ref.g, pairwiseref, pairwise_prepareref, functional_mapintegrand, pairwise_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

FUNCTIONAL_METHOD(Pairwise, total, ref.g, pairwiseref, pairwise_prepareref, functional_sumintegrand, pairwise_integrand, NULL, PAIRWISE_PRP, SYMMETRY_NONE)

value Pairwise_gradient(vm *v, int nargs, value *args) {
    functional_mapinfo info;
    pairwiseref ref;
    value out=MORPHO_NIL;

    if (functional_validateargs(v, nargs, args, &info)) {
        if (pairwise_prepareref(MORPHO_GETINSTANCE(MORPHO_SELF(args)), info.mesh, MESH_GRADE_LINE, info.sel, &ref)) {
            info.g=ref.g;
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
    pairwise_periodicproperty=builtin_internsymbolascstring(PAIRWISE_PERIODIC_PROPERTY);
    pairwise_sigmaproperty=builtin_internsymbolascstring(PAIRWISE_SIGMA_PROPERTY);
    lj_sigmaproperty=builtin_internsymbolascstring(LJ_SIGMA_PROPERTY);

    pairwise_valuemethod=builtin_internsymbolascstring(PAIRWISE_VALUE_METHOD);
    pairwise_derivativemethod=builtin_internsymbolascstring(PAIRWISE_DERIVATIVE_METHOD);

    objectstring objclassname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objclassname));
    builtin_addclass(PAIRWISE_CLASSNAME, MORPHO_GETCLASSDEFINITION(Pairwise), objclass);

    builtin_addclass(COULOMB_CLASSNAME, MORPHO_GETCLASSDEFINITION(Coulomb), objclass);
    builtin_addclass(GRAVITY_CLASSNAME, MORPHO_GETCLASSDEFINITION(Gravity), objclass);
    builtin_addclass(HERTZIAN_CLASSNAME, MORPHO_GETCLASSDEFINITION(Hertzian), objclass);
    builtin_addclass(LENNARDJONES_CLASSNAME, MORPHO_GETCLASSDEFINITION(LennardJones), objclass);

    morpho_defineerror(PAIRWISE_PRP, ERROR_HALT, PAIRWISE_PRP_MSG);
}
