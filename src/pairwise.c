#include <stdio.h>
#include <morpho/morpho.h>
#include <morpho/builtin.h>
#include <morpho/selection.h>
#include <morpho/field.h>
#include <morpho/functional.h>
#include <morpho/integrate.h>

#include "selingernematic.h"

/* ----------------------------------------------
 * Selinger's Nematic Energy
 * ---------------------------------------------- */

bool gradsq_evaluategradient(objectmesh *mesh, objectfield *field, int nv, int *vid, double *out);
bool gradsq_evaluategradient3d(objectmesh *mesh, objectfield *field, int nv, int *vid, double *out);

value selingernematic_k11property;
value selingernematic_k22property;
value selingernematic_k33property;
value selingernematic_k24property;

typedef struct {
    double k11,k22,k33,k24;
    objectfield *field;
    grade grade;
} selingernematicref;

/** Prepares the nematic reference */
bool selingernematic_prepareref(objectinstance *self, objectmesh *mesh, grade g, objectselection *sel, selingernematicref *ref) {
    bool success=false, grdset=false;
    value field=MORPHO_NIL, grd=MORPHO_NIL;
    value val=MORPHO_NIL;
    ref->k11=1.0; ref->k22=1.0; ref->k33=1.0; ref->k24=0.0;

    if (objectinstance_getproperty(self, functional_fieldproperty, &field) &&
        MORPHO_ISFIELD(field)) {
        ref->field=MORPHO_GETFIELD(field);
        success=true;
    }

    if (objectinstance_getproperty(self, selingernematic_k11property, &val) && MORPHO_ISNUMBER(val)) {
        morpho_valuetofloat(val, &ref->k11);
    }
    if (objectinstance_getproperty(self, selingernematic_k22property, &val) && MORPHO_ISNUMBER(val)) {
        morpho_valuetofloat(val, &ref->k22);
    }
    if (objectinstance_getproperty(self, selingernematic_k33property, &val) && MORPHO_ISNUMBER(val)) {
        morpho_valuetofloat(val, &ref->k33);
    }
    if (objectinstance_getproperty(self, selingernematic_k24property, &val) && MORPHO_ISNUMBER(val)) {
        morpho_valuetofloat(val, &ref->k24);
    }

    /*if (objectinstance_getproperty(self, nematic_pitchproperty, &val) && MORPHO_ISNUMBER(val)) {
        morpho_valuetofloat(val, &ref->pitch);
        ref->haspitch=true;
    }*/

    if (objectinstance_getproperty(self, functional_gradeproperty, &grd) &&
        MORPHO_ISINTEGER(grd)) {
        ref->grade=MORPHO_GETINTEGERVALUE(grd);
        if (ref->grade>0) grdset=true;
    }
    if (!grdset) ref->grade=mesh_maxgrade(mesh);

    return success;
}

/** Clones the nematic reference with a given substitute field */
void *selingernematic_cloneref(void *ref, objectfield *field, objectfield *sub) {
    selingernematicref *nref = (selingernematicref *) ref;
    selingernematicref *clone = MORPHO_MALLOC(sizeof(selingernematicref));
    
    if (clone) {
        *clone = *nref;
        if (clone->field==field) clone->field=sub;
    }
    
    return clone;
}

/* Integrates a linear vector function with values at vertices f[0]...f[n]
   Works for dimensions 1-3 at least
   Needs to be multiplied by the volume of the element to get the true value */
double selingernematic_bcintf(unsigned int n, double *f) {
    double sum = 0;
    for (unsigned int i=0; i<n; i++) sum+=f[i];
    return sum/n;
}

/* Integrates a product of two linear functions with values at vertices
   f[0]...f[n] and g[0]...g[n].
   Works for dimensions 1-3 at least
   Needs to be multiplied by the volume of the element to get the true value */
double selingernematic_bcintfg(unsigned int n, double *f, double *g) {
    double sum = 0;
    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) sum+=f[i]*g[j];
        sum+=f[i]*g[i];
    }
    return sum/(n*(n+1));
}

typedef struct {
    objectmatrix *dn;
    double divnn; 
    double *curlnn;
} selingerdeltaref; 

/* Calculates m + outer(u,v) -> m */
void _outer(objectmatrix *m, double scale, double *u, double *v) {
    for (unsigned int i=0; i<m->ncols; i++) matrix_addtocolumn(m, i, scale*v[i], u);
}

/* Integrand function necessary for delta term */
bool deltaintegrand(unsigned int dim, double *t, double *x, unsigned int nquantity, value *quantity, void *data, double *out) {
    selingerdeltaref *ref = (selingerdeltaref *) data;

    // Extract director 
    objectmatrix *director = MORPHO_GETMATRIX(quantity[0]);
    double nn[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<director->nrows*director->ncols; i++) nn[i]=director->elements[i];

    // Create delta matrix 
    double deltast[9], ast[9];
    objectmatrix delta = MORPHO_STATICMATRIX(deltast, 3, 3);
    objectmatrix A = MORPHO_STATICMATRIX(ast, 3, 3); // Temporary matrix 

    matrix_transpose(ref->dn, &A);
    matrix_add(ref->dn, &A, &delta);            // delta = dn + dn^T + other terms 

    double bend[3]; // Bend vector
    functional_veccross(nn, ref->curlnn, bend);

    //printf("bnd->[%g %g %g]\n", bend[0], bend[1], bend[2]);

    _outer(&delta, 1.0, nn, bend);  
    _outer(&delta, 1.0, bend, nn);               // delta += n_i B_j + n_j B_i 
    
    matrix_identity(&A);
    _outer(&A, -1.0, nn, nn); // kronecker_ij - ninj

    matrix_accumulate(&delta, -ref->divnn, &A); // delta += S*(kronecker_ij - ni*nj)

    /*printf("xx->[%g %g]\n", x[0], x[1]);
    matrix_print(&delta);
    printf("\n");
    exit(0);*/

    matrix_mul(&delta, &delta, &A);

    matrix_trace(&A, out);
    *out *= 0.25; // Account for factor of 1/2 

    return true; 
}


/** Calculate the nematic energy */
bool selingernematic_integrand(vm *v, objectmesh *mesh, elementid id, int nv, int *vid, void *ref, double *out) {
    selingernematicref *eref = ref;
    double size=0; // Length area or volume of the element
    double gradnnraw[eref->field->psize*3];
    double gradnn[eref->field->psize*3];
    double divnn, curlnn[3] = { 0.0, 0.0, 0.0 };
    
    for (int i=0; i<eref->field->psize*3; i++) { gradnn[i]=0.0; gradnnraw[i]=0.0; }

    if (!functional_elementsize(v, mesh, eref->grade, id, nv, vid, &size)) return false;

    // Get nematic director components
    double *nn[nv]; // Field value lists
    unsigned int nentries=0;
    for (unsigned int i=0; i<nv; i++) {
        if (!field_getelementaslist(eref->field, MESH_GRADE_VERTEX, vid[i], 0, &nentries, &nn[i])) return false;
    }

    // Evaluate gradients of the director
    if (eref->grade==2) {
        if (!gradsq_evaluategradient(mesh, eref->field, nv, vid, gradnnraw)) return
            false;
    } else if (eref->grade==3) {
        if (!gradsq_evaluategradient3d(mesh, eref->field, nv, vid, gradnnraw)) return
            false;
    }
    
    // Copy into 3x3 matrix
    for (int j=0; j<3; j++) for (int i=0; i<mesh->dim; i++) gradnn[3*j+i] = gradnnraw[mesh->dim*j+i];
    
    // Output of this is the matrix:
    // [ nx,x ny,x nz,x ] [ 0 3 6 ] <- indices
    // [ nx,y ny,y nz,y ] [ 1 4 7 ]
    // [ nx,z ny,z nz,z ] [ 2 5 8 ]
    objectmatrix gradnnmat = MORPHO_STATICMATRIX(gradnn, 3, 3);

    matrix_trace(&gradnnmat, &divnn); // Compute divnn
    curlnn[0]=gradnn[7]-gradnn[5]; // nz,y - ny,z
    curlnn[1]=gradnn[2]-gradnn[6]; // nx,z - nz,x
    curlnn[2]=gradnn[3]-gradnn[1]; // ny,x - nx,y

    //matrix_print(&gradnnmat);
    //printf("\ndivn: %g curln: (%g %g %g)\n", divnn, curlnn[0], curlnn[1], curlnn[2]);

    // Calculate integrals of nx^2, ny^2, nz^2, nx*ny, ny*nz, and nz*nx over the element 
    double nnt[3][nv]; // The transpose of nn
    for (unsigned int i=0; i<nv; i++)
        for (unsigned int j=0; j<3; j++) nnt[j][i]=nn[i][j];

    double integrals[] = {  selingernematic_bcintfg(nv, nnt[0], nnt[0]),
                            selingernematic_bcintfg(nv, nnt[1], nnt[1]),
                            selingernematic_bcintfg(nv, nnt[2], nnt[2]),
                            selingernematic_bcintfg(nv, nnt[0], nnt[1]),
                            selingernematic_bcintfg(nv, nnt[1], nnt[2]),
                            selingernematic_bcintfg(nv, nnt[2], nnt[0])
    };  

    double splay = divnn; // Double splay is div nn 
    double splayen = 0.5*(eref->k11-eref->k24)*size*splay*splay; 

    // Twist pseudoscalar is n . curl(n) and we want to compute integral T^2
    // Multiply out these terms and project onto (nx^2, ny^2, nz^2, nx*ny, ny*nz, nz*nx) polynomial basis
    double ctwst[6] = { curlnn[0]*curlnn[0], curlnn[1]*curlnn[1], curlnn[2]*curlnn[2],
                        2*curlnn[0]*curlnn[1], 2*curlnn[1]*curlnn[2], 2*curlnn[2]*curlnn[0]};

    double twisten = 0.0; 
    for (unsigned int i=0; i<6; i++) twisten+=ctwst[i]*integrals[i];
    twisten *= 0.5*(eref->k22-eref->k24)*size;

    // Bend vector is n x curl(n)
    // As for bend, multiply out and project onto polynomial basis, reusing values from ctwst
    double cbnd[6] = { ctwst[1] + ctwst[2], ctwst[0] + ctwst[2], ctwst[0] + ctwst[1],
                       -ctwst[3], -ctwst[4], -ctwst[5] };

    double benden = 0.0; 
    for (unsigned int i=0; i<6; i++) benden+=cbnd[i]*integrals[i];
    benden *= 0.5*(eref->k33)*size;

    // Compute saddle splay energy via integrate
    selingerdeltaref sref; 
    double ssen;

    sref.curlnn=curlnn;
    sref.divnn=divnn;
    sref.dn=&gradnnmat;

    double *x[3];
    value q0[1], q1[1], q2[1], q3[1];
    value *q[4] = { q0, q1, q2, q3 };

    for (unsigned int i=0; i<nv; i++) {
        mesh_getvertexcoordinatesaslist(mesh, vid[i], &x[i]);
        field_getelement(eref->field, MESH_GRADE_VERTEX, vid[i], 0, &q[i][0]);
    }

    bool success=integrate_integrate(deltaintegrand, mesh->dim, eref->grade, x, 1, q, &sref, &ssen);
    if (success) ssen *=eref->k24*size;

    //printf("%u, [%g %g %g %g]\n", id, splayen, twisten, benden, ssen);

    *out = splayen + twisten + benden + ssen;

    return true;
}

/** Initialize a Nematic object */
value SelingerNematic_init(vm *v, int nargs, value *args) {
    objectinstance *self = MORPHO_GETINSTANCE(MORPHO_SELF(args));

    int nfixed=nargs;
    value k11=MORPHO_FLOAT(2.0),
          k22=MORPHO_FLOAT(2.0),
          k33=MORPHO_FLOAT(2.0),
          k24=MORPHO_FLOAT(1.0);

    if (builtin_options(v, nargs, args, &nfixed, 4,
                        selingernematic_k11property, &k11,
                        selingernematic_k22property, &k22,
                        selingernematic_k33property, &k33,
                        selingernematic_k24property, &k24)) {
        objectinstance_setproperty(self, selingernematic_k11property, k11);
        objectinstance_setproperty(self, selingernematic_k22property, k22);
        objectinstance_setproperty(self, selingernematic_k33property, k33);
        objectinstance_setproperty(self, selingernematic_k24property, k24);
    } else morpho_runtimeerror(v, SELINGERNEMATIC_ARGS);

    if (nfixed==1 && MORPHO_ISFIELD(MORPHO_GETARG(args, 0))) {
        objectinstance_setproperty(self, functional_fieldproperty, MORPHO_GETARG(args, 0));
    } else morpho_runtimeerror(v, SELINGERNEMATIC_ARGS);

    return MORPHO_NIL;
}

FUNCTIONAL_METHOD(SelingerNematic, integrand, (ref.grade), selingernematicref, selingernematic_prepareref, functional_mapintegrand, selingernematic_integrand, NULL, NEMATIC_ARGS, SYMMETRY_NONE);

FUNCTIONAL_METHOD(SelingerNematic, total, (ref.grade), selingernematicref, selingernematic_prepareref, functional_sumintegrand, selingernematic_integrand, NULL, NEMATIC_ARGS, SYMMETRY_NONE);

FUNCTIONAL_METHOD(SelingerNematic, gradient, (ref.grade), selingernematicref, selingernematic_prepareref, functional_mapnumericalgradient, selingernematic_integrand, NULL, NEMATIC_ARGS, SYMMETRY_NONE);

value SelingerNematic_fieldgradient(vm *v, int nargs, value *args) {
    functional_mapinfo info;
    selingernematicref ref;
    value out=MORPHO_NIL;

    if (functional_validateargs(v, nargs, args, &info)) {
        if (selingernematic_prepareref(MORPHO_GETINSTANCE(MORPHO_SELF(args)), info.mesh, MESH_GRADE_AREA, info.sel, &ref)) {
            info.g=ref.grade;
            info.integrand=selingernematic_integrand;
            info.ref=&ref;
            info.cloneref=selingernematic_cloneref;
            functional_mapnumericalfieldgradient(v, &info, &out);
        } else morpho_runtimeerror(v, SELINGERNEMATIC_ARGS);
    }
    if (!MORPHO_ISNIL(out)) morpho_bindobjects(v, 1, &out);
    return out;
}

MORPHO_BEGINCLASS(SelingerNematic)
MORPHO_METHOD(MORPHO_INITIALIZER_METHOD, SelingerNematic_init, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_INTEGRAND_METHOD, SelingerNematic_integrand, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_TOTAL_METHOD, SelingerNematic_total, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_GRADIENT_METHOD, SelingerNematic_gradient, BUILTIN_FLAGSEMPTY),
MORPHO_METHOD(FUNCTIONAL_FIELDGRADIENT_METHOD, SelingerNematic_fieldgradient, BUILTIN_FLAGSEMPTY)
MORPHO_ENDCLASS

/*
bool testintegrand(unsigned int dim, double *t, double *x, unsigned int nquantity, value *quantity, void *data, double *out) {
    *out = MORPHO_GETFLOATVALUE(quantity[0]);
    return true; 
}

void sel_test(void) {
    double x0[3] = { 0,0,0 };
    double x1[3] = { 1,0,0 };
    double x2[3] = { 0,1,0 };
    double x3[3] = { 0,0,1 };
    double *xx[4] = { x0, x1, x2, x3 };

    value v0[1] = { MORPHO_FLOAT(1.0) };
    value v1[1] = { MORPHO_FLOAT(2.0) };
    value v2[1] = { MORPHO_FLOAT(3.0) };
    value v3[1] = { MORPHO_FLOAT(4.0) };
    value *v[4] = { v0, v1, v2, v3 };

    double out;
    integrate_integrate(testintegrand, 3, 2, xx, 1, v, NULL, &out);
    printf("integral value: %g\n", out);
}*/

void selingernematic_initialize(void) { 
    selingernematic_k11property=builtin_internsymbolascstring(SELINGERNEMATIC_K11_PROPERTY);
    selingernematic_k22property=builtin_internsymbolascstring(SELINGERNEMATIC_K22_PROPERTY);
    selingernematic_k33property=builtin_internsymbolascstring(SELINGERNEMATIC_K33_PROPERTY);
    selingernematic_k24property=builtin_internsymbolascstring(SELINGERNEMATIC_K24_PROPERTY);

    objectstring objclassname = MORPHO_STATICSTRING(OBJECT_CLASSNAME);
    value objclass = builtin_findclass(MORPHO_OBJECT(&objclassname));
    builtin_addclass(SELINGERNEMATIC_CLASSNAME, MORPHO_GETCLASSDEFINITION(SelingerNematic), objclass);

    morpho_defineerror(SELINGERNEMATIC_ARGS, ERROR_HALT, SELINGERNEMATIC_ARGS_MSG);
}
