[comment]: # (Pairwise help)

# pairwise
[pairwise]: # (pairwise)

The `pairwise` extension provides fast pairwise potentials and functionals for morpho.

Import it with:

    import pairwise

It supplies:

* potential classes with `value(r)` and `derivative(r)` methods
* a `Pairwise` functional for interactions on a mesh
* a `SpherocylinderOverlap` functional for oriented rod-like particles

[showsubtopics]: # (subtopics)

## GravityPotential
[GravityPotential]: # (GravityPotential)

Gravitational potential of the form `-1/r`.

    var g = GravityPotential()
    print g.value(2) // expect: -0.5
    print g.derivative(2) // expect: 0.25

## CoulombPotential
[CoulombPotential]: # (CoulombPotential)

Coulomb potential of the form `1/r`.

    var c = CoulombPotential()
    print c.value(2) // expect: 0.5
    print c.derivative(2) // expect: -0.25

## HertzianPotential
[HertzianPotential]: # (HertzianPotential)

Soft-sphere Hertzian repulsion with range `sigma`.

For `r < sigma`,

    ((1 - r/sigma))^2.5

and zero otherwise. Construct with a positive `sigma`:

    var h = HertzianPotential(1)
    print h.value(0.5)
    print h.derivative(0.5)

The range can also be updated through the `sigma` property:

    h.sigma = 0.2

## LJPotential
[LJPotential]: # (LJPotential)

Lennard-Jones potential with length scale `sigma`:

    4*((sigma/r)^12 - (sigma/r)^6)

Construct with a positive `sigma`:

    var lj = LJPotential(1)
    print lj.value(1.3)
    print lj.derivative(1.3)

## Pairwise
[Pairwise]: # (Pairwise)

A mesh functional that sums a pairwise potential over pairs of mesh elements.

Construct with a potential object and optional keyword arguments:

    var lp = Pairwise(CoulombPotential())
    var lp = Pairwise(CoulombPotential(), cutoff=2)
    var lp = Pairwise(CoulombPotential(), box=1)
    var lp = Pairwise(CoulombPotential(), grade=2)

Arguments:

* potential — an object providing `value(r)` and `derivative(r)`
* `cutoff` — optional maximum separation; pairs farther than this are skipped
* `box` — optional periodic box side length
* `grade` — optional mesh grade to pair over (default: vertices)

Typical use:

    import pairwise
    import meshtools

    var mb = MeshBuilder()
    mb.addvertex([0,0,0])
    mb.addvertex([1,0,0])
    var m = mb.build()

    var lp = Pairwise(CoulombPotential(), cutoff=2)
    print lp.total(m)
    print lp.integrand(m)
    print lp.gradient(m)

Methods:

* `total(mesh)` — total pairwise energy
* `integrand(mesh)` — per-element contributions
* `gradient(mesh)` — gradient with respect to vertex positions

The `cutoff` property can be updated after construction:

    lp.cutoff = 1.5

## SpherocylinderOverlap
[SpherocylinderOverlap]: # (SpherocylinderOverlap)

A functional for interactions between spherocylinders (capsules).

Each particle is defined by a mesh vertex position and an orientation vector stored in a `Field`. By default the vertex is treated as the center of the spherocylinder (`center=true`).

Construct with a required orientation `Field`, and optionally a potential, a cutoff `sigma`, and the `center` flag:

    var f = Field(m, Matrix([1,0,0]))
    var sc = SpherocylinderOverlap(f)
    var sc = SpherocylinderOverlap(f, HertzianPotential(0.2))
    var sc = SpherocylinderOverlap(f, 0.2)
    var sc = SpherocylinderOverlap(f, HertzianPotential(0.2), center=true)

If no potential is supplied, the functional uses the shortest separation distance itself.

Methods:

* `total(mesh)` — total interaction energy
* `integrand(mesh)` — per-vertex contributions
* `gradient(mesh)` — gradient with respect to vertex positions
* `fieldgradient(mesh)` — gradient with respect to the orientation field

Example:

    import pairwise
    import meshtools

    var mb = MeshBuilder()
    mb.addvertex([-2,0,0])
    mb.addvertex([2,0,0])
    var m = mb.build()

    var f = Field(m, Matrix([1,0,0]))
    var sc = SpherocylinderOverlap(f)

    print sc.total(m)
    print sc.gradient(m)
    print sc.fieldgradient(m)
