Following the successful discussion that we had a few year's ago on
conventions for storing unstructured grid data in netCDF files which
eventually led to the [UGRID
conventions](https://github.com/ugrid-conventions/ugrid-conventions), I
would like to propose a simple convention for documenting staggered data
on structured grids that is consistent with the UGRID conventions. My
proposal, which I will refer to as SGRID convention, is described below.

The link <http://bit.ly/sgrid_cf> points to this page.

Introduction
============

The [CF-conventions](http://cfconventions.org/) are widely used for
storing and distributing environmental / earth sciences / climate data.
The CF-conventions use a data perspective: every data value points to
the latitude and longitude at which that value has been defined; the
combination of latitude and longitude bounds and cell methods attributes
can be used to define spatially averaged rather than point values. This
is all great for the distribution of (interpolated) data for general
visualization and spatial data processing, but it doesn't capture the
relationship of the variables as computed by a numerical model (such as
[Arakawa staggering](http://en.wikipedia.org/wiki/Arakawa_grids)). Many
models use staggered grids (using finite differences, or finite volume
approach) or use a finite element approach of which the correct meaning
may not be captured easily by simple cell methods descriptors. This
becomes a problem if you don't want to just look at the big picture of
the model results, but also at the details at the grid resolution e.g.
what is the exact meaning of a flux on the output file in discrete
terms? Can we verify the mass balance? Can the data be used for
restarting the model? Correctly handling the staggered data has always
been a crucial element of the Delft3D post-processing tools. In the
UGRID conventions, we have defined the (unstructured) grid as a separate
entity on the file which consists of nodes and connections of nodes
defining edges, faces, and volumes. For a structured (staggered) grid we
are currently lacking a consistent convention. Although one could store
structured grid data using UGRID conventions, some fundamental aspects
such as distinction between grid directions would be lost. In this
context I propose the lightweight SGRID conventions to define the core
aspects of a structured staggered grid without trying to capture the
details of finite element formulations. This proposal serves merely the
purpose of getting the conventions for structured grids on par with
those for unstructured grids.

Proposal
========

Consistent with the UGRID conventions we use the following terms for
points, lines, and cells that make up a grid.

Dimensionality

Name

Comments

0

node

A point, a coordinate pair or triplet: the most basic element of the
topology (also known as "vertex").

1

edge

A line or curve bounded by two nodes.

2

face

A plane or surface enclosed by a set of edges. In a Cartesian 2D model,
you might refer to this as a "cell" or "square".

3

volume

A volume enclosed by a set of faces.

In the UGRID conventions the focus is on describing the topology of the
mesh (connectivity of the nodes, edges, faces, and volumes as
appropriate). The topology of a structured grid is inherently defined;
the focus for the SGRID conventions is therefore on the numbering used.
Still we need to distinguish between 2D and 3D grids (1D conventions may
be defined consistently).

### 2D grid

Required topology attributes

Value

cf\_role

mesh\_topology

topology\_dimension

2

node\_dimensions

node\_dimension1 node\_dimension2

face\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*)
face\_dimension2:node\_dimension2 (padding:*type2*)

Optional attributes

Default value

edge1\_dimensions

node\_dimension1 face\_dimension2:node\_dimension2 (padding:*type2*)

edge2\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*) node\_dimension2

node\_coordinates

edge1\_coordinates

edge2\_coordinates

face\_coordinates

vertical\_dimensions

where the padding type may be one of the four literal strings: "none",
"low", "high", or "both" depending on whether the face\_dimension is one
shorter than the corresponding node\_dimension (padding:none), one
longer than the corresponding node\_dimension (padding:both), or of
equal length with one extra value stored on the low or high end of the
dimension (see the figure below). The edge1\_dimensions and
edge2\_dimensions attributes may be used to define separate dimensions
for the edges (see the ROMS example below), but by default the edge
dimensions are assumed to be consistent with the dimensions used by the
edges and faces respectively. The optional vertical\_dimensions
attribute may be used to specify the names of the dimensions for the
layers and layer interfaces respectively using the same syntax:
layer\_dimension:layer\_interface\_dimension (padding:*type*).

Note:

-   The numbering of the edges corresponds to the order of the
    dimensions in the dimensions attributes. The order of the dimensions
    listed here does not change since these are specified as a string.
    The actual order of dimensions in the netCDF API (and hence in the
    dump generated by different tools) depends on the programming
    language used. Hence, the order of the dimensions listed here may
    differ from the order of the dimensions in the data definition.

Figure 1: illustrating the formulations used for expressing the
relationship between face/edge and node dimensions. Please note that the
numbering of the faces and nodes can be adjusted using face(face) and
node(node) coordinate variables.

Example:

dimensions: time = UNLIMITED ; inode = 10 ; jnode = 20 ; icell = 9 ;
jcell = 19 ; variables: char time(time) ; time:standard\_name = "time" ;
time:long\_name = "time" ; time:units = "seconds since 2015-01-01
00:00:00" ; float u(time, jcell, inode) ; u:description = "x-velocity" ;
u:units = "m s-1" ; u:grid = "MyGrid" ; u:location = "edge1" ; float
v(time, jnode, icell) ; u:description = "y-velocity" ; u:units = "m s-1"
; u:grid = "MyGrid" ; u:location = "edge2" ; float c(time, jcell, icell)
; c:description = "some concentration" ; c:grid = "MyGrid" ; c:location
= "face" ; float node\_lat(jnode, inode) ; node\_lat:standard\_name =
"latitude" ; node\_lat:units = "degree\_north" ; float node\_lon(jnode,
inode) ; node\_lon:standard\_name = "longitude" ; node\_lon:units =
"degree\_east" ; int MyGrid ; grid:cf\_role = grid\_topology
grid:topology\_dimension = 2 ; grid:node\_dimensions = "inode jnode" ;
grid:face\_dimensions = "icell: inode (padding: none) jcell: jnode
(padding: none)" ; grid:node\_coordinates = "node\_lon node\_lat" ;

### 3d grid

Required topology attributes

Value

cf\_role

mesh\_topology

topology\_dimension

3

node\_dimensions

node\_dimension1 node\_dimension2 node\_dimension3

volume\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*)
face\_dimension2:node\_dimension2 (padding:*type2*)
face\_dimension3:node\_dimension3 (padding:*type3*)

Optional attributes

Default value

edge1\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*) node\_dimension2
node\_dimension3

edge2\_dimensions

node\_dimension1 face\_dimension2:node\_dimension2 (padding:*type2*)
node\_dimension3

edge3\_dimensions

node\_dimension1 node\_dimension2 face\_dimension3:node\_dimension3
(padding:*type3*)

face1\_dimensions

node\_dimension1 face\_dimension2:node\_dimension2 (padding:*type2*)
face\_dimension3:node\_dimension3 (padding:*type3*)

face2\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*) node\_dimension2
face\_dimension3:node\_dimension3 (padding:*type3*)

face3\_dimensions

face\_dimension1:node\_dimension1 (padding:*type1*)
face\_dimension2:node\_dimension2 (padding:*type2*) node\_dimension3

node\_coordinates

edge<span style="color: rgb(255,0,0);"> *\<i\>* </span>\_coordinates

face *<span style="color: rgb(255,0,0);">\<i\></span>* \_coordinates

volume\_coordinates

Notes:

-   The edge1, edge2, and edge3 are in a 3D grid aligned to the
    dimensions1, 2, and 3 respectively, whereas the edge1 and edge2 are
    in a 2D grid perpendicular to the dimensions 1 and 2 respectively.
    The face1, face2, and face3 play that role in the 3D grid.
-   The 3d grid option should not be used be used for layered grids,
    such as typical ocean and atmosphere models. Use the 2d grid with
    vertical dimensions instead. This allows 2D quantities (such as
    surface quantities) and 3D quantities to be linked to the same mesh.

Example:

dimensions: time = UNLIMITED ; inode = 10 ; jnode = 20 ; knode = 30 ;
iface = 9 ; jface = 19 ; kface = 29 ; variables: char time(time) ;
time:standard\_name = "time" ; time:long\_name = "time" ; time:units =
"seconds since 2015-01-01 00:00:00" ; float u(time, kface, jface, inode)
; u:description = "x-velocity" ; u:units = "m s-1" ; u:grid = "MyGrid3"
; u:location = "face1" ; float v(time, kface, jnode, iface) ;
u:description = "y-velocity" ; u:units = "m s-1" ; u:grid = "MyGrid3" ;
u:location = "face2" ; float w(time, knode, jface, iface) ;
u:description = "z-velocity" ; u:units = "m s-1" ; u:grid = "MyGrid3" ;
u:location = "face3" ; float c(time, kface, jface, iface) ;
c:description = "some concentration" ; c:grid = "MyGrid3" ; c:location =
"volume" ; float node\_lat(knode, jnode, inode) ;
node\_lat:standard\_name = "latitude" ; node\_lat:units =
"degree\_north" ; float node\_lon(knode, jnode, inode) ;
node\_lon:standard\_name = "longitude" ; node\_lon:units =
"degree\_east" ; float node\_elevation(knode, jnode, inode) ;
node\_elevation:description = "elevation" ; node\_elevation:units = "m"
; int MyGrid3 ; grid:cf\_role = grid\_topology grid:topology\_dimension
= 3 ; grid:node\_dimensions = "inode jnode knode" ;
grid:volume\_dimensions = "iface: inode (padding: none) jface: jnode
(padding: none) kface: knode (padding: none)" ; grid:node\_coordinates =
"node\_lon node\_lat node\_elevation" ;

Data variables
==============

The use of the attributes to associate a data variable with a specific
grid and stagger location is copied from the UGRID conventions: To map
the variable onto the topology of the underlying grid, two new
attributes have been introduced. First, the attribute `grid` points to
the grid\_topology variable containing the meta-data attributes of the
grid on which the variable has been defined. Second, the attribute
`location` points to the (stagger) location within the grid at which the
variable is defined.

Example:

double waterlevel(time,j,i) ; waterlevel:standard\_name =
"sea\_surface\_height\_above\_geoid" ; waterlevel:units = "m" ;
waterlevel:grid = "MyGrid" waterlevel:location = "face" ;
waterlevel:coordinates = "lat\_face\_MyGrid lon\_face\_MyGrid" ;

Examples
========

Delft3D
-------

Delft3D uses an Arikawa C-grid with the water level (pressure) computed
in the cell centres, and the normal velocities at the cell edges. This
example shows the use of asymmetric padding (at the low end of the
horizontal coordinate indices there is an extra line of face/mid-point
values). In the vertical there is no padding, so the number of layer
interfaces is one more than the number of layers. The integer coodinate
variables KMAX and KMAX1 are used to indicate that layer interfaces are
numbered 0 to KMAX whereas all other indices use the default numbering
from 1 to the maximum value.

netcdf trim-f34 { dimensions: NMAX = 22 ; NMAXZ = 22 ; MMAX = 15 ; MMAXZ
= 15 ; KMAX = 5 ; KMAX1 = 6 ; time = UNLIMITED ; // (6 currently)
variables: int KMAX(KMAX) ; int KMAX1(KMAX1) ; float XCOR(MMAX, NMAX) ;
XCOR:standard\_name = "projection\_x\_coordinate" ; XCOR:long\_name =
"X-coordinate of grid points" ; XCOR:units = "m" ; float YCOR(MMAX,
NMAX) ; YCOR:standard\_name = "projection\_y\_coordinate" ;
YCOR:long\_name = "Y-coordinate of grid points" ; YCOR:units = "m" ;
float XZ(MMAXZ, NMAXZ) ; XZ:standard\_name = "projection\_x\_coordinate"
; XZ:long\_name = "X-coordinate of cell centres" ; XZ:units = "m" ;
float YZ(MMAXZ, NMAXZ) ; YZ:standard\_name = "projection\_y\_coordinate"
; YZ:long\_name = "Y-coordinate of cell centres" ; YZ:units = "m" ;
float THICK(KMAX) ; THICK:long\_name = "Fraction part of layer thickness
of total water-height" ; THICK:units = "[ .01\*% ]" ; float time(time) ;
time:standard\_name = "time" ; time:long\_name = "time" ; time:units =
"seconds since 1990-08-05 00:00:00" ; float S1(time, MMAXZ, NMAXZ) ;
S1:long\_name = "Water-level in zeta point" ; S1:units = "m" ;
S1:coordinates = "XZ YZ" ; S1:grid = "grid" ; // SGRID attribute
S1:location = "face" ; // SGRID attribute float U1(time, KMAX, MMAX,
NMAXZ) ; U1:long\_name = "U-velocity per layer in U-point (Eulerian)" ;
U1:units = "m/s" ; U1:grid = "grid" ; // SGRID attribute U1:location =
"edge1" ; // SGRID attribute float V1(time, KMAX, MMAXZ, NMAX) ;
V1:long\_name = "V-velocity per layer in V-point (Eulerian)" ; V1:units
= "m/s" ; V1:grid = "grid" ; // SGRID attribute V1:location = "edge2" ;
// SGRID attribute float W(time, KMAX1, MMAXZ, NMAXZ) ; W:long\_name =
"W-omega per layer in zeta point" ; W:units = "m/s" ; W:grid = "grid" ;
// SGRID attribute W:location = "face" ; // SGRID attribute // SGRID
variable int grid ; grid:cf\_role = grid\_topology
grid:topology\_dimension = 2 ; grid:node\_dimensions = "MMAX NMAX" ;
grid:face\_dimensions = "MMAXZ: MMAX (padding: low) NMAXZ: NMAX
(padding: low)" ; grid:node\_coordinates = "XCOR YCOR" ;
grid:face\_coordinates = "XZ YZ" ; grid:vertical\_dimensions = "KMAX:
KMAX1 (padding: none)" ; // global attributes: :title = "Het Friesche
Zeegaatje" ; data: KMAX = 1, 2, 3, 4, 5 ; KMAX1 = 0, 1, 2, 3, 4, 5 ; }

The edge\_dimension attributes are not needed.

ROMS
----

ROMS uses also a C-grid, but it uses on the output file different
dimensions for each staggered location. In this case we need all
attributes defined above including the edge<span
style="color: rgb(255,0,0);"> *\<i\>* </span>\_dimension attributes..

netcdf sed023\_last { dimensions: ocean\_time = UNLIMITED ; // (1
currently) s\_w = 21 ; eta\_rho = 60 ; xi\_rho = 160 ; tracer = 10 ;
s\_rho = 20 ; eta\_u = 60 ; xi\_u = 159 ; eta\_v = 59 ; xi\_v = 160 ;
eta\_psi = 59 ; xi\_psi = 159 ; variables: double lat\_psi(eta\_psi,
xi\_psi) ; lat\_psi:long\_name = "latitude of PSI-points" ;
lat\_psi:units = "degree\_north" ; double lat\_rho(eta\_rho, xi\_rho) ;
lat\_rho:long\_name = "latitude of RHO-points" ; lat\_rho:units =
"degree\_north" ; double lat\_u(eta\_u, xi\_u) ; lat\_u:long\_name =
"latitude of U-points" ; lat\_u:units = "degree\_north" ; double
lat\_v(eta\_v, xi\_v) ; lat\_v:long\_name = "latitude of V-points" ;
lat\_v:units = "degree\_north" ; double lon\_psi(eta\_psi, xi\_psi) ;
lon\_psi:long\_name = "longitude of PSI-points" ; lon\_psi:units =
"degree\_east" ; double lon\_rho(eta\_rho, xi\_rho) ;
lon\_rho:long\_name = "longitude of RHO-points" ; lon\_rho:units =
"degree\_east" ; double lon\_u(eta\_u, xi\_u) ; lon\_u:long\_name =
"longitude of U-points" ; lon\_u:units = "degree\_east" ; double
lon\_v(eta\_v, xi\_v) ; lon\_v:long\_name = "longitude of V-points" ;
lon\_v:units = "degree\_east" ; double ocean\_time(ocean\_time) ;
ocean\_time:long\_name = "time since initialization" ; ocean\_time:units
= "seconds since 1968-05-23 00:00:00 GMT" ; ocean\_time:calendar =
"gregorian" ; double s\_rho(s\_rho) ; s\_rho:long\_name = "S-coordinate
at RHO-points" ; s\_rho:valid\_min = -1. ; s\_rho:valid\_max = 0. ;
s\_rho:standard\_name = "ocean\_s\_coordinate" ; s\_rho:formula\_terms =
"s: s\_rho eta: zeta depth: h a: theta\_s b: theta\_b depth\_c: hc" ;
double s\_w(s\_w) ; s\_w:long\_name = "S-coordinate at W-points" ;
s\_w:valid\_min = -1. ; s\_w:valid\_max = 0. ; s\_w:standard\_name =
"ocean\_s\_coordinate" ; s\_w:formula\_terms = "s: s\_w eta: zeta depth:
h a: theta\_s b: theta\_b depth\_c: hc" ; float u(ocean\_time, s\_rho,
eta\_u, xi\_u) ; u:long\_name = "u-momentum component" ; u:units =
"meter second-1" ; u:coordinates = "lat\_u lon\_u" ; u:grid = "grid" ;
// SGRID attribute u:location = "edge1" ; // SGRID attribute float
v(ocean\_time, s\_rho, eta\_v, xi\_v) ; v:long\_name = "v-momentum
component" ; v:units = "meter second-1" ; v:coordinates = "lat\_v
lon\_v" ; v:grid = "grid" ; // SGRID attribute v:location = "edge2" ; //
SGRID attribute float zeta(ocean\_time, eta\_rho, xi\_rho) ;
zeta:long\_name = "free-surface" ; zeta:units = "meter" ; zeta:time =
"ocean\_time" ; zeta:coordinates = "lat\_rho lon\_rho" ; zeta:grid =
"grid" ; // SGRID attribute zeta:location = "face" ; // SGRID attribute
// SGRID variable int grid ; grid:cf\_role = grid\_topology
grid:topology\_dimension = 2 ; grid:node\_dimensions = "xi\_psi
eta\_psi" ; grid:face\_dimensions = "xi\_rho: xi\_psi (padding: both)
eta\_rho: eta\_psi (padding: both)" ; grid:edge1\_dimensions = "xi\_u:
xi\_psi eta\_u: eta\_psi (padding: both)" ; grid:edge2\_dimensions =
"xi\_v: xi\_psi (padding: both) eta\_v: eta\_psi" ;
grid:node\_coordinates = "lon\_psi lat\_psi" ; grid:face\_coordinates =
"lon\_rho lat\_rho" ; grid:edge1\_coordinates = "lon\_u lat\_u" ;
grid:edge2\_coordinates = "lon\_v lat\_v" ; grid:vertical\_dimensions =
"s\_rho: s\_w (padding: none)" ; // global attributes: :Conventions =
"CF-1.0" ; :title = "ROMS/TOMS 2.2 - Adria02 Uber Run" ; }

WRF (ARW version)
-----------------

The WRF-ARW also uses a C-grid. Again the model results can best be
captured by a 2D grid.

It might be interesting to verify the result for WRF-NMM since that
model uses an E-grid, but I couldn't find an example file.

netcdf wrfout\_v2\_Lambert { dimensions: Time = UNLIMITED ; // (13
currently) DateStrLen = 19 ; west\_east = 73 ; south\_north = 60 ;
west\_east\_stag = 74 ; bottom\_top = 27 ; south\_north\_stag = 61 ;
bottom\_top\_stag = 28 ; variables: char Times(Time, DateStrLen) ; float
U(Time, bottom\_top, south\_north, west\_east\_stag) ; U:description =
"x-wind component" ; U:units = "m s-1" ; U:grid = "grid" ; // SGRID
attribute U:location = "edge1" ; // SGRID attribute float V(Time,
bottom\_top, south\_north\_stag, west\_east) ; V:description = "y-wind
component" ; V:units = "m s-1" ; V:grid = "grid" ; // SGRID attribute
U:location = "edge2" ; // SGRID attribute float W(Time,
bottom\_top\_stag, south\_north, west\_east) ; W:description = "z-wind
component" ; W:units = "m s-1" ; W:grid = "grid" ; // SGRID attribute
W:location = "face" ; // SGRID attribute float T(Time, bottom\_top,
south\_north, west\_east) ; T:description = "perturbation potential
temperature (theta-t0)" ; T:units = "K" ; W:grid = "grid" ; // SGRID
attribute W:location = "face" ; // SGRID attribute float XLAT(Time,
south\_north, west\_east) ; XLAT:description = "LATITUDE, SOUTH IS
NEGATIVE" ; XLAT:units = "degree\_north" ; float XLONG(Time,
south\_north, west\_east) ; XLONG:description = "LONGITUDE, WEST IS
NEGATIVE" ; XLONG:units = "degree\_east" ; float ZNU(Time, bottom\_top)
; ZNU:description = "eta values on half (mass) levels" ; ZNU:units = ""
; float ZNW(Time, bottom\_top\_stag) ; ZNW:description = "eta values on
full (w) levels" ; ZNW:units = "" ; // SGRID variable int grid ;
grid:cf\_role = grid\_topology grid:topology\_dimension = 2 ;
grid:node\_dimensions = "west\_east\_stag south\_north\_stag
bottom\_top\_stag" ; grid:face\_dimensions = "west\_east:
west\_east\_stag (padding: none) south\_north: south\_north\_stag
(padding: none)" ; grid:face\_coordinates = "XLONG XLAT" ; // what to do
with ZNU/ZNW vertical coordinates? grid:vertical\_dimensions =
"bottom\_top: bottom\_top\_stag (padding: none)" ; // global attributes:
:TITLE = "OUTPUT FROM WRF V2.0 MODEL" ; }
