// Magma code to support the calculations in the paper Quadratic.Points on Non-Split Cartan Modular Curves.

// This code carries out the computations of Section 5.1 to obtain a new model for X_ns{13}. 
// The final part of the code verifies non-singularity.

load "eqn_data.m";

X := Curve(ProjectiveSpace(R),old_eqns);  // The curve X_ns(13)

X_plus := Curve(ProjectiveSpace(S),eqn_X_plus); // The curve X_ns^+(13),

rho :=map < X -> X_plus | old_rho_eqns >; 

SvnPts:=PointSearch(XNSplus13,100);   
assert #SvnPts eq 7;

/////////////////////////////////////////////////////////////////////////////////////////////////


// We compute the pullbacks of the seven rational points on XNSplus13

// We first compute the fields of definition and some pullback schemes

Ds := [];
ds := {};
for pt in SvnPts do
    S := Pullback(phi, pt);
    BS := BaseScheme(phi);
    D := Difference(S, BS);
    Ds := Ds cat [D];
    pb, K1 := PointsOverSplittingField(D);
    K2 := NumberField(AbsolutePolynomial(K1));
    d := Squarefree(Discriminant(Integers(K2)));
    K := QuadraticField(d);
    ds := ds join {d};
    
end for;

T<x> := PolynomialRing(Rationals());   // Setting up a more general field that contains all the square roots we need
ds := { -163, -67, -19, -11, -7, -2 };
QQ := ext<Rationals() | [x^2 -d : d in ds] >;
PQ := AmbientSpace(XNS13);
PQQ := BaseChange(PQ,QQ);

quad_pts1 := [ ];
quad_pts2 := [ ];

for D in Ds do
    pair := Points(Intersection(PQ,D),QQ);
    quad_pts1 := quad_pts1 cat [PQQ ! Eltseq(pair[1])];
    quad_pts2 := quad_pts2 cat [PQQ ! Eltseq(pair[2])];
end for;

T1 := TranslationOfSimplex(PQQ,quad_pts1 cat [quad_pts2[1],quad_pts2[2]]);
T2 := TranslationOfSimplex(PQQ,quad_pts2 cat [quad_pts1[1],quad_pts1[2]]);
TofS := T2^(-1)*T1;
EqTS:=DefiningEquations(TofS); 

/////////////////////////////////////////////////////////////////////////////////////////////////

w:=map<XNS13->XNS13 | EqTS>;          // The modular involution on the curve
Mw:=Transpose((Matrix(w)));           // Matrix of the modular involution 
Diag,T:=PrimaryRationalForm(Mw);      
assert T*Mw*(T^-1) eq Diag;
                              // We use T^-1 to find our change of coordinate map
Eqg := [&+[(T^-1)[i][j]*R.j : j in [1..8]] : i in [1..8]];
g:=hom<R->R | Eqg>;                   // Change of coordinate map


/////////////////////////////////////////////////////////////////////////////////////////////////

// Apply our change of coordinates to obtain new equations

Neqns := [g(ee) : ee in old_eqns];
assert Neqns eq new_eqns; // Matches data file

// Apply change of coordinates to obtain new equations for map (to same bottom curve)

Nphis := [g(ee) : ee in eqnsphi];

// We now have the following new data:
NX:=Curve(ProjectiveSpace(R),Neqns);                         // New model of our curve
Nw:=map<NX -> NX | [x_1,x_2,x_3,-x_4,-x_5,-x_6,-x_7,-x_8]>;  // New modular involution
Nphi:=map< NX -> XNSplus13 | Nphis >;                        // New equations for map
 
// Check that this new model is nonsingular at the primes used (rather long).
/*
for p in [3,5,31,43,53,61,73] do 
    print "Starting p =", p;
    NXp:=ChangeRing(NX,GF(p));
    assert IsNonsingular(NXp);
    print "Nonsingular mod", p;
end for;
*/


