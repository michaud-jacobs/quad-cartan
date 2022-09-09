// Magma code to support the calculations in the paper Quadratic.Points on Non-Split Cartan Modular Curves.

// This code carries out the computations of Section 5.1 to obtain a new model for X_ns{13}. 
// The final part of the code verifies non-singularity.

load "eqn_data.m";

XNS13:=Curve(ProjectiveSpace(R),old_eqns);  // The curve X_ns(13)

S<X,Y,Z>:=PolynomialRing(Rationals(),3);   
f:=(-Y-Z)*X^3+(2*Y^2+Z*Y)*X^2+(-Y^3+Z*Y^2-2*(Z^2)*Y+Z^3)*X+(2*Z^2*Y^2-3*Z^3*Y);
XNSplus13:=Curve(ProjectiveSpace(S),f); // The curve X_ns^+(13),

Eqphi1:=-3*x_1+2*x_2;  // Equations for rho
Eqphi2:=-3*x_1+x_2+2*x_4-2*x_5;
Eqphi3:=x_1+x_2+x_4-x_5;
eqnsphi:=[Eqphi1,Eqphi2,Eqphi3];
phi:=map< XNS13->XNSplus13 | eqnsphi >; 

SvnPts:=PointSearch(XNSplus13,100);   

/////////////////////////////////////////////////////////////////////////////////////////////////


// We compute the pullbacks of the seven rational points on XNSplus13

PQ := AmbientSpace(XNS13);

T<x> := PolynomialRing(Rationals());   // Setting up a more general field that contains all the square roots we need
ds := { -163, -67, -19, -11, -7, -2 };
QQ := ext<Rationals() | [x^2 -d : d in ds] >;
PQQ := BaseChange(PQ,QQ);

quad_pts1 := [ ];
quad_pts2 := [ ];
ds := {};
for pt in SvnPts do
    S := Pullback(phi, pt);
    BS := BaseScheme(phi);
    D := Difference(S, BS);
    pb, K1 := PointsOverSplittingField(D);
    K2 := NumberField(AbsolutePolynomial(K1));
    d := Squarefree(Discriminant(Integers(K2)));
    K := QuadraticField(d);
    ds := ds join {d};
    pair := Points(Intersection(PQ,D),QQ);
    quad_pts1 := quad_pts1 cat [PQQ ! Eltseq(pair[1])];
    quad_pts2 := quad_pts2 cat [PQQ ! Eltseq(pair[2])];
end for;

T1 := TranslationOfSimplex(PQQ,quad_pts1 cat [quad_pts2[1],quad_pts2[2]]);
T2 := TranslationOfSimplex(PQQ,quad_pts2 cat [quad_pts1[1],quad_pts1[2]]);
TofS := T2^(-1)*T1;
EqTS:=DefiningEquations(TofS); 
*/

/////////////////////////////////////////////////////////////////////////////////////////////////

w:=map<XNS13->XNS13 | EqTS>;          // The modular involution on the curve
Mw:=Transpose((Matrix(w)));           // Matrix of the modular involution 
Diag,T:=PrimaryRationalForm(Mw);      
assert T*Mw*(T^-1) eq Diag;

Eqg:=[];                              // We use T^-1 to find our change of coordinate map
for i in [1..8] do 
    ri:=&+[(T^-1)[i][j]*R.j : j in [1..8]];   
    Eqg:=Eqg cat [ri];
end for;
g:=hom<R->R | Eqg>;                   // Change of coordinate map


/////////////////////////////////////////////////////////////////////////////////////////////////

// Apply our change of coordinates to each of the 15 equations
Neqns:=[];            
for i in [1..15] do
    Neqn:=g(eqns[i]); 
    Neqns:=Neqns cat [Neqn];
end for;

// Apply change of coordinates to obtain new equations for map (to same bottom curve)
Nphis:=[];
for i in [1..3] do
    Nphi:=g(eqnsphi[i]);      
    Nphis:=Nphis cat [Nphi];
end for;

// We now have the following new data:
NX:=Curve(ProjectiveSpace(R),Neqns);                         // New model of our curve
Nw:=map<NX -> NX | [x_1,x_2,x_3,-x_4,-x_5,-x_6,-x_7,-x_8]>;  // New modular involution
Nphi:=map< NX -> XNSplus13 | Nphis >;                        // New equations for map
 
// Check that this new model is nonsingular at the primes used (rather long).
for p in [3,5,31,43,53,61,73] do 
    print "Starting p =", p;
    NXp:=ChangeRing(NX,GF(p));
    assert IsNonsingular(NXp);
    print "Nonsingular mod", p;
end for;


