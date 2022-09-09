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

T<x> := PolynomialRing(Integers());   // Setting up a more general fieldthat contains all the square roots we need
QQ<d1,d2,d4,d5,d6,d7>:=ext<Rationals() | x^2+11,x^2+67,x^2+2,x^2+19,x^2+163,x^2+7>;
D1:=13*d1; D2:=13*d2; D7:=13*d7; D4:=13*d4; D5:=13*d5; D6:=13*d6;

PQQ<u1,u2,u3,u4,u5,u6,u7,u8> := ProjectiveSpace(QQ,7); // Projective space to be used for translation to simplex afterwards

S1a :=  [44 , 66 , -D1-11 , -D1-55 , -D1+55 , -2*D1+22 , D1+33 , 22]; // The simplified pairs of points
S1b :=  [44 , 66 , +D1-11 , +D1-55 , +D1+55 , +2*D1+22 ,-D1+33 , 22];   
S2a :=  [536 , 804 , -4*D2-134 , D2+201 , D2-201 , 4*D2+268 , 2*D2+402 , 268];
S2b :=  [536 , 804 , +4*D2-134 ,-D2+201 ,-D2-201 ,-4*D2+268 ,-2*D2+402 , 268];
S3a :=  [-5*D7+175 , -D7+35 , -2*D7+70 , -D7+207 , 6*D7-38 , 172 , -3*D7-67 , 172];
S3b :=  [+5*D7+175 , +D7+35 , +2*D7+70 , +D7+207 ,-6*D7-38 , 172 , +3*D7-67 , 172];
S4a :=  [4 , -20 , 12 , 8 , -8 ,-D4+2 , -10 , 2];
S4b :=  [4 , -20 , 12 , 8 , -8 ,+D4+2 , -10 , 2];
S5a :=  [-3*D5+171 , 2*D5-114 , -D5+227 , -D5-113 , -2*D5-56 , 170 , D5-57 , 170];
S5b :=  [+3*D5+171 ,-2*D5-114 , +D5+227 , +D5-113 , +2*D5-56 , 170 ,-D5-57 , 170];
S6a :=  [-24*D6+652 ,-36*D6+978 , 11*D6+4907 , -119*D6+10521 , 133*D6+3675 , -2*D6+10466 , -43*D6-24861 , 12494];
S6b :=  [+24*D6+652 ,+36*D6+978 ,-11*D6+4907 , +119*D6+10521 ,-133*D6+3675 , +2*D6+10466 , +43*D6-24861 , 12494];
S7a :=  [-3*D7+63 , -11*D7+231 , 4*D7-84 , -D7-95 , -2*D7-74 , -6*D7-222 ,-7*D7+31 , 116];
S7b :=  [+3*D7+63 , +11*D7+231 ,-4*D7-84 , +D7-95 , +2*D7-74 , +6*D7-222 ,+7*D7+31 , 116];

S1aP:=PQQ ! S1a; S1bP:=PQQ ! S1b; // Coercing points into projective space
S2aP:=PQQ ! S2a; S2bP:=PQQ ! S2b;
S3aP:=PQQ ! S3a; S3bP:=PQQ ! S3b;
S4aP:=PQQ ! S4a; S4bP:=PQQ ! S4b;
S5aP:=PQQ ! S5a; S5bP:=PQQ ! S5b;
S6aP:=PQQ ! S6a; S6bP:=PQQ ! S6b;
S7aP:=PQQ ! S7a; S7bP:=PQQ ! S7b;

T1 := TranslationOfSimplex(PQQ,[S1aP,S2aP,S3aP,S4aP,S5aP,S6aP,S7aP,S1bP,S2bP]);   // Map to take points to standard simplex in P8     
T2 := TranslationOfSimplex(PQQ,[S1bP,S2bP,S3bP,S4bP,S5bP,S6bP,S7bP,S1aP,S2aP]); 
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


