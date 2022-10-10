// Magma code to support the calculations in the paper Quadratic points on non-split Cartan modular curves.

// This code carries out the computations of Section 5.5 to test saturation.

C := Curve(ProjectiveSpace(S), eqn_X_plus); //  The curve X_ns^+(13)

l := 41; // Prime to test 

Zl3:=AbelianGroup([l,l,l]);          
Kps:=[];          
for p in [2,37,83,101] do                                                  
    Cp:=ChangeRing(C,GF(p));           
    ClGrp,phi,psi:=ClassGroup(Cp);
    Z:=FreeAbelianGroup(1);
    degr:=hom<ClGrp->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(ClGrp)]>;  
    JFp:=Kernel(degr);                 // Jacobian mod p as an abelian group
    JFpmodl,pi:=quo<JFp | l*JFp>;      // J(F_p) / l*J(F_p)
    
    Del_1 := psi(Place(Cp ! [0,1,0])  - Place(Cp ! [1,0,0]));
    Del_2 := psi(Place(Cp ! [0,0,1])  - Place(Cp ! [1,0,0]));
    Del_3 := psi(Place(Cp ! [-1,0,1]) - Place(Cp ! [1,0,0]));

    pi_p:=hom<Zl3->JFpmodl | [pi(Del_1),pi(Del_2),pi(Del_3)]>; 
    Kp:=Kernel(pi_p);
    Kps:=Kps cat [Kp];
    print  "Calculations completed for p =", p;
end for;

IntKp:=&meet(Kps);  // Intersection of kernels
if #IntKp eq 1 then 
print "Index not divisible by", l;          
end if;

