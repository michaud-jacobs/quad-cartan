// Magma code to support the calculations in the paper Quadratic.Points on Non-Split Cartan Modular Curves.

// This code carries out the computations of Section 5.5 to test saturation.

// Defining X_{ns}^+(13).
S<X,Y,Z>:=PolynomialRing(Rationals(),3);   
f:=(-Y-Z)*X^3+(2*Y^2+Z*Y)*X^2+(-Y^3+Z*Y^2-2*(Z^2)*Y+Z^3)*X+(2*Z^2*Y^2-3*Z^3*Y);
C:=Curve(ProjectiveSpace(S),f);

p := 41; // Primes to test 

Zp3:=AbelianGroup([p,p,p]);          
Kls:=[];                            
for i in [1,12,23,26] do //   
    l:=NthPrime(i);                                  
    Cl:=ChangeRing(C,GF(l));           
    ClGrp,phi,psi:=ClassGroup(Cl);
    Z:=FreeAbelianGroup(1);
    degr:=hom<ClGrp->Z | [ Degree(phi(a))*Z.1 : a in OrderedGenerators(ClGrp)]>;  
    JFl:=Kernel(degr);                 // Jacobian mod l as an abelian group
    JFlmodp,pi:=quo<JFl | p*JFl>;   
    
    h1p:= Place(Cl ! [0,1,0]);
    h2p:= Place(Cl ! [0,0,1]);
    h3p:= Place(Cl ! [-1,0,1]);
    bp := Place(Cl ! [1,0,0]);
    Delta1p:=h1p-bp; Delta2p:=h2p-bp; Delta3p:=h3p-bp;  // Generators of H

    pil:=hom<Zp3->JFlmodp | [pi(psi(Delta1p)),pi(psi(Delta2p)),pi(psi(Delta3p))]>; 
    Kl:=Kernel(pil);
    Kls:=Kls cat [Kl];
    print  "Calculations completed for l =", l;
end for;

IntKl:=&meet(Kls);  // Intersection of kernels
if #IntKl eq 1 then 
print "Index not divisible by", p;          
end if;

