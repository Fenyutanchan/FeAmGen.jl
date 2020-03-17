

#########################################
function make_contractor_script()::String
#########################################
  result_str = """

symbol diim, epsilon;
dimension diim;

vector q1, q2, q3, q1C, q2C, q3C;
vector k1,...,k10;
vector r1,...,r10;
vector K1,...,K10;
vector barK1,...,barK10;
vector mom,mom0,mom1,...,mom10;

Set Kn: K1,...,K10;
Set kn: k1,...,k10;
Set rn: r1,...,r10;
Set NULL: k1,...,k10,r1,...,r10,barK1,...,barK10;
Set MASSIVE: K1,...,K10;

set ALLLOOP: q1,q2,q3,q1C,q2C,q3C;
set LOOP: q1,q2,q3;
set LOOPC: q1C,q2C,q3C;
set NonLOOP: k1,...,k10,K1,...,K10,r1,...,r10,barK1,...,barK10;
set ALLMOM: k1,...,k10,K1,...,K10,r1,...,r10,barK1,...,barK10,q1,q2,q3,q1C,q2C,q3C;

symbol pi, im, I, sqrt2, shat;
auto symbol ver, gc;

***Dirac indices
auto symbol spa, spb, spv spw;

***Lorentz indices
auto index mua, mub, muv, muw, rho, dummyMU, epsMU;

Set RHO: rho,rho0,...,rho200;
Set EPSMU: epsMU1,...,epsMU100;
Set EPSMUC: epsMUC1,...,epsMUC100;
Set DUMMYMU: dummyMU1,...,dummyMU100;
Set NonEPSMU: mua,mua0,...,mua100,mub,mub0,...,mub100,muv, muv0,...,muv10000,muw, muw0,...,muw10000,dummyMU1,...,dummyMU100;

Set LOR: mua,mua0,...,mua100,mub,mub0,...,mub100,muv,muv0,...,muv10000,muw,muw0,...,muw10000,
         rho,rho0,...,rho100,dummyMU,dummyMU0,...,dummyMU100,epsMU1,...,epsMU100;

Set LORC: muaC,muaC0,...,muaC100,mubC,mubC0,...,mubC100,muvC,muvC0,...,muvC10000,muwC,muwC0,...,muwC10000,
          rhoC,rhoC0,...,rhoC100,dummyMUC,dummyMUC0,...,dummyMUC100,epsMUC1,...,epsMUC100;

Set ALLLOR: mua,mua0,...,mua100,mub,mub0,...,mub100,muv,muv0,...,muv10000,muw,muw0,...,muw10000,
            rho,rho0,...,rho100,dummyMU,dummyMU0,...,dummyMU100,epsMU1,...,epsMU100,
            muaC,muaC0,...,muaC100,mubC,mubC0,...,mubC100,muvC,muvC0,...,muvC10000,muwC,muwC0,...,muwC10000,
            rhoC,rhoC0,...,rhoC100,dummyMUC,dummyMUC0,...,dummyMUC100,epsMUC1,...,epsMUC100;


***For matching both indices and momenta
auto index var;
***For matching numbers
auto symbol int, setint;

vector ref, ref0, ref1,...,ref4;
symbol mass, mass1,...,mass4, width;

CFunction Den, PowDen, FV, FermionLoopPow, GhostLoopPow;
CFunction VecEpsilon, VecEpsilon1,...,VecEpsilon4, VecEps, VecEpsC;
CFunction Spinor, Spinor1,...,Spinor4, FermionChain;

CFunction SpUB, SpVB, SpU, SpV;
Set SPSET: SpUB, SpVB, SpU, SpV;
Set LSPSET: SpUB, SpVB;
Set RSPSET: SpU, SpV;

CFunction UB, VB, U, V;
Set INSPSET: UB, VB, U, V;
Set ILSPSET: UB, VB;
Set IRSPSET: U, V;

CFunction SP(symmetric), LMT(symmetric), Levi(antisymmetric);
CFunction SPC(symmetric);
CFunction GAij, PLij, PRij, ONEij(symmetric), Trace, Trace5, Trace5sym;
CFunction GA, PL, PR;

CFunction JJ(antisymmetric), FF(antisymmetric);

symbol dZ3x1, dZ3x2; 
symbol dtZ3x1, dtZ3x2; 
symbol dZ2tx1, dZ2tx2; 
symbol dZmtx1, dZmtx2; 
symbol dZ2x1, dZ2x2; 
symbol dZ1x1, dZ1x2; 
symbol dZ4x1, dZ4x2; 
symbol dtZ1x1, dtZ1x2; 
symbol dZ1Fx1, dZ1Fx2; 
symbol dZ1Ftx1, dZ1Ftx2; 
symbol dZgx1, dZgx2;
symbol d2s1, d2s2, d2s3, d2s4;
symbol f2s1, f2s2, f2s3, f2s4;



*----------------------------------------
#procedure Simplification()

id I^2 = -1;
id sqrt2^2 = 2;
id sqrt2^(-2) = 1/2;

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA(5))/2, GA(6)
* PL-->tagF, (1-GA(5))/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id Spinor?SPSET(int?,spa1?,mom?,ref?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,ref,mass);
endrepeat;
.sort
*
* trivial Lorentz indices contraction
*
repeat;
  id LMT(rho1?,rho1?) = diim;
  id LMT(rho1?,rho2?)^2 = diim;
  id LMT(rho1?,rho2?)*LMT(rho2?,rho3?) = LMT(rho1,rho3);

  id LMT(rho1?,rho2?)*FV(mom1?,rho1?) = FV(mom1,rho2);
  id LMT(rho1?,rho2?)*GAij(spa1?,spa2?,rho1?) = GAij(spa1,spa2,rho2);
  id LMT(rho1?,rho2?)*Levi(rho1?,rho3?,rho4?,rho5?) = Levi(rho2,rho3,rho4,rho5);
  id LMT(rho1?,rho2?)*VecEpsilon?{VecEps,VecEpsC}(int1?,rho2?,mom1?,ref1?,mass1?) = VecEpsilon(int1,rho1,mom1,ref1,mass1);
endrepeat;

id GAij(spa1?,spa2?,rho?)*FV(mom?,rho?) = GAij(spa1,spa2,mom);
id FV(mom1?,rho?)*FV(mom2?,rho?) = SP(mom1,mom2);

repeat;
  id FV(mom?,rho?)*Levi(var1?,var2?,var3?,rho?) = Levi(var1,var2,var3,mom);
  id LMT(rho1?,rho2?)*FermionChain(?vars1,GA(rho2?),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
endrepeat;

repeat id LMT(rho?NonEPSMU,rho0?)*VecEpsilon?{VecEps,VecEpsC}(int?,rho?NonEPSMU,?vars)
  = LMT(EPSMU[int],rho0)*VecEpsilon(int,EPSMU[int],?vars);
id GAij(spa1?,spa2?,rho?NonEPSMU)*VecEpsilon?{VecEps,VecEpsC}(int?,rho?NonEPSMU,?vars)
  = GAij(spa1,spa2,EPSMU[int])*VecEpsilon(int,EPSMU[int],?vars);
.sort

*
* vanishing 
*
id FV(mom?,rho?) * VecEpsilon?{VecEps,VecEpsC}(int?,rho?,mom?,ref?,mass?) = 0;
id FV(ref?,rho?) * VecEpsilon?{VecEps,VecEpsC}(int?,rho?,mom?,ref?,mass?) = 0;
.sort

*
* use EPSMU indices
*
id SP( FV(mom1?,0), VecEpsilon?{VecEps,VecEpsC}(int?,0,mom2?,ref?,mass?) )
  = FV(mom1,EPSMU[int]) * VecEpsilon(int,EPSMU[int],mom2,ref,mass);
id FV(mom1?,rho?) * VecEpsilon?{VecEps,VecEpsC}(int?,rho?,mom2?,ref?,mass?)
  = FV(mom1,EPSMU[int]) * VecEpsilon(int,EPSMU[int],mom2,ref,mass);
id VecEpsilon1?{VecEps,VecEpsC}(int1?, rho?, ?vars1) * VecEpsilon2?{VecEps,VecEpsC}(int2?, rho?, ?vars2)
  = SP( VecEpsilon1(int1, 0, ?vars1), VecEpsilon2(int2, 0, ?vars2) );
.sort

*
*Linearly expand momentum polynomial in FV and GA
*
id FV(var?,rho?) = FV(var,rho);
id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
id SP(rho1?,rho2?) = SP(rho1,rho2);
id Levi(rho1?,rho2?,rho3?,rho4?) = Levi(rho1,rho2,rho3,rho4);
.sort

*
* vanishing momentum scalar product also vanishes
*
id SP(mom?NULL,mom?NULL) = 0;

*
* Explain FermionLoopPow and GhostLoopPow
*
id FermionLoopPow(-1,int?) = (-1)^int;
id GhostLoopPow(-1,int?) = (-1)^int;
.sort

repeat id Levi(var1?,var2?,var3?,mom?Kn[int])
  = Levi(var1,var2,var3,kn[int]) + SP(Kn[int],Kn[int])/2/SP(kn[int],rn[int])*Levi(var1,var2,var3,rn[int]);
.sort


#endprocedure




































*----------------------------------------
#procedure ArrangeTrace()

repeat;
  id Trace(?vars1,PL,PL,?vars2) = Trace(?vars1,PL,?vars2);
  id Trace(?vars1,PR,PR,?vars2) = Trace(?vars1,PR,?vars2);
  id Trace(?vars1,PL,PR,?vars2) = 0;
  id Trace(?vars1,PR,PL,?vars2) = 0;

  id Trace(?vars1,GA(rho?ALLLOR),PL,?vars2) = Trace(?vars1,PR,GA(rho),?vars2);
  id Trace(?vars1,GA(rho?ALLLOR),PR,?vars2) = Trace(?vars1,PL,GA(rho),?vars2);

  id Trace(?vars1,GA(mom?),PL,?vars2) = Trace(?vars1,PR,GA(mom),?vars2);
  id Trace(?vars1,GA(mom?),PR,?vars2) = Trace(?vars1,PL,GA(mom),?vars2);
endrepeat;
.sort


repeat;
  id FV(mom?,rho?ALLLOR)*Trace(?vars1,GA(rho?ALLLOR),?vars2) = Trace(?vars1,GA(mom),?vars2);
  id LMT(rho1?ALLLOR,rho2?ALLLOR)*Trace(?vars1,GA(rho2?ALLLOR),?vars2) = Trace(?vars1,GA(rho1),?vars2);

  id Trace(?vars1,GA(mom?ALLMOM),GA(mom?ALLMOM),?vars2) = SP(mom,mom)*Trace(?vars1,?vars2);
  id SP(mom?NULL,mom?NULL) = 0;

  id Trace(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = Trace(?vars1,?vars2)*diim;
endrepeat;
.sort
 
id Trace(PL,?vars) = (1+sign_(nargs_(?vars)))/2*Trace(PL,?vars);
id Trace(PR,?vars) = (1+sign_(nargs_(?vars)))/2*Trace(PR,?vars);
id Trace(GA(var?),?vars) = (1+sign_(1+nargs_(?vars)))/2*Trace(GA(var),?vars);
.sort

id Trace(PL,?vars) = 1/2*Trace(?vars)-1/2*Trace5(?vars);
id Trace(PR,?vars) = 1/2*Trace(?vars)+1/2*Trace5(?vars);
.sort

***
***
repeat;
  id once Trace5(?vars) = 1/24*e_(rho100,rho101,rho102,rho103)*Trace(GA(rho100),GA(rho101),GA(rho102),GA(rho103),?vars);
  sum rho100;
  sum rho101;
  sum rho102;
  sum rho103;
endrepeat;
.sort
***
***

repeat id Trace(?vars1, GA(var?), ?vars2) = Trace(?vars1, var, ?vars2);
.sort

repeat;
  id once Trace(?vars) = g_(1,?vars);
  tracen, 1;
endrepeat;
.sort

contract;
.sort

id VecEpsilon?{VecEps,VecEpsC}(int?,mom0?,mom?,ref?,mass?) = FV(mom0,EPSMU[int])*VecEpsilon(int,EPSMU[int],mom,ref,mass);
id mom?NULL.mom?NULL = 0;
id mom1?.mom2? = SP(mom1,mom2);
id mom?(rho?ALLLOR) = FV(mom,rho);
id e_(rho1?,rho2?,rho3?,rho4?) = I*Levi(rho1,rho2,rho3,rho4);
.sort

id d_(rho1?,rho2?) = LMT(rho1,rho2);
repeat id LMT(rho1?ALLLOR,rho2?ALLLOR)*FermionChain(?vars1,GA(rho2?ALLLOR),?vars2) = FermionChain(?vars1,GA(rho1),?vars2);
.sort

id VecEpsilon?{VecEps,VecEpsC}(int?,rho?NonEPSMU,mom?,ref?,mass?)*FermionChain(?vars1,GA(rho?NonEPSMU),?vars2) 
  = VecEpsilon(int,EPSMU[int],mom,ref,mass)*FermionChain(?vars1,GA(EPSMU[int]),?vars2);
.sort

#endprocedure































































*----------------------------------------
#procedure contractDiracIndices()

* FeynCalc has definition GA[6]=(1+GA[5])/2, GA[7]=(1-GA[5])/2
* FORM has a little different g_(6)=1+g_(5), g_(7)=1-g_(5)
* According to textbooks (Peskin's QFT etc.) PL=GA(7), PR=GA(6)
* PR-->tagJ, (1+GA(5))/2, GA(6)
* PL-->tagF, (1-GA(5))/2, GA(7)
* Same as definition in FeynCalc
*
***id PR(spa1?,spa2?) = GA(spa1,spa2,6);
***id PL(spa1?,spa2?) = GA(spa1,spa2,7);

*
* trivial Dirac indices contraction
*
repeat;
  id ONEij(spa?,spa?) = 4;
  id ONEij(spa1?,spa2?)*ONEij(spa2?,spa3?) = ONEij(spa1,spa3);

  id GAij(spa1?,spa2?,rho?ALLLOR) * GAij(spa2?,spa3?,rho?ALLLOR) = diim*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,mom?ALLMOM) * GAij(spa2?,spa3?,mom?ALLMOM) = SP(mom,mom)*ONEij(spa1,spa3);
  id GAij(spa1?,spa2?,rho?)*ONEij(spa3?,spa2?) = GAij(spa1,spa3,rho);
  id GAij(spa2?,spa1?,rho?)*ONEij(spa3?,spa2?) = GAij(spa3,spa1,rho);

  id PLij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*ONEij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PLij(spa2?,spa3?) = PLij(spa1,spa3);
  id PRij(spa1?,spa2?)*PRij(spa2?,spa3?) = PRij(spa1,spa3);

  id PLij(spa1?,spa2?)*PRij(spa2?,spa3?) = 0;
  id PRij(spa1?,spa2?)*PLij(spa2?,spa3?) = 0;


  id Spinor?SPSET(int?,spa1?,mom?,ref?,mass?)*ONEij(spa1?,spa2?) = Spinor(int,spa2,mom,ref,mass);
endrepeat;
.sort

id GAij(spa1?,spa2?,var?) = GAij(spa1,spa2,var);
.sort



*
* now chainin the Dirac objects in FermionChain according to Dirac indices
*
***Here var? could be mom? or rho?
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,mom?) = FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa1?,spa2?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );

***flip
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,mom?) = -FermionChain( ILSPSET[setint](int,?vars), GA(mom), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ILSPSET[setint](int,?vars), GA(rho), spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PLij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PL, spa2 );
id Spinor?LSPSET[setint](int?,spa1?,?vars)*PRij(spa2?,spa1?) = FermionChain( ILSPSET[setint](int,?vars), PR, spa2 );


id Spinor1?LSPSET[setint1](int1?,spa?,?var1)*Spinor2?RSPSET[setint2](int2?,spa?,?var2) = FermionChain( ILSPSET[setint1](int1,?var1), IRSPSET[setint2](int2,?var2) );
.sort

repeat;
***Here var? could be mom? or rho?
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,mom?) = FermionChain( ?vars, GA(mom), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa1?,spa2?,rho?ALLLOR) = FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa1?,spa2?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa1?,spa2?) = FermionChain( ?vars, PR, spa2 );

***flip
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,mom?) = -FermionChain( ?vars, GA(mom), spa2 );
  id FermionChain(?vars,spa1?)*GAij(spa2?,spa1?,rho?ALLLOR) = -FermionChain( ?vars, GA(rho), spa2 );

  id FermionChain(?vars,spa1?)*PLij(spa2?,spa1?) = FermionChain( ?vars, PL, spa2 );
  id FermionChain(?vars,spa1?)*PRij(spa2?,spa1?) = FermionChain( ?vars, PR, spa2 );
endrepeat;

id FermionChain(?vars1,spa?)*Spinor?RSPSET[setint](int?,spa?,?vars2) = FermionChain( ?vars1, IRSPSET[setint](int,?vars2) );
.sort


*
* Look for Trace
*
repeat;
  id once, GAij(spa1?,spa2?,var?) = Trace(GA(var),spa1,spa2);
  repeat;
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,mom?) = Trace(?vars,GA(mom),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa2?,spa3?,rho?ALLLOR) = Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa2?,spa3?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa2?,spa3?) = Trace(?vars,PR,spa1,spa3);

    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,mom?) = -Trace(?vars,GA(mom),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*GAij(spa3?,spa2?,rho?ALLLOR) = -Trace(?vars,GA(rho),spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PLij(spa3?,spa2?) = Trace(?vars,PL,spa1,spa3);
    id Trace(?vars,spa1?,spa2?)*PRij(spa3?,spa2?) = Trace(?vars,PR,spa1,spa3);
  endrepeat;
  id Trace(?vars,spa?symbol_,spa?symbol_) = Trace(?vars);
  id PLij(spa?symbol_,spa?symbol_) = 2;
  id PRij(spa?symbol_,spa?symbol_) = 2; 
  id ONEij(spa?symbol_,spa?symbol_) = 4;
endrepeat;
.sort

#call ArrangeTrace();

***move PL and PR to left-side of FermionChain, right-side of spinor
repeat;
  id FermionChain( ?vars1, PL, PR, ?vars2 ) = 0;
  id FermionChain( ?vars1, PR, PL, ?vars2 ) = 0;
  id FermionChain( ?vars1, PL, PL, ?vars2 ) = FermionChain( ?vars1, PL, ?vars2 );
  id FermionChain( ?vars1, PR, PR, ?vars2 ) = FermionChain( ?vars1, PR, ?vars2 );

  id FermionChain( ?vars1, GA(mom?), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(mom), ?vars2 );
  id FermionChain( ?vars1, GA(mom?), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(mom), ?vars2 );

  id FermionChain( ?vars1, GA(rho?ALLLOR), PL, ?vars2 ) = FermionChain( ?vars1, PR, GA(rho), ?vars2 );
  id FermionChain( ?vars1, GA(rho?ALLLOR), PR, ?vars2 ) = FermionChain( ?vars1, PL, GA(rho), ?vars2 );
endrepeat;
.sort


repeat;
  id FV(mom?,rho?ALLLOR)*FermionChain(?vars1,GA(rho?ALLLOR),?vars2) = FermionChain(?vars1,GA(mom),?vars2);
  id LMT(rho1?ALLLOR,rho2?ALLLOR)*FermionChain(?vars1,GA(rho1?ALLLOR),?vars2) = FermionChain(?vars1,GA(rho2),?vars2);
endrepeat;
.sort

repeat;
  id FermionChain(?vars1,GA(mom?),GA(mom?),?vars2) = SP(mom,mom)*FermionChain(?vars1,?vars2);
  id SP(mom?NULL,mom?NULL) = 0;
  id FermionChain(?vars1,GA(rho?ALLLOR),GA(rho?ALLLOR),?vars2) = diim*FermionChain(?vars1,?vars2);

*** Dirac equation for U and V (UB and VB) 
  id FermionChain( ?vars, GA(mom?), U(int?,mom?,ref?,0) ) = 0;
  id FermionChain( ?vars, GA(mom?), V(int?,mom?,ref?,0) ) = 0;

  id FermionChain( UB(int?,mom?,ref?,0), GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,ref?,0), GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,ref?,0), PL, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,ref?,0), PL, GA(mom?), ?vars ) = 0;

  id FermionChain( UB(int?,mom?,ref?,0), PR, GA(mom?), ?vars ) = 0;
  id FermionChain( VB(int?,mom?,ref?,0), PR, GA(mom?), ?vars ) = 0;


endrepeat;
.sort

#endprocedure

"""

  return result_str

end # function make_contractor_script









###############################################################################
function make_baseINC_script( graph::GenericGraph )::String
###############################################################################

  result_str = string()

  v0 = vertex_from_label( "graph property", graph )
  n_inc = v0.attributes["n_inc"]
  n_out = v0.attributes["n_out"]
  n_leg = n_inc+n_out

  ext_edge_list = filter( e_ -> ( e_.attributes["style"]=="External" ), edges(graph) )
  sorted_ext_edge_list = sort( ext_edge_list, by=edge_index )

  momN = sorted_ext_edge_list[n_leg].attributes["momentum"]
  momNm1 = sorted_ext_edge_list[n_leg-1].attributes["momentum"]
  momNm2 = sorted_ext_edge_list[n_leg-2].attributes["momentum"]

  #-------------------------------------------------------------
  sorted_notN_ext_edge_list = filter( e_ -> ( edge_index(e_) != n_leg ), sorted_ext_edge_list )

  result_str *= "id FV($(momN),rho?) = ";
  for edge in sorted_notN_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1 : (-1)
    result_str *= "+($(inc_sign))*FV($(mom),rho)"
  end # for edge
  result_str *= ";\n"
 
  result_str *= "id SP($(momN),rho?) = ";
  for edge in sorted_notN_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1 : (-1)
    result_str *= "+($(inc_sign))*SP($(mom),rho)"
  end # for edge
  result_str *= ";\n"

  result_str *= 
    "id FermionChain( Spinor1?ILSPSET(int1?!{,$(n_leg)},mom1?,ref1?,mass1?),"*
    " ?vars1, GA($(momN)), ?vars2,"*
    " Spinor2?IRSPSET(int2?!{,$(n_leg)},mom2?,ref2?,mass2?) ) = \n"
  for edge in sorted_notN_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1 : (-1)
    result_str *= 
    "  +($(inc_sign))*FermionChain( Spinor1(int1,mom1,ref1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2(int2,mom2,ref2,mass2) )\n"
  end # for edge
  result_str *= ";\n"

  #-------------------------------------------------------------
  sorted_notNm1_ext_edge_list = filter( e_ -> ( edge_index(e_) != n_leg-1 ), sorted_ext_edge_list )
  Nm1_sign = n_leg-1 <= n_inc ? (-1) : (+1)

  result_str *=
    "id FermionChain( Spinor?ILSPSET($(n_leg),mom?,ref?,mass?), ?vars1, GA($(momNm1)), ?vars2 ) = \n"
  for edge in sorted_notNm1_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1*Nm1_sign : (-1)*Nm1_sign
    result_str *=
    "  +($(inc_sign))*FermionChain( Spinor($(n_leg),mom,ref,mass), ?vars1, GA($(mom)), ?vars2 )\n"
  end # for edge
  result_str *= ";\n"

  result_str *=
    "id FermionChain( ?vars1, GA($(momNm1)), ?vars2, Spinor?IRSPSET($(n_leg),mom?,ref?,mass?) ) = \n"
  for edge in sorted_notNm1_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1*Nm1_sign : (-1)*Nm1_sign
    result_str *=
    "  +($(inc_sign))*FermionChain( ?vars1, GA($(mom)), ?vars2, Spinor($(n_leg),mom,ref,mass) )\n"
  end # for edge
  result_str *= ";\n"

  #-------------------------------------------------------------
  sorted_notNm2_ext_edge_list = filter( e_ -> ( edge_index(e_) != n_leg-2 ), sorted_ext_edge_list )
  Nm2_sign = n_leg-2 <= n_inc ? (-1) : (+1)
  
  result_str *=
    "id FermionChain( Spinor1?ILSPSET($(n_leg-1),mom1?,ref1?,mass1?), ?vars1, GA($(momNm2)), ?vars2, Spinor2?IRSPSET($(n_leg),mom2?,ref2?,mass2?) ) = \n"
  for edge in sorted_notNm2_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1*Nm2_sign : (-1)*Nm2_sign
    result_str *=
    "  +($(inc_sign))*FermionChain( Spinor1($(n_leg-1),mom1,ref1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2($(n_leg),mom2,ref2,mass2) )\n"
  end # for edge
  result_str *= ";\n"

  result_str *=
    "id FermionChain( Spinor1?ILSPSET($(n_leg),mom1?,ref1?,mass1?), ?vars1, GA($(momNm2)), ?vars2, Spinor2?IRSPSET($(n_leg-1),mom2?,ref2?,mass2?) ) = \n"
  for edge in sorted_notNm2_ext_edge_list
    mom = edge.attributes["momentum"]
    inc_sign = edge_index(edge) <= n_inc ? 1*Nm2_sign : (-1)*Nm2_sign
    result_str *=
    "  +($(inc_sign))*FermionChain( Spinor1($(n_leg),mom1,ref1,mass1), ?vars1, GA($(mom)), ?vars2, Spinor2($(n_leg-1),mom2,ref2,mass2) )\n"
  end # for edge
  result_str *= ";\n"

  #-----------------------------------------------------------------------------------
  result_str *=
    "id FV($(momNm1),rho?)*VecEpsilon?{VecEps,VecEpsC}($(n_leg),rho?,$(momN),r$(n_leg)?,mass?) = \n"
  for index in 1:(n_leg-2)
    edge = sorted_ext_edge_list[index]
    mom = edge.attributes["momentum"]
    inc_sign = index <= n_inc ? (+1) : (-1)
    result_str *=
    "  +($(inc_sign))*FV($(mom),rho)*VecEpsilon($(n_leg),rho,$(momN),r$(n_leg),mass)\n"
  end # for index
  result_str *= 
    ";\n"*
    "id FV(mom?,rho?)*VecEpsilon?{VecEps,VecEpsC}(int?,rho?,mom?,ref?,mass?) = 0;\n"*
    "\n"

  return result_str

end # function make_baseINC_script






##############################################################################
function make_amp_contraction_script( expr::Basic, file_name::String )::String
##############################################################################

  result_str = """
#-

Off Statistics;
Off FinalStats;

#include model_parameters.frm
#include contractor.frm

format nospaces;
format maple;

Local expression = $(expr);

#call Simplification();

#call contractDiracIndices();

#call Simplification();

#include kin_relation.frm
.sort

repeat;
  id once FermionChain(?vars1, GA(mom?), ?vars2 ) = FV(mom,rho100)*FermionChain(?vars1, GA(rho100), ?vars2 );
  sum rho100;
endrepeat;


id FV(mom?,rho?)*VecEpsilon?{VecEps,VecEpsC}(int?,rho?,mom?,ref?,mass?) = 0;
.sort

while( match(FermionChain(?vars1,GA(rho?NonEPSMU\$LORENTZ),?vars2)) );
  sum \$LORENTZ;
endwhile;
.sort
*
* Replace system dummy indices Nm_? by our dummy indices dummyMU in case to read back to GiNaC.
* We assume this should give the canonical form of FermionChain, 
*   since it seems dummy indices Nm_? can make canonical form of an expression automatically.
*

repeat;
if( match( SP(mom1?{q1,q2,q3}\$MOM1,mom2?\$MOM2) ) );
  id once SP(\$MOM1,\$MOM2) = FV(\$MOM1,rho1)*FV(\$MOM2,rho2)*LMT(rho1,rho2);
  sum rho1;
  sum rho2;
endif;
endrepeat; 
.sort


#do MUIDX = 1, 20, 1
  Multiply replace_(N`MUIDX'_?,dummyMU`MUIDX');
#enddo
.sort

id FV(rho1?,rho2?) = FV(rho1,rho2);
id SP(rho1?,rho2?) = SP(rho1,rho2);
.sort

#write <$(file_name).out> "%E", expression
#close <$(file_name).out>
.sort

#system tr -d "[:space:]" < $(file_name).out > $(file_name).out.trim
#system mv $(file_name).out.trim $(file_name).out
.sort

.end

"""

  return result_str

end # function make_amp_contraction_script 



#########################################
function make_color_script()::String
#########################################

  result_str = """
CFunctions DeltaFun(symmetric), DeltaAdj(symmetric), SUNTrace(cyclesymmetric);
CFunctions SUNT, SUNF, SUNTConj, sunTraceConj, sunTrace;
CFunctions SUNTChain, SUNTraceChain;

symbols cla0,...,cla100;
symbols claC0,...,claC100;
symbols claM0,...,claM100;
symbols clb0,...,clb100;
symbols clbC0,...,clbC100;
symbols clbM0,...,clbM100;
symbols clv0,...,clv10000;
symbols clw0,...,clw10000;

Symbols colorX, colorY;

Symbols n0,...,n100;
Symbols m0,...,m100;
Symbols I, im, ca, cf;

*----------------------------------------
#procedure calc1_CF()
**** we use the definition f(a,b,c) = -2*i*Tr(a,b,c)+2*i*Tr(c,b,a) as same as defined in MadGraph etc.

repeat id SUNTConj(?n1,n2?,n3?) = SUNT(reverse_(?n1),n3,n2);
repeat id sunTraceConj(?n0) = sunTrace(reverse_(?n0));
.sort

*
* Trim the cyclesymmetric of SUNTrace
* NB: But we should have replaced SUNTrace by sunTrace in MIRACLE master code.
*multiply replace_(SUNTrace,sunTrace);
*.sort

repeat;
  id SUNF(m1?,m2?,m3?) = -2*I*sunTrace(m1,m2,m3)+2*I*sunTrace(m3,m2,m1);
  id SUNT(m0?,n1?,n2?)*SUNT(m0?,n3?,n4?) = 1/2*( DeltaFun(n1,n4)*DeltaFun(n2,n3)-1/ca*DeltaFun(n1,n2)*DeltaFun(n3,n4) );
endrepeat;
.sort

repeat;
  id DeltaFun(n1?,n0?)*SUNT(m0?,n0?,n2?) = SUNT(m0,n1,n2);

  id DeltaFun(n1?,n0?)*DeltaFun(n0?,n2?) = DeltaFun(n1,n2);
  id DeltaFun(n0?,n0?) = ca;
  id DeltaAdj(n1?,n0?)*DeltaAdj(n0?,n2?) = DeltaAdj(n1,n2);
  id DeltaAdj(n0?,n0?) = 2*ca*cf;

  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?,?m2) = sunTrace(?m1,m3,?m2);
  id DeltaAdj(m0?,m3?)*SUNT(?m1,m0?,?m2,n1?,n2?) = SUNT(?m1,m3,?m2,n1,n2);

  id SUNT(?m3,m0?,?m4,n1?,n2?)*sunTrace(?m1,m0?,?m2) = 1/2*( SUNT(?m3,?m2,?m1,?m4,n1,n2)-1/ca*SUNT(?m3,?m4,n1,n2)*sunTrace(?m2,?m1) );

  id SUNT(?m1,m0?,?m2,n1?,n2?)*SUNT(?m3,m0?,?m4,n3?,n4?) 
    = 1/2*( SUNT(?m1,?m4,n1,n4)*SUNT(?m3,?m2,n3,n2) -1/ca*SUNT(?m1,?m2,n1,n2)*SUNT(?m3,?m4,n3,n4) );
  id SUNT(n1?,n2?) = DeltaFun(n1,n2);

  id SUNT(?m1,n1?,n2?)*SUNT(?m2,n2?,n3?) = SUNT(?m1,?m2,n1,n3);
  id SUNT(?m1,n1?,n1?) = sunTrace(?m1);
  id SUNT(m0?,m0?,n1?,n2?) = cf*DeltaFun(n1,n2);
  id SUNT(?m1,m0?,m0?,?m2,n1?,n2?) = cf*SUNT(?m1,?m2,n1,n2);

  id SUNT(m0?,?m2,m0?,n1?,n2?) = 1/2*( DeltaFun(n1,n2)*sunTrace(?m2) - 1/ca*SUNT(?m2,n1,n2) );

  id SUNT(?m1,m0?,?m2,m0?,?m3,n1?,n2?) = 1/2*( SUNT(?m1,?m3,n1,n2)*sunTrace(?m2) - 1/ca*SUNT(?m1,?m2,?m3,n1,n2) );
  id sunTrace(m0?) = 0;

  id sunTrace(m1?,m1?) = ca*cf;
  id sunTrace(m1?,m2?) = 1/2*DeltaAdj(m1,m2);
  id sunTrace(?m1,m0?,?m2)*sunTrace(?m3,m0?,?m4) = 1/2*( sunTrace(?m1,?m4,?m3,?m2)-1/ca*sunTrace(?m2,?m1)*sunTrace(?m4,?m3) );
endrepeat;
.sort

*
* Recover the cyclesymmetric of SUNTrace
*
multiply replace_(sunTrace,SUNTrace);
.sort


#endprocedure




*----------------------------------------
#procedure calc2_CF()

*
* Trim the cyclesymmetric of SUNTrace
*
multiply replace_(SUNTrace,sunTrace);
.sort

repeat id SUNF(m1?,m2?,m3?) = -2*I*sunTrace(m1,m2,m3)+2*I*sunTrace(m3,m2,m1);
.sort

repeat;
  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?) = sunTrace(?m1,m3);
  id DeltaAdj(m0?,m3?)*sunTrace(?m1,m0?,?m2) = sunTrace(?m1,m3,?m2);
  id DeltaAdj(m0?,m3?)*SUNT(?m1,m0?,?m2) = SUNT(?m1,m3,?m2);
  id DeltaAdj(m3?,m0?)*SUNT(?m1,m0?,?m2) = SUNT(?m1,m3,?m2);

  id SUNT(?m1,n1?,n2?)*DeltaFun(n2?,n3?) = SUNT(?m1,n1,n3);
  id SUNT(?m1,n1?,n2?)*DeltaFun(n3?,n2?) = SUNT(?m1,n1,n3);

  id DeltaFun(n1?,n2?)*DeltaFun(n2?,n3?) = DeltaFun(n1,n3);
  id DeltaFun(n1?,n1?) = ca;
  id DeltaAdj(n0?,n0?) = 2*ca*cf;

  id SUNT(?m1,n1?,n2?)*SUNT(?m2,n2?,n3?) = SUNT(?m1,?m2,n1,n3);
  id SUNT(?m0,m1?,n1?,n1?) = sunTrace(?m0,m1);

  id sunTrace(m1?,m1?) = ca*cf;
  id sunTrace(m1?,m2?) = 1/2*DeltaAdj(m1,m2);

  id sunTrace(?m1,m0?,m0?,?m2) = cf*sunTrace(?m1,?m2);

  id sunTrace(?m1,m0?,?m2,m0?,?m3) = 1/2*( sunTrace(?m1,?m3)*sunTrace(?m2) - 1/ca*sunTrace(?m1,?m2,?m3) );

  id sunTrace(m0?) = 0;

  id DeltaAdj(n1?,n0?)*DeltaAdj(n0?,n2?) = DeltaAdj(n1,n2);
  id sunTrace(?m1,m0?,?m2)*sunTrace(?m3,m0?,?m4) = 1/2*( sunTrace(?m1,?m4,?m3,?m2)-1/ca*sunTrace(?m2,?m1)*sunTrace(?m4,?m3) );
endrepeat;
.sort

id sunTrace() = ca;
.sort

*
* Recover the cyclesymmetric of SUNTrace
*
multiply replace_(sunTrace,SUNTrace);
.sort

#endprocedure

"""

  return result_str

end # function make_color_script 


############################################################################################
function make_simplify_color_factor_script( color_factor::Basic, file_name::String )::String
############################################################################################

  result_str = """
#-
Off Statistics;

format nospaces;
format maple;

#include color.frm

Local colorFactor = $(color_factor);

#call calc1_CF();
.sort 

#call calc2_CF();
.sort 

#write <$(file_name).out> "%E", colorFactor
#close <$(file_name).out>
.sort

#system tr -d "[:space:]" < $(file_name).out > $(file_name).out.trim
#system mv $(file_name).out.trim $(file_name).out
.sort

.end
"""

  return result_str

end # function make_simplify_color_factor_script









