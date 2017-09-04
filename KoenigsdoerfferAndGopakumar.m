
BeginPackage["EccentricIMR`KoenigsdoerfferAndGopakumar`", 
  {"EccentricIMR`EccentricPNSymbols`", "EccentricIMR`Kappas`"}];

(* Notation:

Eps    -2 E
C      -2 E h^2

*)

rInne;
EpsInne;
jInne;
omInne;
lInne;
vInu;
betaphIneph;
hPlus;
hCross;
kappaEIne;
kappaJIne;
nDotInne;
eDotInne;
betaPhi;
kappaE;
kappaJ;

Begin["`Private`"];

rInne = 
  Module[{rPN0, rPN1, rPN2, rPN3, rExpr1, rExpr2},
    rPN0 = n^(-2/3)(1 - e Cos[u]);
    rPN1 = rPN0 xi^(2/3)/(6(1 - e Cos[u]))(-18 + 2eta - (6 - 7eta)e Cos[u]);
    rPN2 = rPN0 xi^(4/
            3)/(72(1 - e^2)(1 - e Cos[u]))(-72(4 - 7eta) + (72 + 30eta + 
              8eta^2 - (72 - 231eta + 35eta^2)e Cos[u])(1 - e^2) - 
        36(5 - 2eta)(2 + e Cos[u])Sqrt[1 - e^2]);
    rPN3 = rPN0 xi^2/(1 - e^2)^2/(1 - e Cos[u])(-70/3 + 56221/840eta - 
        123/64Pi^2eta - 151/36eta^2 + 
        2/81eta^3 + (-2/3 + 87/16eta - 437/144eta^2 + 49/1296eta^3) * 
          e Cos[u] + (-52/3 + 2099/35eta - 41/64Pi^2eta - 341/18eta^2 - 
              4/81eta^3 + (4/3 - 87/8eta + 437/72eta^2 - 49/648eta^3)e Cos[
                  u])e^2 + (2/3 - 1/8eta + 5/36eta^2 + 
              2/81eta^3 + (-2/3 + 87/16eta - 437/144eta^2 + 
                    49/1296eta^3)e Cos[u])e^4 + (-30 + 412/9eta - 
              41/96Pi^2eta - 
              10/3eta^2 + (-45/2 + 1247/36eta - 41/192Pi^2eta - 
                    31/6eta^2)e Cos[
                  u] + (-10 + 29/3eta - 
                    11/3eta^2 + (5/2 - 83/12eta + 5/3eta^2)e Cos[u])e^2)*
          Sqrt[1 - e^2]);

    rExpr1 = rPN0 + rPN1 + rPN2 + rPN3 /. xi->n;
    rExpr2 = rExpr1 + O[n]^(6/3);
    rExpr2];

EpsInne = 
  Module[{x1,x2},
    x1= n^(2/3)(1 +
        xi^(2/3)/12(15 - eta)+
        xi^(4/3)/24(15 - 15eta - eta^2 + 24/Sqrt[1 - e^2](5 - 2eta)) + 
        xi^2/5184(-4995 - 6075eta - 450eta^2 - 35eta^3 + 
                864/Sqrt[1 - e^2](15 + 23eta - 20eta^2) + 
                18/(1 - e^2)^(3/2)(11520 - 15968eta + 123Pi^2eta + 
                      2016eta^2))) /. xi -> n;
    x2 = x1 + O[n]^(10/3)];

(* This is called jInne, but j == C == -2 E L^2 *)

jInne = 
  Module[{x1,x2},

    x1 = (1 - e^2)(1 + xi^(2/3)/(4(1 - e^2))(9 + eta - (17 - 7eta)e^2) + 
          xi^(4/3)/(24(1 - e^2)^2)(189 - 45eta + eta^2 - 
                2(111 + 7eta + 15eta^2)e^2 + (225 - 277eta + 
                      29eta^2)e^4 - (360 - 144eta)e^2Sqrt[1 - e^2]) + 
          xi^2/(6720(1 - e^2)^3)(35(5535 - 9061eta + 246Pi^2eta + 142eta^2 - 
                      eta^3) + (299145 - 1197667eta + 25830Pi^2eta + 
                      173250eta^2 + 2345eta^3)e^2 + 
                35(3549 - 12783eta + 6154eta^2 - 131eta^3)e^4 - 
                35(2271 - 7381eta + 2414eta^2 - 65eta^3)e^6 + 
                70(24(45 - 13eta - 2eta^2) - (17880 - 20120eta + 123Pi^2eta + 
                            2256eta^2)e^2 + 96(55 - 40eta + 3eta^2)e^4)Sqrt[
                    1 - e^2])) /. xi -> n ;

    x2 = x1 + O[n]^(8/3)];

omInne =
 Module[{om0,om1,om2,om3,omFinal1,omFinal2, om3LaTeX},
    om0 = n Sqrt[1 - e^2]/(1 - e Cos[u])^2;
    om1 = om0 xi^(2/3)/(1 - e^2)/(1 - e Cos[u]) (3 - (4-eta)e^2+(1-eta)e Cos[u]);
    om2 = om0 xi^(4/3)/12/(1 - e^2)^2/(1 - e Cos[u])^3(144 - 
          48eta - (162 + 68eta - 2eta^2)e^2 + (60 + 26eta - 
                20eta^2)e^4 + (18eta + 12eta^2)e^6 + (-216 + 125eta + 
                eta^2 + (102 + 188eta + 16eta^2)e^2 - (12 + 97eta - 
                      eta^2)e^4)e Cos[
              u] + (108 - 97eta - 
                5eta^2 + (66 - 136eta + 4eta^2)e^2 - (48 - 17eta + 
                      17eta^2)e^4)(e Cos[u])^2 + (-36 + 2eta - 
                8eta^2 - (6 - 70eta - 14eta^2)e^2)(e Cos[u])^3 + 
          18(1 - e Cos[u])^2(1 - 2 e^2 + e Cos[u])(5 - 2eta)Sqrt[1 - e^2]);

(*    om3 = om0 Get["! extract-expression /usr/center/raid1/hinder/projects/numrel-papers/BBHeccentric/konigs.tex 2411 2695"] /. phiDotN -> 1;*)

om3LaTeX = (phiDotN*xi^2*(75 - (1447*eta)/12 + 2*eta^2 + 
   e^10*((9*eta)/8 - 3*eta^2 - (3*eta^3)/4) + 
   e^8*((-863*eta)/24 + eta^2 + (47*eta^3)/12) + (205*eta*Pi^2)/64 + 
   e^4*(-50 + (175193*eta)/840 - (317*eta^2)/6 + (47*eta^3)/12 - 
     (41*eta*Pi^2)/8) + e^6*(18 + (2987*eta)/210 + (127*eta^2)/6 - 
     (25*eta^3)/4 + (41*eta*Pi^2)/32) + 
   e^2*(1/2 - (57021*eta)/280 + (361*eta^2)/6 - eta^3/3 + 
     (205*eta*Pi^2)/64) + e*(-285 + (488539*eta)/840 - (367*eta^2)/24 + 
     eta^3/24 + e^8*((379*eta)/12 + (25*eta^2)/2 + eta^3/6) - 
     (1025*eta*Pi^2)/64 + e^2*(-121/2 + (94097*eta)/140 - (4571*eta^2)/24 - 
       (73*eta^3)/24 - (451*eta*Pi^2)/64) + 
     e^6*(-54 + (16531*eta)/210 - (769*eta^2)/24 - (17*eta^3)/8 - 
       (41*eta*Pi^2)/16) + e^4*(182 - (114683*eta)/168 + (1987*eta^2)/24 + 
       (59*eta^3)/24 + (205*eta*Pi^2)/16))*Cos[u] + 
   e^2*(411 - (165061*eta)/168 + (121*eta^2)/24 + (3*eta^3)/8 - 
     e^8*((5*eta)/4 - (3*eta^2)/2 + eta^3/3) + (451*eta*Pi^2)/16 + 
     e^4*(-243 + (529223*eta)/840 - (329*eta^2)/24 + (45*eta^3)/8 - 
       (369*eta*Pi^2)/32) + e^6*(54 - (43177*eta)/840 - (227*eta^2)/8 - 
       (7*eta^3)/24 + (41*eta*Pi^2)/32) + 
     e^2*(213 - (268137*eta)/280 + (7693*eta^2)/24 - (3*eta^3)/8 + 
       (123*eta*Pi^2)/16))*Cos[u]^2 + 
   e^3*(-273 + (288269*eta)/420 + (189*eta^2)/8 + (55*eta^3)/24 + 
     e^6*(-18 - (347*eta)/12 + (245*eta^2)/24 - eta^3/24) - 
     (697*eta*Pi^2)/32 + e^4*(137 - (102349*eta)/420 + (53*eta^2)/24 - 
       (55*eta^3)/24 + (41*eta*Pi^2)/8) - 
     e^2*(281 - (79717*eta)/84 + (7705*eta^2)/24 + (119*eta^3)/24 + 
       (287*eta*Pi^2)/32))*Cos[u]^3 + 
   e^4*(78 - (77011*eta)/420 - (281*eta^2)/24 - (13*eta^3)/8 - 
     e^6*((11*eta)/8 + (23*eta^2)/24 - (29*eta^3)/24) + (451*eta*Pi^2)/64 + 
     e^4*(-23 + (5699*eta)/105 + (641*eta^2)/24 - (29*eta^3)/24 - 
       (41*eta*Pi^2)/32) + e^2*(325/2 - (92555*eta)/168 + (3083*eta^2)/24 + 
       (33*eta^3)/8 + (451*eta*Pi^2)/64))*Cos[u]^4 + 
   e^5*(-6 + (139*eta)/8 - (2*eta^2)/3 - eta^3/3 - 
     e^4*(3 - (533*eta)/24 + (91*eta^2)/6 + eta^3) - (41*eta*Pi^2)/64 + 
     e^2*(-69/2 + (13537*eta)/140 - (38*eta^2)/3 + (5*eta^3)/6 - 
       (123*eta*Pi^2)/64))*Cos[u]^5 + 
   (Sqrt[1 - e^2]*(1 - e*Cos[u])^3*(15840 - 16064*eta + 960*eta^2 + 
      e^4*(9600 - 4416*eta - 576*eta^2) + 123*eta*Pi^2 - 
      e^2*(38400 - 38464*eta + 2976*eta^2 + 246*eta*Pi^2) + 
      e*(9600 - 7680*eta + 1536*eta^2 + e^4*(7680 - 6816*eta + 2304*eta^2) + 
        e^2*(8640 - 21472*eta + 1344*eta^2 + 246*eta*Pi^2))*Cos[u] + 
      e^2*(-8160 + 12512*eta - 768*eta^2 - e^2*(4800 - 5472*eta + 
          1824*eta^2) - 123*eta*Pi^2)*Cos[u]^2))/192))/
 ((1 - e^2)^3*(1 - e*Cos[u])^5);

    om3 = om0 om3LaTeX /. phiDotN -> 1;

    omFinal1 = om0 + om1 + om2 + om3 /. {xi -> n};
    omFinal2 = omFinal1 + O[n]^(11/3)
 ];

vInu = 2ArcTan[betaPhi Sin[u]/(1 - betaPhi Cos[u])] + u;

lInne =
  Module[{kEq,kEq2},

    kEq = Module[{cu = Cos[u], su = Sin[u]},
        u - e su +
        xi^(4/3)/(8 Sqrt[1 - e^2](1 - e cu))((15eta - eta^2)e su Sqrt[
                  1 - e^2] + 12(5 - 2eta)(v - u)(1 - e cu)) + 
        xi^2/(6720(1 - e^2)^(3/2)(1 - e cu)^3)((67200 + 143868eta - 
                    4305Pi^2 eta - 62160eta^2 - 
                    280eta^3 - (134400 + 139896eta - 8610Pi^2eta - 
                          67200eta^2 - 3920eta^3)e Cos[
                        u] + (67200 - 752eta - 4305Pi^2eta - 15260eta^2 - 
                          1820eta^3)(e cu)^2 + (-148960eta + 45500eta^2 - 
                          1540eta^3 + (143640eta - 13440eta^2 - 
                                3920eta^3)e Cos[
                              u] - (1120eta + 11620eta^2 - 
                                1820eta^3)*(e cu)^2)e^2 + (3220eta - 
                          10220eta^2 + 1820eta^3)e^4)e su Sqrt[
                  1 - e^2] + (302400 - 461440eta + 4305Pi^2eta + 
                    33600eta^2 + (100800 - 97440eta + 36960eta^2)e^2)(v - 
                    u)(1 - e cu)^3) /. xi -> n] /. {v -> vInu}  (* /. (u -> u0 + du) /. Sin[du] -> du - du^3/6 *);

    kEq2 = kEq + O[n]^(8/3)
];

betaphIneph = (1 - Sqrt[1 - eph^2])/eph;

hPlus = -eta/R((1 + Cos[th]^2)(2 rDot r phiDot Sin[2 (phi - ph)] + 
  (1/r + r^2 phiDot^2 - rDot^2)Cos[2(phi - ph)]) - (* This sign is wrong, but doesn't affect the 2,2 mode *)
  Sin[th]^2(1/r - r^2 phiDot^2 - rDot^2));

hCross = -2eta/R Cos[th]((1/r + r^2 phiDot^2 - rDot^2)Sin[2(phi - ph)] - 
  2 rDot r phiDot Cos[2(phi - ph)]);


nDotInne = 
  Module[{nDot0,nDot1,nDot2,nDotTail,nDot},
    nDot0 = 1/5/(1 - e^2)^(7/2)(96 + 292e^2 + 37e^4);
    nDot1 = xi^(2/3)/280/(1 - e^2)^(9/2)(20368 - 
            14784eta + (219880 - 159600eta)e^2 + (197022 - 
                  141708eta)e^4 + (11717 - 8288eta)e^6) /. xi -> n;

    nDot2 = xi^(4/3)/30240/(1 - e^2)^(11/2)(12592864 - 13677408eta + 
              1903104eta^2 + (133049696 - 185538528eta + 
                    61282032eta^2)e^2 + (284496744 - 411892776eta + 
                    16650606060eta^2)e^4 + (112598442 - 142089066eta + 
                    64828848eta^2)e^6 + (3523113 - 3259980eta + 
                    1964256eta^2)e^8 + 
              3024(96 + 4268e^2 + 175e^6)(5 - 2eta)Sqrt[1 - e^2]) /. xi -> n;

    nDotTail = 384/5n^(14/3)Pi eta kappaE[e];

    nDot = n^(5/3) n^2 eta (nDot0 + nDot1 + nDot2) + nDotTail + O[n]^(16/3)
];

eDotInne = 
  Module[{eDot0,eDot1,eDot2,eSquaredDotTail,eDotTail,eDot},
    eDot0 = 1/15/(1 - e^2)^(5/2)(304 + 121e^2);
    eDot1 = xi^(2/3)/2520/(1 - e^2)^(7/2)(340968 - 
              228704eta + (880632 - 651252eta)e^2 + (125361 - 93184eta)e^4) /. xi -> n;
    eDot2 = xi^(4/3)/30240/(1 - e^2)^(9/2)(20815216 - 25375248eta + 
                4548096eta^2 + (87568332 - 128909916eta + 
                  48711348eta^2)e^2 + (69916862 - 93522570eta + 
                  42810096eta^2)e^4 + (3786543 - 4344852eta + 
                  2758560eta^2)e^6 + 
          1008(2672 + 6963e^2 + 565e^4)(5 - 2eta)Sqrt[1 - e^2]) /. xi -> n;

    eSquaredDotTail = -256/
          5n^(11/3) Pi eta((1 - e^2)kappaE[e] - Sqrt[1 - e^2] kappaJ[e]);

    eDotTail = (* If[e==0, 0, *) 1/(2e) eSquaredDotTail (* ] *);

    eDot1 = (-xi^(5/3)n eta e(eDot0 + eDot1 + eDot2) + eDotTail)  /. {xi -> n};
    eDot2 = eDot1 + O[n]^(13/3)
];

konigsKappaEExprHeld = HoldForm[Sum[
      p^3/4(BesselJ[p, p e]^2(1/e^4 - 1/e^2 + 1/3 + 
                  p^2(1/e^4 - 3/e^2 + 3 - e^2)) + 
            p(-4/e^3 + 7/e - 3e)BesselJ[p, p e] BesselJprime[p, p e] + 
            BesselJprime[p, p e]^2(1/e^2 - 1 + p^2(1/e^2 - 2 + e^2))), {p, 1, 
        pmax}]];

kappaEIne = ReleaseHold[konigsKappaEExprHeld];

konigsKappaJExprHeld = HoldForm[Sum[
        p^2/2 Sqrt[
            1 - e^2](p(3/e^2 - 2/e^4 - 1)BesselJ[p, p e]^2 + (2/e^3 - 1/e + 
                    2p^2(1/e^3 - 2/e + e))BesselJ[p, p e]BesselJprime[p, p e] + 
              2p(1 - 1/e^2)BesselJprime[p, p e]^2), {p, 1, pmax}]]

kappaJIne = ReleaseHold[konigsKappaJExprHeld];

End[];

EndPackage[];
