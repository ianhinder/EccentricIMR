
BeginPackage["EccentricIMR`", 
  {"EccentricPN`"}];

EccentricIMRWaveform;

Begin["`Private`"];

etaOfq[q_] :=
 q/(1 + q)^2;

pnSolnToAssoc[soln_List, q_] :=
    Join[Association@Table[ToString[var[[1]]] -> var[[2]], {var,soln}],Association["q"->q]];

EccentricPNWaveform[params_Association, {t1_, t2_}] :=
  Module[{pnSolnRules, pnSoln, dt=1.0},
    pnSolnRules = EccentricSoln[xeModel, N@etaOfq[params["q"]],
      {N@params["x0"], N@params["e0"], N@params["l0"], N@params["phi0"]},
      N@params["t0"], N/@{t1, t2, dt}];
    pnSoln = pnSolnToAssoc[pnSolnRules, params["q"]];
    pnSoln["h"]];

End[];

EndPackage[];

