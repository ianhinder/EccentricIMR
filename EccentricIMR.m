
BeginPackage["EccentricIMR`", 
  {"EccentricPN`"}];

EccentricIMRWaveform;
EccentricPNSolution;

Begin["`Private`"];

$maxPNX = 0.11;
$blendStartX = 0.12;
$circularisationTime = -30;

etaOfq[q_] :=
 q/(1 + q)^2;

pnSolnToAssoc[soln_List, q_] :=
    Join[Association@Table[ToString[var[[1]]] -> var[[2]], {var,soln}],Association["q"->q]];

EccentricPNSolution[params_Association, {t1_, t2_}] :=
  Module[{pnSolnRules, dt=1.0},
    pnSolnRules = EccentricSoln[xeModel, N@etaOfq[params["q"]],
      {N@params["x0"], N@params["e0"], N@params["l0"], N@params["phi0"]},
      N@params["t0"], N/@{t1, t2, dt}];
    pnSolnToAssoc[pnSolnRules, params["q"]]];

timeToMergerFunction[] :=
  Function[{qq, ee, 
  ll}, 391.1958112997977 + 
  3.133913420161147*ee - 
  2492.9450481148215*ee^2 + 
  2.772121863176699*qq - 
  17.91996245212786*ee*qq + 
  8.118416384384208*qq^2 + 
  76.49439836997179*ee*
   Cos[63.45850634016995 - ll] + 
  1.154966280839567*qq*
   Cos[1.0227481831786127 + ll]];

EccentricIMRWaveform[params_Association, {t1_, t2_}] :=
  Module[{pnSoln, hCirc, ttm},
    pnSoln = EccentricPNSolution[params, {t1, t2}];
    hCirc = CircularWaveform[params["q"], 0.];
    ttm = timeToMergerFunction[];
    EccentricIMRWaveform[pnSoln, hCirc, params["q"]]];

EccentricIMRWaveform[pnSoln_Association, hCirc_List, q_, ttm_] :=
 Module[{xFn, xRef, tRef, t, eRef, ttpCal, tPeakEcc, tPeakCirc, hCircShifted, hFn, 
   hPN, tBlendWindow, ePN, lPN, pnSoln, tBlendStart, xPN, tPN},

  xFn = pnSoln["x"];
  xPN = $maxPNX;

  tPN = t /. FindRoot[xFn[t] == xPN, {t, xFn[[1, 1, 2]] - 100}];
  tBlendStart = t /. FindRoot[xFn[t] == $blendStartX, {t, xFn[[1, 1, 2]] - 100}];
  ePN = pnSoln["e"][tPN];
  lPN = pnSoln["l"][tPN];
  ttpCal = ttm[q,ePN,lPN]; 
  tPeakEcc = tPN + ttpCal;
  tPeakCirc = timeOfPeak[hCirc];
  
  hCircShifted = shifted[hCirc, -tPeakCirc + tPeakEcc];
  hFn = pnSoln["h"];
  hPN = Table[{t, hFn[t]}, {t, hFn[[1, 1, 1]], hFn[[1, 1, 2]], 1}];
  tBlendWindow = {tBlendStart, tPeakEcc + $circularisationTime};
  blendWaveforms[{hPN, hCircShifted}, tBlendWindow]];

End[];

EndPackage[];
