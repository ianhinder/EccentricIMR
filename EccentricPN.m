
BeginPackage["EccentricPN`", 
  {"EccentricPNSymbols`", "KoenigsdoerfferAndGopakumar`",
   "MemmesheimerEtAl`", "Kappas`"}];

ephSquaredInne;
PrecomputeEccentric;
EccentricSoln;

Begin["`Private`"];

(* FIXME *)
t = Global`t;

Col[i_, l_] := Map[{#[[1]], #[[i]]} &, l];

AmpPhase[tb_] :=
  Module[{ampPhaseTb,z,t,previousPhase, i, currentPhase = 0, cycles =
          0, nPoints},
  nPoints = Length[tb];
  ampPhaseTb = Table[i, {i, 1, nPoints}];
  For[i = 1, i <= nPoints, i++,
   z = tb[[i, 2]];
   t = tb[[i, 1]];
   currentPhase = ArcTan[Re[z], Im[z]];
   If[currentPhase - previousPhase > Pi, cycles--];
   If[currentPhase - previousPhase < -Pi, cycles++];
   previousPhase = currentPhase;
   ampPhaseTb[[i]] = {t, Abs[z], 2 Pi cycles + currentPhase}];
  Return[ampPhaseTb]]

Phase[tb_] := Col[3, AmpPhase[tb]];

Amplitude[tb_] := Col[2, AmpPhase[tb]];



ephSquaredInne = 
  Simplify[ComposeSeries[ephSquaredInEpsj /. j -> jInne, EpsInne]];

(* This is computed analytically; we should do that here *)
h22 = (-4*eta*Sqrt[Pi/5]*(1 + r*(phiDot*r + I*rDot)^2))/(E^((2*I)*phi)*r*R);


FirstTerms[SeriesData_[var_, about_, terms_,nmin_,nmax_,den_], n_] :=
  Module[{coeffs,powers},
    coeffs = Take[terms,n];
    powers = Take[Table[var^(i/den), {i, nmin, nmax}],n];
    Dot[coeffs,powers]
  ];

(*

* Make the eccentric solution generator generic

* Allow arbitrary eta

* Make sure all quantities are clearly derived from papers

* Allow arbitrary PN order

*)


PrecomputeEccentric[] :=
  Module[{},
    nDotRHS = FirstTerms[nDotInne, 5] /. {n -> n[t], e -> e[t]};
    eDotRHS = FirstTerms[eDotInne, 5] /. {n -> n[t], e -> e[t]};

    ephExpr = 
      Simplify[Sqrt[Normal[ephSquaredInne]]] /. eta -> 0.25;
    keplerRHS = 
     Simplify[
      FirstTerms[lInne, 7] /. {eta -> 0.25}];
    rExpr = Normal[rInne] /. eta -> 0.25;
    omExpr = Normal[omInne] /. eta -> 0.25;
    konigsh22 = h22 /. {eta -> 0.25, R -> 1};
    hAmpPhaseExpr = hamp[t] Exp[I hphase[t]];
    psi4AmpPhaseExpr = Simplify[D[hAmpPhaseExpr, t, t]];
    psi4DotAmpPhaseExpr = Simplify[D[hAmpPhaseExpr, t, t, t]];
];

EccentricSoln[{n0_?NumberQ, e0_?NumberQ, l0_?NumberQ, phi0_}, t0_?NumberQ, 
              {t1p_?NumberQ, t2p_?NumberQ, dt_?NumberQ}] :=
  Module[{l, neSoln, nFn, eFn, nTb,eTb,
          tTb, ephTb, betaphTb, eph0, betaPhi0, lSoln, lFn, lTb, uTb, uFn, rTb,
          rFn, rDotFn, rDotTb, omTb, omFn, phiSoln, phiFn, phiTb, hTb, hFn, hPhase, 
          hPhaseFn, hAmp, hAmpFn, 
          psi4Tb, psi4DotTb, psi4OmTb, psi4OmFn, psi4phSoln, psi4phFn,
          psi4Phase, psi4PhaseFn, psi4phTb, 
          ord = 5, t2, delta = 3*dt, indeterminate, extend, t1, return},
    indeterminate = {_ -> Indeterminate};
    extend = !(phi0 === None);
    t1 = If[extend, Min[Floor[t0,dt],t1p], Floor[t1p,dt]];
    t2 = Ceiling[t2p,dt];
(*    Print["t1 = ", t1];
    Print["t2 = ", t2];*)
    return = False;
    Quiet[Check[
      neSoln = NDSolve[{D[n[t], t] == nDotRHS, D[e[t], t] == eDotRHS, 
         n[t0] == n0, e[t0] == e0} /. eta -> 0.25,
       {n, e}, {t, Floor[Min[t1,t0]-delta,dt], Ceiling[t2+delta,dt]}][[1]],
       return = True,
       NDSolve::ndsz], NDSolve::ndsz];

    If[return, Return[{psi4Om -> Indeterminate, psi4Phi -> Indeterminate}]];

    nFn = n /. neSoln; eFn = e /. neSoln;
(*    t2 = Min[t2p, nFn[[1]][[1]][[2]]];*)
    If[nFn[[1]][[1]][[2]] < t1, Return[indeterminate]];
    nTb = Table[nFn[t], {t, t1, t2, dt}];
    eTb = Table[eFn[t], {t, t1, t2, dt}];
    tTb = Table[t, {t, t1, t2, dt}];
    ephTb = MapThread[
      ephExpr /. {n -> #1, e -> #2} &, {nTb, eTb}];
    betaphTb = MapThread[betaphIneph /. {eph -> #1} &, {ephTb}];

(*    eph0 = ephExpr /. {n -> n0, e -> e0};
    betaPhi0 = konigsbetaphi /. eph -> eph0;
    l0 = keplerRHS /. {n -> n0, e -> e0, u -> u0, betaPhi -> betaPhi0};*)

    lSoln = NDSolve[{D[l[t], t] == n[t], l[t0] == l0} /. neSoln, {l}, {t, Min[t1,t0], 
        Max[t2,t0]}][[1]];

    lFn = l /. lSoln;
    lTb = Map[lFn, tTb];

    uTb = MapThread[(u /. 
         FindRoot[#1 == (keplerRHS) /. {n -> #4, e -> #3, 
            betaPhi -> #2}, {u, #1, #1+0.1}]) &, {lTb, betaphTb, eTb, nTb}];

    uFn = Interpolation[MapThread[List, {tTb, uTb}], InterpolationOrder->ord];

    rTb = MapThread[
       rExpr /. {u -> #1, n -> #2, e -> #3} &, {uTb, nTb, eTb}];
    rFn = Interpolation[MapThread[List, {tTb, rTb}], InterpolationOrder->ord];
    rDotFn = Derivative[1][rFn];
    rDotTb = Table[rDotFn[t], {t, t1, t2, dt}];
    omTb = MapThread[
       omExpr /. {u -> #1, n -> #2, e -> #3} &, {uTb, nTb, eTb}];
    omFn = Interpolation[MapThread[List, {tTb, omTb}], InterpolationOrder->ord];


(*    Quiet[ *)
    If[phi0 === None,
      phiSoln = 
       NDSolve[{phi'[t] == omFn[t], phi[t1] == 0}, {phi}, {t, t1, t2}][[1]],
      phiSoln = 
       NDSolve[{phi'[t] == omFn[t], phi[t0] == phi0}, {phi}, {t, t1, t2}][[1]]];
    phiFn = phi /. phiSoln;

    phiTb = Table[phiFn[t], {t, t1, t2, dt}];
    hTb = MapThread[
       konigsh22 /. {phi -> #1, r -> #2, rDot -> #3, 
         phiDot -> #4} &, {phiTb, rTb, rDotTb, omTb}];
    hFn = Interpolation[MapThread[List, {tTb, hTb}], InterpolationOrder->ord];
    hPhase = Phase[MapThread[List, {tTb, hTb}]];
    hPhaseFn = Interpolation[hPhase, InterpolationOrder->ord];
    hAmp = Amplitude[MapThread[List, {tTb, hTb}]];
    hAmpFn = Interpolation[hAmp, InterpolationOrder->ord];
    psi4Tb = 
      Table[psi4AmpPhaseExpr /. {hamp -> hAmpFn, hphase -> hPhaseFn,t -> tx}, {tx, t1, t2, 
        dt}];

    psi4Fn = Interpolation[MapThread[List, {tTb, psi4Tb}], InterpolationOrder->ord];
    psi4Phase = Phase[MapThread[List, {tTb, psi4Tb}]];
    psi4PhaseFn = Interpolation[psi4Phase, InterpolationOrder->ord];

    psi4DotTb = 
      Table[psi4DotAmpPhaseExpr /. {hamp -> hAmpFn, hphase -> hPhaseFn, t -> tx}, {tx, t1, 
        t2, dt}];
    psi4OmTb = MapThread[Im[#1/#2] &, {psi4DotTb, psi4Tb}];
    psi4OmFn = Interpolation[MapThread[List, {tTb, psi4OmTb}], InterpolationOrder->ord];
(*    ,
    InterpolatingFunction::dmval];*)

    Return[{r -> rFn, om -> omFn, h -> hFn, psi4 -> psi4Fn, phi -> phiFn, psi4Om -> psi4OmFn, psi4Phi -> psi4PhaseFn,
            n -> nFn, e -> eFn, u -> uFn}];
];

End[];

EndPackage[];
