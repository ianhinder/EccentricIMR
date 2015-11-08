
BeginPackage["EccentricPN`", 
  {"EccentricPNSymbols`", "KoenigsdoerfferAndGopakumar`",
   "MemmesheimerEtAl`", "Kappas`"}];

(* TODO:

* Allow arbitrary eta (easy)

* Allow arbitrary PN order 

* Make sure all quantities are clearly derived from papers

* Rename "Tb" variables as "Vec" variables where appropriate

* Add much more error checking

* Cache all the slow symbolic calculations in a file

* Verify that numerical answers are identical with version from before
  cleanup, and that speed is not significantly worse

* Include cross-check between expressions from Koenigsdoerffer and
  Gopakumar and Memmesheimer et al..  This is currently in a notebook.

* Compute rDot analytically instead of numerically, and omega without
  having to convert from LaTeX or type in very long expressions.

* Add more comments

* Can we express the circular models in this framework?

*)

(*******************************************************************)
(* Public functions *)
(*******************************************************************)

EccentricSoln;
FirstTerms;

(*******************************************************************)
(* Eccentric PN models *)
(*******************************************************************)

neModel;
xeModel;
xeNewtModel;

(*******************************************************************)
(* PN expressions *)
(*******************************************************************)

ephSquaredInne;
nInxe;
xInne;
rInxe;
omInxe;
lInxe;
ephInxe;
xDotInxe;
eDotInxe;
xInEpsj;
eInEpsj;
EpsInxe;
jInxe;

(*******************************************************************)
(* Symbols *)
(*******************************************************************)

{om,h,psi4,psi4Om,hOm,hPhi,psi4Phi,X0,Y0,x0,e0,n0,l0,phi0,XDot,YDot,lInXY, nInXY, lInXY,
ephSquaredInXY, omInXY, rInXY};

Begin["`Private`"];

(*******************************************************************)
(* Utility functions *)
(*******************************************************************)

(* We might want to move these into their own package *)

Col[i_, l_] := Map[{#[[1]], #[[i]]} &, l];

AmpPhase[tb_List] :=
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

Phase[tb_List] := Col[3, AmpPhase[tb]];

Amplitude[tb_List] := Col[2, AmpPhase[tb]];

(* Truncate a power series to its first n terms *)
FirstTerms[SeriesData_[var_, about_, terms_,nmin_,nmax_,den_], n_] :=
  Module[{coeffs,powers},
    coeffs = Take[terms,n];
    powers = Take[Table[var^(i/den), {i, nmin, nmax}],n];
    Dot[coeffs,powers]];

(*******************************************************************)
(* PN computations *)
(*******************************************************************)

ephSquaredInne = 
  Simplify[ComposeSeries[ephSquaredInEpsj /. j -> jInne, EpsInne]];

(* This has been derived from Memmesheimer et al. paper in another
   notebook; we should do this here instead *)
nInxe = SeriesData[x, 0, {1, 0, 3/(-1 + e^2), 0, (-18 + 28*eta + e^2*(-51 + 26*eta))/
   (4*(-1 + e^2)^2), 0, 
  -(16*e^4*(-156 + 240/Sqrt[1 - e^2] + (110 - 96/Sqrt[1 - e^2])*eta - 
       65*eta^2) + 4*(48 - 480/Sqrt[1 - e^2] - 224*eta^2 + 
       eta*(3656 + 192/Sqrt[1 - e^2] - 123*Pi^2)) + 
     e^2*(96*(-89 - 20/Sqrt[1 - e^2]) - 5120*eta^2 + 
       eta*(17856 + 768/Sqrt[1 - e^2] - 123*Pi^2)))/(128*(-1 + e^2)^3)}, 3, 
 11, 2];

xInne = InverseSeries[nInxe, n];
rInxe = Simplify[rInne /. n -> nInxe];
omInxe = Simplify[omInne /. n -> nInxe];
lInxe = Simplify[lInne /. n -> nInxe];
ephInxe = Simplify[PowerExpand[Sqrt[ephSquaredInne /. n -> nInxe]]];
xDotInxe = 
 Simplify[
  D[xInne, n] nDotInne + 
    D[xInne, e] eDotInne /. n -> nInxe];
eDotInxe = Simplify[eDotInne /. n -> nInxe];

(* From Arun et al. *)
xInEpsj = Eps(1+Eps(-5/4+nu/12+2/j)
 +Eps^2(5/2+5/24nu+nu^2/18+(-5+2nu)/Sqrt[j]+1/j(-5+7/6nu)+1/j^2(33/2-5nu))
 +Eps^3(-235/48-25/24nu-25/576nu^2+35/1296nu^3+1/j(35/4-5/3nu+25/36nu^2)
        +1/Sqrt[j](145/8-235/24nu+29/12nu^2)+1/j^(3/2)(-45+(472/9-41/96Pi^2)nu-5nu^2)
        +1/j^2(-565/8+(1903/24-41/64Pi^2)nu-95/12nu^2)+1/j^3(529/3+(-610/3+205/64Pi^2)nu+35/4nu^2))) +O[Eps]^4 /.nu->eta;

etADM = Sqrt[(1-j+Eps/4(-8+8nu-j(-17+7nu))
  +Eps^2/8(8+4nu+20nu^2-j(112-47nu+16nu^2)-24Sqrt[j](-5+2nu)
  +4/j(17-11nu)-24/Sqrt[j](5-2nu))
  +Eps^3/192(24(-2+5nu)(-23+10nu+4nu^2)-15j(-528+200nu-77nu^2+24nu^3)
    -72Sqrt[j](265-193nu+46nu^2)-2/j(6732+117nu Pi^2-12508nu+2004nu^2)
    +2/Sqrt[j](16380-19964nu+123nu Pi^2+3240nu^2)
    -2/j^(3/2)(10080+123nu Pi^2-13952nu+1440nu^2)+96/j^2(134-281nu+5nu Pi^2+16nu^2))) + O[Eps]^4];

etMHminusADM = Eps^2/Sqrt[1-j](1/4+17/4nu)(1-1/j)+
  Eps^3/Sqrt[1-j](-19/32-52/3nu+225/32nu^2+1/j(29/16+(79039/1680-21/16Pi^2)nu-201/16nu^2)
  +1/j^2(-3/2+(-14501/420+21/16Pi^2)nu+5nu^2))+O[Eps]^4;

eInEpsj = etADM + etMHminusADM/.nu->eta;

EpsInxe = Simplify[ComposeSeries[EpsInne, nInxe]];

jInxe = Simplify[ComposeSeries[jInne, nInxe]];



(* FIXME: This is computed analytically; we should do that here instead *)
h22 = (-4*eta*Sqrt[Pi/5]*(1 + r*(phiDot*r + I*rDot)^2))/(E^((2*I)*phi)*r*R);

konigsh22 = h22 /. {eta -> 0.25, R -> 1};
hAmpPhaseExpr = hamp[t] Exp[I hphase[t]];
psi4AmpPhaseExpr = Simplify[D[hAmpPhaseExpr, t, t]];
psi4DotAmpPhaseExpr = Simplify[D[hAmpPhaseExpr, t, t, t]];

etaFixed = 0.25; (*2./9;*)

neModel = {
  X -> n,
  Y -> e,
  X0 -> n0,
  Y0 -> e0,
  XDot -> FirstTerms[nDotInne, 5] /. {n -> n[t], e -> e[t]},
  YDot -> FirstTerms[eDotInne, 5] /. {n -> n[t], e -> e[t]},
  ephSquaredInXY -> (Simplify[Sqrt[Normal[ephSquaredInne]]] /. eta -> etaFixed),
  lInXY -> Simplify[FirstTerms[lInne, 7] /. {eta -> etaFixed}],
  nInXY -> n,
  rInXY -> Normal[rInne] /. eta -> etaFixed,
  omInXY -> Normal[omInne] /. eta -> etaFixed};

xeModel = {
  X -> x,
  Y -> e,
  X0 -> x0,
  Y0 -> e0,
  XDot -> FirstTerms[xDotInxe, 5] /. {x -> x[t], e -> e[t]},
  YDot -> FirstTerms[eDotInxe, 5] /. {x -> x[t], e -> e[t]},
  ephSquaredInXY -> (Simplify[Normal[ephInxe^2]] /. eta -> etaFixed),
  lInXY -> Simplify[FirstTerms[lInxe, 7] /. {eta -> etaFixed}],
  nInXY -> Normal[nInxe] /. eta -> etaFixed,
  rInXY -> Normal[rInxe] /. eta -> etaFixed,
  omInXY -> Normal[omInxe] /. eta -> etaFixed};

xeNewtModel = {
  X -> x,
  Y -> e,
  X0 -> x0,
  Y0 -> e0,
  XDot -> FirstTerms[xDotInxe, 1] /. {x -> x[t], e -> e[t]},
  YDot -> FirstTerms[eDotInxe, 1] /. {x -> x[t], e -> e[t]},
  ephSquaredInXY -> (Simplify[FirstTerms[ephInxe^2,1]] /. eta -> etaFixed),
  lInXY -> Simplify[FirstTerms[lInxe, 1] /. {eta -> etaFixed}],
  nInXY -> FirstTerms[nInxe,1] /. eta -> etaFixed,
  rInXY -> FirstTerms[rInxe,1] /. eta -> etaFixed,
  omInXY -> FirstTerms[omInxe,1] /. eta -> etaFixed};

EccentricSoln[model_, eta0_?NumberQ, {x0_?NumberQ, y0_?NumberQ, l0_?NumberQ, phi0_},
              t0_?NumberQ, {t1p_?NumberQ, t2p_?NumberQ, 
              dt_?NumberQ}] :=
  Module[{x,y, xySoln, xFn, yFn, xTb,yTb,
          tTb, ephTb, betaphTb, eph0, betaPhi0, lSoln, lFn, lTb, uTb, uFn, rTb,
          rFn, rDotFn, rDotTb, omTb, omFn, phiSoln, phiFn, phiTb, 
          ord = 5, t2, delta = 3*dt, indeterminate, extend, t1, return, 
          adiabatic, lEqs},
    x = (X/.model);
    y = (Y/.model);
    indeterminate = {_ -> Indeterminate};
    extend = !(phi0 === None);
    (* t1 = If[extend, Min[Floor[t0,dt],t1p], Floor[t1p,dt]]; *)
    (* t2 = Ceiling[t2p,dt]; *)

    (* The above rounds t1, which I think isn't necessary.  What is
       necessary is that dt divides t2-t1. *)

    t1 = If[extend, Min[t0,t1p], t1p];
    t2 = t1 + Ceiling[t2p-t1,dt];

    return = False;
    adiabatic = {D[x[t], t] == (XDot /. model), 
                        D[y[t], t] == (YDot /. model), 
         x[t0] == x0, y[t0] == y0} /. eta -> eta0;

    Quiet[xySoln = NDSolve[adiabatic,
      {x, y}, {t, Floor[Min[t1,t0]-delta,dt], Ceiling[t2+delta,dt]}][[1]],NDSolve::ndsz];

    (* Print[xySoln]; *)
    (* Print["old t2 = ",t2]; *)
    t2 = t1 + Floor[(x/.xySoln)[[1,1,2]]-t1,dt];
    (* Print["new t2 = ",t2]; *)

    (* Quiet[Check[ *)
    (*   xySoln = NDSolve[adiabatic, *)
    (*    {x, y}, {t, Floor[Min[t1,t0]-delta,dt], Ceiling[t2+delta,dt]}][[1]], *)
    (*    return = True, *)
    (*    NDSolve::ndsz], NDSolve::ndsz]; *)

    (* If[return, Return[{psi4Om -> Indeterminate, psi4Phi -> Indeterminate, *)
    (*   hOm -> Indeterminate}]]; *)

    xFn = x /. xySoln; yFn = y /. xySoln;
    If[xFn[[1]][[1]][[2]] < t1, Return[indeterminate]];
    xTb = Table[xFn[t], {t, t1, t2, dt}];
    yTb = Table[yFn[t], {t, t1, t2, dt}];
    tTb = Table[t, {t, t1, t2, dt}];
    ephTb = MapThread[
      (ephSquaredInXY /. model) /. {x -> #1, y -> #2} &, {xTb, yTb}];
    betaphTb = MapThread[betaphIneph /. {eph -> #1} &, {ephTb}];

    lEqs = {D[l[t], t] == ((nInXY/.model)/.{x->x[t],y->y[t]}), l[t0] == l0} /. xySoln;
    lSoln = NDSolve[lEqs, {l}, 
      {t, Min[t1,t0], Max[t2,t0]}][[1]];

    lFn = l /. lSoln;
    lTb = Map[lFn, tTb];

    uTb = MapThread[(u /. 
         FindRoot[#1 == (lInXY/.model) /. {x -> #4, y -> #3, 
            betaPhi -> #2}, {u, #1, #1+0.1}]) &, {lTb, betaphTb, yTb, xTb}];

    (* FIXME: at this point, we might have a very small number number
       of points.  Is the calculation even meaningful?  Shouldn't we
       have some sort of automatic error check? *)

    uFn = Interpolation[MapThread[List, {tTb, uTb}], InterpolationOrder->ord];
    rTb = MapThread[
       (rInXY/.model) /. {u -> #1, x -> #2, y -> #3} &, {uTb, xTb, yTb}];
    rFn = Interpolation[MapThread[List, {tTb, rTb}], InterpolationOrder->ord];
    rDotFn = Derivative[1][rFn];
    rDotTb = Table[rDotFn[t], {t, t1, t2, dt}];
    omTb = MapThread[
       (omInXY/.model) /. {u -> #1, x -> #2, y -> #3} &, {uTb, xTb, yTb}];
    omFn = Interpolation[MapThread[List, {tTb, omTb}], InterpolationOrder->ord];

    If[phi0 === None,
      phiSoln = 
       NDSolve[{phi'[t] == omFn[t], phi[t1] == 0}, {phi}, {t, t1, t2}][[1]],
      phiSoln = 
       NDSolve[{phi'[t] == omFn[t], phi[t0] == phi0}, {phi}, {t, t1, t2}][[1]]];
    phiFn = phi /. phiSoln;
    phiTb = Table[phiFn[t], {t, t1, t2, dt}];

    coords = {r -> rFn, om -> omFn, phi -> phiFn};
    vars = {x -> xFn, y -> yFn, u -> uFn, l -> lFn};
    waveform = EccentricWaveform[tTb, phiTb, rTb, rDotTb, omTb, ord];

    Return[Join[coords, vars, waveform]];
];

EccentricWaveform[tTb_List, phiTb_List, rTb_List, rDotTb_List, omTb_List, 
                  ord_Integer] :=
  Module[{t1, t2, dt, hTb, hFn, hPhase, hPhaseFn, hOmFn, hAmp, hAmpFn, psi4Tb,
          psi4Fn, psi4Phase, psi4PhaseFn, psi4DotTb, psi4OmTb, psi4OmFn},

    t1 = First[tTb];
    t2 = Last[tTb];
    dt = tTb[[2]] - tTb[[1]];
(* Print["Creating waveform with dt = ", dt]; *)
    hTb = MapThread[
       konigsh22 /. {phi -> #1, r -> #2, rDot -> #3, 
         phiDot -> #4} &, {phiTb, rTb, rDotTb, omTb}];
    hFn = Interpolation[MapThread[List, {tTb, hTb}], InterpolationOrder->ord];
    hPhase = Phase[MapThread[List, {tTb, hTb}]];
    hPhaseFn = Interpolation[hPhase, InterpolationOrder->ord];
    hOmFn = Derivative[1][hPhaseFn];
    hAmp = Amplitude[MapThread[List, {tTb, hTb}]];
    hAmpFn = Interpolation[hAmp, InterpolationOrder->ord];
    psi4Tb = 
      Table[psi4AmpPhaseExpr /. {hamp -> hAmpFn, hphase -> hPhaseFn,t -> tx}, {tx, t1, t2, 
        dt}];

    psi4Fn = Interpolation[MapThread[List, {tTb, psi4Tb}], InterpolationOrder->ord];
    psi4Phase = Phase[MapThread[List, {tTb, psi4Tb}]];
    psi4PhaseFn = Interpolation[psi4Phase, InterpolationOrder->ord];

    psi4DotTb = 
      Table[psi4DotAmpPhaseExpr /. {hamp -> hAmpFn, hphase -> hPhaseFn, 
        t -> tx}, {tx, t1, t2, dt}];
    psi4OmTb = MapThread[Im[#1/#2] &, {psi4DotTb, psi4Tb}];
    psi4OmFn = Interpolation[MapThread[List, {tTb, psi4OmTb}], 
      InterpolationOrder->ord];

    Return[{h -> hFn, hPhi -> hPhaseFn, psi4 -> psi4Fn, psi4Om -> psi4OmFn, 
            psi4Phi -> psi4PhaseFn, hOm -> hOmFn}];
  ];


End[];

EndPackage[];
