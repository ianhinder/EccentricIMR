
BeginPackage["EccentricIMR`Kappas`", 
  {"EccentricIMR`EccentricPNSymbols`"}];

kappaE;
kappaJ;

$KappaECacheFile;
$KappaJCacheFile;

$KappasFromRepository; (* Set this to False to regenerate the cache;
                          you also need to set the cache file
                          variables *)

Begin["`Private`"];

$dataDirectory = FileNameJoin[{FileNameDrop[FindFile["EccentricIMR`Kappas`"],-1],"Data"}];
$kappaERepoCacheFile = FileNameJoin[{$dataDirectory,"kappaE.m"}];
$kappaJRepoCacheFile = FileNameJoin[{$dataDirectory,"kappaJ.m"}];

If[FileExistsQ[$kappaERepoCacheFile] =!= True, Print["ERROR: ", $kappaERepoCacheFile, " not found"]; Abort[]];
If[FileExistsQ[$kappaJRepoCacheFile] =!= True, Print["ERROR: ", $kappaERepoCacheFile, " not found"]; Abort[]];

BesselJprime[p_, x_] := Module[{y, d}, d = D[BesselJ[p, y], y]; d /. y -> x];


(* I have checked that the error made by using the interpolated values
   of kappaE or kappaJ from the table instead of the computed values
   is of the order of 10^-12 when the table spacing is 0.001.  This
   was measured by using a spacing ten times smaller and measuring the
   largest difference over the range [0,0.4]. *)

computeKappaE[ep_] :=
  Module[{p, kappaE, e, term, sumPrev, sum, n, ser},

(*    If[ep < 0 || ep >= 1,
      Throw["computeKappaE: range error in e"]]; *)

    term[e_, p_] := p^3/4(BesselJ[p, p e]^2(1/e^4 - 1/e^2 + 1/3 + 
                    p^2(1/e^4 - 3/e^2 + 3 - e^2)) + 
              p(-4/e^3 + 7/e - 3e)BesselJ[p, p e] BesselJprime[p, p e] + 
              BesselJprime[p, p e]^2(1/e^2 - 1 + p^2(1/e^2 - 2 + e^2)));

    If[Abs[ep] > 10^-8,

      sumPrev = 0;
      sum = term[ep, 1];
      n = 2;

      While[Abs[sum - sumPrev] > 10^-15,
        sumPrev = sum;
        sum = sum + term[ep, n];
        n = n + 1;
(*        Print["n = ", n, "  sum = ", sum // FullForm];*)];
      Return[sum],

(*      Print["Using Series"];*)
      kappaE = Sum[term[e,p],{p,1,5}];
(*      Print["kappaE = ", kappaE];*)
      ser = Normal[Series[kappaE,{e,0,4}]];
(*      Print[ser];*)
      Return[(ser /. e->ep)//N]]];

generateKappaETable[] :=
  Table[{e, computeKappaE[e]}, {e, -0.8, 0.8, 0.001}];

saveKappaETable[file_] :=
  Export[file, generateKappaETable[], "TSV"];

loadKappaETable[file_] :=
  kappaEFn = Interpolation[ReadList[file, Real, RecordLists->True], InterpolationOrder->4];

generateKappaEFn[] :=
  Module[{},
    Print["Computing kappaE"];
    Print[Timing[kappaEFn = Interpolation[generateKappaETable[], InterpolationOrder->4]][[1]]]];

(*generateKappaEFn[];*)

kappaE[e_?NumberQ] := 
(
  If[$KappasFromRepository =!= False,
    kappaEFn = Get[$kappaERepoCacheFile],
    (* else *)
    If[StringQ[$KappaECacheFile] && FileExistsQ[$KappaECacheFile],
      kappaEFn = Get[$KappaECacheFile],
      (* else *)
      generateKappaEFn[];
      If[StringQ[$KappaECacheFile],
        Put[kappaEFn, $KappaECacheFile]]]];
  
  kappaE[ee_?NumberQ] :=
  kappaEFn[ee];
  kappaE[e]);




computeKappaJ[ep_] :=
  Module[{p, kappaJ, e, term, sumPrev, sum, n, ser},
(*    If[ep < 0 || ep >= 1,
      Throw["computeKappaJ: range error in e"]];*)

    term[e_, p_] := p^2/2 Sqrt[
            1 - e^2](p(3/e^2 - 2/e^4 - 1)BesselJ[p, p e]^2 + (2/e^3 - 1/e + 
                    2p^2(1/e^3 - 2/e + e))BesselJ[p, p e]BesselJprime[p, p e] + 
              2p(1 - 1/e^2)BesselJprime[p, p e]^2);

    If[Abs[ep] > 10^-8,

      sumPrev = 0;
      sum = term[ep, 1];
      n = 2;

      While[Abs[sum - sumPrev] > 10^-15,
        sumPrev = sum;
        sum = sum + term[ep, n];
        n = n + 1;
(*        Print["n = ", n, "  sum = ", sum // FullForm];*)];
      Return[sum],
(*      Print["Using Series"];*)
      kappaJ = Sum[term[e,p],{p,1,5}];
(*      Print["kappaJ = ", kappaJ];*)
      ser = Normal[Series[kappaJ,{e,0,4}]];
      Return[(ser /. e->ep)//N]]];

generateKappaJTable[] :=
  Table[{e, computeKappaJ[e]}, {e, -0.8, 0.8, 0.001}];

generateKappaJFn[] :=
  Module[{},
    Print["Computing kappaJ"];
    Print[Timing[kappaJFn = Interpolation[generateKappaJTable[], InterpolationOrder->4]][[1]]]];

(*generateKappaJFn[];*)

kappaJ[e_?NumberQ] := 
(
  (* TODO: eliminate duplication with kappaE above *)
  If[$KappasFromRepository =!= False,
    kappaJFn = Get[$kappaJRepoCacheFile],
    (* else *)
    If[StringQ[$KappaJCacheFile] && FileExistsQ[$KappaJCacheFile],
      kappaJFn = Get[$KappaJCacheFile],
      (* else *)
      generateKappaJFn[];
      If[StringQ[$KappaJCacheFile],
        Put[kappaJFn, $KappaJCacheFile]]]];
  
  kappaJ[ee_?NumberQ] :=
  kappaJFn[ee];
  kappaJ[e]);

(* generateKappaEFn[]; *)

End[];

EndPackage[];
