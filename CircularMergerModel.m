
BeginPackage["EccentricIMR`CircularMergerModel`"];

CircularWaveformModel;
CircularWaveform;
EccentricIMR`$EccentricIMRCalibrationWarnings;

Begin["`Private`"];

$CircularMergerModelDirectory = FileNameDrop[FindFile["EccentricIMR`CircularMergerModel`"],-1];

ReadCircularWaveformModel[] :=
  Get[FileNameJoin[{$CircularMergerModelDirectory,"Data","circular-model.m"}]];

CircularWaveformModel[] :=
  CircularWaveformModel[] =
  ReadCircularWaveformModel[];

etaOfq[q_] :=
 q/(1 + q)^2;

qOfeta[eta_] :=
  (1 + Sqrt[1 - 4 eta] - 2 eta)/(2 eta);

antiDerivative[d_List, {tbc_, fbc_}] :=
  Module[{tMin, tMax, dFn, gFn, g, t, dt, gTb},
    {tMin, tMax} = {d[[1,1]], d[[-1,1]]};
    If[tbc < tMin || tbc > tMax,
      Error["antiDerivative: boundary condition is not within range of data"]];
    dFn = Interpolation[d];
    gFn = g /.
    NDSolve[{D[g[t], t] == dFn[t], g[tbc] == fbc}, {g}, {t, tMin, tMax}, MaxSteps -> 1000000][[
      1]];

    gTb = Thread[{d[[All,1]], gFn/@d[[All,1]]}]];

waveformOfAOm[A_List, om_List, phi0_?NumericQ] :=
   Module[{phi},
     phi = antiDerivative[om, {A[[1,1]], phi0}];
     Thread[{A[[All,1]], A[[All,2]] Exp[-I phi[[All,2]]]}]];

waveformOfAOm[model_Association, q_, phi0_] :=
   waveformOfAOm[model["AmplitudeFunction"][q], 
     model["FrequencyFunction"][q], phi0];

CircularWaveform[q_?NumericQ, phi0: (_?NumericQ) : 0.] :=
   CircularWaveform[CircularWaveformModel[], q, phi0];

CircularWaveform[model_Association, q_?NumericQ, phi0_?NumericQ] :=
   Module[{AOfqIneta, omOfqIneta, etaRange, qRange, eps=10.^-3},
     etaRange = model["AmplitudeFunctions"][[1]][[1,1]];
     qRange = Reverse[qOfeta/@etaRange];

     If[EccentricIMR`$EccentricIMRCalibrationWarnings =!= False &&
       (q < qRange[[1]]-eps || q > qRange[[2]]+eps),
       Print["WARNING: Circular waveform model evaluated at q=", q, " which is outside its calibration range ", qRange, ".  Set $EccentricIMRCalibrationWarnings to False to disable this warning."]];

     AOfqIneta =
       MapThread[{#1, 
         Quiet[#2[etaOfq[q]], InterpolatingFunction::dmval]} &, {model[
           "Times"], model["AmplitudeFunctions"]}];
     omOfqIneta =
       MapThread[{#1, 
         Quiet[#2[etaOfq[q]], InterpolatingFunction::dmval]} &, {model[
           "Times"], model["FrequencyFunctions"]}];
     waveformOfAOm[AOfqIneta, omOfqIneta, phi0]];

End[];

EndPackage[];

