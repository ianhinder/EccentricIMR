BeginPackage["EccentricIMR`"];

$EccentricIMRCalibrationWarnings;

EndPackage[];

If[$VersionNumber < 10.,
  Print["The EccentricIMR package requires Mathematica version 10 or greater"];
  Abort[]];

Get["EccentricIMR`EccentricIMR`"];
