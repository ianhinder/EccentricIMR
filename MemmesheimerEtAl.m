
BeginPackage["EccentricIMR`MemmesheimerEtAl`", {"EccentricIMR`EccentricPNSymbols`"}];

ephSquaredInEpsj;
eSquaredInEpsj;

Begin["`Private`"];

ephSquaredInEpsj =
  Module[{ephSquared},
    ephSquared = 1 - j + (-2En)/4(24 + j(-15 + eta)) + (-2En)^2/
              16(-40 + 34eta + 18eta^2 - j(160 - 31eta + 3eta^2) - 
                1/j(-416 + 91eta + 15eta^2)) + (-2En)^3/
              13440(-584640 - 17482eta - 4305Pi^2eta - 7350eta^2 + 
                8190eta^3 - 420(j)(744 - 248eta + 31eta^2 + 3eta^3) - 
                14/(j)(36960 - 341012eta + 4305Pi^2eta - 225eta^2 + 
                      150eta^3) - 
                1/(j)^2(-2956800 + 5627206eta - 81795Pi^2eta + 
                      14490eta^2 + 7350eta^3)) + O[En]^4;

    ComposeSeries[ephSquared, -Eps/2 + O[Eps]^100]
];

(* This does not come directly from some paper.  We should derive it again. *)

eSquaredInEpsj = ComposeSeries[1 - j + (-2 En)/4 (-8+8eta-j(-17+7eta))+
(  (-2En)^2/8(12+72eta+20eta^2-24Sqrt[j](-5+2eta)-j(112-47eta+16eta^2)-16/j(-4+7eta)+24/Sqrt[j](-5+2eta)) + (-2En)^3/6720 (23520-464800eta+179760eta^2+16800eta^3-2520Sqrt[j](265-193eta+46eta^2)-525 j (-528+200eta-77eta^2+24eta^3)-6/j(73920-260272eta+4305Pi^2eta+61040eta^2)+70/Sqrt[j](16380-19964eta+123Pi^2eta+3240eta^2)+8/j^2(53760-176024eta+4305Pi^2eta+15120eta^2)-70/j^(3/2)(10080-13952eta+123Pi^2eta+1440eta^2)) +O[En]^4), -Eps/2 + O[Eps]^100];

End[];

EndPackage[];
