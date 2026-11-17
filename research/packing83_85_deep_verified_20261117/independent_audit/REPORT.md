# Independent 83D and 84D audit

Both stronger alternative candidates PASS. No publication was performed by this audit.

## Frozen inputs and results
83D `candidate_alt.json` SHA256 `975621f33449add5f757f5a75d8493bc38419cb842e98d5026b7fccca69e6205`.
84D `candidate84_alt.json` SHA256 `038a6a387e1b01d092948b9b38fb12114c8dd961956344e45bd88acf97fef5a4`.

83D vertical Gram I+J (3x3), determinant4; center density squared2^81. Safe ordinary lower bound3.08087078338732e-18.
84D vertical Gram I+J (4x4), determinant5; center density squared2^84/5. Safe ordinary lower bound1.06265246407636e-18.
256-bit Arb intervals appear in densities.json; density.py independently computes them.

## Exact transfer and geometry
check_alt.py and check84.py independently reconstruct each full80-coordinate isometry from its E8 action and block permutation. Integral unimodularity and A G A^t=G hold exactly. Every singleton/pair binary class in the alternative candidates is thereby identified with source case101, already independently published with depth>=6. No case504 premise is needed for either alternative.

The new remaining classes have integral half-coset norm: r G r is divisible4. The decoder inputs agree byte-for-byte in all fields with the audited parent decoder input except the target, whose 80 Euclidean numerators are independently rebuilt from r,M,E. A fresh local E8 shell was recursively generated and compared as a set with all9121 input vectors.

The vertical forms have quadratic value sum(z_i^2)+(sum z_i)^2. Exact enumeration of [-4,4]^k gives minima8 for nonzero even parity,2 for weight1/2 parities,4 for weight3/4 parities. Outside that box, V>=I implies value>=25, so no lower values are missed. The checker uses the weaker positive definite V-(1/2)I and bound25/2>8. Combining vertical minima with horizontal depths6 or4 proves total squared minimum>=8. Parent horizontal squared minimum>=8 remains the published Bachoc–Nebe premise tied to the frozen reconstruction.

## Fresh exhaustive depth replays and algorithm review
Compiled fresh replay_alt7 from producer depth4.cpp, reviewed nearest-E8 rounding and least-cost parity correction, translated complete local shell coverage, code-state recursion, optimistic branch bounds and final XOR carry DP. Nearest local squared distance<=1 is checked for each target. With cutoff384, required local radius²<=1.5, hence translated local shell norm<=6 is sufficient. Unknown-column bit masks overapproximate all completion choices; consequently pruning is only by a valid lower bound. No numerical arithmetic is used in the decoder. The final sums' denominator128 is bound by the independent input-map normalization audit.

83D mask7: fresh replay exhausted849936 nodes, zero leaves, all16 glue classes. Therefore no norm<=3; integral norms imply depth>=4.
84D masks7,11,13,14,15: fresh replays each exhausted all16 glue classes, respectively849936,849936,800784,849936,849936 nodes, zero leaves (see result84 files for binding). Again depth>=4 by integrality. Replays use the reviewed production decoder, not an independently reimplemented global algorithm; exact transfer, geometry, local shell generation and density are independently implemented.

For completeness the superseded first candidate's source504 was also freshly replayed:5877535660 nodes,100 leaves, all16glues, no vector<=5. Masks6/7 were replayed too. These are not dependencies of the stronger alternative results.

## Inherited dependencies
Frozen Bachoc–Nebe80 reconstruction/minimum theorem, and previously independently checked case101 depth>=6 and parent decoder input mapping. These were not silently replaced by a catalogue assumption. Worldwide priority is not established by this audit; compare current repository and external tables separately.

## 85D transfer added
85D also PASS. Frozen candidate SHA256 d2766829119c7632dd6aae0dd84e967ff7495faf25716c60173016eb43eb29b3. All15 singleton/pair parity classes transfer by independently checked full80 integer isometries to case101. Remaining16 parity classes: fresh exhaustive depth4 decoder replays total19,515,180 nodes, all16glues perclass completed, no vector<=3. Target half-coset norms are integral, so all have depth>=4. Weight5 vertical parity minimum6 (even stronger than required); weights3/4 minimum4. Exact vertical Gram I+J is5x5,det6. Same geometric tail argument and exact input checks apply. Density safe floor3.71893675292617e-19, verified by256bit Arb. Frozen input, check85.py/check85_result.json and allresult85/replay85 files included. The remaining sources and audit limitations are identical to83/84 above.

## Stronger 84D D4 vertical replacement
Also PASS. FrozenSHA32561174612588e520eb92903439295a69b126a8312febf394a3e58651cbe9a5. All12 masks except5,9,12 transfer by exact full80 isometry to source101. Fresh depth4 decoders for5,9,12 each exhausted849936nodes/all16glues/zero leaves; half-coset norms integral. Vertical D4 Gramdet4 checked positive definite; V−I/4 positive definite ensures outside[-5,5]^4 norm>8; finite exact search confirms minimum2 for12deepclasses,4for5/9/12,8fornonzeroevenclass. Centerdensity2^41, safe occupieddensity>=1.18808157306620e-18, computed with256bitArb. This supersedes earlier84A4det5 candidate. Checkcheck84D4.py and resultD4 files.

## Stronger 85D D5 vertical replacement
Also PASS. FrozenSHA57c865b443726e914de1d608958712e294a9920e3a16dd9dc4e312b130dbec5a. All20 advertised depth6 classes transfer by full exact isometry to case101. Fresh independent replays for11remainingmasks5,9,12,18,20,21,22,24,25,26,30 exhaust9,255,088nodes acrossall16glues perclass, no vector<=3; allnorms integral. VerticalD5Gramdet4; V−I/8positive definite and box[-8,8]^5exactly enumerated, outside norm>81/8>8. Minima2 on20deepclasses,4onremaining11classes,8onzero nonzerovectorclass. Safe density>=4.55474871517601e-19 by256bitArb. This supersedes85A5det6. Independentcheckercheck85D5.py and resultD5 files, frozeninput included.
