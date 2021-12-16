(* ::Package:: *)

(* ::Section:: *)
(*Code Constants, Libraries, Coordinate Grid and Integrals*)


(* ::Subsection:: *)
(*Code Constants: s0=-1, precgrid=350, prec=300. *)


(* 
s0: is the center of the conformal map; 
precgrid: is the intermediate precision when performing Mathematica operations;
prec: is the final precision in which the output file is passed to sdpb;
*)


Unprotect[s0];ClearAll[s0];s0=-1;Protect[s0];
Unprotect[precgrid,prec]; precgrid= 350; prec=300;Protect[precgrid,prec];


(* ::Subsection:: *)
(*Load SDPB Library (UPDATE THE FOLDER)*)


SetDirectory["/gpfs/hmurali/string_bootstrap"];
<< "/home/hmurali/string_bootstrap/sdpbPkg.wl";


(* ::Subsection:: *)
(*Coordinates*)


(* ::Subsubsection:: *)
(*Basic Variables*)


(* 
t,u: Mandelstam invariants in terms of energy and angles;
\[Rho]v[s]: is the conformal map; 
\[ScriptS][\[Rho]]: is the inverse map;
\[Rho]vExp[x]: is the expansion of the conformal map around x=0 up to order x^10;
*)


ClearAll[t,u]
t[s_,x_]:=-(s/2)(1-x); u[s_,x_]:=-(s/2)(1+x);

ClearAll[\[Rho]v,\[ScriptS]];
 \[Rho]v[s_]:=(Sqrt[-s0  ]-Sqrt[-I] Sqrt[-s I])/(Sqrt[-s0 ]+Sqrt[-I] Sqrt[-s I]);
 \[ScriptS][\[Rho]_]:=(s0 (-1+\[Rho])^2)/(1+\[Rho])^2 ;
 
 ClearAll[\[Rho]vExp];
\[Rho]vExp[x_]:=1+2 I Sqrt[x]-2 x-2 I x^(3/2)+2 x^2+2 I x^(5/2)-2 x^3-2 I x^(7/2)+2 x^4+2 I x^(9/2)-2 x^5-2 I x^(11/2)+2 x^6+2 I x^(13/2)-2 x^7-2 I x^(15/2)+2 x^8+2 I x^(17/2)-2 x^9-2 I x^(19/2)+2 x^10;


(* ::Subsubsection:: *)
(*Basis of monomials for the crossing symmetric ansatz*)


(* 
BasisMandelstamAll: lists all the elements in the totally symmetric triple series ansatz in terms of monomials of rho variables mon[a,b,c]\[Equal]\[Rho]s^a \[Rho]t^b \[Rho]u^c;
BasisMandelstam: contains the gauged fixed ansatz, where we removed all the redundant monomials related by the onshellness condition s+t+u=0; 
*)


Unprotect[mon];ClearAll[mon];Protect[mon];

ClearAll[BasisMandelstamAll]
BasisMandelstamAll[nmax_]:=BasisMandelstamAll[nmax]=Table[(mon@@@Permutations[{a,b,c}]//Total)Boole[a+b+c<= nmax],{a,0,nmax},{b,a,nmax},{c,b,nmax}]//Flatten//DeleteCases[#,0]&;

ClearAll[BasisMandelstam];
BasisMandelstam[nmax_]:=BasisMandelstam[nmax]=BasisMandelstamAll[nmax]/.mon[x_,y_,z_]/;x-2>=0&&y-2>= 0&&z-2>=0->0//DeleteCases[#,0]&;


(* ::Subsubsection:: *)
(*Unitarity Grid*)


(* 
numPTS=1000: number of MAX unitarity points >> the actual points used;
\[Phi]gridCheb: Chebyschev grid in the angle \[Rho]=Exp[\[ImaginaryI] \[Phi]]; 
\[Rho]gridCheb: \[Phi]gridCheb grid in \[Rho];
gridIndexCheb: position of an element of \[Rho]gridCheb in the list;
\[Rho]gridCheb2: 200 grid points contained in \[Rho]gridCheb actually used to impose unitarity;
*)


ClearAll[\[Phi]gridCheb,\[Rho]gridCheb,gridIndexCheb];

Unprotect[numPTS];ClearAll[numPTS]; 
numPTS=1000;
Protect[numPTS];

\[Phi]gridCheb=SetPrecision[\[Pi]/2 (1+Cos[\[Pi] Range[numPTS+1]/(numPTS+1)])//Drop[#,-1]&,precgrid]//Reverse;
\[Rho]gridCheb=SetPrecision[(Exp[I \[Phi]gridCheb]),precgrid];
sgridCheb=Tan[\[Phi]gridCheb/2]^2;
dsgridCheb=SetPrecision[Tan[\[Phi]gridCheb/2]Sec[\[Phi]gridCheb/2]^2 \[Pi]^2/(2numPTS) (Sin[\[Pi] Range[numPTS+1]/(numPTS+1)]//Drop[#,-1]&//Reverse), precgrid];
gridIndexCheb={SetPrecision[\[Rho]gridCheb,precgrid],Range[1,Length[\[Rho]gridCheb]]}//Transpose;
constraintsgrid=Join[Range[1,numPTS,4],Range[numPTS-100,numPTS,1]]//DeleteDuplicates//Sort;

ClearAll[\[Rho]gridCheb2];
\[Rho]gridCheb2=SetPrecision[(Exp[I \[Phi]gridCheb]),precgrid]//Table[#[[i]],{i,1,Length[#],5}]&;


(* ::Section:: *)
(*Building the Ansatz *)


(* ::Subsection:: *)
(*Ansatz Definition*)


(* 
\[Rho]mon[s,t,u][a,b,c]: monomials \[Rho]s^a \[Rho]t^b \[Rho]u^c in terms of s,t,u variables;
\[Alpha]abcAll[Nmax]: list of free variables in the ansatz ac[a,b,c] multiplying the corrsponding symmetrized monomial;
\[Alpha]abc[Nmax]: gauge fixed free variables;
RhoAmplitude[s,t,u][Nmax]: ANSATZ as a function of s,t,u for the stripped amplitude;
Amplitude[s,t,u][Nmax]: ANSATZ + low energy term 1/stu for the axidilaton component;
*)


ClearAll[\[Rho]mon];
\[Rho]mon[s_,t_,u_][a_List]:=\[Rho]v[s]^a[[1]] \[Rho]v[t]^a[[2]] \[Rho]v[u]^a[[3]];

Unprotect[ac];ClearAll[ac];Protect[ac];
ClearAll[\[Alpha]abcAll];
\[Alpha]abcAll[maxdegree_]:=
Table[ac[a,b,c]Boole[a+b+c<= maxdegree],{a,0,maxdegree},{b,a,maxdegree},{c,b,maxdegree}]//Flatten//DeleteCases[#,0]&;

ClearAll[\[Alpha]abc];
\[Alpha]abc[maxdegree_]:=
{\[DoubleStruckCapitalG]N}~Join~(\[Alpha]abcAll[maxdegree]/.ac[x_,y_,z_]/;x-2>=0&&y-2>=0&&z-2>=0->0//DeleteCases[#,0]&);

ClearAll[RhoAmplitude];
RhoAmplitude[s_,t_,u_][maxdegree_]:=(1+\[Rho]v[s])^2 (1+\[Rho]v[t])^2 (1+\[Rho]v[u])^2 (BasisMandelstam[maxdegree]/.{mon[x_,y_,z_]:> \[Rho]mon[s,t,u][{x,y,z}]});

ClearAll[Amplitude];
Amplitude[s_,t_,u_][maxdegree_]:=s^4 ({1/(s t u)}~Join~RhoAmplitude[s,t,u][maxdegree]);


(* 
THIS FUNCTIONS MAKE THE EVALUATION OF THE AMPLITUDE FASTER

\[Rho]term[s,x][a,b,c]: PRE-SIMPLIFIED MONOMIALS \[Rho]s^a \[Rho]t^b \[Rho]u^c in terms of s and x variables;
RhoAmplitudeOpt[s,x][Nmax]: Optimized ANSATZ as a function of s and x;
AmplitudeOpt[s,x][Nmax]: Optimized ANSATZ + low energy term 1/stu;
*)


ClearAll[\[Rho]term];
\[Rho]term[s_,x_][a_,b_,c_]:=-((512 I ((1+2 (-1)^(3/4) Sqrt[-I s]-s)/(1+s))^a ((2+s-(2-2 I) Sqrt[-I s (-1+x)]-s x)/(2+s (-1+x)))^b (((1+I)-Sqrt[I s (1+x)])/((1+I)+Sqrt[I s (1+x)]))^c)/(((1+I)+Sqrt[2] Sqrt[-I s])^2 ((1+I)+Sqrt[-I s (-1+x)])^2 ((1+I)+Sqrt[I s (1+x)])^2));
ClearAll[RhoAmplitudeOpt];
RhoAmplitudeOpt[s_,x_][maxdegree_]:= (BasisMandelstam[maxdegree]/.{mon[a_,b_,c_]:> \[Rho]term[s,x][a,b,c]});

ClearAll[AmplitudeOpt];
AmplitudeOpt[s_,x_][maxdegree_]:=s^4 ({1/(s t[s,x] u[s,x])}~Join~RhoAmplitudeOpt[s,x][maxdegree]);


(* ::Subsection:: *)
(*Linear Constraints for high energy convergence (UPDATE THE FOLDER)*)


(* 
LARGE s EXPANSION OF THE PARTIAL WAVE PROJECTIONS OF THE MONOMIAL \[Rho]s^a \[Rho]t^b \[Rho]u^c

\[Mu]10: 10-dimensional measure proejecting onto partial waves;
int[k,n][b,c]: large s expansion of the integrals Integrate[\[Rho]t^b \[Rho]u^c (1-z)^n,{z,0,1}] at order s^{-k} for the terms that require renormalization;
\[Rho]the,\[Rho]uhe: large s, fixed angle expansion of the conformal maps for t,u;
\[Rho]bcexp: lists of the expansion coefficients of \[Rho]t^b \[Rho]u^c at large s up to s^{-7};
\[ScriptA]int[j,n][b,c]: regular large s expansion at order s^{-j} of \[Rho]t^b \[Rho]u^c (1-z)^n;
d[k,a]: large s expansion of \[Rho]s^a at order s^{-k};
\[Mu]: coefficient list of \[Mu]10 in powers of (1-z)^n; 

\[ScriptE][j,n][b,c]: Full analytic part of the expansion of Integrate[\[Rho]t^b \[Rho]u^c (1-z)^n,{z,0,1}] at order s^{-j} for s large;
\[ScriptG][a,b,c][\[ScriptL],j]: Full analytic part of the expansion of the projection of \[Rho]s^a \[Rho]t^b \[Rho]u^c for spin \[ScriptL] at order s^{-j};
\[ScriptF][j,n][b,c]: Full non-analytic part of the expansion of Integrate[\[Rho]t^b \[Rho]u^c (1-z)^n,{z,0,1}] at order Log[s]s^{-j};
\[ScriptH][a,b,c][\[ScriptL],j]: Full non-analytic part of the expansion of the projection of \[Rho]s^a \[Rho]t^b \[Rho]u^c for spin \[ScriptL] at order Log[s] s^{-j};
*)


$MaxExtraPrecision=10000;

ClearAll[\[Mu]10];
\[Mu]10[\[ScriptL]_,x_]:=(1-x^2)^3 GegenbauerC[\[ScriptL],7/2,x]/GegenbauerC[\[ScriptL],7/2,1];


ClearAll[int];
Do[Do[int[k,n][b_,c_]=ToExpression[Import["/gpfs/aguirrieri/integral_storage_SUGRA2/LargeEnergyIntegrals2/aintegrals_"<>ToString[2k]<>"_"<>ToString[n]<>".txt"]],{k,0,7,1/2}],{n,3,40,1}];


ClearAll[\[Rho]the,\[Rho]uhe];
\[Rho]the=\[Rho]v[t[s,x]]/.x->1-\[ScriptK]//Series[#,{s,\[Infinity],10},Assumptions->s>0&&\[ScriptK]>0]&//Simplify[#,Assumptions->s>0&&\[ScriptK]>0]&//FullSimplify[#,Assumptions->s>0&&\[ScriptK]>0]&//Normal;
\[Rho]uhe=\[Rho]v[u[s,x]]/.x->-1+\[ScriptJ]//Series[#,{s,\[Infinity],10},Assumptions->s>0&&\[ScriptJ]>0]&//Simplify[#,Assumptions->s>0&&\[ScriptJ]>0]&//FullSimplify[#,Assumptions->s>0&&\[ScriptJ]>0]&//Normal;

ClearAll[\[Rho]bcexp];
\[Rho]bcexp=Table[(\[Rho]the^b \[Rho]uhe^c//SeriesCoefficient[#,{s,\[Infinity],kk},Assumptions->s>0&&b\[Element] Integers&&b>0]&) 1/s^kk//Simplify,{kk,0,7,1/2}];


ClearAll[\[ScriptA]int];
Do[\[ScriptA]int[j,n_][b_,c_]=\[Rho]bcexp[[2 j+1]]\[ScriptK]^n//Expand//#/.\[ScriptK]^x_ \[ScriptJ]^y_:>2^(x+1+y) Beta[1/2,x+1,y+1]&//#/.\[ScriptK]^x_:> 1/(1+x)&,{j,0,7,1/2}];


ClearAll[d];
d[k_,a_]:=d[k,a]=SeriesCoefficient[\[Rho]v[s]^a,{s,\[Infinity],k/2},Assumptions->s>0]//Simplify[#,Assumptions->s>0]&;


ClearAll[\[Mu]];
\[Mu][n_,\[ScriptL]_]:=\[Mu][n,\[ScriptL]]=SeriesCoefficient[\[Mu]10[\[ScriptL],1-\[ScriptK]],{\[ScriptK],0,n}];


ClearAll[\[ScriptE]];
\[ScriptE][j_,n_][b_,c_]/;(n<=6 &&b<=c):=\[ScriptE][j,n][b,c]=int[j/2,n][b,c]+int[j/2,n][c,b]//Coefficient[#,Log[s],0]&//#/.s->1&//N[#,{350,350}]&//Rationalize;
\[ScriptE][j_,n_][b_,c_]/;(n<=6 &&b>c):=\[ScriptE][j,n][c,b];
\[ScriptE][j_,n_][b_,c_]/;(n>=7&&b<=c):=\[ScriptE][j,n][b,c]=\[ScriptA]int[j/2,n][b,c]+\[ScriptA]int[j/2,n][c,b]//#/.s->1&//N[#,{350,350}]&//Rationalize;
\[ScriptE][j_,n_][b_,c_]/;(n>=7 &&b>c):=\[ScriptE][j,n][c,b];

ClearAll[\[ScriptG]];
\[ScriptG][a_,b_,c_][\[ScriptL]_,j_]:=\[ScriptG][a,b,c][\[ScriptL],j]=Sum[Sum[d[k,a]\[Mu][n,\[ScriptL]]\[ScriptE][j-k,n][b,c],{k,0,j}],{n,3,6+\[ScriptL]}];


ClearAll[\[ScriptF]]
\[ScriptF][j_,n_][b_,c_]/;(n<=6 &&b<=c):=\[ScriptF][j,n][b,c]=int[j/2,n][b,c]+int[j/2,n][c,b]//Coefficient[#,Log[s],1]&//#/.s->1&//N[#,{350,350}]&//Rationalize;
\[ScriptF][j_,n_][b_,c_]/;(n<=6 &&b>c):=\[ScriptF][j,n][c,b];
\[ScriptF][j_,n_][b_,c_]/;n>=7:=0;

ClearAll[\[ScriptH]];
\[ScriptH][a_,b_,c_][\[ScriptL]_,j_]:=\[ScriptH][a,b,c][\[ScriptL],j]=Sum[Sum[d[k,a]\[Mu][n,\[ScriptL]]\[ScriptF][j-k,n][b,c],{k,0,j}],{n,3,6+\[ScriptL]}];


(* 
SOLVING THE HIGH ENERGY LINEAR CONSTRAINTS

PreExpansion[Nmax]: actual monomials rrr[a,b,c]\[Equal]\[Rho]s^a \[Rho]t^b \[Rho]u^c contained in each term of the ansatz
HEConstraintOpt[k][Nmax,Spin]: real and imaginary part of the Spin projection of the ansatz at order s^{-k};

(THIS CAN BE WRITTEN MORE COMPACTLY, I LEFT IT IN THIS WAY IN MY VERSION OF THE CODE 
SINCE I CHECKED THE EXPANSION ORDER BY ORDER TO MAKE SURE THERE WERE NO BUGS)

HEConstraintModuleOpt[Nmax]: solves the linear constraints, the resulting amplitude decays as 1/s^7 for any spin;

(IT CAN BE WRITTEN IN A CLEANER WAY)

HEConditionOpt[7][Nmax,Spin]: leading order of the ansatz for large s, it will be used to impose unitarity at infinity;
*)


ClearAll[PreExpansion];
PreExpansion[nmax_]:=PreExpansion[nmax]=\[Rho]s^\[Alpha] \[Rho]t^\[Beta] \[Rho]u^\[Gamma] (1+\[Rho]s)^2 (1+\[Rho]t)^2 (1+\[Rho]u)^2 (BasisMandelstam[nmax]/.{mon[x_,y_,z_]:> \[Rho]s^x \[Rho]t^y \[Rho]u^z})//ExpandAll//#/.\[Rho]s^a_ \[Rho]t^b_ \[Rho]u^c_-> rrr[a,b,c]&//#/.\[Alpha]->0/.\[Beta]->0/.\[Gamma]->0&;


ClearAll[HEConstraintOpt];
HEConstraintOpt[3][nmax_,\[ScriptL]_]:=
HEConstraintOpt[3][nmax,\[ScriptL]]=Join[{s^3 Integrate[1/(s t[s,x]u[s,x]) (1-x^2)^3 GegenbauerC[\[ScriptL],7/2,x]/GegenbauerC[\[ScriptL],7/2,1],{x,-1,1}]},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],6]]//Rationalize;

HEConstraintOpt[7/2][nmax_,\[ScriptL]_]:=
HEConstraintOpt[7/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],7]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[4][nmax_,\[ScriptL]_]:=
HEConstraintOpt[4][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],8]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[4][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[4][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],8]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[9/2][nmax_,\[ScriptL]_]:=
HEConstraintOpt[9/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],9]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[9/2][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[9/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],9]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[5][nmax_,\[ScriptL]_]:=
HEConstraintOpt[5][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],10]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[5][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[5][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],10]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[11/2][nmax_,\[ScriptL]_]:=
HEConstraintOpt[11/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],11]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[11/2][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[11/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],11]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[6][nmax_,\[ScriptL]_]:=
HEConstraintOpt[6][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],12]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[6][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[6][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],12]]//Rationalize//{#//Re,#//Im}&;

HEConstraintOpt[13/2][nmax_,\[ScriptL]_]:=
HEConstraintOpt[13/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],13]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[13/2][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[13/2][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],13]]//Rationalize//{#//Re,#//Im}&;

LogHEConstraintOpt[7][nmax_,\[ScriptL]_]:=
LogHEConstraintOpt[7][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptH][a,b,c][\[ScriptL],14]]//Rationalize//{#//Re,#//Im}&;


ClearAll[HEConstraintModuleOpt];
HEConstraintModuleOpt[nnmax_]:=HEConstraintModuleOpt[nnmax]=Block[{vrs,soll1,soll2,soll3,soll4,soll5,soll6,soll7,soll8,soll9},
vrs=\[Alpha]abc[nnmax]//Drop[#,1]&;
soll1=HEConstraintOpt[3][nnmax,0] . \[Alpha]abc[nnmax]==0//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll2=(#==0)&/@(HEConstraintOpt[7/2][nnmax,0] . \[Alpha]abc[nnmax])/.soll1//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll3=(#==0)&/@(HEConstraintOpt[4][nnmax,0] . \[Alpha]abc[nnmax])/.soll1/.soll2//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll4={(#==0)&/@(HEConstraintOpt[9/2][nnmax,0] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[9/2][nnmax,2] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3&//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll5={(#==0)&/@(HEConstraintOpt[5][nnmax,0] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[5][nnmax,2] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3/.soll4&//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll6={(#==0)&/@(HEConstraintOpt[11/2][nnmax,0] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[11/2][nnmax,2] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3/.soll4/.soll5&//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll7={(#==0)&/@(HEConstraintOpt[6][nnmax,0] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[6][nnmax,2] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[6][nnmax,4] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3/.soll4/.soll5/.soll6&//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll8={(#==0)&/@(HEConstraintOpt[13/2][nnmax,0] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[13/2][nnmax,2] . \[Alpha]abc[nnmax]),(#==0)&/@(HEConstraintOpt[13/2][nnmax,4] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3/.soll4/.soll5/.soll6/.soll7&//Expand//NSolve[#,vrs,WorkingPrecision -> 300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
soll9={(#==0)&/@(LogHEConstraintOpt[7][nnmax,0] . \[Alpha]abc[nnmax])}//Flatten//#/.soll1/.soll2/.soll3/.soll4/.soll5/.soll6/.soll7/.soll8&//Expand//NSolve[#,vrs,WorkingPrecision->300]&//First//Chop[#,10^-250]&//Rationalize//Quiet;
Table[Rule[\[Alpha]abc[nnmax][[i]],\[Alpha]abc[nnmax][[i]]/.soll1/.soll2/.soll3/.soll4/.soll5/.soll6/.soll7/.soll8/.soll9//Expand],{i,1,Length[\[Alpha]abc[nnmax]]}]];


ClearAll[HEConditionOpt];
HEConditionOpt[7][nmax_,\[ScriptL]_]:=
HEConditionOpt[7][nmax,\[ScriptL]]=Join[{0},PreExpansion[nmax]/.rrr[a_,b_,c_]:>\[ScriptG][a,b,c][\[ScriptL],14]]//Rationalize;


(* ::Subsection:: *)
(*High Spin Constraints *)


(* 
numHSC=4000: number of points where we impose higher spin constraints;
\[Phi]gridHSC: numHSC Chebyschev points in \[Phi] between 0 and \[Pi];
\[Rho]gridHSC: same grid in \[Rho];

\[Rho]gridHSCReal: \[Rho][-s] for s evaluated \[Rho] in \[Rho]gridHSC;
premon: prefactor (1+\[Rho]s)^2 (1+\[Rho]t)^2(1+\[Rho]u)^2 in terms of monomials pmon[a,b,c] of \[Rho]s,\[Rho]t,\[Rho]u;
RhoAmplitudeNR2[Nmax]: Full ansatz elements in terms of monomilas of \[Rho]s^a \[Rho]t^b \[Rho]u^c \[Equal] \[Rho]mon2[a,b,c];

\[Rho]mon2Num[x,y,z]: large spin Im [\[Rho]s^a \[Rho]t^b \[Rho]u^c] numerically evaluated for \[Rho]gridHSC;
HSConstraintOPT[Nmax]: Large Spin Imaginary part of the full amplitude;
*)


ClearAll[numHSC,\[Phi]gridHSC,\[Rho]gridHSC,\[Rho]gridHSCReal]
numHSC=4000;
\[Phi]gridHSC=SetPrecision[\[Pi]/2 (1+Cos[\[Pi] Range[numHSC+1]/(numHSC+1)])//Drop[#,-1]&,precgrid]//Reverse;
\[Rho]gridHSC=SetPrecision[(Exp[I \[Phi]gridHSC]),precgrid];
dsgridHSC=SetPrecision[Tan[\[Phi]gridHSC/2]Sec[\[Phi]gridHSC/2]^2 \[Pi]^2/(2numHSC) (Sin[\[Pi] Range[numHSC+1]/(numHSC+1)]//Drop[#,-1]&//Reverse), precgrid];
sgridHSC=Tan[\[Phi]gridHSC/2]^2;
\[Rho]gridHSCReal=\[Rho]v[-\[ScriptS][\[Rho]gridHSC]]//Re;


ClearAll[premon];
premon=(1+\[Rho]s)^2 (1+\[Rho]t)^2 (1+\[Rho]u)^2 \[Rho]s^a \[Rho]t^b \[Rho]u^c//ExpandAll//#/.\[Rho]s^a_ \[Rho]t^b_ \[Rho]u^c_:> pmon[a,b,c]&//#/.a->0/.b->0/.c->0&;

ClearAll[RhoAmplitudeNR2];
RhoAmplitudeNR2[maxdegree_]:=Join[{\[Rho]mon2[{0,0,0}]},premon BasisMandelstam[maxdegree]//ExpandAll//#/.mon[x_,y_,z_]pmon[a_,b_,c_]:> \[Rho]mon2[{a+x,b+y,c+z}]&,{\[Rho]mon2[{0,0,0}]}];

ClearAll[\[Rho]mon2Num];
\[Rho]mon2Num[x_,y_,z_]:=\[Rho]mon2Num[x,y,z]=y \[Rho]gridHSC^x \[Rho]gridHSCReal^z;

ClearAll[HSConstraintOPT];
HSConstraintOPT[maxdegree_]:=HSConstraintOPT[maxdegree]=RhoAmplitudeNR2[maxdegree]/.\[Rho]mon2[{x_,y_,z_}] :> \[Rho]mon2Num[x,y,z]//Im//Transpose;


(* ::Subsection:: *)
(*Positivity of the Imaginary Part *)


(* 
ImAmplitude[Nmax]: Imaginary part of the amplitude at t=0 on the grid \[Rho]gridHSC;
*)


ClearAll[ImAmplitude];
ImAmplitude[nmax_]:=ImAmplitude[nmax]=Table[({0}~Join~RhoAmplitudeOpt[\[ScriptS][\[Rho]],1][nmax]),{\[Rho],\[Rho]gridHSC}]//Im;


(* ::Section:: *)
(*Building SDPB Input (UPDATE THE FOLDER)*)


(* ::Subsection:: *)
(*Defining the S-matrix partial wave decomposition: the high energy constraints will be imposed in the SDP code part*)
(*+ Objective Function*)


ClearAll[\[Alpha]abcReduced];
\[Alpha]abcReduced[nmax_]:=\[Alpha]abcReduced[nmax]=\[Alpha]abc[nmax]/.HEConstraintModuleOpt[nmax]//Variables;

ClearAll[matrix];
matrix[nmax_]:=matrix[nmax]=\[Alpha]abc[nmax]/.HEConstraintModuleOpt[nmax]//Table[Table[Coefficient[#[[kk]],x],{x,\[Alpha]abcReduced[nmax]}],{kk,1,Length[#]}]&;

ClearAll[SpwOPT];
SpwOPT[\[ScriptL]_][maxdegree_]:=Block[{tpwReduced,tpw,int1Reg,int2RegOPT},

int1Reg[l_][b_,c_]/;b<=c:=int1Reg[l][b,c]=ToExpression[Import["/gpfs/aguirrieri/integral_storage_SUGRA2/integral_storage_10d/int1Reg_10d_grid1_"<>ToString[l]<>"_"<>ToString[b]<>"_"<>ToString[c]<>".txt"]]//Chop[#,10^-prec]&;
int1Reg[l_][b_,c_]/;(b>c):=int1Reg[l][b,c]=(-1)^(l) int1Reg[l][c,b];

int2RegOPT[l_][a_,b_,c_]:=int2RegOPT[l][a,b,c]= \[Rho]gridCheb^a int1Reg[l][b,c];

tpw=( (\[ScriptS][\[Rho]gridCheb]^7)/(2^18 3 \[Pi]^4))({1/\[ScriptS][\[Rho]gridCheb]^3 (2^10 3)/Pochhammer[\[ScriptL]+1,6]}~Join~(BasisMandelstam[maxdegree]premon//ExpandAll//#/.mon[x_,y_,z_]pmon[a_,b_,c_]:> int2RegOPT[\[ScriptL]][x+a,y+b,z+c]&) //Transpose);
tpwReduced=tpw . matrix[maxdegree];

Print[\[ScriptL]];

Join[I #,{1}]&/@(tpwReduced)//Table[#[[x]],{x,constraintsgrid}]&];


ClearAll[target];
target[nmax_]:=target[nmax]=Block[{targetp,s\[ScriptS],t\[ScriptT],u\[ScriptU]},targetp=(Amplitude[\[Epsilon] s\[ScriptS],\[Epsilon] t\[ScriptT],\[Epsilon] u\[ScriptU]][nmax]//SeriesCoefficient[#,{\[Epsilon],0,4},Assumptions->s\[ScriptS]>0&&t\[ScriptT]<0&&u\[ScriptU]<0]&//#/s\[ScriptS]^4 &) . \[Alpha]abc[nmax]/.HEConstraintModuleOpt[nmax]//Expand//Table[Coefficient[#,x],{x,\[Alpha]abcReduced[nmax]}]&;Join[targetp,{0}]];


ClearAll[SetScale];
SetScale[nmax_,\[CapitalLambda]_]:=SetScale[nmax]=Join[\[Alpha]abcReduced[nmax]/.ac[___]->0/.\[DoubleStruckCapitalG]N->1,{-\[CapitalLambda]}];


ClearAll[HSConstraintReduced];
HSConstraintReduced[nmax_]:=HSConstraintReduced[nmax]=Table[Drop[HSConstraintOPT[nmax][[i]],-1] . \[Alpha]abc[nmax]/.HEConstraintModuleOpt[nmax]//Expand//Table[Coefficient[#,x],{x,\[Alpha]abcReduced[nmax]}]&//Join[#,{0}]&,{i,1,Length[HSConstraintOPT[nmax]]}];


ClearAll[ImAmplitudeReduced];
ImAmplitudeReduced[nmax_]:=ImAmplitudeReduced[nmax]=Table[ImAmplitude[nmax][[i]] . \[Alpha]abc[nmax]/.HEConstraintModuleOpt[nmax]//Expand//Table[Coefficient[#,x],{x,\[Alpha]abcReduced[nmax]}]&//Join[#,{0}]&,{i,1,2}];


ClearAll[SpwInfty];
SpwInfty[\[ScriptL]_][maxdegree_]:=SpwInfty[\[ScriptL]][maxdegree]=
Block[{tpwReduced},
If[FileExistsQ["/gpfs/hmurali/string_bootstrap/HEConditions/spwinfty_"<>ToString[\[ScriptL]]<>"_"<>ToString[maxdegree]<>".txt"],
tpwReduced=ToExpression[Import["/gpfs/hmurali/string_bootstrap/HEConditions/spwinfty_"<>ToString[\[ScriptL]]<>"_"<>ToString[maxdegree]<>".txt"]],
tpwReduced=1/(2^18 3\[Pi]^4) HEConditionOpt[7][maxdegree,\[ScriptL]] . \[Alpha]abc[maxdegree]/.HEConstraintModuleOpt[maxdegree]//Expand//Table[Coefficient[#,x],{x,\[Alpha]abcReduced[maxdegree]}]&;
Export["/gpfs/hmurali/string_bootstrap/HEConditions/spwinfty_"<>ToString[\[ScriptL]]<>"_"<>ToString[maxdegree]<>".txt",ToString[FullForm[tpwReduced]]];
];
Join[I (tpwReduced),{1}]];


(* ::Subsection:: *)
(*Building the Input*)


ClearAll[almostzeroarray];
almostzeroarray[len_,repl_] :=ReplacePart[ConstantArray[0,len],repl];


writeQuadraticProgram[file_,obj_,qc[\[Alpha]vecs_,\[Beta]vecs_]] :=
Block[{objnew, norm, tomatrix, pseudoidmatrix, quadraticconstraints,linearconstraints},
tomatrix[a_]:={{-2 Re[a],Im[a]-Re[a]},{Im[a]-Re[a],Im[a]}};
pseudoidmatrix:={{2,1},{1,1}};
quadraticconstraints =
Table[
PositiveMatrixWithPrefactor[1,(Transpose[Append[tomatrix /@ (\[Alpha]vecs[[k]])//Drop[#,-1]&,pseudoidmatrix+tomatrix[\[Alpha]vecs[[k]][[-1]]]],{3,2,1}])],{k,Length[\[Alpha]vecs]}];
linearconstraints=Table[
PositiveMatrixWithPrefactor[1,{{\[Beta]}}],{\[Beta],\[Beta]vecs}];

norm = almostzeroarray[Length[\[Alpha]vecs[[1]]],-1 -> 1];
WriteBootstrapSDP[file,SDP[obj,norm,quadraticconstraints~Join~linearconstraints]];
]


ClearAll[buildSmatrixprogram]
buildSmatrixprogram[file_][\[ScriptL]max_,nmax_,Linfty_]:=
Module[
{objective,\[ScriptL]grid,\[ScriptL]gridinfty,\[Alpha]vecsinfty,cvals,\[Beta]vecsHE,\[Beta]vecsHS,\[Beta]scale,\[Beta]impart},

\[ScriptL]grid = Range[0,\[ScriptL]max,2];
\[ScriptL]gridinfty=Range[0,Linfty,2];

Print["Start building quadratic constraints..."];

\[Alpha]vecs = ParallelTable[SpwOPT[\[ScriptL]][nmax],{\[ScriptL],\[ScriptL]grid}]//Flatten[#,1]&;

ClearAll[SpwOPT];

\[Alpha]vecsinfty=Table[SpwInfty[\[ScriptL]2][nmax],{\[ScriptL]2,\[ScriptL]gridinfty}];

ClearAll[SpwInfty];

Print["Start building linear constraints..."];


If[FileExistsQ["/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/linear_constraints_"<>ToString[nmax]<>".txt"],
\[Beta]vecsHS=ToExpression[Import["/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/linear_constraints_"<>ToString[nmax]<>".txt"]];,
\[Beta]vecsHS=HSConstraintReduced[nmax];
Export["/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/linear_constraints_"<>ToString[nmax]<>".txt", ToString[FullForm[HSConstraintReduced[nmax]]]]
];

ClearAll[HSConstraintReduced];

\[Beta]scale=Join[{SetScale[nmax, 512 \[Pi]^4]},{-SetScale[nmax, 512 \[Pi]^4]}];

\[Beta]impart=ImAmplitudeReduced[nmax];

gammal = 225 \[ScriptL]grid!(\[ScriptL]grid+7/2)/Gamma[\[ScriptL]grid+7];
Cl1vals = GegenbauerC[\[ScriptL]grid, 7/2, 1];
Imfls = Im[ReplacePart[\[Alpha]vecs, {_,-1} -> 0]/I];
alphaobj = Re[(2^19 3 \[Pi]^3 Table[dsgridCheb[[i]]/\[ScriptS][\[Rho]gridCheb[[i]]]^8 ((2^6 Gamma[7/2]^2)/\[Pi] (l!(l+7/2))/Gamma[l+7]) GegenbauerC[l,7/2,1]^2, {l,\[ScriptL]grid},{i, constraintsgrid}]//Flatten) . Imfls];
(*objective=-target[nmax];*)
objective = -alphaobj;
Print["writing xml input file..."];

(*writeQuadraticProgram[file,objective,qc[\[Alpha]vecs~Join~\[Alpha]vecsinfty,\[Beta]vecsHS~Join~\[Beta]scale~Join~\[Beta]impart]];*)
]


(* ::Section:: *)
(*Sdpb*)


ClearAll[maxabk];
maxabk[\[ScriptL]max_, nmax_,linfty_][fname_]:=
Block[{},

buildSmatrixprogram[fname<>".xml"][\[ScriptL]max,nmax,linfty];

Print["End xml input file ..."];

(*ClearAll[int,HSConstraintOPT];
Print["pvm2sdp start ..."];

Run["mpirun -v --oversubscribe -np 160 pvm2sdp 1024 "<>fname<>".xml ./"<>fname ]; 
DeleteFile[fname<>".xml"];

Print["sdpb start ..."];

Run["mpirun -v --oversubscribe -np 160 sdpb -s ./"<>fname<>" --noFinalCheckpoint --checkpointInterval 5000  --precision 1000 --procsPerNode 40  --maxRuntime 1072800 --maxIterations 100000 --initialMatrixScalePrimal 1e30 --initialMatrixScaleDual 1e30 --dualityGapThreshold 1e-12 --maxComplementarity 1e100 >> "<>fname<>".log"];
Run["rm -r "<>fname];*)
];


namerun=test;
Lmax=6;
Nmax=20;
Linfty=100;
(*namerun=ToExpression[$CommandLine[[4]]];
Lmax=ToExpression[$CommandLine[[5]]];
Nmax=ToExpression[$CommandLine[[6]]];
Linfty=ToExpression[$CommandLine[[7]]];*)
(*b3label=ToExpression[$CommandLine[[7]]];*)



maxabk[Lmax,Nmax,Linfty][ToString[namerun]<>"_"<>ToString[Lmax]<>"_"<>ToString[Nmax]<>"_"<>ToString[Linfty]];


ClearAll[importresultparallel,writeresultparallel]
importresultparallel[fname_,Nmax_,Lmax_,Linfty_] :=
importresultparallel[fname,Nmax,Lmax,Linfty]=Module[{res},
res=If[FileExistsQ["/gpfs/hmurali/string_bootstrap/"<>fname<>"_"<>ToString[Lmax]<>"_"<>ToString[Nmax]<>"_"<>ToString[Linfty]<>"_out/y.txt"],
Import["/gpfs/hmurali/string_bootstrap/"<>fname<>"_"<>ToString[Lmax]<>"_"<>ToString[Nmax]<>"_"<>ToString[Linfty]<>"_out/y.txt","Table"],
Print["Chirp!_"<>ToString[Lmax]<>"_"<>ToString[Nmax]]];
Export["/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/"<>fname<>"Out_"<>ToString[Lmax]<>"_"<>ToString[Nmax]<>"_"<>ToString[Linfty]<>".txt", ToString[FullForm[Table[res[[i]][[1]],{i,1,Length[res]}]]]]
]


importresultparallel[ToString[namerun],Nmax,Lmax,Linfty];
Block[{filename},
filename = ToString[namerun]<>"_"<>ToString[Lmax]<>"_"<>ToString[Nmax]<>"_"<>ToString[Linfty];
Run["mv "<>filename<>"_out/out.txt res_"<>filename<>".txt"];
Run["rm -r "<>filename<>"_out"];
Run["rm -r "<>filename<>".ck"];]


(*res = Import["/gpfs/hmurali/string_bootstrap/tmp_6_18_100_out/y.txt", "Table"][[All, 1]];
2/\[Pi] Sum[dsgridCheb[[i]]/sgridCheb[[i]]^5 Amplitude[sgridCheb[[i]]+I 10^-6,0,-sgridCheb[[i]]-I 10^-6][Nmax], {i, 1, Length@sgridCheb}]
ellp = (512 Pi^4 / (64 Pi^7))^(1/8);
Sum[dsgridHSC[[i]] / (32 Pi^8 ellp^14 sgridHSC[[i]]^5) (ImAmplitude[Nmax][[i]] . matrix[Nmax] . res), {i,1,4000}]
imar = Sum[(2^18 3 Pi^4)/ssmallgrid^7 Im@fls2d[[li]] gammal[[li]] Cl1vals[[li]]^2, {li,1,4}];
fls2d = ArrayReshape[fls, {4, Length@constraintsgrid}];
blah = Transpose[Im@RhoAmplitudeOpt[ssmallgrid, 1][Nmax]] . matres[[2;;]];
2/Pi Sum[dsgridCheb[[i]]/sgridCheb[[i]] blahbig[[i]], {i,1,Length@sgridCheb}];
2/Pi Sum[dsgridCheb[[constraintsgrid[[i]]]]/ssmallgrid[[i]] blah[[i]], {i,1,Length@ssmallgrid}];*)


(*ClearAll[SpwOPT2];
int1Reg[l_][b_,c_]/;b<=c:=int1Reg[l][b,c]=ToExpression[Import["/gpfs/aguirrieri/integral_storage_SUGRA2/integral_storage_10d/int1Reg_10d_grid1_"<>ToString[l]<>"_"<>ToString[b]<>"_"<>ToString[c]<>".txt"]]//Chop[#,10^-prec]&;
int1Reg[l_][b_,c_]/;(b>c):=int1Reg[l][b,c]=(-1)^(l) int1Reg[l][c,b];

int2RegOPT[l_][a_,b_,c_]:=int2RegOPT[l][a,b,c]= \[Rho]gridCheb^a int1Reg[l][b,c];
SpwOPT2[\[ScriptL]_][maxdegree_]:=Block[{tpwReduced,tpw},

tpw=({1/\[ScriptS][\[Rho]gridCheb]^3 (2^10 3)/Pochhammer[\[ScriptL]+1,6]}~Join~(BasisMandelstam[maxdegree]premon//ExpandAll//#/.mon[x_,y_,z_]pmon[a_,b_,c_]:> int2RegOPT[\[ScriptL]][x+a,y+b,z+c]&) //Transpose);
tpwReduced=tpw . matrix[maxdegree];

Print[\[ScriptL]];

tpwReduced];*)
