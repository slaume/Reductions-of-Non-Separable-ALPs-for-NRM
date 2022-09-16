

param counter; #auxiliary variable for iterations


#-----------------------------------------------------------------------------------
#SOLVER OPTIONS
#-----------------------------------------------------------------------------------
option solver cplex;
option show_stats 1;
option cplex_options 'baropt bardisplay=1 comptol=1e-8 crossover=0';


#-----------------------------------------------------------------------------------
#SETTING PARAMETERS, see 000Main.run
#-----------------------------------------------------------------------------------
param task; 
param datatype; 
param alp symbolic; 
param primaldual symbolic;


#-----------------------------------------------------------------------------------
#NETWORK PARAMETERS
#-----------------------------------------------------------------------------------
#general
param I; #number of legs
param J; #number of products
param T; #number of time steps
param S; #number of simulations
set Ical := 1..I; #legs
set Jcal := 0..J; #products
set Tcal := 1..T; #time steps
set Scal := 1..S; #simulations
param C; #initial capacity
param c {i in Ical}; #initial capacity per leg
param A {i in Ical, j in Jcal}; #incidence matrix 
param fconst {j in Jcal}; #fares if constant in time
param f {t in Tcal, j in Jcal}; #fares
param p {t in Tcal, j in Jcal}; #arrival probabilities
param pconst {j in Jcal}; #arrival probabilities if constant in time
param p0; #non-arrival probability per time step
param d {j in Jcal}; #total demand
#loads (per leg and total)
param load {i in Ical} := sum{t in Tcal, j in Jcal} p[t,j]*A[i,j]/c[i];
param loadtotal := 1/I*sum{i in Ical} load[i];
#for H&S setting
param L; #number of non-hub destinations
#for SBL setting
param loverline; #maximal length of an OD-pair
param lunderline; #minimal length of an OD-pair
param m; #number of SBL stringed together for CBL
param demandTwoLegs; #demand for OD-pair BD for small running example
param demandThreeLegs; #demand for OD-pair AD for small running example
param rexact {i in Ical}; #for solving the exact Bellman equation


#-----------------------------------------------------------------------------------
#PARTITIONS OF {1,...I} AND CORRESPONDING METHODS
#-----------------------------------------------------------------------------------
param N; #number of subnetworks Icaln
param qBar; #size of subnetworks
#sets Icaln etc.
set Icaln {n in 1..I}; #partition (Icaln)_n
#all resources i that are not covered by non-separable subnetworks 
#and are therefore treated affine-linear (AL)
set IcalAL := Ical diff (union {n in 1..N} Icaln[n]);
#sizes and indices of subnetworks
param qn {n in 1..I}; #sizes of Icaln
param inq {n in 1..I, q in 1..I}; #indices i of Icaln
#translating vector rn and decision u into a single scalar Rn and U,
#maximal Rn and maximal U
param Rnmax {n in 1..N} := sum{q in 1..qn[n]}c[inq[n,q]]*prod{qq in 1..q-1}(c[inq[n,qq]]+1);
param Umax := sum{j in 0..J} 2^j;


#-----------------------------------------------------------------------------------
#GO THROUGH ALL PARTITIONS
#-----------------------------------------------------------------------------------
param qCompPcal;
param ICompPcal;
set Pcal {qbar in 0..I, Ibar in 0..I}; 
param idiff {q in 1..I}; 
param ptilde {n in 1..I}; 
param Pbar {n in 1..I};
set PcalTilde := 0..sum{n in 1..floor(I/qBar)}(Pbar[n]-1)*(prod{nn in 1..n-1}Pbar[nn]);
param IcalnPsiPart {n in 1..I, ptil in 1..Pbar[n]};
set IcalReverseFromScalar;
param qReverseFromScalar;
param IcalnPsiReverseFromScalar;
param IReverseFromScalar;
param iAllReverseFromScalar {i in 1..I};
param iReverseFromScalar {q in 1..I};
set IcalnReverseFromScalar;
param iAux;
param inqAux {q in 1..qBar};
set IcalnAux;
set IcalRest; 
set IcalLeft; 
set IcalnAllPartitions {IcalPsiRest in Pcal[I mod qBar,I], ptildePsi in PcalTilde, n in 1..N};
set IcalnAllSubsets {IcalPsi in Pcal[qBar,I]};
set IcalnMaxM {n in 1..N};
param maxM;


#-----------------------------------------------------------------------------------
#SIMULATION
#-----------------------------------------------------------------------------------
#general
param rCapacity {i in Ical}; #remaining capacity
param available; #1 or 0 whether product is available or not
param jArrived; #requested product j
param gain {s in Scal}; #gain per simulation s
param revenue; #average revenue of all simulations s
param standarderror; #standard error of average revenue
param marginalvalue; #marginal value using the value function approximation
#dices for randomness
param dice;
param dice1;
param dice2;
param upperbound; #optimal value of linear program, Z
param computingtime; #time used for computing linear program
param tol := 0.0000001; #tolerance for accepting customer request


#-----------------------------------------------------------------------------------
#VALUE FUNCTION APPROXIMATION (VFA)
#-----------------------------------------------------------------------------------
#general NSEP including AL, SPL, NSEP
param RnTop := (C+1)^qBar;
param Vnonlin {t in 1..T+1, n in 1..I, Rn in 0..RnTop}; #variable W in dual LP
param Vlin {t in 1..T+1, i in Ical}; #variable V in dual LP


#-----------------------------------------------------------------------------------
#DETERMINE NETWORK MEASURE
#-----------------------------------------------------------------------------------
param Woverline {i in 1..I};
param betaoverline {t in Tcal, i in 1..I, j in Jcal, u in 0..1};
param Wstar {t in 1..T+1, n in 1..N, Rn in 0..Rnmax[n]};
param rNWM {n in 1..N, q in 1..qn[n]};
param Msubnetwork {n in 1..N};
param Mpartition;


#-----------------------------------------------------------------------------------
#VARIABLES
#-----------------------------------------------------------------------------------

#primal
var varrho {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]} >= 0;
var mu {t in Tcal, j in Jcal, u in 0..1} >= 0;
var xi {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n], u in 0..1} >= 0;
var rho {t in Tcal, i in Ical} >= 0;
var Upsilon {t in Tcal, i in Ical, n in 1..N, Rn in 0..Rnmax[n]} >= 0;
var WNM {t in 1..T+1, i in 1..I, r in 0..c[i]};
var alphaNM {t in Tcal, i in 1..I, j in Jcal, r in 0..c[i]};
var betaNM {t in Tcal, i in 1..I, j in Jcal, u in 0..1};
var xiHalf {t in Tcal, n in 1..N, Rn in 0..Rnmax[n], U in 0..Umax} >= 0;
var muHalf {t in Tcal, U in 0..Umax} >= 0;

#dual
var W {t in 1..T+1, n in 1..N, Rn in 0..Rnmax[n]};
var V {t in 1..T+1, i in Ical};
var alpha {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n]};
var beta {t in Tcal, n in 1..I, j in Jcal, u in 0..1};
var gamma {t in Tcal, n in 1..N, i in Ical, j in Jcal, Rn in 0..Rnmax[n]} >= 0;
var delta {t in Tcal, n in 1..I, i in Ical} >= 0;
var epsilon {t in Tcal, n in 1..I, nn in 1..I, i in Ical} >= 0;
var theta {t in Tcal, j in Jcal};
var phi {t in Tcal, i in Ical, j in Jcal} >= 0;
var alphaHalf {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]};
var betaHalf {t in Tcal, n in 1..N, U in 0..Umax};
var gammaHalf {t in Tcal, i in Ical} >= 0;
var thetaHalf {t in Tcal};


#-----------------------------------------------------------------------------------
#LINEAR PROGRAMS
#-----------------------------------------------------------------------------------

#D' for Network Measure, SPL (N = I)
#-----------------------------------------
minimize ZNM: sum{i in 1..I} WNM[1,i,c[i]];
subject to conNMa {t in Tcal, i in 1..I, r in 0..c[i]}:
	WNM[t,i,r] - sum{j in Jcal} alphaNM[t,i,j,r] >= 0;
subject to conNMb {t in Tcal, i in 1..I, j in Jcal, r in 0..c[i], 
	u in 0..1 diff {2+min(0,r - A[i,j])}}:
	alphaNM[t,i,j,r]+betaNM[t,i,j,u] 
	- p[t,j]*WNM[t+1,i,r-u*A[i,j]] >= 0;
subject to conNMc {t in Tcal, j in Jcal, u in 0..1}:
	-sum{i in 1..I} betaNM[t,i,j,u] >= p[t,j]*f[t,j]*u;
subject to conNMd {i in 1..I, r in 0..c[i]}:
	WNM[T+1,i,r] = 0;
	
#P' and H
#-----------------------------------------
maximize Z: sum{t in Tcal, j in Jcal, u in 0..1} p[t,j]*f[t,j]*u*mu[t,j,u];
#flow constraints for non-sep parts
subject to conP1a {n in 1..N}: varrho[1,n,Rnmax[n]] = 1;
subject to conP1b {t in 2..T, n in 1..N, Rn in 0..Rnmax[n]}:
	varrho[t,n,Rn] = sum{j in Jcal} 
	sum{u in 0..min(1,min{q in 1..qn[n]} (1+c[inq[n,q]]-A[inq[n,q],j]-
		floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
		/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
	))}
	p[t-1,j]*xi[t-1,n,j,
		sum{q in 1..qn[n]} (prod{qq in 1..q-1} (c[inq[n,qq]]+1))*
			(u*A[inq[n,q],j] +
				floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
				/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
			),
	u];
#probability-distribution consistency constraints
subject to conP2a {t in Tcal, n in 1..N}:
	sum{Rn in 0..Rnmax[n]} varrho[t,n,Rn] = 1;
subject to conP2b {t in Tcal, j in Jcal}: sum{u in 0..1} mu[t,j,u] = 1;
subject to conP2c {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n]}:
	sum{u in 0..1} xi[t,n,j,Rn,u] = varrho[t,n,Rn];
subject to conP2d {t in Tcal, n in 1..N, j in Jcal, u in 0..1}:
	sum{Rn in 0..Rnmax[n]} xi[t,n,j,Rn,u] = mu[t,j,u];
#AL constraints
subject to conP3a {i in IcalAL}: rho[1,i] = c[i];
subject to conP3b {t in 2..T, i in IcalAL}: rho[t,i] = 
	rho[t-1,i] - sum{j in Jcal} p[t-1,j]*A[i,j]*mu[t-1,j,1];
subject to conP3c {t in Tcal, j in Jcal, i in IcalAL}:
	mu[t,j,1]*A[i,j] <= rho[t,i]; 
#additional Upsilon-constraints for post-decision heuristic (H)
subject to conP4a {t in Tcal, n in 1..N, j in Jcal, i in Ical diff Icaln[n], Rn in 0..Rnmax[n]}:
	A[i,j]*xi[t,n,j,Rn,1] <= Upsilon[t,i,n,Rn];
subject to conP4b {t in Tcal, n in 1..N, i in IcalAL}:
	sum {Rn in 0..Rnmax[n]} Upsilon[t,i,n,Rn] <= rho[t,i];
subject to conP4c {t in Tcal, n in 1..N, nn in 1..N diff {n}, q in 1..qn[nn]}:
	sum {Rn in 0..Rnmax[n]} Upsilon[t,inq[nn,q],n,Rn] <= sum {Rnn in 0..Rnmax[nn]}
		varrho[t,nn,Rnn]*min(1,
			floor((Rnn mod prod{qq in 1..q}(c[inq[nn,qq]]+1))
			/(prod{qq in 1..q-1}(c[inq[nn,qq]]+1)))
		);
#make sure that only valid xi-variables are non-zero
subject to conPX {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n]}:
	sum{u in 
		2+min(0,min{q in 1..qn[n]}(-A[inq[n,q],j]+
			floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
			/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
		))
		..1} xi[t,n,j,Rn,u] = 0;
		
#D^H and D'
#-----------------------------------------
minimize ZD: sum{n in 1..N} W[1,n,Rnmax[n]] + sum{i in IcalAL} V[1,i]*c[i] 
	+ sum{t in Tcal, j in Jcal} theta[t,j];
subject to conD1 {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]}:
	W[t,n,Rn] - sum{j in Jcal} alpha[t,n,j,Rn] 
	- (if alp = 'H' then 1 else 0)*sum{q in 1..qn[n], nn in 1..N diff {n}}
	epsilon[t,nn,n,inq[n,q]]*min(1,
		floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
		/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
	) >= 0;
subject to conD2 {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n], 
	u in 0..1 diff {2+min(0,min{q in 1..qn[n]} (
		floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
		/(prod{qq in 1..q-1}(c[inq[n,qq]]+1))) - A[inq[n,q],j]
	))}}:
	alpha[t,n,j,Rn]+beta[t,n,j,u] 
	+ (if alp = 'H' then 1 else 0)*sum{i in Ical diff Icaln[n]} A[i,j]*u*gamma[t,n,i,j,Rn]
	- p[t,j]*W[t+1,n,
		sum{q in 1..qn[n]} (prod{qq in 1..q-1} (c[inq[n,qq]]+1))*
			(-u*A[inq[n,q],j] +
				floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
				/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
			)
	] >= 0;
subject to conD3 {t in Tcal, j in Jcal, u in 0..1}:
	-sum{n in 1..N} beta[t,n,j,u] + sum{i in IcalAL}p[t,j]*A[i,j]*u*V[t+1,i]
	+ sum{i in IcalAL} u*A[i,j]*phi[t,i,j] + theta[t,j] >= p[t,j]*f[t,j]*u;
subject to conD4 {t in Tcal, i in IcalAL}:
	V[t,i]-V[t+1,i] - (if alp = 'H' then 1 else 0)*sum{n in 1..N} delta[t,n,i]
	- sum{j in Jcal} phi[t,i,j] >= 0;
subject to conD5a {t in Tcal, n in 1..N, Rn in 0..Rnmax[n], i in IcalAL}:
	-(if alp = 'H' then 1 else 0)*sum{j in Jcal} gamma[t,n,i,j,Rn] 
	+ (if alp = 'H' then 1 else 0)*delta[t,n,i] >= 0;
subject to conD5b {t in Tcal, n in 1..N, Rn in 0..Rnmax[n], nn in 1..N diff {n}, i in Icaln[nn]}:
	-(if alp = 'H' then 1 else 0)*sum{j in Jcal} gamma[t,n,i,j,Rn] 
	+ (if alp = 'H' then 1 else 0)*epsilon[t,n,nn,i] >= 0;
subject to conD6a {n in 1..N, Rn in 0..Rnmax[n]}:
	W[T+1,n,Rn] = 0;
subject to conD6b {i in IcalAL}:
	V[T+1,i] = 0;
	
#P only reduced in the state space
#-----------------------------------------
maximize Zhalf: sum{t in Tcal, j in Jcal, U in 0..Umax} p[t,j]*f[t,j]*floor((U mod 2^(j+1))/(2^j))*muHalf[t,U];
#flow constraints for non-sep parts
subject to conPhalf1a {n in 1..N}: varrho[1,n,Rnmax[n]] = 1;
subject to conPhalf1b {t in 2..T, n in 1..N, Rn in 0..Rnmax[n]}:
	varrho[t,n,Rn] = sum{j in Jcal} 
	sum{U in 0..Umax} min(1,min{q in 1..qn[n]} (
		1+c[inq[n,q]]-
			floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
			/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
		- A[inq[n,q],j]*
			floor((U mod 2^(j+1))/(2^j))
	))*
	p[t-1,j]*xiHalf[t-1,n,min(Rnmax[n],
		sum{q in 1..qn[n]} (prod{qq in 1..q-1} (c[inq[n,qq]]+1))*
			(floor((U mod 2^(j+1))/(2^j))*A[inq[n,q],j] +
				floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
				/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
			)),
	U];
#probability-distribution consistency constraints
subject to conPhalf2a {t in Tcal, n in 1..N}:
	sum{Rn in 0..Rnmax[n]} varrho[t,n,Rn] = 1;
subject to conPhalf2b {t in Tcal}:
	sum{U in 0..Umax} muHalf[t,U] = 1;
subject to conPhalf2c {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]}:
	sum{U in 0..Umax} xiHalf[t,n,Rn,U] = varrho[t,n,Rn];
subject to conPhalf2d {t in Tcal, n in 1..N, U in 0..Umax}:
	sum{Rn in 0..Rnmax[n]} xiHalf[t,n,Rn,U] = muHalf[t,U];
#AL constraints
subject to conPhalf3a {i in IcalAL}: rho[1,i] = c[i];
subject to conPhalf3b {t in 2..T, i in IcalAL}: rho[t,i] = 
	rho[t-1,i] - sum{j in Jcal, U in 0..Umax} p[t-1,j]*A[i,j]
	*floor((U mod 2^(j+1))/(2^j))
	*muHalf[t-1,U];
subject to conPhalf3c {t in Tcal, i in IcalAL}:
	sum{U in 0..Umax} (max{j in Jcal}(floor((U mod 2^(j+1))/(2^j))*
		A[i,j])) * muHalf[t,U] <= rho[t,i]; 
#make sure that only valid xiHalf-variables are non-zero	
subject to conPhalfX {t in Tcal, n in 1..N, Rn in 0..Rnmax[n], U in 0..Umax}:
	max(0,- min{j in Jcal, q in 1..qn[n]}(
		floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
		/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
		- A[inq[n,q],j]*floor((U mod 2^(j+1))/(2^j))
	))*xiHalf[t,n,Rn,U] = 0;

#D only reduced in the state space
#-----------------------------------------
minimize ZDhalf: sum{n in 1..N} W[1,n,Rnmax[n]] + sum{i in IcalAL}V[1,i]*c[i] + sum{t in Tcal} thetaHalf[t];
subject to conDhalf1 {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]}:
	W[t,n,Rn] - alphaHalf[t,n,Rn] >= 0;
subject to conDhalf2 {t in Tcal, n in 1..N, Rn in 0..Rnmax[n], U in 0..Umax}:
	min(1,1+min{j in Jcal, q in 1..qn[n]}(
		floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
		/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
		- A[inq[n,q],j]*floor((U mod 2^(j+1))/(2^j))
	))*
	(alphaHalf[t,n,Rn] + betaHalf[t,n,U] - sum{j in Jcal} p[t,j]*W[t+1,n,
		max(0,
			sum{q in 1..qn[n]} (prod{qq in 1..q-1} (c[inq[n,qq]]+1))*
			(-floor((U mod 2^(j+1))/(2^j))*A[inq[n,q],j] +
				floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
				/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
			))
	]) >= 0;
subject to conDhalf3 {t in Tcal, U in 0..Umax}:
	thetaHalf[t]
	-sum{n in 1..N} betaHalf[t,n,U] 
	+ sum{i in IcalAL,j in Jcal} A[i,j]*
		floor((U mod 2^(j+1))/(2^j))*p[t,j]*V[t+1,i]
	+ sum{i in IcalAL} (max{j in Jcal} A[i,j]*floor((U mod 2^(j+1))/(2^j)))*gammaHalf[t,i]
	>= sum{j in Jcal} p[t,j]*f[t,j]*floor((U mod 2^(j+1))/(2^j));
subject to conDhalf4 {t in Tcal, i in IcalAL}:
	V[t,i] - V[t+1,i] - gammaHalf[t,i] >= 0;
subject to conDhalf5a {n in 1..N, Rn in 0..Rnmax[n]}:
	W[T+1,n,Rn] = 0;
subject to conDhalf5b {i in IcalAL}:
	V[T+1,i] = 0;
		


#D' for SPL, used for determining network measure
problem DNM: WNM, alphaNM, betaNM, conNMa, conNMb, conNMc, conNMd, ZNM;	
#H 	
problem Hnonsep: varrho, mu, xi, rho, Upsilon, 
	conP1a, conP1b, conP2a, conP2b, conP2c, conP2d, 
	conP3a, conP3b, conP3c, conP4a, conP4b, conP4c, 
	conPX,
	Z;
#P' 
problem Pnonsepprime: varrho, mu, xi, rho, 
	conP1a, conP1b, conP2a, conP2b, conP2c, conP2d, 
	conP3a, conP3b, conP3c, 
	conPX,
	Z;
#D^H and D'
problem DHnonsep: W,V,alpha,beta,gamma,delta,epsilon,theta,phi,
	conD1,conD2,conD3,conD4,conD5a,conD5b,conD6a,conD6b,
	ZD;
#P only reduced in the state space
problem PnonsepHalf: varrho, rho, xiHalf, muHalf,
	conPhalf1a, conPhalf1b, conPhalf2a, conPhalf2b, conPhalf2c, conPhalf2d, 
	conPhalf3a, conPhalf3b, conPhalf3c,
	conPhalfX,
	Zhalf;
#D only reduced in the state space
problem DnonsepHalf: W,V,alphaHalf,betaHalf,gammaHalf,thetaHalf,
	conDhalf1, conDhalf2, conDhalf3, conDhalf4, conDhalf5a, conDhalf5b, 
	ZDhalf;



