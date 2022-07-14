


#-----------------------------------------------------------------------------------
#SOLVER OPTIONS
#-----------------------------------------------------------------------------------
option solver cplex;
option show_stats 0;
option cplex_options 'baropt bardisplay=0 comptol=1e-8 crossover=0';



param counter; #auxiliary
param solveH; #1 or 0, solve H or P'
param solveP := 1-solveH;
param tol := 0.0000001;
param task;
param datatype;


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



#-----------------------------------------------------------------------------------
#PARTITIONS OF {1,...I} AND CORRESPONDING METHODS
#-----------------------------------------------------------------------------------
param N; #number of subnetworks Icaln
param qBar; #size of subnetworks for NSEP
#sets Icaln etc.
set Icaln {n in 1..I}; #partition (Icaln)_n
#all resources i that are not covered by non-sep parts 
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
param jArrived; #requested product
param gain {s in Scal}; #gain per simulation s
param revenue; #total revenue, sum of gain over all simulations s
param standarderror; #standard error of gain over all simulations s
param marginalvalue; #approximate difference of value function, marginal value
#dices for randomness
param dice;
param dice1;
param dice2;
param upperbound;



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

#for Primal of NSEP
var varrho {t in Tcal, n in 1..N, Rn in 0..Rnmax[n]} >= 0;
var mu {t in Tcal, j in Jcal, u in 0..1} >= 0;
var xi {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n], u in 0..1} >= 0;
var rho {t in Tcal, i in Ical} >= 0;
var Upsilon {t in Tcal, i in Ical, n in 1..N, Rn in 0..Rnmax[n]} >= 0;

#for Dual of SPL
var W {t in 1..T+1, i in 1..I, r in 0..c[i]};
var alpha {t in Tcal, i in 1..I, j in Jcal, r in 0..c[i]};
var beta {t in Tcal, i in 1..I, j in Jcal, u in 0..1};


#(P') and (H)
#-----------------------------------------
maximize Z: sum{t in Tcal, j in Jcal, u in 0..1} p[t,j]*f[t,j]*u*mu[t,j,u];
#flow constraints for non-sep parts
subject to con1a {n in 1..N}: varrho[1,n,Rnmax[n]] = 1;
subject to con1b {t in 2..T, n in 1..N, Rn in 0..Rnmax[n]}:
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
subject to con2a {t in Tcal, n in 1..N}:
	sum{Rn in 0..Rnmax[n]} varrho[t,n,Rn] = 1;
subject to con2b {t in Tcal, j in Jcal}: sum{u in 0..1} mu[t,j,u] = 1;
subject to con2c {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n]}:
	sum{u in 0..1} xi[t,n,j,Rn,u] = varrho[t,n,Rn];
subject to con2d {t in Tcal, n in 1..N, j in Jcal, u in 0..1}:
	sum{Rn in 0..Rnmax[n]} xi[t,n,j,Rn,u] = mu[t,j,u];
#AL constraints
subject to con3a {i in IcalAL}: rho[1,i] = c[i];
subject to con3b {t in 2..T, i in IcalAL}: rho[t,i] = 
	rho[t-1,i] - sum{j in Jcal} p[t-1,j]*A[i,j]*mu[t-1,j,1];
subject to con3c {t in Tcal, j in Jcal, i in IcalAL}:
	mu[t,j,1]*A[i,j] <= rho[t,i]; 
#additional Upsilon-constraints for post-decision heuristic (H)
subject to con4a {t in Tcal, n in 1..N, j in Jcal, i in Ical diff Icaln[n], Rn in 0..Rnmax[n]}:
	A[i,j]*xi[t,n,j,Rn,1] <= Upsilon[t,i,n,Rn];
subject to con4b {t in Tcal, n in 1..N, i in IcalAL}:
	sum {Rn in 0..Rnmax[n]} Upsilon[t,i,n,Rn] <= rho[t,i];
subject to con4c {t in Tcal, n in 1..N, nn in 1..N diff {n}, q in 1..qn[nn]}:
	sum {Rn in 0..Rnmax[n]} Upsilon[t,inq[nn,q],n,Rn] <= sum {Rnn in 0..Rnmax[nn]}
		varrho[t,nn,Rnn]*min(1,
			floor((Rnn mod prod{qq in 1..q}(c[inq[nn,qq]]+1))
			/(prod{qq in 1..q-1}(c[inq[nn,qq]]+1)))
		);
#make sure that only valid xi-variables are non-zero
subject to conX {t in Tcal, n in 1..N, j in Jcal, Rn in 0..Rnmax[n]}:
	sum{u in 
		2+min(0,min{q in 1..qn[n]}(-A[inq[n,q],j]+
			floor((Rn mod prod{qq in 1..q}(c[inq[n,qq]]+1))
			/(prod{qq in 1..q-1}(c[inq[n,qq]]+1)))
		))
		..1} xi[t,n,j,Rn,u] = 0;
	
		
#(D') for Network Measure, SPL (N = I)
#-----------------------------------------
minimize ZSPL: sum{i in 1..I} W[1,i,c[i]];
subject to con5a {t in Tcal, i in 1..I, r in 0..c[i]}:
	W[t,i,r] - sum{j in Jcal} alpha[t,i,j,r] >= 0;
subject to con5b {t in Tcal, i in 1..I, j in Jcal, r in 0..c[i], 
	u in 0..1 diff {2+min(0,r - A[i,j])}}:
	alpha[t,i,j,r]+beta[t,i,j,u] 
	- p[t,j]*W[t+1,i,r-u*A[i,j]] >= 0;
subject to con5c {t in Tcal, j in Jcal, u in 0..1}:
	-sum{i in 1..I} beta[t,i,j,u] >= p[t,j]*f[t,j]*u;
subject to con5d {i in 1..I, r in 0..c[i]}:
	W[T+1,i,r] = 0;
		
		
		
#(H) 	
problem Hnonsep: varrho, mu, xi, rho, Upsilon, 
	con1a, con1b, con2a, con2b, con2c, con2d, con3a, con3b, con3c, 
	con4a, con4b, con4c, 
	conX,
	Z;
#(P') 
problem Pnonsepprime: varrho, mu, xi, rho, 
	con1a, con1b, 
	con2a, 
	con2b, 
	con2c, con2d, 
	con3a, con3b, con3c,  
	conX,
	Z;
#(D') for SPL
problem DSPL: W, alpha, beta, con5a, con5b, con5c, con5d, ZSPL;



