# Reductions-of-Non-Separable-ALPs-for-NRM

The files in this depository can be used for reproducing the results presented in the paper "Reductions of Non-Separable Approximate Linear Programs for Network Revenue Management" (Simon Laumer and Christiane Barz) submitted to EJOR. 

For every experiment, you have to execute the file "001Main.run" in AMPL. In this file, you can specify the following parameters:
- The parameter "task": task = 1: Solve an ALP and simulate the policy for a manually chosen partition; task = 2: Determine the partition for which the network measure is largest, which can then be used for NSEP(I/q,q,gl.); task = 3: Determine the subset for which the network measure is largest, which can then be used for NSEP(1,q,gr.); task = 4: Solve ALP and simulate the policy for all partitions, the result is printed to "000correlation.txt".
- The parameter "datatype": datatype = 1: Small Running Example (SRE, Sections 2 to 7), parameters demandTwoLegs and demandThreeLegs have to be specified; datatype = 2: Hub-and-spoke network (H&S, Section 8.1.1), parameters L,T,C have to be specified in "002datafileHS.dat"; datatype = 3: Simple bus line (SBL, Section 8.1.2), parameters I,C,T,l^,l_ have to be specified in "002datafileSBL.dat"; datatype = 4: Consecutive bus line (CBL, Section 8.1.3), parameters m,I,C,T have to be specified in "002datafileCBL.dat"; datatype = 5: Real bus line (RBL, Section 8.1.4 and Appendix H); datatype = 6: Example to show non-equivalance of ALPs (NE, see Appendix C). 
- The parameter "qBar": This is the size of the subnetwork(s) or subset(s) of {1,...,I}.
- The parameter "S": This is the number of simulation runs done for determining the average revenue of the policy.
- The parameter "alp": If alp = 'H', then the ALP (H_NSEP) is solved; If alp = 'Pprime', then the ALP (P_NSEP') is solved; If alp = 'P', then the unreduced ALP P_NSEP is solved.
- The parameter "primaldual": If primaldual = 'primal', then the primal ALP is solved; If primaldual = 'dual', then the dual ALP is solved. In the experiment in the paper, we have always used 'primal'.
- The parameter "N" (only if the chosen task is 1): The number of subnetworks.
- The sets "Icaln[1]",...,"Icaln[N]" (only if the chosen task is 1): The explicit choice of subnetworks.
