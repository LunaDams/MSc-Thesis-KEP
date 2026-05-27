# MSc-Thesis-KEP
WIP - This Repository contains all code used for my master's thesis in kidney exchange programmes, where highly-insentised and long-waiting patients are prioritised simultaneously. A banch-and-price algorithm is implemented for exact solutions.  Codes are written using Java in IntelliJ Idea, and library Gurobi version 12.0.3.

# Edge.java
A java class used to define edges used in compatibility graphs. An edge contains a origin node (int), destination node (int) and a weight (double)

# Graph.java
A java class used to define a compatibility graph of the KEP problem. A graph contains a list of pair nodes (Integer), a list of NDD (Non-directed-donor) nodes (Integers), and a list of edges (Edge). 

# Instance.java
A java class used to define an instance pool of the KEP. An instance constians the number of pairs (int), the number of NDDs (int), the list of node indices (Integer), the list of patient blood types (String), the list of donor blood types (String), the list of patient vPRA levels (Double), the binary list of NDD indications (Integer), and a list of patient waiting times (Double). 

# WaitingTimes.java
A class used for waiting time simulation of each patient using a point score system. Waiting times depend on blood-types and vPRA levels. Functions in the class are: Score, returns a list of waiting time scores for each patient. Pmatrix, returns a probability waiting time distribution for each patient. waitingtimes, returns waiting time of each patient using sampling from Pmatrix and exponential distributions. 

# Read.java
A class used to read the instance files containing a list of patient data and files containing the graph data. Data files are obtained from ohn P. Dickerson, Kidney Data (PrefLib), saved as .txt files. Functions in the class are: readInstance, returns the instance from the File. readGraph, returns the compatibility graph from the File. ExtractIndex, returns the index from a pair or NDD. isABOcompatible, returns the boolean determining if a patient and donor are blood-type compatible. 

# PIEF.java
A class used to find an LP-solution to PIEF for a given instance, in a single iteration. Functions in the class are: initialisePIEF, initialises the RMP of PIEF using only loop-edges. solvePIEF, solves the current RMP, obtains duals and prints results. subgraph, makes subgraphs for each given node including only nodes with higher or equal index as the given node. getAlphas, returns the alpha (capacity) duals. getBetas, returns the beta (flow) duals. PricingProblemPIEF, solves the pricing problem of PIEF for each subgraph. addColumnsPIEF, adds new columns to the RMP of PIEF. 

# PICEF.java
A class used to find an LP-solution to PICEF for a given instance, in a single iteration. Functions in the class are: initialisePICEF, initialises the RMP of PICEF using only loop-edges. solvePIEF, solves the current RMP, obtains duals and prints results. subgraph, makes subgraphs for each given node including only nodes with higher or equal index as the given node. getAlphas, returns the alpha (capacity) duals for pair nodes. getGammas, returns the gammas (capacity) duals for NDD nodes. getBetas, returns the beta (flow) duals. PricingProblemChains, solves the pricing problem of PICEF that finds chains. PricingProblemCycles, solves the pricing problem of PICEF that finds cycles. addchainColumns, adds new columns representing chains to the RMP of PICEF. addcycleColumns, adds new columns representing cycles to the RMP of PICEF. 

# CG.java
A class used to find an LP-solution of the KEP problem using PIEF or PICEF, and column generation. It iteratively solves the RMP and Pricing Problem until no improvements are found. Functions in the class are: CGPIEF: runs the column generation algorithm for PIEF. CGPICEF: runs the column generation algorithm for PICEF. 


