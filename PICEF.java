// Packages used.
import java.util.*;
import com.gurobi.gurobi.*;

public class PICEF {

    // This class is for the PICEF (combined chain and cycle) algorithm
    // It produces solutions for the KEP when NDDs are available

    // An environment so that we can solve the RMP
    public GRBEnv env;
    // A GRB model so that we can solve the RMP
    public GRBModel model;

    // A multidimensional map of edge variables (one for each edge and position)
    Map<Edge, List<GRBVar>> edgevars;
    // A map of cycle variables (one for each cycle) and their corresponding lengths (number of edges)
    Map<String, GRBVar> cyclevars;
    Map<String, Integer> cyclelengths;
    // A map of loop variables (one for each node) and corresponding loop weights
    Map<Integer, GRBVar> loopvars;
    private Map<Integer, Double> loopweights;

    // A map of capacity constraints of the pair nodes (one for each node)
    private Map<Integer, GRBConstr> paircapconstrs;
    // A list of capacity constraints of the NDD nodes (one for each node)
    private Map<Integer, GRBConstr> NDDcapconstrs;
    // A multidimensional map of flow constraints (one for each node and position)
    // Structure = Map<i, Map<c, constr>>
    private Map<Integer, Map<Integer, GRBConstr>> flowconstrs;

    // A map of capacity duals for the pairs (one for each node)
    private Map<Integer, Double> alphas;
    // A map of capacity duals for the NDDs (one for each node)
    private Map<Integer, Double>  gammas;
    // A multidimensional map of flow duals (one for each node and position)
    private Map<Integer, Map<Integer, Double>> betas;

    // A set of edge variables that have been added to the RMP
    private Set<String> addedvars = new HashSet<>();
    // A set of cycles that have been added to the RMP
    Set<String> addedcycles = new HashSet<>();
    // A set of chains that have been added to the RMP
    private Set<String> addedchains = new HashSet<>();

    // A function to initialise the PIEF RMP
    public void initialisePICEF(Graph G, int C) throws GRBException{

        // Initialise the environment
        env =new GRBEnv(true);
        env.set("logfile", "rmp.log");
        env.start();

        // Create empty model
        model = new GRBModel(env);

        // Initialise the list of edge variables, but do not add them yet
        edgevars = new HashMap<>();

        // Initialise the list of cycle variables and their lengths, but do not add them yet
        cyclevars = new HashMap<>();
        cyclelengths = new HashMap<>();

        // Initialise the list of loop variables and their weights (for the objective);
        loopvars = new HashMap<>();
        loopweights = new HashMap<>();
        for (Edge e: G.edges){
            // A loop exists if the edge from a node goes to itself
            if (e.from == e.to){
                GRBVar loopvar = model.addVar(0, 1, 0, GRB.CONTINUOUS, "loop" + e.from);
                loopvars.put(e.from, loopvar);
                loopweights.put(e.from, e.weight);
            }
        }

        // update the model
        model.update();

        // Initialise the capacity constraint for the pair nodes
        paircapconstrs = new HashMap<>();
        // (only add the loop variables since there are no edge variables yet)
        for (int p : G.pairnodes){
            GRBLinExpr expr = new GRBLinExpr();

            // Add the loop to the constraint if it exists
            // Add the loop to the constraint if it exists
            if ( loopvars.containsKey(p)){
                expr.addTerm(1.0 ,loopvars.get(p));
            }

            GRBConstr capconstr = model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "pair capacity" + p);
            paircapconstrs.put(p, capconstr);
        }

        // Initialise the capacity constraint for the NDD nodes
        NDDcapconstrs = new HashMap<>();
        for (int d : G.NDDnodes){
            GRBLinExpr expr = new GRBLinExpr();
            GRBConstr capconstr = model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "NDD capacity" + d);
            NDDcapconstrs.put(d, capconstr);
        }

        // Initialise the flow constraints
        flowconstrs = new HashMap<>();
        for (int i : G.pairnodes){
            Map<Integer, GRBConstr> nodeflowconstrs = new HashMap<>();
                for (int c=1; c < C; c++){
                    GRBLinExpr expr = new GRBLinExpr();
                    // Do not add any edges yet
                    // When sums are added, use the form Sum_in - Sum_out = 0
                    GRBConstr flowconstr = model.addConstr(expr, GRB.EQUAL, 0.0,
                            "flow^" + i + "_" + c);
                    nodeflowconstrs.put(c, flowconstr);
                }
            flowconstrs.put(i, nodeflowconstrs);
        }

        // Set objective
        // Only the loops are in the initial objective, other variables will be added iteratively
        GRBLinExpr obj =new GRBLinExpr();
        for (int p : loopvars.keySet()){
            obj.addTerm(loopweights.get(p), loopvars.get(p));
        }
        model.setObjective(obj, GRB.MAXIMIZE);
        model.update();
    }

    // Solve the master problem
    public void solvePICEF(Graph G,int C) throws GRBException{

        // Use interior point solving instead of simplex
        // Because of degeneracy possibility
        model.set(GRB.IntParam.Method, 2);

        // Optimise the model
        model.optimize();

        // Initialise the duals
        alphas = new HashMap<>();
        gammas = new HashMap<>();
        betas = new HashMap<>();

        // Obtain the pair capacity duals for each node
        for (int p : G.pairnodes){
            GRBConstr capconstr = paircapconstrs.get(p);
            double alpha = capconstr.get(GRB.DoubleAttr.Pi);
            alphas.put(p, alpha);
        }

        // Obtain the NDD capacity duals for each node
        for (int d : G.NDDnodes){
            GRBConstr capconstr = NDDcapconstrs.get(d);
            double gamma = capconstr.get(GRB.DoubleAttr.Pi);
            gammas.put(d, gamma);
        }

        // Obtain the flow duals for each node and position
        for (int i : G.pairnodes) {
            Map<Integer, GRBConstr> nodeflowconstrs = flowconstrs.get(i);
            Map<Integer, Double> nodebetas = new HashMap<>();
            for (int c = 1; c < C; c++) {
                GRBConstr flowconstr = nodeflowconstrs.get(c);
                double beta = flowconstr.get(GRB.DoubleAttr.Pi);
                nodebetas.put(c, beta);
            }
            // Store the duals in a map with the node ID as key
            betas.put(i, nodebetas);
        }

        // Check the status
        int status = model.get(GRB.IntAttr.Status);
        if (status == GRB.OPTIMAL || status == GRB.TIME_LIMIT || status == GRB.SUBOPTIMAL){
            double obj = model.get(GRB.DoubleAttr.ObjVal);
            System.out.println("Solution found!");
            System.out.println("Objective = " + obj);
        } else if (status == GRB.INFEASIBLE){
            System.out.println("RMP is infeasible!");
        } else {
            System.out.println("Model not solved to optimality. Status = " + status);
        }

        // Check the used variables in the solution
        System.out.println("Solution values: ");

        // Initialise the number of matches
        double edgecount = 0;

        // Check the loops first
        for (int p : loopvars.keySet()){
            GRBVar var = loopvars.get(p);
            String name = var.get(GRB.StringAttr.VarName);
            double val = var.get(GRB.DoubleAttr.X);
            if (val > 1E-6){
                System.out.println("Used loop variable " + name + " with value " + val);
                edgecount += val;
            }
        }

        // Check the edges
        if (edgevars.isEmpty()){
            System.out.println("No edges are used.");
        } else {
            for ( Edge e : edgevars.keySet()){
                List<GRBVar> nodelist = edgevars.get(e);

                for (int c=0; c < nodelist.size(); c++){
                    GRBVar var = nodelist.get(c);
                    String name = var.get(GRB.StringAttr.VarName);
                    double val = var.get(GRB.DoubleAttr.X);

                    if (val > 1E-6){
                        System.out.println("Used edge variable " + name + " with value " + val);
                        edgecount += val;
                    }
                }
            }
        }

        // Check the cycles
        if (cyclevars.isEmpty()){
            System.out.println("No cycles are used.");
        } else {
            for (String cycle : cyclevars.keySet()){
                GRBVar var = cyclevars.get(cycle);
                double val = var.get(GRB.DoubleAttr.X);
                if (val > 1E-6){
                    System.out.println("Used cycle variable " + cycle + " with value " + val);
                    edgecount+= cyclelengths.get(cycle)*val;
                }
            }
        }

        System.out.println("Number of matches: " + edgecount);

        // DEBUG PRINT: total number of variables in the model
        System.out.println("Total variables in model: " + model.get(GRB.IntAttr.NumVars));
    }

    // A function that returns the pair constraint duals (alphas)
    public Map<Integer, Double> getAlphas(){
        return alphas;
    }

    // A function that returns the NDD constraint duals (gammas)
    public Map<Integer, Double> getGammas(){
        return gammas;
    }

    // A function that returns the flow duals (betas)
    public Map<Integer, Map<Integer, Double>> getBetas(){
        return betas;
    }

    // A function to find chains with positive reduced costs
    // So that they can be added as columns in the RMP
    // Based on the polynomial-time algorithm of Glorie et. al. (2014)
    public Map<Integer, List<Integer>> PricingProblemChains(Graph G, int C,
                                                           Map<Integer, Double> alphas,
                                                           Map<Integer, Double> gammas,
                                                          Map<Integer, Map<Integer, Double>> betas) {

        // Initialise a numerical infinity (divide by 4 for numerical stability)
        double INF = Double.MAX_VALUE / 4;

        // Initialise reduced costs for each edge variable
        // The reduced costs of a chain is the sum of reduced costs of edges in their respective positions
        // Structure: Map<e, Map< c, RC>>>
        Map<Edge, Map<Integer, Double>> varRCs = new HashMap<>();

        // Compute reduced costs from edge weights and duals
        for (Edge e : G.edges) {
            Map<Integer, Double> edgeRCs = new HashMap<>();
            Map<Integer, Double> ibetas = betas.get(e.from);
            Map<Integer, Double> jbetas = betas.get(e.to);

            // Case c=1: then i \in D
            // Then beta_j_1 exists but beta_i_0 does not
            // RC = w_{i,j} - alpha_j - gamma_i + beta_j_1
            if (G.NDDnodes.contains(e.from) && jbetas != null) {
                double beta_j_1 = jbetas.get(1);
                double RC1 = e.weight - alphas.get(e.to) - gammas.get(e.from) + beta_j_1;
                edgeRCs.put(1, RC1);
            }

            // Case c=2 to C-1: then i \in P
            // Then both beta_j_c and beta_i_c-1 exist
            // RC = w_{i,j} - alpha_j + beta_j_c - beta_i_c-1
            if (ibetas != null && jbetas != null) {
                for (int c = 2; c < C; c++) {
                    double beta_j_c = jbetas.get(c);
                    double beta_i_cminusone = ibetas.get(c - 1);
                    double RC = e.weight - alphas.get(e.to) + beta_j_c - beta_i_cminusone;
                    edgeRCs.put(c, RC);
                }
            }

            // Case c=C: then i \in P
            // Then beta_j_C does not exist since it is the last node in the chain (no out-flow)
            // Then beta_i_C-1 does exist
            // RC = w_{i,j} - alpha_j - beta_i_C-1
            if (ibetas != null) {
                double beta_i_Cminusone = ibetas.get(C - 1);
                double RC = e.weight - alphas.get(e.to) - beta_i_Cminusone;
                edgeRCs.put(C, RC);
            }

            varRCs.put(e, edgeRCs);
        }

        // Initialise a map of best chains for each NDD
        // Structure: Map<d, Chainstar_d>
        Map<Integer, List<Integer>> Chainstarmap = new HashMap<>();

        // Find the maximum cost (= minimum negative cost) chain for each NDD
        for (int d : G.NDDnodes){

            // Initialise the optimal chain
            List<Integer> chainstar = new ArrayList<>();

            // Step 1: Initialise variables
            // f_c^d(p) = weight of the shortest path in between d and p \in P using at most c edges
            // Structure = Map<c, Map<p, weight>>
            Map<Integer, Map<Integer, Double>> f_d = new HashMap<>();
            // g_c^d(p) = predecessor of node p \in P in such a shortest path from d to p
            // Structure = Map<c, Map<p, predecessor>>
            Map<Integer, Map<Integer, Integer>> g_d = new HashMap<>();

            //f_0^d(d) = 0; f_0^d(p) = INF for all p \in P
            Map<Integer, Double> f_d_0 = new HashMap<>();
            f_d_0.put(d, 0.0);
            for (int p : G.pairnodes){
                f_d_0.put(p, INF);
            }
            f_d.put(0, f_d_0);

            //g_0^d(p) = emptyset for d, and for all p \in P
            Map<Integer, Integer> g_d_0 = new HashMap<>();
            g_d_0.put(d, null);
            for (int p: G.pairnodes){
                g_d_0.put(p, null);
            }
            g_d.put(0, g_d_0);

            //Step 2: Recurrence loop to find all values f and g
            for (int c=0; c < C; c++){
                int position = c+1;
                Map<Integer, Double> f_d_cplus1 = new HashMap<>();
                Map<Integer, Integer> g_d_cplus1 = new HashMap<>();

                for (int p : G.pairnodes) {
                    // Edge (p', p) = arg min_{(u,p) \in E} {f_c^d(u) - RC(z_u_p_c+1)}
                    Integer pprime = null;
                    double bestweight = INF;

                    for (Edge e : G.edges) {
                        if (e.to == p) {
                            int u = e.from;

                            // Skip if no RC exists for this edge at position c+1
                            Map<Integer, Double> edgeRCs = varRCs.get(e);
                            if (edgeRCs == null || !edgeRCs.containsKey(position)) {
                                continue;
                            }

                            //TEST
                            // Skip if u is not reachable
                            Double f_u = f_d.get(c).get(u);
                            if (f_u == null || f_u >= INF/2){
                                continue;
                            }
                            // Compute the weight of the path of length c+1
                            double weight = f_d.get(c).get(u) - varRCs.get(e).get(position);

                            // Check if e.to is not already in the predecessor list (so loops are prevented)
                            boolean inpath = false;
                            int current = p;
                            for (int step = c; step >= 0; step--){
                                Integer predecessor = g_d.get(step).get(current);

                                // Stop if all nodes are checked
                                if (predecessor == null){
                                    break;
                                }
                                // Stop if u is in the path, then inpath = true
                                if (predecessor == p){
                                    inpath = true;
                                    break;
                                }
                                current = predecessor;
                            }

                            if (inpath){
                                continue;
                            }

                            if (weight < bestweight - 1e-6) {
                                pprime = u;
                                bestweight = weight;
                            }
                        }
                    }

                    // f_c+1^d(p) = min {f_c^d(p), f_c^d(p') - RC(z_p'_p_c+1)}
                    double f;
                    double leftside = f_d.get(c).get(p);
                    double rightside = (pprime != null)? bestweight : INF;
                    if ( leftside < rightside - 1e-6 ){
                        f = leftside;
                    } else{
                        f = rightside;
                    }
                    f_d_cplus1.put(p, f);

                    // g_c+1^d(p) = p' if f_c^d(p') - RC(x_p'_p_c+1) < f_c^d(p) and g_c^d(p) otherwise
                    Integer g;
                    if ( pprime != null &&  rightside < leftside - 1e-6 ){
                        g = pprime;
                    } else{
                        g = g_d.get(c).get(p);
                    }
                    g_d_cplus1.put(p, g);
                }
                f_d.put(c+1, f_d_cplus1);
                g_d.put(c+1, g_d_cplus1);
            }

            // Step 3: Retrieve the maximum cost chain
            double chainstarcost = 0;

            for (int p: G.pairnodes){

                // If f_C^d(p) < 0, then there is a positive cost chain for NDD node d
                // By initialising chainstarcost =0, we automatically find negative minimum costs if they exist
                double chaincost = f_d.get(C-1).get(p);
                if (chaincost < chainstarcost -1e-4){

                    // Backtrack the chain from node p using predecessors g
                    List<Integer> chain = new ArrayList<>();
                    chain.add(d);

                    int current = p;
                    for (int step = C-1; step >=1; step--){
                        // Add the current node to the chain in the second position
                        // So that the chain is not reversed
                        chain.add(1,current);
                        Integer predecessor = g_d.get(step).get(current);
                        if (predecessor == null){
                            // If there are no predecessors, then the first node is reached
                            break;
                        }
                        current = predecessor;
                    }

                    // Remove the duplicate nodes
                    // If a node is its own predecessor then the chain is shorter than length C
                    Set<Integer> seen = new LinkedHashSet<>();
                    seen.add(d);
                    List<Integer> trimmedchain = new ArrayList<>();
                    trimmedchain.add(d);

                    for (int i =1; i < chain.size(); i++){
                        int node = chain.get(i);
                        if (!seen.contains(node)){
                            trimmedchain.add(node);
                            seen.add(node);
                        }
                    }

                    // Update the best found chain
                    chainstar = trimmedchain;
                    chainstarcost = -chaincost;
                }
            }
            // Print the best found chain if it is nonempty and not found earlier
            if (!chainstar.isEmpty() && !addedchains.contains(chainstar.toString()) && chainstarcost > 1e-4 ){
                System.out.println("Found maximum reduced cost chain: " + chainstar);
                System.out.println("With cost " + chainstarcost);
            }

            // Add best chain to the map only if reduced costs are positive (double-check)
            // And if it is not already added
            // Sometimes RC is postive multiple times if there are multiple optimal solutions
            if (chainstarcost > 1e-4 && !addedchains.contains(chainstar.toString())){
                Chainstarmap.put(d, chainstar);
            }
        }
        return Chainstarmap;
    }

    // A function to find cycles with positive reduced costs
    // In a similar manner as for the chains
    public Map<Integer, List<Integer>> PricingProblemCycles(Graph G, int K,
                                                            Map<Integer, Double> alphas) {

        // Initialise a numerical infinity (divide by 4 for numerical stability)
        double INF = Double.MAX_VALUE / 4;

        // Initialise reduced costs for each edge, such that reduced costs of cycles can be computed
        // The reduced costs of a cycle is the sum of reduced costs of its edges
        Map<Edge, Double> edgeRCs = new HashMap<>();

        // Compute reduced costs for all edges
        // RC((i,j)) = w_{i,j} - alpha_j
        for (Edge e : G.edges){
            // Skip the RCs for self-cycles (loops)
            if (e.from == e.to){
                continue;
            }
            double RC = e.weight - alphas.get(e.to);
            edgeRCs.put(e, RC);
        }

        // Initialise a map of best cycles in each subgraph
        // Structure: Map<l, Cstar_l>
        Map<Integer, List<Integer>> Cstarmap = new HashMap<>();

        // Find the maximum cost (= minimum negative cost) cycle for each node (only looking in their subgraph)
        for (int l : G.pairnodes){
            Graph subgraph = PIEF.subgraph(G,l);
            List<Integer> nodes = subgraph.pairnodes;
            List<Edge> edges = subgraph.edges;

            // Initialise the optimal cycle
            List<Integer> Cstar = new ArrayList<>();
            double Cstarcost = 0;

            // Bellman-ford algorithm with path length restriction

            // Step 1: initialise the distances and predecessors
            Map<Integer, Double> dist = new HashMap<>();
            Map<Integer, Map<Integer, Integer>> g_l = new HashMap<>();

            // For all nodes, initialise the distance as INF (except for node l itself)
            for (int p : nodes){
                if (p != l){
                    dist.put(p, INF);
                }
            }
            // The distance from node l to itself is zero
            dist.put(l, 0.0);

            //g_0^l(p) = emptyset for all p \in V^l
            Map<Integer, Integer> g_l_0 = new HashMap<>();
            for (int p: nodes){
                g_l_0.put(p, null);
            }
            g_l.put(0, g_l_0);

            // Step 2: Relax the edges until paths of length K-1 are reached
            for (int k=0; k < K-1; k++){
                // Make a copy of the dist map of the previous iteration
                // Such that a path of length k is not extended in the same iteration
                Map<Integer, Double> distcopy = new HashMap<>(dist);
                // Initialise the predecessor map for the next iteration
                Map<Integer, Integer> g_l_kplus1 = new HashMap<>();

                // Update the distances (costs) and predecessors
                for (Edge e : edges){
                    // Skip self-cycles (loops)
                    if (e.from == e.to){
                        continue;
                    }
                    // Skip unreachable nodes
                    if (distcopy.get(e.from) >= INF/2){
                        continue;
                    }
                    // Skip if the reduced costs are null (for loops)
                    if (edgeRCs.get(e) == null){
                        continue;
                    }

                    // Compute the weight (cost) of the path of length k
                    double newcost = distcopy.get(e.from) - edgeRCs.get(e);

                    // Check if e.to is not already in the predecessor list (so loops are prevented)
                    boolean inpath = false;
                    int current = e.to;
                    for (int step = k; step >= 0; step--){
                        Integer predecessor = g_l.get(step).get(current);

                        // Stop if all nodes are checked
                        if (predecessor == null){
                            break;
                        }
                        // Stop if p is in the path, then inpath = true
                        if (predecessor == e.to){
                            inpath = true;
                            break;
                        }
                        current = predecessor;
                    }

                    if (inpath){
                        continue;
                    }

                    if (newcost < dist.get(e.to) - 1e-6){
                        dist.put(e.to, newcost);
                        g_l_kplus1.put(e.to, e.from);
                    }
                }
                // Add the predecessor map for k+1 to g_l
                g_l.put(k+1, g_l_kplus1);
            }

            // Step 3: Retrieve the maximum cost cycle for node l
            for (Edge e : edges){
                // if (p, l) in edges and dist(p) - RC((p,l)) <0, then there is a positive cost cycle for node l
                // By initialising Cstarcost =0, we automatically find negative minimum costs if they exist
                if (e.to == l && dist.get(e.from) < INF/2){

                    // Skip if no RC exists for this edge (for loops)
                    if (edgeRCs.get(e) == null){
                        continue;
                    }

                    // Cost of the cycle = cost of path to p + RC((p,l))
                    double cyclecost = dist.get(e.from) - edgeRCs.get(e);
                    if (cyclecost < Cstarcost - 1e-4){

                        // Backtrack the cycle from p to l using pred
                        List<Integer> cycle = new ArrayList<>();
                        cycle.add(l);

                        // Initialise current as p
                        int current = e.from;
                        for (int step = K-1; step >=1; step--){
                            // Add the current node to the cycle in the second position
                            // So that the cycle is not reversed
                            cycle.add(1,current);
                            Integer predecessor = g_l.get(step).get(current);
                            if (predecessor == null){
                                // If there are no predecessors, then the first node is reached
                                break;
                            }
                            current = predecessor;
                        }

                        // Remove the duplicate nodes
                        // If a node is its own predecessor then the cycle is shorter than length K
                        // Initialise a set of seen nodes and add node l
                        Set<Integer> seen = new LinkedHashSet<>();
                        seen.add(l);

                        // Initialise the cycle without duplicates and add node l
                        List<Integer> trimmedcycle = new ArrayList<>();
                        trimmedcycle.add(l);

                        for (int i =1; i < cycle.size(); i++){
                            int node = cycle.get(i);
                            if (!seen.contains(node)){
                                // Only add the node if it has not been added yet, and add it to the seen set
                                trimmedcycle.add(node);
                                seen.add(node);
                            }
                        }

                        // Close the cycle by returning to node l
                        trimmedcycle.add(l);
                        Cstar = trimmedcycle;
                        Cstarcost = - cyclecost;
                    }
                }
            }
            // Print the best found cycle if it is nonempty and not found earlier
            if (!Cstar.isEmpty() && !addedcycles.contains(Cstar.toString()) && Cstarcost > 1e-4){
                System.out.println("Found maximum reduced cost cycle: " + Cstar);
                System.out.println("With cost " + Cstarcost);
            }

            // Add best cycle to the map only if reduced costs are positive (double-check)
            // And if it is not already added
            // Sometimes RC is postive multiple times if there are multiple optimal solutions
            if (Cstarcost > 1e-4 && !addedcycles.contains(Cstar.toString())){
                Cstarmap.put(l, Cstar);
            }
        }
        return Cstarmap;
    }

    public void addchainColumns (Graph G, Map<Integer, List<Integer>> chainstarmap, int C)
            throws GRBException {
        //EXAMPLE: for a chain c=(d,a,b,c,e) of NDD d \in D:
        // Add to variables: z_d_a_1, z_a_b_2, z_b_c_3, z_c_e_4
        // Add to pair capacity constraints: +z_d_a_1 for a, +  z_a_b_2 for b, z_b_c_3 for c, + z_c_e_4 for e
        // Add to NDD capacity constraints: + z_d_a_1 for d
        // Add to flow constraints:
        // + z_d_a_1 - z_a_b_2 for node a, and k=1
        // + z_a_b_2 - z_b_c_3 for node b, and k=2
        // + z_b_c_3 - z_c_e_4 for node c, and k=3
        // Add to objective: + w_d_a*z_d_a_1 + w_a_b*z_a_b_2 + w_b_c*z_b_c_3 + w_c_e*z_c_e_4

        for (int d : G.NDDnodes){
            // Obtain the chain from this NDD
            List<Integer> chain = chainstarmap.get(d);

            // Skip if no chain was found for this NDD
            if (chain == null || chain.size() <=1 ){
                continue;
            }

            // DEBUG PRINT
            System.out.println("Already added chains: " + addedchains);
            System.out.println("Checking chain  " + chain);

            // Skip if the chain is already in the RMP
            String chainKey = chain.toString();
            if (addedchains.contains(chainKey)){
                System.out.println("Skipping chain " + chain + " since it is already added");
                continue;
            }
            addedchains.add(chainKey);

            // Loop over edges in the chain
            for (int c =0; c < chain.size()-1; c++){
                int i = chain.get(c);
                int j = chain.get(c+1);
                int position = c+1;

                /*// DEBUG PRINT: flow constraints for this edge
                System.out.println("  Capacity constraint for j=" + j + ": " + (paircapconstrs.get(j) != null));
                System.out.println("  Flow inflow j=" + j + " pos=" + position + ": " +
                        (j != d && position < C && flowconstrs.get(j) != null &&
                                flowconstrs.get(j).containsKey(position)));
                System.out.println("  Flow outflow i=" + i + " pos=" + (position-1) + ": " +
                        (i != d && position > 1 && flowconstrs.get(i) != null &&
                                flowconstrs.get(i).containsKey(position-1)));
                 */

                // Find the corresponding edge and it's weight
                Edge edge = null;
                double weight = 0;
                for (Edge e : G.edges){
                    if (e.from == i && e.to == j){
                        edge = e;
                        weight = e.weight;
                        break;
                    }
                }

                // Check if the edge variable already exists and skip if it is a duplicate
                String varName = "z_" + i + "_" + j + "_" + position;
                if (addedvars.contains(varName)){
                    System.out.println("Skipping duplicate variable: " + varName);
                    continue;
                } else {
                    System.out.println("Adding edge z_" + i + "_" + j + "_" + position);
                }
                addedvars.add(varName);

                // Initialise gurobi column
                GRBColumn column = new GRBColumn();

                // Add the edge variable +z_i_j_k to the pair capacity constraint of j
                column.addTerm(1.0, paircapconstrs.get(j));

                // Add the edge variable +z_d_j_1 to the NDD capacity constraint of d
                if (G.NDDnodes.contains(i) && position == 1){
                    column.addTerm(1.0, NDDcapconstrs.get(i));
                }

                // Add the edge to the flow constraints:
                // Inflow to j at position c: + z_i_j_c in flow constr of (j,c) and j != d (double-check)
                if (position < C && j!= d){
                    Map<Integer, GRBConstr> jconstrs = flowconstrs.get(j);
                    if (jconstrs != null && jconstrs.containsKey(position)){
                        column.addTerm(1.0, jconstrs.get(position));
                    }
                }
                // Outflow from i at position c+1: -z_i_j_c+1 at flow constr of (i,c) and i != d (double-check)
                if (position > 1  && i!=d){
                    Map<Integer, GRBConstr> iconstrs = flowconstrs.get(i);
                    if (iconstrs != null && iconstrs.containsKey(position-1)){
                        column.addTerm(-1.0, iconstrs.get(position-1));
                    }
                }

                // Add the variable to the model with objective coefficient
                GRBVar edgevar = model.addVar(0, 1, weight, GRB.CONTINUOUS, column,
                        "z_"+ i + "_" + j + "_" + position);
                //Store the edge variable
                edgevars.putIfAbsent(edge, new ArrayList<>());
                edgevars.get(edge).add(edgevar);

                // Update the model
                model.update();
            }

        }
    }

    public void addcycleColumns (Graph G, Map<Integer, List<Integer>> Cstarmap)
            throws GRBException {
        //EXAMPLE: for a cycle c=(l,a,b,l) of G^l:
        // Add to variables: x_c
        // Add to capacity constraints: +x_c for all nodes i in c (so for l, a, and b)
        // Add to objective: + u_c*x_c

        for (int l : G.pairnodes){
            // Obtain the cycle from this subgraph
            List<Integer> cycle = Cstarmap.get(l);

            // Skip if no cycle was found for this subgraph
            if (cycle == null || cycle.size() <=1 ){
                continue;
            }

            // DEBUG PRINT
            System.out.println("Already added cycles: " + addedcycles);
            System.out.println("Checking cycle: " + cycle);

            // Skip if the cycle is already in the RMP
            String cycleKey = cycle.toString();
            if (addedcycles.contains(cycleKey)){
                System.out.println("Skipping cycle " + cycle + " since it is already added");
                continue;
            }
            addedcycles.add(cycleKey);

            // Add the variable x_c for this cycle
            System.out.println("Adding cycle: " + cycleKey);

            // Initialise the Gurobi column
            GRBColumn column = new GRBColumn();

            // Initialise the weight of the cycle
            double cycleweight = 0;

            // Initialise the length of the cycle
            int cyclelength = 0;

            // Loop over nodes in the cycle
            for (int k =0; k < cycle.size()-1; k++) {
                int i = cycle.get(k);
                int j = cycle.get(k + 1);

                // Add the cycle variable to the capacity constraint of j
                column.addTerm(1.0, paircapconstrs.get(j));

                // Update the cycle length
                cyclelength++;

                // Find the corresponding edge (i,j) and it's weight
                double weight = 0;
                for (Edge e : G.edges) {
                    if (e.from == i && e.to == j) {
                        weight = e.weight;
                        break;
                    }
                }
                // Add the weight of the edge to the cycle weight
                cycleweight += weight;
            }

            // Add the variable to the model with objective coefficient
            GRBVar cyclevar = model.addVar(0, 1, cycleweight, GRB.CONTINUOUS, column,
                    cycleKey);

            //Store the cycle variable and its length
            cyclevars.put(cycleKey, cyclevar);
            cyclelengths.put(cycleKey, cyclelength);

            // Update the model
            model.update();

        }
    }
}
