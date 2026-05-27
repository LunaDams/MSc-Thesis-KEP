// Packages used.
import java.util.*;
import com.gurobi.gurobi.*;

public class PIEF {

    // This class is for the PIEF (cycle-only) algorithm
    // It produces solutions for the KEP when no NDDs are available

    // An environment so that we can solve the RMP
    public GRBEnv env;
    // A GRB model so that we can solve the RMP
    public GRBModel model;

    // A multidimensional map of edge variables (one for each subgraph, edge within the subgraph and position)
    private Map<Integer, Map<Edge, List<GRBVar>>>  edgevars;
    // Maps of loop variables (one for each node) and corresponding loop weights, with indexing
    private Map<Integer, GRBVar> loopvars;
    private Map<Integer, Double> loopweights;

    // A map of capacity constraints (one for each node), with node indices
    private Map<Integer, GRBConstr> capconstrs;
    // A multidimensional map of flow constraints (one for each subgraph, node within the subgraph and position)
    // Structure: Map<l, Map< i, Map<k, constr>>>
    private Map<Integer, Map<Integer, Map<Integer, GRBConstr>>> flowconstrs;

    // A map of capacity duals (one for each node), with corresponding index
    private Map<Integer, Double> alphas;
    // A multidimensional map of flow duals (one for each subgraph, node within the subgraph and position)
    private Map<Integer, Map<Integer, Map<Integer, Double>>> betas;

    // A set of edge variables that have been added to the RMP
    private Set<String> addedvars = new HashSet<>();
    // A set of cycles that have been added to the RMP
    private Set<String> addedcycles = new HashSet<>();

    // A function to initialise the PIEF RMP
    public void initialisePIEF(Graph G, int K) throws GRBException{

        // Initialise the environment
        env =new GRBEnv(true);
        env.set("logfile", "rmp.log");
        env.start();

        // Create empty model
        model = new GRBModel(env);

        // Initialise the list of edge variables, but do not add them yet
        edgevars = new HashMap<>();

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

        // Update the model
        model.update();

        // Initialise the capacity constraint
        capconstrs = new HashMap<>();
        // (only add the loop variables since there are no edge variables yet)
        for (int p : G.pairnodes){
            GRBLinExpr expr = new GRBLinExpr();

            // Add the loop to the constraint if it exists
            if ( loopvars.containsKey(p)){
                    expr.addTerm(1.0 ,loopvars.get(p));
            }

            // Add the constraint to the model
            GRBConstr capconstr = model.addConstr(expr, GRB.LESS_EQUAL, 1.0, "capacity" + p);
            capconstrs.put(p, capconstr);
        }

        // Initialise the flow constraints
        flowconstrs = new HashMap<>();
        for (int l : G.pairnodes){
            Map<Integer, Map<Integer, GRBConstr>> subgraphconstrs = new HashMap<>();
            Graph subgraph = subgraph(G, l);

            for (int i : subgraph.pairnodes){
                // Only add flow constraints for i >= l+1
                if (i == l){
                    continue;
                }
                Map<Integer, GRBConstr> nodeflowconstrs = new HashMap<>();
                for (int k=1; k < K; k++){
                    GRBLinExpr expr = new GRBLinExpr();
                    // Do not add any edges yet
                    // When sums are added, use the form Sum_in - Sum_out = 0
                    GRBConstr flowconstr = model.addConstr(expr, GRB.EQUAL, 0.0,
                            "flow^" + l + "_" + i +  "_" + k);
                    nodeflowconstrs.put(k, flowconstr);
                }
                subgraphconstrs.put(i, nodeflowconstrs);
            }
            // Store the list of subgraph constraints with the node ID as the key
            flowconstrs.put (l, subgraphconstrs);
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
    public void solvePIEF(Graph G, int K) throws GRBException{

        // Use interior point solving instead of simplex
        // Because of degeneracy possibility
        model.set(GRB.IntParam.Method, 2);

        // Optimise the model
        model.optimize();

        // Initialise the duals
        alphas = new HashMap<>();
        betas = new HashMap<>();

       // Obtain the capacity duals for each node
        for (int p : G.pairnodes){
            GRBConstr capconstr = capconstrs.get(p);
            double alpha = capconstr.get(GRB.DoubleAttr.Pi);
            alphas.put(p, alpha);
        }

        // Obtain the flow duals for each subgraph, node in the subgraph and position
        for (int l : G.pairnodes){
            Map<Integer, Map<Integer, GRBConstr>> subgraphconstrs = flowconstrs.get(l);
            Graph subgraph = subgraph(G, l);
            Map<Integer, Map<Integer, Double>> subgraphbetas = new HashMap<>();
            for (int i : subgraph.pairnodes){
                if (i == l){
                    continue;
                }
                Map<Integer, GRBConstr> nodeflowconstrs = subgraphconstrs.get(i);
                Map<Integer, Double> nodebetas = new HashMap<>();
                for (int k=1; k < K; k++){
                    GRBConstr flowconstr = nodeflowconstrs.get(k);
                    double beta = flowconstr.get(GRB.DoubleAttr.Pi);
                    nodebetas.put(k, beta);
                }
                subgraphbetas.put(i, nodebetas);
            }
            // Store the duals in a map with the node ID as key
            betas.put(l, subgraphbetas);
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

        /*// DEBUG PRINT of flow constraint, for subgraph l=1, node i=6
        GRBConstr testConstr = flowconstrs.get(1).get(6).get(1);
        System.out.println("Flow constraint l=1, i=6, k=1: " + testConstr);
        System.out.println("  RHS: " + testConstr.get(GRB.DoubleAttr.RHS));
        System.out.println("  Dual: " + testConstr.get(GRB.DoubleAttr.Pi));
        System.out.println("  Slack: " + testConstr.get(GRB.DoubleAttr.Slack));

        // Check if duals are actually non-zero
        System.out.println("Alpha for node 1: " + alphas.get(1));
        System.out.println("Alpha for node 6: " + alphas.get(6));
        System.out.println("Beta for l=1, i=6, k=1: " + betas.get(1).get(6).get(1));
        System.out.println("Beta for l=1, i=8, k=3: " + betas.get(1).get(8).get(3));
        */

        // Check the used variables in the solution
        System.out.println("Solution values: ");

        // Initialise the number of matches
        double edgecount = 0;

        // Check the loops first
        for (int p : loopvars.keySet()){
            GRBVar var = loopvars.get(p);
            String name = var.get(GRB.StringAttr.VarName);
            double val = var.get(GRB.DoubleAttr.X);
            if (val > 1e-6){
                System.out.println("Used loop variable " + name + " with value " + val);
                edgecount += val;
            }
        }

        // Check the edges
        if (edgevars.isEmpty()){
            System.out.println("No edges are used.");
        } else {
            for (int l : edgevars.keySet()){
                Map<Edge, List<GRBVar>> subgraph = edgevars.get(l);

                for (Edge e : subgraph.keySet()){
                    List<GRBVar> nodelist = subgraph.get(e);

                    for (int k=0; k < nodelist.size(); k++){
                        GRBVar var = nodelist.get(k);
                        String name = var.get(GRB.StringAttr.VarName);
                        double val = var.get(GRB.DoubleAttr.X);

                        if (val > 1e-6){
                            System.out.println("Used edge variable " + name + " with value " + val);
                            edgecount += val;
                        }
                    }
                }
            }
        }
        System.out.println("Number of matches: " + edgecount);

        // DEBUG PRINT: total number of variables in the model
        System.out.println("Total variables in model: " + model.get(GRB.IntAttr.NumVars));
    }


    // A function that makes subgraphs of the compatibility graph G for every node l in P
    public static Graph subgraph(Graph G, int l){

        // Only use the nodes with index >=l
        List<Integer> subnodes = new ArrayList<>();
        for (int p : G.pairnodes){
            if (p >= l){
                subnodes.add(p);
            }
        }

        // Only use the edges from the nodes in the subgraph
        // Exclude NDDs since those will not be included in the cycle
        List<Edge> subedges = new ArrayList<>();
        for (Edge e : G.edges){
            if ((e.from >= l && e.to >= l) && !G.NDDnodes.contains(e.from)){
                subedges.add(e);
            }
        }

        return new Graph(subnodes, null, subedges);

    }

    // A function that returns the constraint duals (alphas)
    public Map<Integer, Double> getAlphas(){
        return alphas;
    }

    // A function that returns the flow duals (betas)
    public Map<Integer, Map<Integer, Map<Integer, Double>>> getBetas(){
        return betas;
    }

    // A function to find cycles with positive reduced costs
    // So that they can be added as columns in the RMP
    // Based on the polynomial-time algorithm of Glorie et. al. (2014)
    public Map<Integer, List<Integer>> PricingProblemPIEF(Graph G, int K, Map<Integer, Double> alphas,
                                                          Map<Integer, Map<Integer, Map<Integer, Double>>> betas){

        // Initialise a numerical infinity (divide by 4 for numerical stability)
        double INF = Double.MAX_VALUE/4;

        // Initialise reduced costs for each variable
        // The reduced costs of a cycle is the sum of reduced costs of edges in their respective positions
        // Structure: Map< l, Map<e, Map< k, RC>>>
        Map<Integer, Map<Edge, Map<Integer, Double>>> varRCs = new HashMap<>();

        // Obtain reduced costs from edge weights and duals
        for (int l : G.pairnodes){
            Map<Edge, Map<Integer, Double>> subgraphRCs = new HashMap<>();
            Map<Integer, Map<Integer, Double>> subgraphbetas = betas.get(l);
            Graph subgraph = subgraph(G,l);

            for (Edge e: subgraph.edges){
                Map<Integer, Double> edgeRCs = new HashMap<>();
                Map<Integer, Double> ibetas = subgraphbetas.get(e.from);
                Map<Integer, Double> jbetas = subgraphbetas.get(e.to);

                // Case k=1: beta_j_1^l exists but beta_i_0^l does not
                // RC = w_{i,j} - alpha_j + beta_j_1^l
                if (jbetas != null){
                    double beta_j_1 = jbetas.get(1);
                    double RC1 = e.weight - alphas.get(e.to) + beta_j_1;
                    edgeRCs.put(1, RC1);
                }

                // Case k=2 to K-1:
                // RC = w_{i,j} - alpha_j + beta_j_k^l - beta_i_k-1^l
                if (ibetas != null && jbetas != null){
                    for (int k= 2; k < K; k++){
                        double beta_j_k = jbetas.get(k);
                        double beta_i_kminusone = ibetas.get(k-1);
                        double RC = e.weight - alphas.get(e.to) +beta_j_k - beta_i_kminusone;
                        edgeRCs.put(k, RC);
                    }
                }

                // Case k=K: the edge must return to node l
                // beta_j_K^l does not exist since betas only go up to node K-1
                // But beta_i_K-1^l does exist
                // RC = w_{i,j} - alpha_j - beta_i_K-1^l
                if (ibetas != null){
                    double beta_i_Kminusone = ibetas.get(K-1);
                    double RCK = e.weight -alphas.get(l) - beta_i_Kminusone;
                    edgeRCs.put(K, RCK);
                }

                subgraphRCs.put(e, edgeRCs);
           }

           varRCs.put(l, subgraphRCs);
        }

        // Initialise a map of best cycles in each subgraph
        // Structure: Map<l, Cstar_l>
        Map<Integer, List<Integer>> Cstarmap = new HashMap<>();

        // Find the maximum cost (= minimum negative cost) cycle in each subgraph
        for (int l : G.pairnodes){
            Graph subgraph = subgraph(G,l);
            List<Integer> nodes = subgraph.pairnodes;
            List<Edge> edges = subgraph.edges;

            // Initialise the optimal cycle
            List<Integer> Cstar = new ArrayList<>();

            // Step 1: Initialise variables
            // f_k^l(p) = weight of the shortest path in between l and p \in V^l using at most k edges
            // Structure = Map<k, Map<p, weight>>
            Map<Integer, Map<Integer, Double>> f_l = new HashMap<>();
            // g_k^l(p) = predecessor of node p \in V^l in such a shortest path from l to p
            // Structure = Map<k, Map<p, predecessor>>
            Map<Integer, Map<Integer, Integer>> g_l = new HashMap<>();

            //f_0^l(l) = 0; f_0^l(p) = INF for all p \in V^l\{l}
            Map<Integer, Double> f_l_0 = new HashMap<>();
            f_l_0.put(l, 0.0);
            for (int p : nodes){
                if (p != l){
                    f_l_0.put(p, INF);
                }
            }
            f_l.put(0, f_l_0);

            //g_0^l(p) = emptyset for all p \in V^l
            Map<Integer, Integer> g_l_0 = new HashMap<>();
            for (int p: nodes){
                g_l_0.put(p, null);
            }
            g_l.put(0, g_l_0);

            // Obtain the RCs for this subgraph
            Map<Edge, Map<Integer, Double>> subgraphRCs = varRCs.get(l);

            //Step 2: Recurrence loop to find all values f and g
            for (int k =0; k < K-1; k++){
                int position = k+1;
                Map<Integer, Double> f_l_kplus1 = new HashMap<>();
                Map<Integer, Integer> g_l_kplus1 = new HashMap<>();

                for (int p : nodes){
                    // Edge (p', p) = arg min_{(u,p) \in E^l} {f_k^l(u) - RC(x_u_p_k+1^l)}
                    Integer pprime = null;
                    double bestweight = INF;

                    for (Edge e : edges){
                        if (e.to == p){
                            int u =e.from;

                            // Skip if no RC exists for this edge at position k+1
                            Map<Integer, Double> edgeRCs = subgraphRCs.get(e);
                            if (edgeRCs == null || !edgeRCs.containsKey(position)){
                                continue;
                            }

                            // Compute the weight of the path to node p
                            double weight = f_l.get(k).get(u) - subgraphRCs.get(e).get(position);

                            // Check if u is not already in the predecessor list (so loops are prevented)
                            boolean inpath = false;
                            int current = p;
                            for (int step = k; step >= 0; step--){
                                Integer predecessor = g_l.get(step).get(current);

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


                            if (weight < bestweight - 1e-6){
                                pprime = u;
                                bestweight = weight;
                            }
                        }
                    }

                    // f_k+1^l(p) = min {f_k^l(p), f_k^l(p') - RC(x_p'_p_k+1^l)}
                    double f;
                    double leftside = f_l.get(k).get(p);
                    double rightside = (pprime != null)? bestweight : INF;
                    if ( leftside < rightside - 1e-6 ){
                        f = leftside;
                    } else{
                        f = rightside;
                    }
                    f_l_kplus1.put(p, f);

                    // g_k+1(p) = p' if f_k^l(p') - RC(x_p'_p_k+1^l) < f_k^l(p) and g_k^l(p) otherwise
                    Integer g;
                    if ( pprime != null &&  rightside < leftside - 1e-6 ){
                        g = pprime;
                    } else{
                        g = g_l.get(k).get(p);
                    }
                    g_l_kplus1.put(p, g);
                }
                f_l.put(k+1, f_l_kplus1);
                g_l.put(k+1, g_l_kplus1);
            }

            // Step 3: Retrieve the maximum cost cycle
            double Cstarcost = 0;
            for (int p : nodes){

                // Exlude self-cycles (loops)
                if (p==l){
                    continue;
                }

                // If (p,l) in edges and f_K-1^l(p) - RC(x_p_l^l(K) < 0, then there is a positive cost cycle
                // By initialising Cstarcost =0, we automatically find negative minimum costs if they exist
                for (Edge e : edges){
                    if (e.from == p && e.to ==l){

                        // Skip if no RC exists for this edge at position k+1
                        Map<Integer, Double> edgeRCs = subgraphRCs.get(e);
                        if (edgeRCs == null || !edgeRCs.containsKey(K)){
                            continue;
                        }

                        double cyclecost = f_l.get(K-1).get(p) - subgraphRCs.get(e).get(K);
                        if (cyclecost < Cstarcost - 1e-4){

                            // Backtrack the cycle from p using predecessors g
                            List<Integer> cycle = new ArrayList<>();
                            cycle.add(l);

                            int current = p;
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
                            Cstarcost = -cyclecost;
                        }
                    }
                }
            }
            // Print the best found cycle if it is nonempty and not found earlier
            if (!Cstar.isEmpty() && !addedcycles.contains(Cstar.toString()) && Cstarcost > 1e-4 ){
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

    public void addColumnsPIEF (Graph G, Map<Integer, List<Integer>> Cstarmap, int K)
            throws GRBException {
        //EXAMPLE: for a cycle c=(l,a,b,l) of G^l:
        // Add to variables: x_l_a_1^l, x_a_b_2^l, x_b_l_3^l
        // Add to capacity constraints: +x_l_a_1^l for a, +  x_a_b_2^l for b, x_b_l_3^l for l
        // Add to flow constraints:
        // + x_l_a_1^l - x_a_b_2^l for l, node a, and k=1
        // + x_a_b_2^l - x_b_l_3^l for l, node b, and k=2
        // Add to objective: + w_l_a*x_l_a_1^l + w_a_b*x_a_b_2^l + w_b_l*x_b_l_3^l

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

            // Add all edges to the edge variables
            // Initialise the edge variable if it does not exist in the edgevars map yet
            edgevars.putIfAbsent(l, new HashMap<>());

            // Loop over edges in the cycle
            for (int k =0; k < cycle.size()-1; k++){
                int i = cycle.get(k);
                int j = cycle.get(k+1);
                int position = k+1;

                /*// DEBUG PRINT: flow constraints for this edge
                System.out.println("  Capacity constraint for j=" + j + ": " + (capconstrs.get(j) != null));
                System.out.println("  Flow inflow j=" + j + " pos=" + position + ": " +
                        (j != l && position < K && flowconstrs.get(l).get(j) != null &&
                                flowconstrs.get(l).get(j).containsKey(position)));
                System.out.println("  Flow outflow i=" + i + " pos=" + (position-1) + ": " +
                        (i != l && position > 1 && flowconstrs.get(l).get(i) != null &&
                                flowconstrs.get(l).get(i).containsKey(position-1)));
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
                String varName = "x_" + i + "_" + j + "_" + position + "^" + l;
                if (addedvars.contains(varName)){
                    System.out.println("Skipping duplicate variable: " + varName);
                    continue;
                } else {
                    System.out.println("Adding edge x_" + i + "_" + j + "_" + position + "^" + l);
                }
                addedvars.add(varName);

                // Initialise gurobi column
                GRBColumn column = new GRBColumn();

                // Add the edge variable +x_i_j_k^l to the capacity constraint of j
                column.addTerm(1.0, capconstrs.get(j));

                // Add the edge to the flow constraints:
                // Inflow to j at position k: + x_i_j_k^l in flow constr of (l,j,k) and j != l
                if (position < K && j!= l){
                    Map<Integer, GRBConstr> jconstrs = flowconstrs.get(l).get(j);
                    if (jconstrs != null && jconstrs.containsKey(position)){
                        column.addTerm(1.0, jconstrs.get(position));
                    }
                }
                // Outflow from i at position k+1: -x_i_j_k+1^l at flow constr of (l,i,k) and i != l
                if (position > 1  && i!=l){
                    Map<Integer, GRBConstr> iconstrs = flowconstrs.get(l).get(i);
                    if (iconstrs != null && iconstrs.containsKey(position-1)){
                        column.addTerm(-1.0, iconstrs.get(position-1));
                    }
                }

                // Add the variable to the model with objective coefficient
                GRBVar edgevar = model.addVar(0, 1, weight, GRB.CONTINUOUS, column,
                        "x_"+ i + "_" + j + "_" + position +"^" + l);
                //Store the edge variable
                edgevars.get(l).putIfAbsent(edge, new ArrayList<>());
                edgevars.get(l).get(edge).add(edgevar);

                // Update the model
                model.update();
            }

        }
    }
}
