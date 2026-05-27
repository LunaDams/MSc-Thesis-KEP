// Current version - to be possibly updated

// Packages used.
import java.util.*;
import com.gurobi.gurobi.*;

public class CG {

    // This class is used for finding LP solutions using a column generation algorithm
    // Functions used are obtained from the PIEF and PICEF class


    // A function that finds CG LP solutions for PIEF
    public static void CGPIEF(Graph G, int K){
        PIEF pief = new PIEF();

        try {
            // Initialise the running time
            long start = System.nanoTime();

            // Initialise the RMP
            pief.initialisePIEF(G, K);

            // Initialise the number of iterations
            int iterations = 0;

            // Set a boolean to induce the CG loop
            boolean improving = true;

            // Column generation loop
            while (improving){
                iterations++;

                // Solve the RMP
                pief.solvePIEF(G, K);

                // Obtain the duals
                Map<Integer, Double> alphas = pief.getAlphas();
                Map<Integer, Map<Integer, Map<Integer, Double>>> betas = pief.getBetas();

                // Solve the pricing problem and save positive cost cycles
                Map<Integer, List<Integer>> Cstarmap =
                        pief.PricingProblemPIEF(G, K, alphas, betas);

                // Stop if there are no more positive cost cycles
                if (Cstarmap.isEmpty()) {
                    improving = false;

                    // Save the final time and compute the runtime in seconds
                    long end = System.nanoTime();
                    double runtime = (end - start)/1e9;

                    // Print termination statement
                    System.out.println("No more positive reduced cost cycles, final solution found!");
                    System.out.println("Number of iterations: " + iterations);
                    System.out.println("Total runtime = " + runtime + " seconds");

                } else{
                    // Add the improving cycles as columns to the RMP
                    pief.addColumnsPIEF(G, Cstarmap, K);
                }
            }


        } catch (GRBException e) {
            throw new RuntimeException(e);
        }
    }

    public static void CGPICEF(Graph G, int K, int C){
        PICEF picef = new PICEF();
        try {
            // Initialise the running time
            long start = System.nanoTime();

            // Initialise the RMP
            picef.initialisePICEF(G, C);

            // Initialise the number of iterations
            int iterations = 0;

            // Set a boolean to induce the CG loop
            boolean improving = true;

            // Column generation loop
            while (improving){
                iterations++;

                // Solve the RMP
                System.out.println("Iteration " + iterations);
                picef.solvePICEF(G, C);

                // Obtain the duals
                Map<Integer, Double> alphas = picef.getAlphas();
                Map<Integer, Double> gammas = picef.getGammas();
                Map<Integer, Map<Integer, Double>> betas = picef.getBetas();

                // Solve the pricing problem for chains and save the positive cost chains
                Map<Integer, List<Integer>> chainstarmap =
                        picef.PricingProblemChains(G, C, alphas, gammas, betas);

                // Solve the pricing problem for cycles and save the positive cost cycles
                Map<Integer, List<Integer>> cyclestarmap =
                        picef.PricingProblemCycles(G, K, alphas);

                // Stop if there are no more positive cost chains or cycles
                if (chainstarmap.isEmpty() && cyclestarmap.isEmpty()) {
                    improving = false;

                    // Save the final time and compute the runtime in seconds
                    long end = System.nanoTime();
                    double runtime = (end - start)/1e9;

                    // Print termination statement
                    System.out.println("No more positive reduced cost chains or cycles, final solution found!");
                    System.out.println("Number of iterations: " + iterations);
                    System.out.println("Total runtime = " + runtime + " seconds");

                } else{
                    // Add the improving cycles as columns to the RMP
                    if (!chainstarmap.isEmpty()){
                        picef.addchainColumns(G, chainstarmap, C);
                    }
                    if (!cyclestarmap.isEmpty()){
                        picef.addcycleColumns(G, cyclestarmap);
                    }
                }
            }


            /*// TEST: Find IP-solution with columns found in the loop
            // Set the model back to simplex if it is better than interior point
            // Set to automatic and let Gurobi choose
            picef.model.set(GRB.IntParam.Method, -1);

            // Set all variables to binary
            for (Edge e: picef.edgevars.keySet()){
                for (GRBVar var : picef.edgevars.get(e)){
                    var.set(GRB.CharAttr.VType, GRB.BINARY);
                }
            }
            for (String cycle: picef.cyclevars.keySet()){
                picef.cyclevars.get(cycle).set(GRB.CharAttr.VType, GRB.BINARY);
            }
            for (int p: picef.loopvars.keySet()){
                picef.loopvars.get(p).set(GRB.CharAttr.VType, GRB.BINARY);
            }

            // Solve one more time to find the integer solution
            picef.model.update();
            picef.model.optimize();
            // Check the status
            int status = picef.model.get(GRB.IntAttr.Status);
            if (status == GRB.OPTIMAL || status == GRB.TIME_LIMIT || status == GRB.SUBOPTIMAL){
                double obj = picef.model.get(GRB.DoubleAttr.ObjVal);
                System.out.println("Solution found!");
                System.out.println("Objective = " + obj);
            } else if (status == GRB.INFEASIBLE){
                System.out.println("RMP is infeasible!");
            } else {
                System.out.println("Model not solved to optimality. Status = " + status);
            }


            System.out.println("Solution values: ");
            // Check the used variables in the solution
            int edgecount = 0;

            // Check the loops first
            for (int p : picef.loopvars.keySet()){
                GRBVar var = picef.loopvars.get(p);
                String name = var.get(GRB.StringAttr.VarName);
                double val = var.get(GRB.DoubleAttr.X);
                if (val > 1E-6){
                    System.out.println("Used loop variable " + name + " with value " + val);
                    edgecount++;
                }
            }

            // Check the edges
            if (picef.edgevars.isEmpty()){
                System.out.println("No edges are used.");
            } else {
                for ( Edge e : picef.edgevars.keySet()){
                    List<GRBVar> nodelist = picef.edgevars.get(e);

                    for (int c=0; c < nodelist.size(); c++){
                        GRBVar var = nodelist.get(c);
                        String name = var.get(GRB.StringAttr.VarName);
                        double val = var.get(GRB.DoubleAttr.X);

                        if (val > 1E-6){
                            System.out.println("Used edge variable " + name + " with value " + val);
                            edgecount++;
                        }
                    }
                }
            }

            // Check the cycles
            if (picef.cyclevars.isEmpty()){
                System.out.println("No cycles are used.");
            } else {
                for (String cycle : picef.cyclevars.keySet()){
                    GRBVar var = picef.cyclevars.get(cycle);
                    double val = var.get(GRB.DoubleAttr.X);
                    if (val > 1E-6){
                        System.out.println("Used cycle variable " + cycle + " with value " + val);
                        edgecount+= picef.cyclelengths.get(cycle);
                    }
                }
            }

            System.out.println("Number of matches: " + edgecount);

             */


        } catch (GRBException e) {
            throw new RuntimeException(e);
        }
    }
}
