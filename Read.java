// Packages used.
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.SplittableRandom;


public class Read {

    // This class is used for reading instances
    // Instances obtained from John P. Dickerson, Kidney Data (PrefLib)

    // Define the baseline probabilities of each waiting time category:
    // Preemptive: 36%, 0-1 years: 16%, 1-3 years: 32%, 3+ years: 16%.
    // This data is received from van de Laar et.al. (2025)
    private static final List<Double> baseprobs = new ArrayList<>(List.of(0.36, 0.16, 0.32, 0.16));

    // Define the reweighting factors to sample waiting times
    // Preemptive: -1.0, 0-1 years: -0.5, 1-3 years: 0.5, 3+ years: 1.0.
    private static final List<Double> weights = new ArrayList<>(List.of(-1.0, -0.5, 0.5, 1.0));




    // Read the instance with patients, donors, blood types and vPRA levels from the .txt file
    public static Instance readInstance(File filename) throws FileNotFoundException{

        //Try to open the file and initialise the use of a buffered reader
        try (BufferedReader s = new BufferedReader(new FileReader(filename))){

            // Ignore the first line of text (headers)
            String headers = s.readLine();

            // For each line, obtain the number of pairs, number of NDDs, the pair index,
            // patient and donor blood types, vPRA level, and NDD if applicable
            // Note that we ignore the columns "Wife-P?" and "Out-Deg"
            int NrPairs = 0;
            int NrNDDs = 0;
            List<Integer> indices = new ArrayList<>();
            List<String> patients = new ArrayList<>();
            List<String> donors = new ArrayList<>();
            List<Double> vPRAs = new ArrayList<>();
            List<Integer> NDDs = new ArrayList<>();

            String line;

            while ((line = s.readLine()) != null){

                // Skip the empty lines
                if (line.isEmpty()){
                    continue;
                }

                String[] tokens = line.split(",");


                // Collect the index and add to the list
                int index = Integer.parseInt(tokens[0]);
                // Sanity check
                if (index < 1) {
                    throw new RuntimeException("Invalid pair index (must be >= 1): " + index);
                }
                indices.add(index);

                // Collect the donor blood type and add to the list
                String donor = tokens[2];
                donors.add(donor);


                // Collect the NDD identity (1 if NDD, 0 otherwise) and add to the list and number of NDDs
                int NDDbinary= Integer.parseInt(tokens[6]);
                NDDs.add(NDDbinary);
                if (NDDbinary == 1){
                    NrNDDs++;
                } else{
                    // Count the number of pairs
                    NrPairs++;
                }

                // Collect the patient blood type and vPRA level and add to the list
                String patient;
                Double vPRA;

                if (NDDbinary == 1){
                    patient = "-";
                    vPRA = 0.0;
                } else {
                    patient = tokens[1];
                    vPRA = Double.parseDouble(tokens[4]);
                }

                patients.add(patient);
                vPRAs.add(vPRA);
            }

            // Simulate the waiting times for the given instance
            List<Double> waitingtimes = WaitingTimes.waitingtimes(indices, patients, donors,
                    vPRAs, baseprobs, weights);

            // Create the instance from the read file
            return new Instance(NrPairs, NrNDDs, indices, patients, donors, vPRAs, NDDs, waitingtimes);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // Build the corresponding compatability graph from the additional .txt file
    public static Graph readGraph(File instancefile, File graphfile,
                                  double thresholdHI, double thresholdLW, double M_HI, double M_LW)
            throws FileNotFoundException{

        // Read the instance file to extract the waiting times and vRPA levels of patients
        // These are used to adjust edge weights in the graph
        Instance instance = readInstance(instancefile);
        List<Double> waitingtimes = instance.waitingtimes;
        List<Double> vPRAs = instance.vPRAs;

        // Extract patients, donors and indices for loops in the graph
        List<String> patients = instance.patients;
        List<String> donors = instance.donors;
        List<Integer> indices = instance.indices;


        //Try to open the graph file and initialise the graph using a buffered reader
        try (BufferedReader s = new BufferedReader(new FileReader(graphfile))){

            // Initialise the pair nodes, NDD nodes, and edges lists
            List<Integer> pairnodes = new ArrayList<>();
            List<Integer> NDDnodes = new ArrayList<>();
            List<Edge> edges = new ArrayList<>();

            String line;

            while ((line = s.readLine()) != null){
                line = line.trim();

                // Skip the empty lines
                if (line.isEmpty()){
                    continue;
                }

                // Header lines: only count the pairs and NDDs
                if (line.startsWith("#")){

                    // Line contains "Pair": add its index to the pair list
                    if (line.contains("Pair")){
                        int index = extractindex(line);
                        pairnodes.add(index);
                    } else if (line.contains("Alturist")){
                        int index = extractindex((line));
                        NDDnodes.add(index);
                    }
                }

                // Edge lines: they do not contain a "#"
                else{

                    // Split the line on the comma to get the start node, end node, and weight (default = 1.0)
                    String[] parts = line.split(",");

                    // Obtain  start node, end node, and weight
                    int from = Integer.parseInt(parts[0]);
                    int to = Integer.parseInt(parts[1]);
                    double weight = Double.parseDouble(parts[2]);

                    // Sanity check
                    if (from < 1 || from > vPRAs.size()) {
                        throw new RuntimeException("Invalid 'from' index: " + from);
                    }

                    if (to < 1 || to > vPRAs.size()) {
                        throw new RuntimeException("Invalid 'to' index: " + to);
                    }

                    // If patient (the "to" node) is HI, multiply the weight by M_HI
                    double vPRA = vPRAs.get(to-1);
                    if (vPRA > thresholdHI - 1E-6){
                        weight = weight*M_HI;
                    }

                    // If patient is LW, multiply the weight by M_LW
                    double waitingtime = waitingtimes.get(to-1);
                    if (waitingtime > thresholdLW - 1E-6){
                        weight = weight*M_LW;
                    }

                    // Add the edge to the edge list, but only if weight is non-zero
                    // If edge.weight = 0, then the edge goes to an NDD so it does nothing
                    Edge edge = new Edge(from, to, weight);
                    if (edge.weight > 1E-6){
                        edges.add(edge);
                    }
                }
            }

            // Add additional loop edges with probability 1-vPRA for blood-type compatible pairs
            for (int i : indices){

                String patient = patients.get(i-1);
                String donor = donors.get(i-1);
                double vPRA = vPRAs.get(i-1);



                if (isABOcompatible(patient,donor)){

                    double p = 1-vPRA;

                    // Random draw for patient-donor pair i
                    SplittableRandom rand = new SplittableRandom(i);
                    double r = rand.nextDouble();

                    // Add the loop with probability p
                    if (r <= p){

                        double weight = 1.0;

                        // If patient is HI, multiply the weight by M_HI
                        if (vPRA > thresholdHI - 1E-6){
                            weight = weight*M_HI;
                        }

                        // If patient is LW, multiply the weight by M_LW
                        double waitingtime = waitingtimes.get(i-1);
                        if (waitingtime > thresholdLW - 1E-6){
                            weight = weight*M_LW;
                        }

                        Edge loop = new Edge(i,i, weight);
                        edges.add(loop);

                    }
                }
            }

            // Sanity check for edge indexing
            for (Edge e : edges) {
                if (e.from < 1 || e.to < 1) {
                    throw new RuntimeException("Edge has invalid node: " + e.from + " -> " + e.to);
                }
            }

            // Create graph
            Graph g = new Graph(pairnodes, NDDnodes, edges);
            return g;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // A function to extract an index from the line of a .txt file
    // Example: "# ALTERNATIVE NAME 1: Pair 1", then the index is 1
    private static int extractindex(String line){

        String[] parts = line.split(":")[0].split(" ");

        // Get the last token before the ":", this is the index
        return Integer.parseInt(parts[parts.length-1]);
    }


    // A boolean function to determine whether a donor-patient pair is blood-type compatible
    // Automatically returns false for NDDs
    private static boolean isABOcompatible(String patient, String donor){

        // Patient blood type A: compatible with donor type A and O
        if (Objects.equals(patient, "A")){
            return Objects.equals(donor, "A") || Objects.equals(donor, "O");
        }
        // Patient blood type B: compatible with donor type B and O
        else if (Objects.equals(patient, "B")) {
            return Objects.equals(donor, "B") || Objects.equals(donor, "O");
        }
        // Patient blood type AB: always compatible
        else if (Objects.equals(patient, "AB")) {
            return true;
        }
        // Patient blood type O: compatible with donor type O
        else if (Objects.equals(patient, "O")) {
            return Objects.equals(donor, "O");
        }
        else{
            return false;
        }
    }
}
