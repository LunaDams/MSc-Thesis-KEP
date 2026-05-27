// Packages used.
import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

public class WaitingTimes {

    // Define infinity as a large number
    private static final double INF = 1e6;

    // A function that gives each patient a score based on blood-type compatibility and vRPA levels
    // Used to determine expected waiting time
    public static List<Integer> score(List<Integer> indices, List<String> patients, List<String> donors, List<Double> vPRAs){

        // Count blood types in the donor list
        int donorA = 0;
        int donorB = 0;
        int donorAB = 0;
        int donorO = 0;

        for (int i : indices){
            String donor = donors.get(i-1);
            if (donor.equals("A")){
                donorA++;
            }
            if (donor.equals("B")){
                donorB++;
            }
            if (donor.equals("AB")){
                donorAB++;
            }
            if (donor.equals("O")) {
                donorO++;
            }
        }

        // Count the number of blood-type compatible donors for each blood type
        int compA = donorA + donorO;
        int compB = donorB + donorO;
        int compAB = donors.toArray().length;
        int compO = donorO;

        // Print the number of donors for each blood type and the number of compatible donors
        System.out.println("Blood type A donors available: " + donorA);
        System.out.println("Blood type B donors available: " + donorB);
        System.out.println("Blood type AB donors available: " + donorAB);
        System.out.println("Blood type O donors available: " + donorO);

        System.out.println("Donors compatible with A-patients: " + compA);
        System.out.println("Donors compatible with B-patients: " + compB);
        System.out.println("Donors compatible with AB-patients: " + compAB);
        System.out.println("Donors compatible with O-patients: " + compO);


        // Rank the blood-types by most to least number of donors
        // Since AB compatible with all: always first
        // Since O only compatible with O: always last

        List<String> RankedBloodTypes = new ArrayList<>();
        RankedBloodTypes.add("AB");

        if (compA >= compB){
            RankedBloodTypes.add("A");
            RankedBloodTypes.add("B");
        } else {
            RankedBloodTypes.add("B");
            RankedBloodTypes.add("A");
        }

        RankedBloodTypes.add("O");

        // For each patient, give a waiting time score form 1 to 7 based on their "match difficulty"
        // 1 being the shortest to match and 7 the hardest

        // Initialise the score list
        List<Integer> scores = new ArrayList<>();

        // Assign scores
        for (int i : indices){
            String patient = patients.get(i-1);
            int score = 1;

            // Add to score based on blood-type ranking
            if (patient.equals(RankedBloodTypes.get(1))){
                score++;
            }
            if (patient.equals(RankedBloodTypes.get(2))){
                score = score+2;
            }
            if (patient.equals(RankedBloodTypes.get(3))){
                score = score+3;
            }

            // Add to score based on vPRA levels
            // vRPA 0-9% : easy
            // vRPA 10-79% : medium
            // vRPA 80-100% : hard

            double vPRA = vPRAs.get(i-1);
            if (vPRA > 0.09){
                score++;
            }
            if (vPRA > 0.79){
                score++;
            }
            scores.add(score);
        }

    return scores;
    }

    // A function that makes a probability matrix U based on the baseline probabilities and scores
    // P_{ik} = probability that patient with score i has waiting time in range k
    public static List<List<Double>> Pmatrix(List<Double> baseprobs, List<Double> weights){

        // Initialise the reweighting matrix U and probability matrix P
        List<List<Double>> U = new ArrayList<>();
        List<List<Double>> P = new ArrayList<>();

        // For each score and time range: reweight the baseline probabilities based on score
        //U_{ik} = log(p_k) + w_k*s_i
        for (int i=0; i < 7; i++){
            List<Double> ulist = new ArrayList<>();
            for (int k =0; k < baseprobs.size(); k++){
                // Center weights such that baseline probabilities are used at score 4
               double u = Math.log(baseprobs.get(k)) + weights.get(k)*((i+1)-4);
               ulist.add(u);
            }
            U.add(ulist);
        }

        // For each score: determine the probability of being in the time range
        //P_{ik} = exp(U_ik)/(Sum_j(exp(U_{ij})))
        for (int i=0; i < 7; i++){
            List<Double> plist = new ArrayList<>();
            List<Double> ulist = U.get(i);
            double ulistexpsum = 0.0;
            for (int j=0; j<baseprobs.size(); j++){
                ulistexpsum =  ulistexpsum + Math.exp(ulist.get(j));
            }

            for (int k=0; k<baseprobs.size(); k++){
                double p = Math.exp(ulist.get(k))/ulistexpsum;
                plist.add(p);
            }
            P.add(plist);
        }
        return P;
    }


    // A function that simulates the waiting time for each patient in the pool
    // Modelled based on waiting time probabilities for each score
    // Times are given in years
    public static List<Double> waitingtimes(List<Integer> indices, List<String> patients,
                                            List<String> donors, List<Double> vPRAs,
                                            List<Double> baseprobs, List<Double> weights){

        // Initialise the waiting time list
        List<Double> times = new ArrayList<>();

        // Get scores and probability matrix
        List<Integer> scores = score(indices,patients,donors,vPRAs);
        System.out.println("Patient scores: " + scores);
        List<List<Double>> P = Pmatrix(baseprobs, weights);

        for (int i : indices){

            int score = scores.get(i-1);

            // Retrieve the waiting time distribution for this score
            List<Double> plist = P.get(score-1);

            // Random draw for patient i
            SplittableRandom rand = new SplittableRandom(i);
            double r = rand.nextDouble();

            // Initialise waiting time for patient i
            double time = 0;

            // Initialise the cumulative probability
            double cumulative = 0;

            // Determine the sampled waiting time category for patient i
            for (int k =0; k < plist.size(); k++){
                cumulative = cumulative+plist.get(k);

                if (r <= cumulative){

                    // Sample the waiting time based on the category
                    switch (k){
                        case 0:
                            // Preemetive, waiting time = 0
                            time = 0.0;
                            break;

                        case 1:
                            // waiting time ~ Un[0,1]
                            time = rand.nextDouble(0.0, 1.0);
                            break;

                        case 2:
                            // waiting time ~ Un[1,3]
                            time = rand.nextDouble(1.0, 3.0);
                            break;

                        case 3:
                            // waiting time ~ Un[3,5]
                            time = rand.nextDouble(3.0, 5.0);
                            break;
                    }

                    break;
                }
            }
            times.add(time);
        }
        return times;
    }
}



