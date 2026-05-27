import java.util.List;

public class Instance {
    // This class is used to define an instance, obtained from a .txt file
    // .txt files are read in the Read.Java class

    // The number of pairs in the pool
    public int NrPairs;

    // The number of non-directed donors in the pool
    public int NrNDDs;

    // The list of indices of the pairs and NDDs
    public List<Integer> indices;

    // The list of patient blood types
    public List<String> patients;

    // The list of donor blood types
    public List<String> donors;

    // The list of vPRA levels of the patients
    public List<Double> vPRAs;

    // The list on whether a pair is donor-only (NDD)
    public List<Integer> NDDs;

    // A list of the waiting times (years on dialysis) for each patient
    public List<Double> waitingtimes;

    public Instance(int nrPairs, int nrNDDs, List<Integer> indices,
                    List<String> patients, List<String> donors, List<Double> vPRAs, List<Integer> NDDs,
                    List<Double> waitingtimes) {
    this.NrPairs = nrPairs;
    this.NrNDDs = nrNDDs;
    this.indices = indices;
    this.patients = patients;
    this.donors = donors;
    this.vPRAs = vPRAs;
    this.NDDs = NDDs;
    this.waitingtimes = waitingtimes;
    }
}
