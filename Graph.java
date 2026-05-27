// Packages used.
import java.util.List;

public class Graph {
    // This class is used to build a graph, obtained from a .txt file
    // And adjusted with weights and possibly extra edges for self-compatible pairs
    // .txt files are read in the Read.Java class

    // The list of donor-patient pair nodes in the graph
    public List<Integer> pairnodes;

    // The list of NDD nodes in the graph
    public List<Integer> NDDnodes;

    // The list of edges in the graph
    public List<Edge> edges;


    public Graph(List<Integer> pairnodes, List<Integer> NDDnodes, List<Edge> edges) {
        this.pairnodes = pairnodes;
        this.NDDnodes = NDDnodes;
        this.edges = edges;
    }
}
