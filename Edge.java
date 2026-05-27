public class Edge {
    // A class used to define edges in a graph

    // The starting node of the edge
    public int from;

    // The end node of the edge
    public int to;

    // The weight on the edge
    public double weight;

    public Edge(int from, int to, double weight){
        this.from = from;
        this.to = to;
        this.weight = weight;
    }
}
