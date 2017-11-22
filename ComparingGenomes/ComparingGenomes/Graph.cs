using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//*************************************************************
//  Graph implementation
class Edge<TLabel, TLoad> where TLabel : IComparable
{
    public Edge(Node<TLabel, TLoad> _source, Node<TLabel, TLoad> _target, double _weight) : this(_source, _target, _weight, "")
    {

    }
    public Edge(Node<TLabel, TLoad> _source, Node<TLabel, TLoad> _target, double _weight, string _label)
    {
        source = _source;
        target = _target;
        weight = _weight;
        label = _label;

        source.outgoing.Add(this);
        target.incoming.Add(this);
    }
    public void Remove()
    {
        source.outgoing.Remove(this);
        target.incoming.Remove(this);
    }
    public Node<TLabel, TLoad> source;
    public Node<TLabel, TLoad> target;
    public double weight;
    public string label;
}
class Node<TLabel, TLoad> where TLabel : IComparable
{
    public Node()
        : this(default(TLabel))
    {

    }
    public void ShallowCopy(Node<TLabel, TLoad> o)
    {
        incoming = o.incoming;
        outgoing = o.outgoing;
    }
    public Node(TLabel _label)
        : this(_label, default(TLoad), 0.0)
    {

    }
    public Node(TLabel _label, TLoad _load, double _weight)
    {
        incoming = new List<Edge<TLabel, TLoad>>();
        outgoing = new List<Edge<TLabel, TLoad>>();
        label = _label;
        load = _load;
        weight = _weight;
    }
    public List<Edge<TLabel, TLoad>> incoming;
    public List<Edge<TLabel, TLoad>> outgoing;
    public TLabel label;
    public TLoad load;
    public double weight;
}
class Graph<TLabel, TLoad> where TLabel : IComparable
{
    public Graph()
    {
        nodes = new List<Node<TLabel, TLoad>>();
    }
    public Graph(Graph<TLabel, TLoad> other)
        : this()
    {
        foreach (Edge<TLabel, TLoad> edge in other.Edges())
            AddEdge(edge.source.label, edge.target.label, edge.weight);
    }
    public List<Edge<TLabel, TLoad>> Edges()
    {
        List<Edge<TLabel, TLoad>> edges = new List<Edge<TLabel, TLoad>>();

        foreach (Node<TLabel, TLoad> node in nodes)
        {
            foreach (Edge<TLabel, TLoad> edge in node.outgoing)
            {
                if (!edges.Contains(edge)) edges.Add(edge);
            }
        }

        return edges;
    }
    public void AddEdge(TLabel source, TLabel target, double weight)
    {
        new Edge<TLabel, TLoad>(Node(source), Node(target), weight);
    }
    public void AddEdge(TLabel source, TLabel target, double weight, string label)
    {
        new Edge<TLabel, TLoad>(Node(source), Node(target), weight, label);
    }
    public Node<TLabel, TLoad> Node(TLabel label)
    {
        foreach (Node<TLabel, TLoad> node in nodes)
        {
            if (node.label.CompareTo(label) == 0) return node;
        }
        nodes.Add(new Node<TLabel, TLoad>(label));
        return nodes.Last();
    }
    public List<Node<TLabel, TLoad>> TopologicalOrdering()
    {
        List<Node<TLabel, TLoad>> ordering = new List<Node<TLabel, TLoad>>();

        Graph<TLabel, TLoad> reduced = new Graph<TLabel, TLoad>(this);

        List<Node<TLabel, TLoad>> candidates = new List<Node<TLabel, TLoad>>();

        foreach (Node<TLabel, TLoad> node in reduced.nodes)
            if (node.incoming.Count == 0)
                candidates.Add(node);

        while (candidates.Count > 0)
        {
            Node<TLabel, TLoad> a = candidates[0];
            ordering.Add(Node(a.label));
            candidates.RemoveAt(0);
            while (a.outgoing.Count > 0)
            {
                if (a.outgoing[0].target.incoming.Count == 1)
                {
                    candidates.Add(a.outgoing[0].target);
                }
                a.outgoing[0].Remove();
            }
        }
        if (reduced.Edges().Count > 0)
            throw new Exception("The graph is not DAG!");

        return ordering;
    }
    public List<Edge<TLabel, TLoad>> Path(TLabel source, TLabel sink)
    {
        List<Edge<TLabel, TLoad>> path = new List<Edge<TLabel, TLoad>>();

        if (!RecursivePath(Node(source), Node(sink), ref path)) ;
        //         throw new Exception("There is no path!");

        return path;
    }
    bool RecursivePath(Node<TLabel, TLoad> actual, Node<TLabel, TLoad> sink, ref List<Edge<TLabel, TLoad>> path)
    {
        if (actual == sink) return true;
        foreach (Edge<TLabel, TLoad> edge in actual.outgoing)
        {
            bool cycle = false;
            foreach (Edge<TLabel, TLoad> visited in path)
            {
                if (edge.target == visited.source)
                {
                    cycle = true;
                    break;
                }
            }
            if (cycle) continue;
            path.Add(edge);
            if (RecursivePath(edge.target, sink, ref path))
            {
                return true;
            }
            path.Remove(edge);
        }
        return false;
    }

    public List<List<Edge<TLabel, TLoad>>> AllPaths(TLabel source, TLabel sink)
    {
        List<List<Edge<TLabel, TLoad>>> allPaths = new List<List<Edge<TLabel, TLoad>>>();
        List<Edge<TLabel, TLoad>> path = new List<Edge<TLabel, TLoad>>();
        RecursiveAllPath(Node(source), Node(sink), path, allPaths);
        return allPaths;
    }
    void RecursiveAllPath(Node<TLabel, TLoad> actual, Node<TLabel, TLoad> sink, List<Edge<TLabel, TLoad>> path, List<List<Edge<TLabel, TLoad>>> allPaths)
    {
        if (actual == sink)
        {
            allPaths.Add(new List<Edge<TLabel, TLoad>>(path));
        }

        foreach (Edge<TLabel, TLoad> edge in actual.outgoing)
        {
            path.Add(edge);
            RecursiveAllPath(edge.target, sink, path, allPaths);
            path.Remove(edge);
        }
    }
    public List<List<Edge<TLabel, TLoad>>> Cycles()
    {
        List<List<Edge<TLabel, TLoad>>> cycles = new List<List<Edge<TLabel, TLoad>>>();

        List<Edge<TLabel, TLoad>> edges = Edges();
        while (edges.Count > 0)
        {
            List<Edge<TLabel, TLoad>> cycle = new List<Edge<TLabel, TLoad>>();

            cycle.Add(edges[0]);
            edges.Remove(edges[0]);
            while (cycle[0].source != cycle.Last().target)
            {
                cycle.Add(cycle.Last().target.outgoing[0]);
                edges.Remove(cycle.Last());
            }

            cycles.Add(cycle);
        }

        return cycles;
    }
    public List<Edge<TLabel, TLoad>> LongestPath(TLabel source, TLabel sink)
    {
        List<Node<TLabel, TLoad>> ordering = TopologicalOrdering();

        while (ordering[0].label.CompareTo(source) != 0) ordering.RemoveAt(0);
        while (ordering.Last().label.CompareTo(sink) != 0) ordering.RemoveAt(ordering.Count - 1);

        Dictionary<Node<TLabel, TLoad>, Edge<TLabel, TLoad>> backtrack = new Dictionary<Node<TLabel, TLoad>, Edge<TLabel, TLoad>>();
        Dictionary<Node<TLabel, TLoad>, double> s = new Dictionary<Node<TLabel, TLoad>, double>();
        foreach (Node<TLabel, TLoad> node in nodes) s[node] = -1000000;
        s[ordering[0]] = 0;
        foreach (Node<TLabel, TLoad> node in ordering)
        {
            if (node.incoming.Count == 0) continue;
            Edge<TLabel, TLoad> maxEdge = null;
            double max = Int32.MinValue;
            foreach (Edge<TLabel, TLoad> edge in node.incoming)
            {
                if (s[edge.source] + edge.weight > max)
                {
                    max = s[edge.source] + edge.weight;
                    maxEdge = edge;
                }
            }
            s[node] = max;
            backtrack[node] = maxEdge;
        }

        List<Edge<TLabel, TLoad>> path = new List<Edge<TLabel, TLoad>>();

        Edge<TLabel, TLoad> back = backtrack[ordering.Last()];
        while (back.source.label.CompareTo(source) != 0)
        {
            path.Add(back);
            back = backtrack[back.source];
        }
        path.Add(back);

        path.Reverse();

        return path;
    }
    static public double Length(List<Edge<TLabel, TLoad>> path)
    {
        double l = 0;
        foreach (Edge<TLabel, TLoad> edge in path) l += edge.weight;
        return l;
    }
    public double Weight()
    {
        double weight = 0;
        foreach (Edge<TLabel, TLoad> edge in Edges()) weight += edge.weight;
        return weight;
    }
    public List<Node<TLabel, TLoad>> nodes;
}
//  End of Graph
//****************************************************************