using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//*************************************************************
//  Graph implementation
class Edge<TLoad>
{
    public Edge(Node<TLoad> _source, Node<TLoad> _target, double _weight) : this(_source, _target, _weight, "")
    {

    }
    public Edge(Node<TLoad> _source, Node<TLoad> _target, double _weight, string _label)
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
    public Node<TLoad> source;
    public Node<TLoad> target;
    public double weight;
    public string label;
}
class Node<TLoad>
{
    public Node()
        : this(null)
    {

    }
    public void ShallowCopy(Node<TLoad> o)
    {
        incoming = o.incoming;
        outgoing = o.outgoing;
    }
    public Node(string _label)
        : this(_label, default(TLoad), 0.0)
    {

    }
    public Node(string _label, TLoad _load, double _weight)
    {
        incoming = new List<Edge<TLoad>>();
        outgoing = new List<Edge<TLoad>>();
        label = _label;
        load = _load;
        weight = _weight;
    }
    public List<Edge< TLoad>> incoming;
    public List<Edge<TLoad>> outgoing;
    public string label;
    public TLoad load;
    public double weight;
}
class Graph<TLoad>
{
    public Graph()
    {
        nodes = new List<Node<TLoad>>();
    }
    public Graph(Graph<TLoad> other)
        : this()
    {
        foreach (Edge<TLoad> edge in other.Edges())
            AddEdge(edge.source.label, edge.target.label, edge.weight);
    }
    public List<Edge<TLoad>> Edges()
    {
        List<Edge< TLoad>> edges = new List<Edge< TLoad>>();

        foreach (Node< TLoad> node in nodes)
        {
            foreach (Edge<TLoad> edge in node.outgoing)
            {
                if (!edges.Contains(edge)) edges.Add(edge);
            }
        }

        return edges;
    }
    public void AddEdge(string source, string target, double weight)
    {
        new Edge< TLoad>(Node(source), Node(target), weight);
    }
    public void AddEdge(string source, string target, double weight, string label)
    {
        new Edge<TLoad>(Node(source), Node(target), weight, label);
    }
    public Node<TLoad> Node(string label)
    {
        foreach (Node<TLoad> node in nodes)
        {
            if (node.label.CompareTo(label) == 0) return node;
        }
        nodes.Add(new Node<TLoad>(label));
        return nodes.Last();
    }
    public List<Node<TLoad>> TopologicalOrdering()
    {
        List<Node<TLoad>> ordering = new List<Node<TLoad>>();

        Graph<TLoad> reduced = new Graph<TLoad>(this);

        List<Node<TLoad>> candidates = new List<Node<TLoad>>();

        foreach (Node<TLoad> node in reduced.nodes)
            if (node.incoming.Count == 0)
                candidates.Add(node);

        while (candidates.Count > 0)
        {
            Node<TLoad> a = candidates[0];
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
    public List<Edge<TLoad>> Path(string source, string sink)
    {
        List<Edge<TLoad>> path = new List<Edge<TLoad>>();

        if (!RecursivePath(Node(source), Node(sink), ref path)) ;
        //         throw new Exception("There is no path!");

        return path;
    }
    bool RecursivePath(Node<TLoad> actual, Node<TLoad> sink, ref List<Edge<TLoad>> path)
    {
        if (actual == sink) return true;
        foreach (Edge<TLoad> edge in actual.outgoing)
        {
            bool cycle = false;
            foreach (Edge<TLoad> visited in path)
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

    public List<List<Edge<TLoad>>> AllPaths(string source, string sink)
    {
        List<List<Edge<TLoad>>> allPaths = new List<List<Edge<TLoad>>>();
        List<Edge<TLoad>> path = new List<Edge<TLoad>>();
        RecursiveAllPath(Node(source), Node(sink), path, allPaths);
        return allPaths;
    }
    void RecursiveAllPath(Node<TLoad> actual, Node<TLoad> sink, List<Edge<TLoad>> path, List<List<Edge<TLoad>>> allPaths)
    {
        if (actual == sink)
        {
            allPaths.Add(new List<Edge<TLoad>>(path));
        }

        foreach (Edge<TLoad> edge in actual.outgoing)
        {
            path.Add(edge);
            RecursiveAllPath(edge.target, sink, path, allPaths);
            path.Remove(edge);
        }
    }
    public List<List<Edge<TLoad>>> Cycles()
    {
        List<List<Edge<TLoad>>> cycles = new List<List<Edge<TLoad>>>();

        List<Edge<TLoad>> edges = Edges();
        while (edges.Count > 0)
        {
            List<Edge<TLoad>> cycle = new List<Edge<TLoad>>();

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
    public List<Edge<TLoad>> LongestPath(string source, string sink)
    {
        List<Node<TLoad>> ordering = TopologicalOrdering();

        while (ordering[0].label.CompareTo(source) != 0) ordering.RemoveAt(0);
        while (ordering.Last().label.CompareTo(sink) != 0) ordering.RemoveAt(ordering.Count - 1);

        Dictionary<Node<TLoad>, Edge<TLoad>> backtrack = new Dictionary<Node<TLoad>, Edge<TLoad>>();
        Dictionary<Node<TLoad>, double> s = new Dictionary<Node<TLoad>, double>();
        foreach (Node<TLoad> node in nodes) s[node] = -1000000;
        s[ordering[0]] = 0;
        foreach (Node<TLoad> node in ordering)
        {
            if (node.incoming.Count == 0) continue;
            Edge<TLoad> maxEdge = null;
            double max = Int32.MinValue;
            foreach (Edge<TLoad> edge in node.incoming)
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

        List<Edge<TLoad>> path = new List<Edge<TLoad>>();

        Edge<TLoad> back = backtrack[ordering.Last()];
        while (back.source.label.CompareTo(source) != 0)
        {
            path.Add(back);
            back = backtrack[back.source];
        }
        path.Add(back);

        path.Reverse();

        return path;
    }
    static public double Length(List<Edge<TLoad>> path)
    {
        double l = 0;
        foreach (Edge<TLoad> edge in path) l += edge.weight;
        return l;
    }
    public double Weight()
    {
        double weight = 0;
        foreach (Edge<TLoad> edge in Edges()) weight += edge.weight;
        return weight;
    }
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        List<Edge<TLoad>> edges = Edges();
        edges.Sort((a, b) => (a.source.label.CompareTo(b.source.label)));

        //  write result
        foreach (Edge<TLoad> edge in edges)
            sb.AppendFormat("{0}->{1}:{2}\n", edge.source.label, edge.target.label, edge.label);

        return sb.ToString();
    }
    /*
        AddToTrie
        currentNode ← root
            for i ← 1 to |Pattern|
                if there is an outgoing edge from currentNode with label currentSymbol
                    currentNode ← ending node of this edge
                else
                    add a new node newNode to Trie
                    add a new edge from currentNode to newNode with label currentSymbol
                    currentNode ← newNode
    */
    public void AddToTrie(string Pattern)
    {
        Node<TLoad> currentNode = Node("0");
        for(int i=0;i<Pattern.Length; ++i)
        {
            bool found = false;
            foreach(Edge<TLoad> edge in currentNode.outgoing)
            {
                if (edge.label[0] == Pattern[i])
                {
                    currentNode = edge.target;
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                //  AddEdge(currentNode.label, nodes.Count.ToString(), 0, Pattern[i].ToString());
                nodes.Add(new Node<TLoad>(nodes.Count.ToString(),default(TLoad),currentNode.weight + 1));
                new Edge<TLoad>(currentNode, nodes.Last(), 0, Pattern[i].ToString());

                currentNode = nodes.Last();
            }
        }
    }
    /*
    PREFIXTRIEMATCHING(Text, Trie)
        symbol ← first letter of Text
        v ← root of Trie
        while forever
            if v is a leaf in Trie
                return the pattern spelled by the path from the root to v
            else if there is an edge (v, w) in Trie labeled by symbol
                symbol ← next letter of Text
                v ← w
            else
                output "no matches found"
                return
    */
    public bool PrefixTrieMatching(string Text)
    {
        Node<TLoad> currentNode = Node("0");
        int i = 0;
        for(;i<Text.Length;++i)
        {
            if (currentNode.outgoing.Count == 0) return true;
            bool found = false;
            foreach (Edge<TLoad> edge in currentNode.outgoing)
            {
                if (edge.label[0] == Text[i])
                {
                    currentNode = edge.target;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                break;
            }
            if (currentNode.outgoing.Count == 0) return true;
        }
        return false;
    }
    public void Collapse()
    {
        bool collapsed = false;
        while(!collapsed)
        {
            collapsed = true;
            for(int i = 0; i < nodes.Count;++i)
            {
                if(nodes[i].outgoing.Count == 1 && nodes[i].incoming.Count == 1)
                {
                    nodes[i].incoming[0].label += nodes[i].outgoing[0].label;
                    nodes[i].incoming[0].target = nodes[i].outgoing[0].target;
                    nodes[i].outgoing[0].target.incoming[0] = nodes[i].incoming[0];
                    nodes.RemoveAt(i);
                    i--;
                    //int incoming_index = nodes[i].incoming[0].source.outgoing.IndexOf(nodes[i].incoming[0]);
                    collapsed = false;
              //      break;
                }
            }
        }
    }
    public string LongestRepeat()
    {
        string longest = "";

        foreach (Node<TLoad> node in nodes)
        {
            if (node.outgoing.Count == 0)
                {
                    if(node.incoming[0].source.weight > longest.Length)
                    {
                        longest = "";
                        Node<TLoad> current = node.incoming[0].source;
                        while(current.incoming.Count != 0)
                        {
                            longest = current.incoming[0].label + longest;
                            current = current.incoming[0].source;
                        }
                    }
                }
        }

        return longest;
    }
    public bool TrieMatching(string word)
    {
        Node<TLoad> currentNode = Node("0");
        for (int i = 0; i < word.Length; ++i)
        {
            if (currentNode.outgoing.Count == 0) return false;
            bool found = false;
            foreach (Edge<TLoad> edge in currentNode.outgoing)
            {
                if (edge.label[0] == word[i])
                {
                    currentNode = edge.target;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                return false;
            }
        }
        return true;
    }
    public List<string> Substring(int length)
    {
        List<string> subStrings = new List<string>();
        RecursiveSubstring(nodes[0], ref subStrings, "", length);
        return subStrings;
    }
    public void RecursiveSubstring( Node<TLoad> node, ref List<string> subStrings, string way, int remainingDepth )
    {
        if (remainingDepth == 0)
        {
            subStrings.Add(way);
            return;
        }
        foreach(Edge<TLoad> edge in node.outgoing)
        {
            RecursiveSubstring(edge.target, ref subStrings, way + edge.label, remainingDepth - 1);
        }
    }
    public int Leaves()
    {
        int leaves = 0;
        foreach (Node<TLoad> node in nodes)
        {
            if (node.outgoing.Count == 0)
            {
                leaves++;
            }
        }
        return leaves;
    }
    public List<Node<TLoad>> nodes;
}
//  End of Graph
//****************************************************************