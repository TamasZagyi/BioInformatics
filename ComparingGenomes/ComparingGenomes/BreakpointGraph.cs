using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

class Chromosome : List<int>
{
    /*
        ChromosomeToCycle(Chromosome)
             for j ← 1 to |Chromosome|
                  i ← Chromosomej
                  if i > 0
                       Node2j−1 ←2i−1
                       Node2j ← 2i
                  else
                       Node2j−1 ← -2i
                       Node2j ←-2i−1
             return Nodes
        */
    public Nodes ToCycle()
    {
        Nodes nodes = new Nodes();

        for (int j = 0; j < Count; ++j)
        {
            int i = this[j];
            if (i > 0)
            {
                nodes.Add(2 * i - 1);
                nodes.Add(2 * i);
            }
            else
            {
                nodes.Add(-2 * i);
                nodes.Add(-2 * i - 1);
            }
        }

        return nodes;
    }
    public void Parse(string line)
    {
        line = line.Replace("+", "");
        line = line.Replace("(", " ");
        line = line.Replace(")", " ");

        char[] delimiterChars = { ' ', ',', '.', ':', '\t' };

        foreach (string item in line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries))
            Add(Int32.Parse(item));
    }
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();

        sb.Append("(");
        for (int i = 0; i < Count; ++i)
        {
            sb.AppendFormat(this[i] < 0 ? "{0}" : "+{0}", this[i]);
            if (i != Count - 1) sb.Append(" ");
        }
        sb.Append(")");

        return sb.ToString();
    }
    public void Sort()
    {
        //Sort((a, b) => (a.Item1.CompareTo(b.Item1)));
    }
    public Edges ColoredEdges()
    {
        Edges edges = new Edges();

        Nodes nodes = ToCycle();
        for (int j = 0; j < Count; ++j)
            edges.Add(nodes[2 * j + 1], j < Count - 1 ? nodes[2 * j + 2] : nodes[0]);

        return edges;
    }
    public Edges BlackEdges()
    {
        Edges edges = new Edges();
        for (int j = 0; j < Count; ++j)
        {
            if( this[j] > 0 )
                edges.Add(this[j]*2-1, this[j] * 2);
            else
                edges.Add(this[j] * -2, this[j] * -2-1);
        }

        return edges;
    }
}

class Genome : List<Chromosome>
{
    /*
        ColoredEdges(P)
         Edges ← an empty set
         for each chromosome Chromosome in P
              Nodes ← ChromosomeToCycle(Chromosome)
              for j ← 1 to |Chromosome|
                   add the edge (Nodes2j, Nodes2j +1) to Edges
         return Edges
    */
    public Edges ColoredEdges()
    {
        Edges edges = new Edges();

        foreach (Chromosome chromosome in this)
        {
            edges.Add(chromosome.ColoredEdges());
        }

        return edges;
    }
    public Edges BlackEdges()
    {
        Edges edges = new Edges();

        foreach (Chromosome chromosome in this)
        {
            edges.Add(chromosome.BlackEdges());
        }

        return edges;
    }
    public void Parse(string line)
    {
        char[] chromosomeDelimiter = { '(' };
        string[] chromosomes = line.Split(chromosomeDelimiter, StringSplitOptions.RemoveEmptyEntries);

        foreach (string item in chromosomes)
        {
            Chromosome chromosome = new Chromosome();
            chromosome.Parse(item);
            Add(chromosome);
        }
    }
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        foreach(Chromosome chromosome in this)
        {
            sb.Append(chromosome.ToString());
        }
        return sb.ToString();
    }
    public int BlockNumber()
    {
        int blocks = 0;
        foreach(Chromosome chromosome in this)
        {
            blocks += chromosome.Count;
        }
        return blocks;
    }
    /*
           2-BreakOnGenome(P, i, i′, j, j′)
             GenomeGraph ← BlackEdges(P) and ColoredEdges(P)
             GenomeGraph ← 2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′)
             P ← GraphToGenome(GenomeGraph)
             return P
        */
    public void TwoBreak(int i, int im, int j, int jm)
    {
        Edges genomeGraph = ColoredEdges();
  //      genomeGraph.Add(BlackEdges());
        genomeGraph.TwoBreak(i, im, j, jm);       
        Clear();
        AddRange(genomeGraph.ToGenome());
    }
}

class Nodes : List<int>
{
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        sb.Append("(");
        for (int i = 0; i < Count; ++i)
        {
            sb.AppendFormat("{0}", this[i]);
            if (i != Count - 1) sb.Append(" ");
        }
        sb.Append(")");

        return sb.ToString();
    }
    public void Parse(string line)
    {
        line = line.Replace("(", "");
        line = line.Replace(")", "");

        char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
        foreach (string item in line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries))
            Add(Int32.Parse(item));
    }
    /*
        CycleToChromosome(Nodes)
         for j ← 1 to |Nodes|/2
              if Node2j−1 < Node2j
                   Chromosomej ← Node2j /2
              else
                   Chromosomej ← −Node2j−1/2
         return Chromosome
        */
    public Chromosome ToChromosome()
    {
        Chromosome chromosome = new Chromosome();

        for (int j = 0; j < Count / 2; ++j)
        {
            if (this[2 * j] < this[2 * j + 1])
            {
                chromosome.Add(this[2 * j + 1] / 2);
            }
            else
            {
                chromosome.Add(-this[2 * j] / 2);
            }
        }

        return chromosome;
    }
    public void TakeOver(Edges edges, int i)
    {
        Add(edges[i].Item1);
        Add(edges[i].Item2);
        edges.RemoveAt(i);
    }
}

class Edges
{
    public Edges()
    {
        data = new List<Tuple<int, int>>();
    }
    public Edges(Edges other)
    {
        data = new List<Tuple<int, int>>(other.data);
    }
    public void Add(int s, int t)
    {
        data.Add(new Tuple<int, int>(s, t));
    }
    public void Add(Edges other)
    {
        data.AddRange(other.data);
    }
    public void Sort()
    {
        data.Sort((a, b) => (a.Item1.CompareTo(b.Item1)));
    }
    public List<Nodes> CycleNodes()
    {
        List<Nodes> cycles = new List<Nodes>();
        Edges reduced = new Edges(this);
        while (reduced.data.Count > 0)
        {
            Nodes cycle = new Nodes();
            cycle.TakeOver(reduced, 0);
            while (true)
            {
                int next = (cycle.Last() % 2) == 0 ? cycle.Last() - 1 : cycle.Last() + 1;
                int index = reduced.IndexOf(next);
                if (index == -1) break;
                cycle.TakeOver(reduced, index);
            }
            cycle.Insert(0, cycle.Last());
            cycle.RemoveAt(cycle.Count - 1);
            cycles.Add(cycle);
        }

        return cycles;
    }
    public List<Edges> Cycles()
    {
        List<Edges> cycles = new List<Edges>();
        Edges reduced = new Edges(this);

        while(reduced.data.Count>0)
        {
            Edges edges = new Edges();
            edges.Add(reduced[0].Item1, reduced[0].Item2);
            reduced.RemoveAt(0);
            while(true)
            {
                int index = Math.Max(reduced.IndexOf(edges.Last().Item1), reduced.IndexOf(edges.Last().Item2));
                if (index == -1) break;
                edges.Add(reduced[index].Item1, reduced[index].Item2);
                reduced.RemoveAt(index);
            }
            cycles.Add(edges);
        }

        return cycles;
    }
    public int IndexOf(int node)
    {
        for (int i = 0; i < data.Count; ++i)
        {
            if (data[i].Item1 == node || data[i].Item2 == node) return i;
        }
        return -1;
    }
    /*
        GraphToGenome(GenomeGraph)
         P ← an empty set of chromosomes
         for each cycle Nodes in GenomeGraph
              Chromosome ← CycleToChromosome(Nodes)
              add Chromosome to P
         return P
        */
    public Genome ToGenome()
    {
        Genome P = new Genome();
        foreach(Nodes cycle in CycleNodes())
        {
            Chromosome chromosome = cycle.ToChromosome();
            P.Add(chromosome);
        }
        return P;
    }
    public Genome ToGenomeFromCompleteGraph()
    {
        Genome P = new Genome();
        foreach (Edges cycle in Cycles())
        {
            Chromosome chromosome = cycle.ToChromosome();
            P.Add(chromosome);
        }
        return P;
    }
    public Chromosome ToChromosome()
    {
        Chromosome chromosome = new Chromosome();

        for (int j = 0; j < data.Count; ++j)
        {
            //  if not black edge
            if (Math.Abs(data[j].Item1 - data[j].Item2) != 1) continue;
            if (data[j].Item1 % 2 == data[j].Item2 % 2) continue;

            if (data[j].Item1 < data[j].Item2)
            {
                chromosome.Add(data[j].Item2 / 2);
            }
            else
            {
                chromosome.Add(-data[j].Item1 / 2);
            }
        }

        return chromosome;
    }
    /*
        2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′)
             remove colored edges (i, i') and (j, j′) from GenomeGraph
             add colored edges (i, j) and (i′, j') to GenomeGraph
             return GenomeGraph 
        */
    public void TwoBreak(int i, int im, int j, int jm)
    {
        int edge1 = IndexOf(i);
        int edge2 = IndexOf(j);

        data.RemoveAt(IndexOf(i));
        data.RemoveAt(IndexOf(j));
        data.Add( new Tuple<int, int>(i, j));
        data.Add( new Tuple<int, int>(im, jm));

        return;

        if (data[edge1].Item1 == i)
        {
            data.RemoveAt(edge1);
            data.Insert(edge1, new Tuple<int, int>(i, j));
        }
        else
        {
            data.RemoveAt(edge1);
            data.Insert(edge1, new Tuple<int, int>(j, i));
        }

        if (data[edge2].Item1 == j)
        {
            data.RemoveAt(edge2);
            data.Insert(edge2, new Tuple<int, int>(im, jm));
        }
        else
        {
            data.RemoveAt(edge2);
            data.Insert(edge2, new Tuple<int, int>(jm, im));
        }
    }
    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.Count; ++i)
        {
            sb.AppendFormat("({0}, {1})", data[i].Item1, data[i].Item2);
            if (i != data.Count - 1) sb.Append(", ");
        }
        return sb.ToString();
    }
    public void RemoveAt(int i)
    {
        data.RemoveAt(i);
    }
    public Tuple<int, int> this[int key]
    {
        get
        {
            return data[key];
        }
        set
        {
            data[key] = value;
        }
    }
    public Tuple<int, int> Last()
    {
        return data.Last();
    }
    public void Parse(string line)
    {
        char[] chromosomeDelimiter = { '(' };
        line = line.Replace(")", " ");
        string[] nodes = line.Split(chromosomeDelimiter, StringSplitOptions.RemoveEmptyEntries);

        foreach (string node in nodes)
        {
            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] node_string = node.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            Add(Int32.Parse(node_string[0]), Int32.Parse(node_string[1]));
        }
    }
    List<Tuple<int, int>> data;
}