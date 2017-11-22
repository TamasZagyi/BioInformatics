﻿using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GenomeSequencing
{
    class Program
    {
        static List<string> Composition(string Text, int k)
        {
            List<string> kmers = new List<string>();
            for(int i=0;i<=Text.Length - k; ++i)
            {
                kmers.Add(Text.Substring(i, k));
            }

            kmers.Sort();
            return kmers;
        }
        static void DisplayStringList(List<string> list)
        {
            foreach(string s in list)
            {
                Console.WriteLine(s);
            }
        }
        static string StringSpelledByGenomePath(List<string> GenomePath)
        {
            StringBuilder sb = new StringBuilder();
           
            foreach(string s in GenomePath)
            {
                sb.Append(s[0]);
            }

            string last = GenomePath[GenomePath.Count - 1];
            sb.Append(last.Substring(1, last.Length - 1));

            return sb.ToString();
        }
        static Dictionary<string, List<string>> Overlap(List<string> Patterns)
        {
            Dictionary<string, List<string>> Graph = new Dictionary<string, List<string>>() ;
            for(int i=0; i < Patterns.Count;++i)
            {
                for(int j=0;j<Patterns.Count;++j)
                {
                    if(i!=j)
                    {
                        if(Patterns[i].Substring(1) == Patterns[j].Substring(0,Patterns[j].Length-1))
                        {
                            if (!Graph.ContainsKey(Patterns[i])) Graph.Add(Patterns[i], new List<string>());
                            Graph[Patterns[i]].Add(Patterns[j]);
                        }
                    }
                }
            }

            return Graph;
        }
        static Dictionary<string, List<string>> DeBruijn(string Text, int k)
        {
            int l = k - 1;
            List<string> lmers = new List<string>();
            for (int i = 0; i <= Text.Length - l; ++i)
            {
                lmers.Add(Text.Substring(i, l));
            }

            Dictionary<string, List<string>> graph = new Dictionary<string, List<string>>();

            for (int i = 0; i < lmers.Count - 1;++i)
            {
                if (!graph.ContainsKey(lmers[i]))
                {
                    graph.Add(lmers[i], new List<string>());
                    graph[lmers[i]].Add(lmers[i + 1]);
                }
                else
                {
                    graph[lmers[i]].Add(lmers[i + 1]);
                    graph[lmers[i]].Sort();
                }
            }

            return graph;
        }
        static Dictionary<string, List<string>> DeBruijn(List<string> Patterns)
        {
            Dictionary<string, List<string>> graph = new Dictionary<string, List<string>>();

            foreach (string pattern in Patterns)
            {
                string key = pattern.Substring(0, pattern.Length - 1);
                if (!graph.ContainsKey(key))
                {
                    graph.Add(key, new List<string>());
                    graph[key].Add(pattern.Substring(1));
                }
                else
                {
                    graph[key].Add(pattern.Substring(1));
                    graph[key].Sort();
                }
            }

            return graph;
        }
        /*
        EULERIANCYCLE(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
            Cycle ← Cycle’
        return Cycle
            */
        static List<string> EulerianCycle(Dictionary<string, List<string>> Graph)
        {
            List<string> Cycle = new List<string>();

            //  deep copy of graph
            Dictionary<string, List<string>> ReducedGraph = new Dictionary<string, List<string>>();
            foreach(KeyValuePair<string, List<string>> edge in Graph)
            {
                ReducedGraph.Add(edge.Key , new List<string>(edge.Value));
            }

            Cycle = RandomWalk(ReducedGraph, ReducedGraph.ElementAt(0).Key);

            while(ReducedGraph.Count > 0)
            {
                for(int i = 0; i < Cycle.Count; ++i)
                {
                    if(ReducedGraph.ContainsKey(Cycle[i]))
                    {
                        Cycle.AddRange(Cycle.GetRange(1, i-1));
                        Cycle.RemoveRange(0, i);
                        Cycle.AddRange(RandomWalk(ReducedGraph, Cycle[0]));
                        break;
                    }
                }
            }

            return Cycle;
        }
        static List<string> RandomWalk(Dictionary<string, List<string>> ReducedGraph, string Position)
        {
            List<string> Walk = new List<string>();

            if( ReducedGraph.ContainsKey(Position) )
                Walk.Add(Position);

            while( ReducedGraph.ContainsKey(Position) )
            {
                Walk.Add(ReducedGraph[Position][0]);
                ReducedGraph[Position].RemoveAt(0);
                if(ReducedGraph[Position].Count == 0)
                {
                    ReducedGraph.Remove(Position);
                }
                Position = Walk[Walk.Count-1];
            }

            return Walk;
        }
        static List<string> EulerianPath(Dictionary<string, List<string>> Graph)
        {
            string first = "";
            string last = "";
            //  deep copy of graph & determine first / last node
            Dictionary<string, List<string>> ExtendedGraph = new Dictionary<string, List<string>>();
            Dictionary<string,int> balance = new Dictionary<string, int>();
            foreach (KeyValuePair<string, List<string>> edge in Graph)
            {
                ExtendedGraph.Add(edge.Key, new List<string>(edge.Value));
                if (!balance.ContainsKey(edge.Key))
                    balance.Add(edge.Key, edge.Value.Count);
                else balance[edge.Key]+= edge.Value.Count;
                foreach (string target in edge.Value)
                {
                    if (!balance.ContainsKey(target))
                        balance.Add(target, -1);
                    else balance[target]--;
                }
            }
            foreach(KeyValuePair<string, int> edge in balance)
            {
                if(edge.Value < 0)
                {
                    last = edge.Key;
                }
                if (edge.Value > 0)
                {
                    first = edge.Key;
                }
            }

            //  connect the last and first nodes
            ExtendedGraph.Add(last, new List<string>());
            ExtendedGraph[last].Add(first);

            //  calculate Eulerian cycle
            List<string> Cycle = EulerianCycle(ExtendedGraph);
            Cycle.RemoveAt(Cycle.Count - 1);
            
            //  re-arrange the cycle and remove the last - first connection
            if(!(Cycle[0] == first && Cycle[Cycle.Count-1] == last))
            {
                int firstPosition = Cycle.IndexOf(first);
                if(firstPosition == 0) firstPosition = Cycle.IndexOf(first,firstPosition+1);
                while ( Cycle[firstPosition - 1] != last)
                {
                    firstPosition = Cycle.IndexOf(first, firstPosition + 1);
                }
                Cycle.AddRange(Cycle.GetRange(0, firstPosition));
                Cycle.RemoveRange(0, firstPosition);
            }
            
            return Cycle;
        }
        static List<string> UniversalString(int k)
        {
            List<string> Patterns = new List<string>();

            for (int i=0;i<Math.Pow(2,k);++i)
            {
                Patterns.Add(Convert.ToString(i, 2).PadLeft(k,'0'));
            }

            Dictionary<string, List<string>> Graph = DeBruijn(Patterns);

            List<string>  Path = EulerianCycle(Graph);

            return Path;
        }
        static List<string> PairedComposition(string Text, int k, int d)
        {
            List<string> Composition = new List<string>();
            for(int i=0;i<=Text.Length - (2*k + d); ++i)
            {
                Composition.Add(Text.Substring(i, k) + Text.Substring(i + k + d, k));
            }
            //Composition.Sort();
            return Composition;
        }
        static string PrefixString(List<string> patterns, int k)
        {
            StringBuilder sb = new StringBuilder();
            
            foreach(string pattern in patterns)
            {
                sb.Append(pattern[0]);
            }
            sb.Append(patterns[patterns.Count - 1].Substring(1, k - 1));
            return sb.ToString();
        }
        static string SuffixString(List<string> patterns, int k)
        {
            StringBuilder sb = new StringBuilder();

            foreach (string pattern in patterns)
            {
                sb.Append(pattern[k]);
            }
            sb.Append(patterns[patterns.Count - 1].Substring(k + 1));
            return sb.ToString();
        }
        static string StringSpelledByPatterns(List<string> patterns, int k)
        {
            return PrefixString(patterns, k) + SuffixString(patterns, k).Substring(patterns.Count-k);
        }
        /*
        StringSpelledByGappedPatterns(GappedPatterns, k, d)
        FirstPatterns ← the sequence of initial k-mers from GappedPatterns
        SecondPatterns ← the sequence of terminal k-mers from GappedPatterns
        PrefixString ← StringSpelledByPatterns(FirstPatterns, k)
        SuffixString ← StringSpelledByPatterns(SecondPatterns, k)
        for i = k + d + 1 to |PrefixString|
            if the i-th symbol in PrefixString does not equal the(i - k - d)-th symbol in SuffixString
                return "there is no string spelled by the gapped patterns"
        return PrefixString concatenated with the last k + d symbols of SuffixString
        */
        static string StringSpelledByGappedPatterns(List<string> GappedPatterns, int k, int d)
        {
            string prefix = PrefixString(GappedPatterns, k);
            string suffix = SuffixString(GappedPatterns, k);

            string pre_match = prefix.Substring(k + d);
            string suff_match = suffix.Substring(0, prefix.Length - (k + d));

           if ( pre_match != suff_match  )
            {
                return "there is no string spelled by the gapped patterns";
            }
            return prefix + suffix.Substring(suffix.Length - (k + d));
        }
        /*
        MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
        */
        static List<List<string>> MaximalNonBranchingPaths( Dictionary<string, List<string>> Graph )
        {
            //  deep copy of graph
            Dictionary<string, List<string>> ReducedGraph = new Dictionary<string, List<string>>();
            Dictionary<string, int> balance = new Dictionary<string, int>();
            foreach (KeyValuePair<string, List<string>> edge in Graph)
            {
                ReducedGraph.Add(edge.Key, new List<string>(edge.Value));
                if (!balance.ContainsKey(edge.Key))
                    balance.Add(edge.Key, edge.Value.Count);
                else balance[edge.Key] += edge.Value.Count;
                foreach (string target in edge.Value)
                {
                    if (!balance.ContainsKey(target))
                        balance.Add(target, -1);
                    else balance[target]--;
                }
            }

            List<List<string>> Paths = new List<List<string>>();
            List<string> Keys = Graph.Keys.ToList();
            foreach(string key in Keys)
            {
                if(!(balance[key] == 0 && Graph[key].Count == 1))
        //        if (Graph[key].Count > 1)
                {
                    foreach(string edge in Graph[key])
                    {
                        List<string> NonBranchingPath = new List<string>();
                        string w = edge;
                        NonBranchingPath.Add(key);
                        NonBranchingPath.Add(w);
                        while (Graph.ContainsKey(w) && balance[w] == 0)
                        {
                            if (Graph[w].Count != 1) break;
                            ReducedGraph.Remove(w);
                            w = Graph[w][0];
                            NonBranchingPath.Add(w);   
                        }

                        Paths.Add(NonBranchingPath);
                    }
                    ReducedGraph.Remove(key);
                }
            }
            while (ReducedGraph.Count > 0)
            {
                KeyValuePair<string, List<string>> node = ReducedGraph.ElementAt(0);

                List<string> Cycle = new List<string>();
                string w = node.Key;
                Cycle.Add(w);
                do
                {
                    if(ReducedGraph[w].Count != 1)
                    {
                        throw new Exception("Illegal cycle found");
                    }
                    string to_be_removed = w;  
                    w = ReducedGraph[w][0];
                    ReducedGraph.Remove(to_be_removed);
                    Cycle.Add(w);
                    
                } while (ReducedGraph.ContainsKey(w));
                Paths.Add(Cycle);
            }
            return Paths;
        }
        static string[,] RNA_codon_table = new string[,]{
{"AAA", "K"},{"AAC", "N"},{"AAG", "K"},{"AAU", "N"},{"ACA", "T"},{"ACC", "T"},{"ACG", "T"},{"ACU", "T"},{"AGA", "R"},{"AGC", "S"},{"AGG", "R"},{"AGU", "S"},{"AUA", "I"},
{"AUC", "I"},{"AUG", "M"},{"AUU", "I"},{"CAA", "Q"},{"CAC", "H"},{"CAG", "Q"},{"CAU", "H"},{"CCA", "P"},{"CCC", "P"},{"CCG", "P"},{"CCU", "P"},{"CGA", "R"},{"CGC", "R"},
{"CGG", "R"},{"CGU", "R"},{"CUA", "L"},{"CUC", "L"},{"CUG", "L"},{"CUU", "L"},{"GAA", "E"},{"GAC", "D"},{"GAG", "E"},{"GAU", "D"},{"GCA", "A"},{"GCC", "A"},{"GCG", "A"},
{"GCU", "A"},{"GGA", "G"},{"GGC", "G"},{"GGG", "G"},{"GGU", "G"},{"GUA", "V"},{"GUC", "V"},{"GUG", "V"},{"GUU", "V"},{"UAA", ""},{"UAC", "Y"},{"UAG", ""},{"UAU", "Y"},
{"UCA", "S"},{"UCC", "S"},{"UCG", "S"},{"UCU", "S"},{"UGA", ""},{"UGC", "C"},{"UGG", "W"},{"UGU", "C"},{"UUA", "L"},{"UUC", "F"},{"UUG", "L"},{"UUU", "F"}};
        static string RNAToAminoAcid(string RNA)
        {
            StringBuilder sb = new StringBuilder();
            
            for (int i = 0; i <= RNA.Length - 3;i+=3 )
            {
                string codon = RNA.Substring(i, 3);
                for(int j=0;j<RNA_codon_table.GetLength(0);++j )
                { 
                    if(RNA_codon_table[j,0] == codon)
                    {
                        sb.Append(RNA_codon_table[j, 1]);
                        break;
                    }
                }
            }

            return sb.ToString();
        }
        //  A - T G - C
        static string ReverseComplement(string Text)
        {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < Text.Length; ++i)
            {
                char complement = Text[Text.Length - (i + 1)];
                switch (complement)
                {
                    case 'A': complement = 'T'; break;
                    case 'T': complement = 'A'; break;
                    case 'G': complement = 'C'; break;
                    case 'C': complement = 'G'; break;
                }
                builder.Append(complement);
            }

            return builder.ToString();
        }
        static List<string> SubstringsEncodingAminoAcidForward(string DNA, string Peptide)
        {
            List<string> sub = new List<string>();
            for(int i=0;i<=DNA.Length - Peptide.Length*3;++i)
            {
                string RNA = DNA.Substring(i, Peptide.Length * 3).Replace('T', 'U');
                if(RNAToAminoAcid(RNA)==Peptide)
                {
                    sub.Add(DNA.Substring(i, Peptide.Length * 3));
                }
            }

            return sub;
        }
        static List<string> SubstringsEncodingAminoAcid(string DNA, string Peptide)
        {
            List<string> forward = SubstringsEncodingAminoAcidForward(DNA, Peptide);

            string DNAReverse = ReverseComplement(DNA);
            
            foreach (string sub in SubstringsEncodingAminoAcidForward(DNAReverse, Peptide))
            {
                forward.Add(ReverseComplement(sub));
            }

            return forward;
        }
        static void Main(string[] args)
        {
            Console.WriteLine("Where to read from? Console or File? (C/F)");
            string input = Console.ReadLine().ToUpper().Trim();

          //  Console.WriteLine(RNAToAminoAcid(Console.ReadLine()));

             StreamReader sr = new StreamReader("dataset.txt");
            string DNA="";
            string Peptide = "";
         /*   while (!sr.EndOfStream)
            {
                string line = sr.ReadLine();
                if (line.Length > 0) RNA+=line;
            }
            */

            DNA = sr.ReadLine();
            Peptide = sr.ReadLine();
            
            StreamWriter sw = new StreamWriter("result.txt");
            foreach(string sub in SubstringsEncodingAminoAcid(DNA, Peptide))
                sw.WriteLine(sub);
            sw.Flush();

            /*
            List<string> Patterns = new List<string>();
            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();
                    if (line.Length > 0) Patterns.Add(line);
                }
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");

            Dictionary<string, List<string>> Graph = DeBruijn(Patterns);

            List<List<string>> paths = MaximalNonBranchingPaths(Graph);

            foreach (List<string> path in paths)
            {
                foreach (string s in path)
                {
                    sw.Write("{0}", s[0]);
                }
                sw.Write("{0}", path[path.Count - 1].Substring(1));

                sw.Write("\n");
            }

            

            sw.Flush();
            */
            /*
            Dictionary<string, List<string>> Graph = new Dictionary<string, List<string>>();

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();
                    if (line.Length > 0)
                    {
                        string[] tokens = { "->", "," };
                        string[] values = line.Replace(" ", "").Split(tokens, StringSplitOptions.RemoveEmptyEntries);
                        Graph.Add(values[0], new List<string>());
                        for (int i = 1; i < values.GetLength(0); ++i)
                            Graph[values[0]].Add(values[i]);
                    }
                }
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");

            List<List<string>> paths = MaximalNonBranchingPaths(Graph);

            foreach(List<string> path in paths)
            {
                sw.Write(path[0]);

                for (int i = 1; i < path.Count; ++i)
                {
                    sw.Write(" -> {0}", path[i]);
                }

                sw.Write("\n");
            }
            
            sw.Flush();
            */
            
 /*           int k = 1;
            int d = 1;
            
            Dictionary<string, List<string>> Graph = new Dictionary<string, List<string>>();
            if (input == "C")
            {
 
                k = 3;
                d = 1;
                string Text = Console.ReadLine();
                List<string> pairs = PairedComposition("TAATGCCATGGGATGTT",k,d);
                foreach(string line in pairs)
                {
                    string read1 = line.Substring(0, k);
                    string read2 = line.Substring(k);
                    string key = read1.Substring(0, k - 1) + read2.Substring(0, k - 1);
                    string value = read1.Substring(1) + read2.Substring(1);
                    if (!Graph.ContainsKey(key)) Graph.Add(key, new List<string>());
                    Graph[key].Add(value);
                }
            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                k = Int32.Parse(sr.ReadLine());
                d = Int32.Parse(sr.ReadLine());
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();
                    
                    if (line.Length > 0)
                    {
                        string read1 = line.Replace("|", "").Substring(0, k);
                        string read2 = line.Replace("|", "").Substring(k);
                        string key = read1.Substring(0, k - 1) + read2.Substring(0, k - 1);
                        string value = read1.Substring(1) + read2.Substring(1);
                        if (!Graph.ContainsKey(key)) Graph.Add(key, new List<string>());
                            Graph[key].Add(value);
                    }
                }
            }
            else
            {
                return;
            }

            List<string> path = EulerianPath(Graph);
            StreamWriter sw = new StreamWriter("result.txt");
        
            sw.Write(StringSpelledByGappedPatterns(path, k-1, d+1));
            sw.Flush();
    */        
            //StreamWriter sw = new StreamWriter("result.txt");
            /*
                        List<string> composition = PairedComposition(Text,k, d);
                        foreach (string pair in composition)
                        {
                            Console.Write("({0}|{1}) ", pair.Substring(0,k),pair.Substring(k));
                            //Console.Write("{0} {1} ", pair.Substring(0,k),pair.Substring(k));
                        }
                        Console.WriteLine(" ");
                        Console.WriteLine("Prefix: {0}", PrefixString(composition,k));
                        Console.WriteLine("Suffix: {0}", SuffixString(composition, k));
                        Console.WriteLine("Genom: {0}", StringSpelledByPatterns(composition, k));
              */          //sw.Write("{0}", path[path.Count - 1].Substring(1));

            //sw.Flush();
            /*
            int k = 1;
            if (input == "C")
            {
                k = Int32.Parse(Console.ReadLine());
            }
            else if (input == "F")
            {
                
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");
            
            List<string> path = UniversalString(k);
            for (int i = 0; i < path.Count - 1; ++i)
            {
                sw.Write("{0}", path[i][0]);
            }
            //sw.Write("{0}", path[path.Count - 1].Substring(1));

            sw.Flush();
            */
            /*
            List<string> Patterns = new List<string>();
            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                //sr.ReadLine();  //  read k
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();
                    if (line.Length > 0) Patterns.Add(line);
                }
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");

            Dictionary<string, List<string>> Graph = DeBruijn(Patterns);
            List<string> path = EulerianPath(Graph);
            
            foreach (string s in path)
            {
                sw.Write("{0}", s[0]);
            }
            sw.Write("{0}", path[path.Count - 1].Substring(1));

            sw.Flush();
            */
            /*
            Dictionary<string, List<string>> Graph = new Dictionary<string, List<string>>();

           if (input == "C")
           {

           }
           else if (input == "F")
           {
               StreamReader sr = new StreamReader("dataset.txt");
               while (!sr.EndOfStream)
               {
                   string line = sr.ReadLine();
                    if (line.Length > 0)
                    {
                        string[] tokens = { "->", "," };
                        string[] values = line.Replace(" ", "").Split(tokens,StringSplitOptions.RemoveEmptyEntries);
                        Graph.Add(values[0], new List<string>());
                        for (int i = 1; i < values.GetLength(0); ++i)
                            Graph[values[0]].Add(values[i]);
                    }
               }
           }
           else
           {
               return;
           }

           StreamWriter sw = new StreamWriter("result.txt");

           List<string> Cycle = EulerianPath(Graph);

           sw.Write(Cycle[0]);

           for (int i = 1;i < Cycle.Count;++i)
           {
               sw.Write("->{0}", Cycle[i]);
           }

           sw.Flush();
           */

            /*          if (input == "C")
                      {

                      }
                      else if (input == "F")
                      {
                          StreamReader sr = new StreamReader("dataset.txt");
                          while (!sr.EndOfStream)
                          {
                              string line = sr.ReadLine();
                              if (line.Length > 0) Patterns.Add(line);
                          }
                      }
                      else
                      {
                          return;
                      }

                      StreamWriter sw = new StreamWriter("result.txt");

                      Dictionary<string, List<string>> Graph = DeBruijn(Patterns);
                      List<KeyValuePair<string, List<string>>> list = Graph.ToList();
                      list.Sort((a, b) =>
                      {
                          return a.Key.CompareTo(b.Key);
                      });
                      foreach (KeyValuePair<string, List<string>> source in list)
                      {
                          sw.Write("{0} -> ", source.Key);
                          for (int i = 0; i < source.Value.Count; ++i)
                          {
                              sw.Write("{0}", source.Value[i]);
                              if (i != source.Value.Count - 1) sw.Write(",");
                              else sw.Write("\n");
                          }
                      }

                      sw.Flush();
          */
            /*           string Text;
                       int k;

                       if (input == "C")
                       {
                           Console.WriteLine("Enter k: ");
                           k = Int32.Parse(Console.ReadLine());
                           Console.WriteLine("Enter Text: ");
                           Text = Console.ReadLine();
                       }
                       else if (input == "F")
                       {
                           StreamReader sr = new StreamReader("dataset.txt");
                           k = Int32.Parse(sr.ReadLine());
                           Text = sr.ReadLine();
                       }
                       else
                       {
                           return;
                       }

                       StreamWriter sw = new StreamWriter("result.txt");

                       Dictionary<string, List<string>> Graph = DeBruijn(Text,k);
                       List<KeyValuePair<string, List<string>>> list = Graph.ToList();
                       list.Sort((a, b) =>
                       {
                           return a.Key.CompareTo(b.Key);
                       });
                       foreach (KeyValuePair<string, List<string>> source in list)
                       {
                           sw.Write("{0} -> ", source.Key);
                           for (int i = 0; i < source.Value.Count; ++i)
                           {
                               sw.Write("{0}", source.Value[i]);
                               if (i != source.Value.Count - 1) sw.Write(",");
                               else sw.Write("\n");
                           }
                       }

                       sw.Flush();
           */
            /*
                        List<string> Patterns = new List<string>();

                        if (input == "C")
                        {

                        }
                        else if (input == "F")
                        {
                            StreamReader sr = new StreamReader("dataset.txt");
                            while (!sr.EndOfStream)
                            {
                                string line = sr.ReadLine();
                                if (line.Length > 0) Patterns.Add(line);
                            }
                        }
                        else
                        {
                            return;
                        }

                        StreamWriter sw = new StreamWriter("result.txt");

                        Dictionary<string, List<string>> Graph = Overlap(Patterns);
                        List<KeyValuePair<string, List<string>>> list = Graph.ToList();
                        list.Sort((a, b) =>
                        {
                            return a.Key.CompareTo(b.Key);
                        });
                        foreach (KeyValuePair<string, List<string>> source in list)
                        {
                            foreach(string target in source.Value )
                            {
                                sw.WriteLine("{0} -> {1}", source.Key, target);
                            }
                        }

                        sw.Flush();
            */
            /*            List<string> GenomePath = new List<string>();

                        if(input == "C")
                        {

                        }
                        else if(input == "F")
                        {
                            StreamReader sr = new StreamReader("dataset.txt");
                            while(!sr.EndOfStream)
                            {
                                string line = sr.ReadLine();
                                if(line.Length > 0 ) GenomePath.Add(line);                    
                            }
                        }
                        else
                        {
                            return;
                        }

                        StreamWriter sw = new StreamWriter("result.txt");

                        sw.WriteLine(StringSpelledByGenomePath(GenomePath));

                        sw.Flush();
            */

            /*            string Text;
                        int k;

                        if(input == "C")
                        {
                            Console.WriteLine("Enter k: ");
                            k = Int32.Parse(Console.ReadLine());
                            Console.WriteLine("Enter Text: ");
                            Text = Console.ReadLine();
                        }
                        else if(input == "F")
                        {
                            StreamReader sr = new StreamReader("dataset.txt");
                            k = Int32.Parse(sr.ReadLine());
                            Text = sr.ReadLine();
                        }
                        else
                        {
                            return;
                        }

                        DisplayStringList(Composition(Text, k));

                        StreamWriter sw = new StreamWriter("result.txt");
                        foreach(string s in Composition(Text, k))
                        {
                            sw.WriteLine(s);
                        }
                        sw.Flush();
            */
            Console.ReadKey();
        }
    }
}
