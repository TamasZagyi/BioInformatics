using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ComparingGenomes
{
    class Program
    {
        /*
        DPCHANGE(money, Coins)
         MinNumCoins(0) ← 0
         for m ← 1 to money
            MinNumCoins(m) ← ∞
                for i ← 1 to |Coins|
                    if m ≥ coini
                        if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                            MinNumCoins(m) ← MinNumCoins(m - coini) + 1
        output MinNumCoins(money)
        */
        static int DPChange(int money, List<int> coins)
        {
            List<int> minNumCoins = new List<int>();
            minNumCoins.Add(0);

            for (int m = 1; m <= money; ++m)
            {
                minNumCoins.Add(Int32.MaxValue);
                foreach (int coini in coins)
                {
                    if (m >= coini)
                        if (minNumCoins[m - coini] + 1 < minNumCoins[m])
                            minNumCoins[m] = minNumCoins[m - coini] + 1;
                }
            }

            return minNumCoins[money];
        }
        static void SolveChange()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int money = Int32.Parse(input.ReadLine());
            List<int> coins = new List<int>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    foreach (string word in words)
                        coins.Add(Int32.Parse(word));
                }
            }

            //  compute

            int minNumCoins = DPChange(money, coins);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(minNumCoins.ToString());
            output.Flush();
        }
        /*
         MANHATTANTOURIST(n, m, Down, Right)
            s0, 0 ← 0
            for i ← 1 to n
                si, 0 ← si-1, 0 + downi, 0
            for j ← 1 to m
                s0, j ← s0, j−1 + right0, j
            for i ← 1 to n
                for j ← 1 to m
                    si, j ← max{si - 1, j + downi, j, si, j - 1 + righti, j}
            return sn, m
        */
        static int ManhattanTourist(int n, int m, int[,] down, int[,] right)
        {
            int[,] s = new int[n + 1, m + 1];

            s[0, 0] = 0;
            for (int i = 1; i < n; ++i) s[i, 0] = s[i - 1, 0] + down[i, 0];
            for (int j = 1; j < m; ++j) s[0, j] = s[0, j - 1] + right[0, j];
            for (int i = 1; i <= n; ++i)
                for (int j = 1; j <= m; ++j)
                {
                    s[i, j] = Math.Max(s[i - 1, j] + down[i, j],
                                        s[i, j - 1] + right[i, j]);
                }

            return s[n, m];
        }
        static void SolveManhattan()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            int m = Int32.Parse(input.ReadLine());

            // n × (m + 1) matrix Down and an (n + 1) × m matrix Right.

            int[,] down = new int[n + 1, m + 1];
            int[,] right = new int[n + 1, m + 1];

            for (int i = 1; i < down.GetLength(0); ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    int j = 0;
                    foreach (string word in words)
                        down[i, j++] = (Int32.Parse(word));
                }
            }

            input.ReadLine();   //  read -

            for (int i = 0; i < right.GetLength(0); ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    int j = 1;
                    foreach (string word in words)
                        right[i, j++] = (Int32.Parse(word));
                }
            }

            //  compute

            int maxPath = ManhattanTourist(n, m, down, right);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(maxPath.ToString());
            output.Flush();
        }
        /*
        LCSBACKTRACK(v, w)
        for i ← 0 to |v|
            si, 0 ← 0
        for j ← 0 to |w| 
            s0, j ← 0
        for i ← 1 to |v|
            for j ← 1 to |w|
                si, j ← max{si-1, j, si,j-1, si-1, j-1 + 1 (if vi = wj)}
                if si,j = si-1,j
                    Backtracki, j ← "↓"
                else if si, j = si, j-1
                    Backtracki, j ← "→"
                else if si, j = si-1, j-1 + 1 and vi = wj
                    Backtracki, j ← "↘"
        return Backtrack
        */
        static int[,] LCSBacktrack(string v, string w)
        {
            int[,] s = new int[v.Length + 1, w.Length + 1];
            int[,] backtrack = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i) s[i, 0] = 0;
            for (int j = 0; j <= w.Length; ++j) s[0, j] = 0;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    int down = s[i - 1, j];
                    int right = s[i, j - 1];
                    int diag = s[i - 1, j - 1] + (v[i - 1] == w[j - 1] ? 1 : 0);
                    s[i, j] = Math.Max(down, Math.Max(right, diag));
                    if (s[i, j] == down) backtrack[i, j] = -1;
                    else if (s[i, j] == right) backtrack[i, j] = 1;
                    else if (s[i, j] == diag) backtrack[i, j] = 0;
                }

            return backtrack;
        }
        static int[,] LCSBacktrack(string v, string w, ref int score)
        {
            int[,] s = new int[v.Length + 1, w.Length + 1];
            int[,] backtrack = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i) s[i, 0] = -2 * i;
            for (int j = 0; j <= w.Length; ++j) s[0, j] = -2 * j;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    int down = s[i - 1, j] - 2;
                    int right = s[i, j - 1] - 2;
                    int diag = s[i - 1, j - 1] + (v[i - 1] == w[j - 1] ? 1 : 0);
                    s[i, j] = Math.Max(down, Math.Max(right, diag));
                    if (s[i, j] == down) backtrack[i, j] = -1;
                    else if (s[i, j] == right) backtrack[i, j] = 1;
                    else if (s[i, j] == diag) backtrack[i, j] = 0;
                }

            score = s[v.Length, w.Length];

            return backtrack;
        }
        static int[,] LCSBacktrack(string v, string w, ScoringMatrix m, ref int score)
        {
            int[,] s = new int[v.Length + 1, w.Length + 1];
            int[,] backtrack = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i) s[i, 0] = i * -5;
            for (int j = 0; j <= w.Length; ++j) s[0, j] = j * -5;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    int down = s[i - 1, j] - 5;
                    int right = s[i, j - 1] - 5;
                    int diag = s[i - 1, j - 1] + m.Weight(v[i - 1], w[j - 1]);

                    s[i, j] = Math.Max(down, Math.Max(right, diag));

                    if (s[i, j] == down) backtrack[i, j] = -1;
                    else if (s[i, j] == right) backtrack[i, j] = 1;
                    else if (s[i, j] == diag) backtrack[i, j] = 0;
                }

            score = s[v.Length, w.Length];

            return backtrack;
        }
        static int[,] LCSBacktrackLocal(string v, string w, ScoringMatrix m, ref int score, ref int si, ref int sj)
        {
            int[,] s = new int[v.Length + 1, w.Length + 1];
            int[,] backtrack = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i) s[i, 0] = i * -5;
            for (int j = 0; j <= w.Length; ++j) s[0, j] = j * -5;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    int down = s[i - 1, j] - 5;
                    int right = s[i, j - 1] - 5;
                    int diag = s[i - 1, j - 1] + m.Weight(v[i - 1], w[j - 1]);

                    s[i, j] = Math.Max(0, Math.Max(down, Math.Max(right, diag)));

                    if (s[i, j] == 0) backtrack[i, j] = -2;
                    else if (s[i, j] == down) backtrack[i, j] = -1;
                    else if (s[i, j] == right) backtrack[i, j] = 1;
                    else if (s[i, j] == diag) backtrack[i, j] = 0;
                }

            score = Int32.MinValue;
            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    if (s[i, j] > score)
                    {
                        si = i;
                        sj = j;
                        score = s[i, j];
                    }
                }

            return backtrack;
        }
        /*
        OUTPUTLCS(Backtrack, v, i, j)
        if i = 0 or j = 0
            return
        if backtracki, j = "↓"
            OUTPUTLCS(Backtrack, v, i - 1, j)
        else if Backtracki, j = "→"
            OUTPUTLCS(Backtrack, v, i, j - 1)
        else
            OUTPUTLCS(Backtrack v, i - 1, j - 1)
            output vi
        */
        static void OutputLCS(int[,] backtrack, string v, int i, int j, ref List<char> output)
        {
            if (i == 0 || j == 0) return;
            if (backtrack[i, j] == -1) OutputLCS(backtrack, v, i - 1, j, ref output);
            else if (backtrack[i, j] == 1) OutputLCS(backtrack, v, i, j - 1, ref output);
            else
            {
                OutputLCS(backtrack, v, i - 1, j - 1, ref output);
                output.Add(v[i - 1]);
            }
        }
        static void OutputLCS(int[,] backtrack, string v, string w, int i, int j, ref List<char> ov, ref List<char> ow)
        {
            while (i != 0 && j != 0)
            {
                if (backtrack[i, j] == 0)
                {
                    ov.Add(v[i - 1]);
                    ow.Add(w[j - 1]);
                    --i; --j;
                }
                else if (backtrack[i, j] == 1)
                {
                    ov.Add('-');
                    ow.Add(w[j - 1]);
                    --j;
                }
                else if (backtrack[i, j] == -1)
                {
                    ov.Add(v[i - 1]);
                    ow.Add('-');
                    --i;
                }
                else if (backtrack[i, j] == -2)
                {
                    i = j = 0;
                }
            }

            while (i > 0)
            {
                ov.Add(v[i - 1]);
                ow.Add('-');
                i--;
            }

            while (j > 0)
            {
                ow.Add(w[j - 1]);
                ov.Add('-');
                i--;
            }

            ov.Reverse();
            ow.Reverse();
        }
        static void SolveLCS()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            //  compute

            List<char> output = new List<char>();
            OutputLCS(LCSBacktrack(v, w), v, v.Length, w.Length, ref output);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            foreach (char c in output) result.Write(c);
            result.Flush();
        }
        static void SolveGlobalAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            ScoringMatrix m = new ScoringMatrix();
            m.Load("BLOSUM62.txt");

            //  compute
            int score = 0;
            int[,] backtrack = LCSBacktrack(v, w, m, ref score);

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            OutputLCS(backtrack, v, w, v.Length, w.Length, ref ov, ref ow);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(score.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }
        static void SolveLocalAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            ScoringMatrix m = new ScoringMatrix();
            m.Load("PAM250_1.txt");

            //  compute
            int score = 0;
            int i = 0, j = 0;
            int[,] backtrack = LCSBacktrackLocal(v, w, m, ref score, ref i, ref j);

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            OutputLCS(backtrack, v, w, i, j, ref ov, ref ow);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(score.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }

        static void SolveLDAG()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int source = Int32.Parse(input.ReadLine());
            int sink = Int32.Parse(input.ReadLine());
            Graph<int, int> graph = new Graph<int, int>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");
                line = line.Replace(":", " ");

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    int s = Int32.Parse(words[0]);
                    int t = Int32.Parse(words[1]);

                    graph.AddEdge(s,
                                    t,
                                    Int32.Parse(words[2]));
                }
            }

            //  compute

            List<Edge<int, int>> path = graph.LongestPath(source, sink);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            output.WriteLine("{0}", Graph<int, int>.Length(path));

            foreach (Edge<int, int> edge in path)
                output.Write("{0}->", edge.source.label);
            output.Write("{0}", path.Last().target.label);

            output.Flush();
        }
        //  Distance matrix
        static double[,] DistanceMatrix(Graph<int, int> graph, int n)
        {
            double[,] matrix = new double[n, n];
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < i; ++j)
                {
                    if (i == j) matrix[i, j] = 0;
                    else
                    {
                        matrix[i, j] = matrix[j, i] = Graph<int, int>.Length(graph.Path(i, j));
                    }
                }
            return matrix;
        }
        static void SolveDistanceMatrix()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            Graph<int, int> graph = new Graph<int, int>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");
                line = line.Replace(":", " ");

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    int source = Int32.Parse(words[0]);
                    int target = Int32.Parse(words[1]);

                    graph.AddEdge(source,
                                    target,
                                    Int32.Parse(words[2]));
                }
            }

            //  compute

            double[,] matrix = DistanceMatrix(graph, n);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    output.Write("{0} ", matrix[i, j]);
                }
                output.Write("\n");
            }

            output.Flush();
        }
        //  LimbLength
        static int LimbLength(int[,] D, int n, int j)
        {
            if (j < 0 || j > n - 1)
                throw new Exception("j out of range");
            int min = Int32.MaxValue;
            for (int i = 0; i < n - 1; ++i)
            {
                if (i == j) continue;
                for (int k = i + 1; k < n; ++k)
                {
                    if (k == j) continue;
                    if (D[i, j] + D[j, k] - D[i, k] < min)
                        min = D[i, j] + D[j, k] - D[i, k];
                }
            }
            return min / 2;
        }
        static void SolveLimbLength()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            int j = Int32.Parse(input.ReadLine());
            int[,] D = new int[n, n];
            for (int i = 0; i < n; ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < n; ++k)
                        D[i, k] = Int32.Parse(words[k]);
                }
            }

            //  compute

            int length = LimbLength(D, n, j);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(length.ToString());
            output.Flush();
        }
        /*
        AdditivePhylogeny(D, n)
        if n = 2
            return the tree consisting of a single edge of length D1,2
        limbLength ← Limb(D, n)
        for j ← 1 to n - 1
            Dj,n ← Dj,n - limbLength
            Dn,j ← Dj,n
        (i,n,k) ← three leaves such that Di,k = Di,n + Dn,k
        x ← Di,n
        remove row n and column n from D
        T ← AdditivePhylogeny(D, n - 1)
        v ← the (potentially new) node in T at distance x from i on the path between i and k
        add leaf n back to T by creating a limb (v, n) of length limbLength
        return T
        */
        static Graph<int, int> AdditivePhylogeny(int[,] D, int n)
        {
            if (n == 2)
            {
                Graph<int, int> T = new Graph<int, int>();
                T.AddEdge(0, 1, D[0, 1]);
                T.AddEdge(1, 0, D[0, 1]);
                return T;
            }
            int limbLength = LimbLength(D, n, n - 1);
            for (int j = 0; j < n - 1; ++j)
            {
                D[j, n - 1] -= limbLength;
                D[n - 1, j] = D[j, n - 1];
            }
            bool found = false;
            int i, k = 0;
            for (i = 0; i < n - 2; ++i)
            {
                for (k = i + 1; k < n - 1; ++k)
                    if (D[i, k] == D[i, n - 1] + D[n - 1, k]) { found = true; break; }
                if (found) break;
            }
            if (!found) throw new Exception("(i,n,k) not found");
            int x = D[i, n - 1];
            Graph<int, int> tree = AdditivePhylogeny(D, n - 1);
            List<Edge<int, int>> path = tree.Path(i, k);
            double l = 0;
            Node<int, int> v = null;
            if (path.Count == 0) v = tree.Node(D.GetLength(0));
            else
                foreach (Edge<int, int> edge in path)
                {
                    l += edge.weight;
                    if (l == x)
                    {
                        v = edge.target;
                        break;
                    }
                    if (l > x)
                    {
                        int inter = D.GetLength(0);
                        foreach (Node<int, int> node in tree.nodes)
                            if (node.label >= inter)
                                inter = node.label + 1;
                        v = tree.Node(inter);
                        //  link edge.source, edge.target 
                        Node<int, int> left = edge.source;
                        Node<int, int> right = edge.target;

                        foreach (Edge<int, int> tbr in right.outgoing)
                            if (tbr.target == left)
                            {
                                tbr.Remove();
                                break;
                            }

                        new Edge<int, int>(left, v, x - (l - edge.weight));
                        new Edge<int, int>(v, left, x - (l - edge.weight));
                        new Edge<int, int>(right, v, l - x);
                        new Edge<int, int>(v, right, l - x);

                        edge.Remove();

                        break;
                    }
                }

            new Edge<int, int>(tree.Node(n - 1), v, limbLength);
            new Edge<int, int>(v, tree.Node(n - 1), limbLength);
            return tree;
        }
        static void SolveAdditivePhylogeny()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            int[,] D = new int[n, n];
            for (int i = 0; i < n; ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < n; ++k)
                        D[i, k] = Int32.Parse(words[k]);
                }
            }

            //  compute

            Graph<int, int> T = AdditivePhylogeny(D, n);
            List<Edge<int, int>> edges = T.Edges();
            edges.Sort((a, b) => (a.source.label.CompareTo(b.source.label)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            foreach (Edge<int, int> edge in edges)
                output.WriteLine("{0}->{1}:{2}", edge.source.label, edge.target.label, edge.weight);
            output.Flush();
        }
        //**************************************************************************
        //  Scoring matrix : L: label type, V: value type
        class ScoringMatrix
        {
            public void Load(string fileName)
            {
                StreamReader input = new StreamReader(fileName);
                string header = input.ReadLine();
                char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                string[] headLines = header.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                int n = headLines.Length;
                labels = new Dictionary<char, int>();
                for (int i = 0; i < n; ++i)
                    labels.Add(headLines[i][0], i);

                weights = new int[n, n];
                for (int i = 0; i < n; ++i)
                {
                    string line = input.ReadLine();

                    if (line.Length > 0)
                    {
                        string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                        for (int k = 0; k < n; ++k)
                            weights[i, k] = Int32.Parse(words[k + 1]);
                    }
                }
            }
            public int Weight(char a, char b)
            {
                return weights[labels[a], labels[b]];
            }
            public Dictionary<char, int> labels;
            public int[,] weights;
        }

        // Scoring matrix
        //**************************************************************************

        //  Edit distance
        static int EditDistance(string v, string w)
        {
            int[,] s = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i) s[i, 0] = i;
            for (int j = 0; j <= w.Length; ++j) s[0, j] = j;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    int down = s[i - 1, j] + 1;
                    int right = s[i, j - 1] + 1;
                    int diag = s[i - 1, j - 1] + (v[i - 1] == w[j - 1] ? 0 : 1);

                    s[i, j] = Math.Min(down, Math.Min(right, diag));
                }
            return s[v.Length, w.Length];
        }
        static void SolveEditDistance()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            //  compute

            int distance = EditDistance(v, w);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");


            result.Write(distance.ToString());
            result.Flush();
        }
        //  Fitting alignment
        static void SolveFittingAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            //  Position, length and score of best global alignment
            int max_p = 0, max_l = 0, max_s = 0;
            int[,] max_backtrack = null;
            for (int l = 1; l < v.Length; ++l)
                for (int p = 0; p < v.Length - l; ++p)
                {
                    int s = 0;
                    int[,] backtrack = LCSBacktrack(v.Substring(p, l), w, ref s);
                    if (s > max_s)
                    {
                        max_s = s;
                        max_p = p;
                        max_l = l;
                        max_backtrack = backtrack;
                    }
                }

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            OutputLCS(max_backtrack, v.Substring(max_p, max_l), w, max_l, w.Length, ref ov, ref ow);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(max_s.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }
        //  Fitting alignment
        static void SolveOverlapAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            //  Position, length and score of best global alignment
            int max_lv = 0, max_lw = 0, max_s = Int32.MinValue;
            int[,] max_backtrack = null;
            for (int lv = (v.Length / 3) * 2; lv < v.Length - 1; ++lv)
                for (int lw = 1; lw < w.Length * 0.33; ++lw)
                {
                    int s = 0;
                    int[,] backtrack = LCSBacktrack(v.Substring(lv), w.Substring(0, lw), ref s);
                    if (s > max_s)
                    {
                        max_s = s;
                        max_lv = lv;
                        max_lw = lw;
                        max_backtrack = backtrack;
                    }
                }

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            OutputLCS(max_backtrack, v.Substring(max_lv), w.Substring(0, max_lw), max_backtrack.GetLength(0) - 1, max_backtrack.GetLength(1) - 1, ref ov, ref ow);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(max_s.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }
        /*
            UPGMA(D, n) 
            Clusters ← n single-element clusters labeled 1, ... , n 

             construct a graph T with n isolated nodes labeled by single elements 1, ... , n 
             for every node v in T 
              Age(v) ← 0
             while there is more than one cluster 
              find the two closest clusters Ci and Cj
 
              merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
              add a new node labeled by cluster Cnew to T
              connect node Cnew to Ci and Cj by directed edges

              remove the rows and columns of D corresponding to Ci and Cj

              remove Ci and Cj from Clusters

              add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters 
              add Cnew to Clusters 
             root ← the node in T corresponding to the remaining cluster
             for each edge (v, w) in T
              length of (v, w) ← Age(v) - Age(w)
             return T
         */
        static Graph<int, int> UPGMA(int[,] D_in, int n)
        {
            Graph<int, int> tree = new Graph<int, int>();
            Dictionary<int, double> age = new Dictionary<int, double>();
            List<List<int>> clusters = new List<List<int>>();
            List<int> clusterToNode = new List<int>();
            for (int i = 0; i < n; ++i)
            {
                tree.nodes.Add(new Node<int, int>(i));
                age.Add(i, 0);
                clusters.Add(new List<int>());
                clusters.Last().Add(i);
                clusterToNode.Add(i);
            }

            //       double [,] D = new double[n,n];
            //     for(int i=0;i<n;++i) for(int j=0;j<n;++j) D[i,j] = D_in[i,j];

            int cnew = n - 1;

            while (clusters.Count > 1)
            {
                int closest_ci = 0, closest_cj = 0;
                double closest_distance = double.MaxValue;
                for (int ci = 0; ci < clusters.Count - 1; ++ci)
                {
                    for (int cj = ci + 1; cj < clusters.Count; ++cj)
                    {
                        double distance = 0;
                        foreach (int c1 in clusters[ci])
                            foreach (int c2 in clusters[cj])
                                distance += D_in[c1, c2];
                        distance /= (clusters[ci].Count * clusters[cj].Count);
                        if (distance < closest_distance)
                        {
                            closest_distance = distance;
                            closest_ci = ci;
                            closest_cj = cj;
                        }
                    }
                }

                List<int> cluster_new = new List<int>();
                cluster_new.AddRange(clusters[closest_ci]);
                cluster_new.AddRange(clusters[closest_cj]);

                ++cnew;

                age.Add(cnew, closest_distance / 2);

                Node<int, int> node_ci = tree.nodes[clusterToNode[closest_ci]];
                Node<int, int> node_cj = tree.nodes[clusterToNode[closest_cj]];
                Node<int, int> node_cnew = new Node<int, int>(cnew);
                tree.nodes.Add(node_cnew);

                new Edge<int, int>(node_ci, node_cnew, closest_distance);
                new Edge<int, int>(node_cnew, node_ci, closest_distance);

                new Edge<int, int>(node_cj, node_cnew, closest_distance);
                new Edge<int, int>(node_cnew, node_cj, closest_distance);

                clusters.RemoveAt(closest_cj);
                clusters.RemoveAt(closest_ci);
                clusters.Add(cluster_new);

                clusterToNode.RemoveAt(closest_cj);
                clusterToNode.RemoveAt(closest_ci);
                clusterToNode.Add(cnew);
            }

            foreach (Edge<int, int> edge in tree.Edges())
            {
                edge.weight = Math.Abs(age[edge.source.label] - age[edge.target.label]);
            }

            return tree;
        }
        static void SolveUPGMA()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            int[,] D = new int[n, n];
            for (int i = 0; i < n; ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < n; ++k)
                        D[i, k] = Int32.Parse(words[k]);
                }
            }

            //  compute

            Graph<int, int> T = UPGMA(D, n);
            List<Edge<int, int>> edges = T.Edges();
            edges.Sort((a, b) => (a.source.label.CompareTo(b.source.label)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            foreach (Edge<int, int> edge in edges)
                output.WriteLine("{0}->{1}:{2}", edge.source.label, edge.target.label, edge.weight);
            output.Flush();
        }
        /// <summary>
        /// ////////////////////////////////////////////////////////////////////////////////////////////
        /*
        UPGMA(D, n) 
            Clusters ← n single-element clusters labeled 1, ... , n 

             construct a graph T with n isolated nodes labeled by single elements 1, ... , n 
             for every node v in T 
              Age(v) ← 0
             while there is more than one cluster 
              find the two closest clusters Ci and Cj
 
              merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
              add a new node labeled by cluster Cnew to T
              connect node Cnew to Ci and Cj by directed edges

              remove the rows and columns of D corresponding to Ci and Cj

              remove Ci and Cj from Clusters

              add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters 
              add Cnew to Clusters 
             root ← the node in T corresponding to the remaining cluster
             for each edge (v, w) in T
              length of (v, w) ← Age(v) - Age(w)
             return T
        HierarchicalClustering(D, n)
          Clusters ← n single-element clusters labeled 1, ... , n 
          construct a graph T with n isolated nodes labeled by single elements 1, ... , n 
         while there is more than one cluster 
          find the two closest clusters Ci and Cj  
          merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
          add a new node labeled by cluster Cnew to T
          connect node Cnew to Ci and Cj by directed edges 
          remove the rows and columns of D corresponding to Ci and Cj 
          remove Ci and Cj from Clusters 
          add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters 
          add Cnew to Clusters 
            assign root in T as a node with no incoming edges
            return T
        */
        static List<List<int>> UPGMA2(double[,] D_in, int n)
        {
            Graph<int, int> tree = new Graph<int, int>();
            List<List<int>> clusters = new List<List<int>>();
            List<List<int>> legend = new List<List<int>>();
            List<int> clusterToNode = new List<int>();
            for (int i = 0; i < n; ++i)
            {
                tree.nodes.Add(new Node<int, int>(i));
                clusters.Add(new List<int>());
                clusters.Last().Add(i);
                clusterToNode.Add(i);
            }

            //       double [,] D = new double[n,n];
            //     for(int i=0;i<n;++i) for(int j=0;j<n;++j) D[i,j] = D_in[i,j];

            int cnew = n - 1;

            while (clusters.Count > 1)
            {
                int closest_ci = 0, closest_cj = 0;
                double closest_distance = double.MaxValue;
                for (int ci = 0; ci < clusters.Count - 1; ++ci)
                {
                    for (int cj = ci + 1; cj < clusters.Count; ++cj)
                    {
                        double distance = 0;
                        foreach (int c1 in clusters[ci])
                            foreach (int c2 in clusters[cj])
                                distance += D_in[c1, c2];
                        distance /= (clusters[ci].Count * clusters[cj].Count);
                        if (distance < closest_distance)
                        {
                            closest_distance = distance;
                            closest_ci = ci;
                            closest_cj = cj;
                        }
                    }
                }

                List<int> cluster_new = new List<int>();
                cluster_new.AddRange(clusters[closest_ci]);
                cluster_new.AddRange(clusters[closest_cj]);

                ++cnew;

                Node<int, int> node_ci = tree.nodes[clusterToNode[closest_ci]];
                Node<int, int> node_cj = tree.nodes[clusterToNode[closest_cj]];
                Node<int, int> node_cnew = new Node<int, int>(cnew);
                tree.nodes.Add(node_cnew);

                new Edge<int, int>(node_ci, node_cnew, closest_distance);
                new Edge<int, int>(node_cnew, node_ci, closest_distance);

                new Edge<int, int>(node_cj, node_cnew, closest_distance);
                new Edge<int, int>(node_cnew, node_cj, closest_distance);

                clusters.RemoveAt(closest_cj);
                clusters.RemoveAt(closest_ci);
                clusters.Add(cluster_new);
                legend.Add(cluster_new);

                clusterToNode.RemoveAt(closest_cj);
                clusterToNode.RemoveAt(closest_ci);
                clusterToNode.Add(cnew);
            }

            return legend;
        }
        static void SolveUPGMA2()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            double[,] D = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < n; ++k)
                        D[i, k] = Double.Parse(words[k]);
                }
            }

            //  compute

            List<List<int>> legend = UPGMA2(D, n);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            foreach (List<int> step in legend)
            {
                foreach (int e in step)
                    output.Write("{0} ", e+1);
                output.WriteLine();
            }
            output.Flush();
        }
        
        /// /////////////////////////////////////////////////////////////////////////////////////////////////////////
      
        /*
        NeighborJoining(D,n)
         if n = 2
          T ← tree consisting of a single edge of length D1,2
          return T
         D* ← neighbor-joining matrix constructed from the distance matrix D
         find elements i and j such that D*i,j is a minimum non-diagonal element of D*
         Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
         limbLengthi ← (1/2)(Di,j + Δ)
         limbLengthj ← (1/2)(Di,j - Δ)
         add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
         remove rows i and j from D
         remove columns i and j from D
         T ← NeighborJoining(D, n - 1)
         add two new limbs (connecting node m with leaves i and j) to the tree T
         assign length limbLengthi to Limb(i)
         assign length limbLengthj to Limb(j)
         return T
        */
        static Graph<int, int> NeighborJoining(double[,] D, int n, ref int ma, List<int> index)
        {
            if (n == 2)
            {
                Graph<int, int> T = new Graph<int, int>(); T.AddEdge(index[0], index[1], D[0, 1]); T.AddEdge(index[1], index[0], D[0, 1]);
                return T;
            }
            double min_nj = double.MaxValue, delta = 0;
            int min_i = 0, min_j = 0;
            double[,] NJ = new double[n, n];
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    if (i == j) { NJ[i, j] = 0; continue; }
                    double total_i = 0; for (int ti = 0; ti < n; ++ti) total_i += D[i, ti];
                    double total_j = 0; for (int tj = 0; tj < n; ++tj) total_j += D[j, tj];
                    NJ[i, j] = (n - 2) * D[i, j] - total_i - total_j;
                    if (NJ[i, j] < min_nj) { min_nj = NJ[i, j]; min_i = i; min_j = j; delta = (total_i - total_j) / (n - 2); }
                }
            double limbLengthi = (D[min_i, min_j] + delta) / 2, limbLengthj = (D[min_i, min_j] - delta) / 2;
            double[,] D_new = new double[n - 1, n - 1]; int m = n - 2;
            for (int k = 0, rk = 0; k < n; ++k)
            {
                if (k == min_i || k == min_j) continue;
                D_new[rk, m] = D_new[m, rk] = (D[k, min_i] + D[k, min_j] - D[min_i, min_j]) / 2;
                ++rk;
            }
            for (int i = 0, ri = 0; i < n; ++i)
            {
                if (i == min_i || i == min_j) continue;
                for (int j = 0, rj = 0; j < n; ++j)
                {
                    if (j == min_i || j == min_j) continue;
                    D_new[ri, rj] = D[i, j];
                    ++rj;
                }
                ++ri;
            }
            List<int> index_new = new List<int>();
            index_new.AddRange(index);

            if (min_i < min_j)
            {
                index_new.RemoveAt(min_j);
                index_new.RemoveAt(min_i);
            }
            else
            {
                index_new.RemoveAt(min_i);
                index_new.RemoveAt(min_j);
            }
            index_new.Add(ma++);
            Graph<int, int> tree = NeighborJoining(D_new, n - 1, ref ma, index_new);
            Node<int, int> node_m = tree.Node(index_new.Last());
            Node<int, int> node_i = tree.Node(index[min_i]);
            Node<int, int> node_j = tree.Node(index[min_j]);
            new Edge<int, int>(node_m, node_i, limbLengthi);
            new Edge<int, int>(node_i, node_m, limbLengthi);
            new Edge<int, int>(node_m, node_j, limbLengthj);
            new Edge<int, int>(node_j, node_m, limbLengthj);
            return tree;
        }
        static void SolveNeighborJoining()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            double[,] D = new double[n, n];
            for (int i = 0; i < n; ++i)
            {
                string line = input.ReadLine();

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                    for (int k = 0; k < n; ++k)
                        D[i, k] = Int32.Parse(words[k]);
                }
            }

            //  compute
            int ma = n;
            List<int> index = new List<int>();
            for (int i = 0; i < n; ++i) index.Add(i);
            Graph<int, int> T = NeighborJoining(D, n, ref ma, index);
            List<Edge<int, int>> edges = T.Edges();
            edges.Sort((a, b) => (a.source.label.CompareTo(b.source.label)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            foreach (Edge<int, int> edge in edges)
                output.WriteLine("{0}->{1}:{2}", edge.source.label, edge.target.label, edge.weight);
            output.Flush();
        }
        //************************************************************
        //  Affine Gap Penalties
        static int[,] AGPBacktrack(string v, string w, ScoringMatrix m, ref int score)
        {
            int epsilon = 1;
            int sigma = 11;

            int[,] lower = new int[v.Length + 1, w.Length + 1];
            int[,] middle = new int[v.Length + 1, w.Length + 1];
            int[,] upper = new int[v.Length + 1, w.Length + 1];
            int[,] backtrack = new int[v.Length + 1, w.Length + 1];

            for (int i = 0; i <= v.Length; ++i)
            {
                lower[i, 0] = i * -epsilon;
                middle[i, 0] = 0;
                upper[i, 0] = 0;
            }

            for (int j = 0; j <= w.Length; ++j)
            {
                lower[0, j] = 0;
                middle[0, j] = 0;
                upper[0, j] = j * -epsilon;
            }

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                {
                    lower[i, j] = Math.Max(lower[i - 1, j] - epsilon, middle[i - 1, j] - sigma);
                    upper[i, j] = Math.Max(upper[i, j - 1] - epsilon, middle[i, j - 1] - sigma);
                    middle[i, j] = Math.Max(Math.Max(lower[i, j], upper[i, j]),
                                            middle[i - 1, j - 1] + m.Weight(v[i - 1], w[j - 1]));

                    if (middle[i, j] == lower[i, j]) backtrack[i, j] = -1;
                    else if (middle[i, j] == upper[i, j]) backtrack[i, j] = 1;
                    else backtrack[i, j] = 0;
                }

            score = middle[v.Length, w.Length];
            return backtrack;
        }

        static void OutputAGP(int[,] backtrack, string v, string w, int i, int j, ref List<char> ov, ref List<char> ow)
        {
            while (i != 0 && j != 0)
            {
                if (backtrack[i, j] == 0)
                {
                    ov.Add(v[i - 1]);
                    ow.Add(w[j - 1]);
                    --i; --j;
                }
                else if (backtrack[i, j] == 1)
                {
                    ov.Add('-');
                    ow.Add(w[j - 1]);
                    --j;
                }
                else if (backtrack[i, j] == -1)
                {
                    ov.Add(v[i - 1]);
                    ow.Add('-');
                    --i;
                }
                else if (backtrack[i, j] == -2)
                {
                    i = j = 0;
                }
            }

            while (i > 0)
            {
                ov.Add(v[i - 1]);
                ow.Add('-');
                i--;
            }

            while (j > 0)
            {
                ow.Add(w[j - 1]);
                ov.Add('-');
                i--;
            }

            ov.Reverse();
            ow.Reverse();
        }
        static void SolveAGP()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            ScoringMatrix m = new ScoringMatrix();
            m.Load("BLOSUM62.txt");

            //  compute
            int score = 0;
            int[,] backtrack = AGPBacktrack(v, w, m, ref score);

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            OutputAGP(backtrack, v, w, v.Length, w.Length, ref ov, ref ow);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(score.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }
        //**************************************************************************************
        //  Efficient alignment
        static void SolveMiddleEdge()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            ScoringMatrix m = new ScoringMatrix();
            m.Load("BLOSUM62.txt");

            //  compute
            int[,] edge = MiddleEdge(v, w, m, 5);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine("({0}, {1}) ({2}, {3})", edge[0, 0], edge[0, 1], edge[1, 0], edge[1, 1]);
            result.Flush();
        }
        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }
        static void AlignmentLastColumn(string v, string w, ScoringMatrix m, int indel, out int[] backtrack, out int[] b)
        {
            b = new int[v.Length + 1];
            int[] a = new int[v.Length + 1];

            backtrack = new int[v.Length + 1];

            for (int i = 0; i <= v.Length; ++i) a[i] = -i * indel;

            for (int j = 1; j <= w.Length; ++j)
            {
                b[0] = (j) * -indel;
                backtrack[0] = 1;
                for (int i = 1; i <= v.Length; ++i)
                {
                    int down = b[i - 1] - indel;
                    int right = a[i] - indel;
                    int diag = a[i - 1] + m.Weight(v[i - 1], w[j - 1]);

                    b[i] = Math.Max(down, Math.Max(right, diag));

                    if (b[i] == down) backtrack[i] = -1;
                    else if (b[i] == right) backtrack[i] = 1;
                    else if (b[i] == diag) backtrack[i] = 0;
                }
                if (j < w.Length)
                {
                    int[] tmp = a; a = b; b = tmp;
                }
            }
        }
        static int[,] MiddleEdge(string v, string w, ScoringMatrix m, int indel)
        {
            if (v.Length == 1 && w.Length == 1)
            {
                return new int[2, 2] { { 0, 0 }, { 1, 1 } };
            }

            int middle = w.Length / 2;

            int[] lastColumnLeft;
            int[] backtrackLeft;
            AlignmentLastColumn(v, w.Substring(0, middle), m, indel, out backtrackLeft, out lastColumnLeft);

            int[] lastColumnRight;
            int[] backtrackRight;
            AlignmentLastColumn(Reverse(v), Reverse(w.Substring(middle)), m, indel, out backtrackRight, out lastColumnRight);
            Array.Reverse(lastColumnRight);
            Array.Reverse(backtrackRight);

            int max = lastColumnLeft[0] + lastColumnRight[0];
            int maxPos = 0;
            for (int i = 1; i <= v.Length; ++i)
            {
                if (lastColumnLeft[i] + lastColumnRight[i] > max)
                {
                    max = lastColumnLeft[i] + lastColumnRight[i];
                    maxPos = i;
                }
            }

            int[,] edge = new int[2, 2];
            edge[0, 1] = middle;
            edge[0, 0] = maxPos;
            edge[1, 1] = middle + (backtrackRight[maxPos] != -1 ? 1 : 0);
            edge[1, 0] = maxPos + (backtrackRight[maxPos] != 1 ? 1 : 0);
            return edge;
        }
        /*
        LinearSpaceAlignment(top, bottom, left, right)
        if left = right
            return alignment formed by bottom − top vertical edges
        if top = bottom
            return alignment formed by right − left horizontal edges
        middle ← ⌊ (left + right)/2⌋
        midNode ← MiddleNode(top, bottom, left, right)
        midEdge ← MiddleEdge(top, bottom, left, right)
        LinearSpaceAlignment(top, midNode, left, middle)
        output midEdge
        if midEdge = "→" or midEdge = "↘"
            middle ← middle + 1
        if midEdge = "↓" or midEdge = "↘"
            midNode ← midNode + 1 
        LinearSpaceAlignment(midNode, bottom, middle, right)
        */
        static void LinearSpaceAlignment(int top, int bottom, int left, int right,
                                            string v, string w, ScoringMatrix m, int indel,
                                            ref List<char> ov, ref List<char> ow)
        {
            if (left == right)
            {
                for (int i = top; i < bottom; ++i)
                {
                    ov.Add(v[i]);
                    ow.Add('-');
                }
                return;
                //  return alignment formed by bottom − top vertical edges
            }
            if (top == bottom)
            {
                for (int i = left; i < right; ++i)
                {
                    ov.Add('-');
                    ow.Add(w[i]);
                }
                return;
                //  return alignment formed by right − left horizontal edges
            }

            int[,] edge = MiddleEdge(v.Substring(top, bottom - top), w.Substring(left, right - left), m, indel);
            edge[0, 0] += top;
            edge[0, 1] += left;
            edge[1, 0] += top;
            edge[1, 1] += left;

            LinearSpaceAlignment(top, edge[0, 0], left, edge[0, 1],
                                    v, w, m, indel,
                                    ref ov, ref ow);

            bool d = (edge[0, 0] != edge[1, 0]);
            bool r = (edge[0, 1] != edge[1, 1]);

            if (d && r)  //  "↘"
            {
                ov.Add(v[edge[0, 0]]);
                ow.Add(w[edge[0, 1]]);
            }
            else if (d)  //  "↓"
            {
                ov.Add(v[edge[0, 0]]);
                ow.Add('-');
            }
            else if (r)  //  "→"
            {
                ov.Add('-');
                ow.Add(w[edge[0, 1]]);
            }

            LinearSpaceAlignment(edge[1, 0], bottom, edge[1, 1], right,
                                    v, w, m, indel,
                                    ref ov, ref ow);
        }
        static void SolveLinearSpaceAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();

            ScoringMatrix m = new ScoringMatrix();
            m.Load("BLOSUM62.txt");

            //  compute

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            LinearSpaceAlignment(0, v.Length, 0, w.Length, v, w, m, 5, ref ov, ref ow);
            int score = 0;
            for (int i = 0; i < ov.Count; ++i)
            {
                if (ov[i] == '-' || ow[i] == '-') score -= 5;
                else score += m.Weight(ov[i], ow[i]);
            }

            //  write result
            StreamWriter result = new StreamWriter("result.txt");


            result.WriteLine(score.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.Flush();
        }
        //***********************************************************************
        const int ip = 1;
        const int jp = 2;
        const int kp = 4;
        static int[,,] MultipleLCSBacktrack(string v, string w, string u, ref int score)
        {
            int[,,] s = new int[v.Length + 1, w.Length + 1, u.Length + 1];
            int[,,] backtrack = new int[v.Length + 1, w.Length + 1, u.Length + 1];

            for (int i = 0; i <= v.Length; ++i)
                for (int j = 0; j <= w.Length; ++j)
                    for (int k = 0; k <= u.Length; ++k)
                        s[i, j, k] = 0;

            for (int i = 1; i <= v.Length; ++i)
                for (int j = 1; j <= w.Length; ++j)
                    for (int k = 1; k <= u.Length; ++k)
                    {
                        int diag = s[i - 1, j - 1, k - 1] + (((v[i - 1] == w[j - 1]) && (v[i - 1] == u[k - 1])) ? 1 : 0);

                        s[i, j, k] = Math.Max(diag,
                                        Math.Max(s[i - 1, j - 1, k - 0],
                                        Math.Max(s[i - 1, j - 0, k - 1],
                                        Math.Max(s[i - 0, j - 1, k - 1],
                                        Math.Max(s[i - 1, j - 0, k - 0],
                                        Math.Max(s[i - 0, j - 1, k - 0],
                                                    s[i - 0, j - 0, k - 1]
                                        ))))));

                        if (s[i, j, k] == s[i - 1, j - 1, k - 0]) backtrack[i, j, k] = ip + jp;
                        if (s[i, j, k] == s[i - 1, j - 0, k - 1]) backtrack[i, j, k] = ip + kp;
                        if (s[i, j, k] == s[i - 0, j - 1, k - 1]) backtrack[i, j, k] = jp + kp;
                        if (s[i, j, k] == s[i - 1, j - 0, k - 0]) backtrack[i, j, k] = ip;
                        if (s[i, j, k] == s[i - 0, j - 1, k - 0]) backtrack[i, j, k] = jp;
                        if (s[i, j, k] == s[i - 0, j - 0, k - 1]) backtrack[i, j, k] = kp;
                        if (s[i, j, k] == diag) backtrack[i, j, k] = ip + jp + kp;

                    }

            score = s[v.Length, w.Length, u.Length];

            return backtrack;
        }
        static void OutputMultipleLCS(int[,,] backtrack, string v, string w, string u, int i, int j, int k, ref List<char> ov, ref List<char> ow, ref List<char> ou)
        {
            while (i != 0 && j != 0 && k != 0)
            {
                BitArray b = new BitArray(new int[] { backtrack[i, j, k] });
                bool[] bits = new bool[b.Count];
                b.CopyTo(bits, 0);

                if (b[0])
                {
                    ov.Add(v[i - 1]);
                    --i;
                }
                else
                {
                    ov.Add('-');
                }

                if (b[1])
                {
                    ow.Add(w[j - 1]);
                    --j;
                }
                else
                {
                    ow.Add('-');
                }

                if (b[2])
                {
                    ou.Add(u[k - 1]);
                    --k;
                }
                else
                {
                    ou.Add('-');
                }
            }

            while (i > 0) ov.Add(v[--i]);
            while (j > 0) ow.Add(w[--j]);
            while (k > 0) ou.Add(u[--k]);

            int max = Math.Max(ou.Count, Math.Max(ow.Count, ov.Count));

            while (ov.Count < max) ov.Add('-');
            while (ow.Count < max) ow.Add('-');
            while (ou.Count < max) ou.Add('-');

            ov.Reverse();
            ow.Reverse();
            ou.Reverse();
        }
        static void SolveMultipleLCS()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            string v = input.ReadLine();
            string w = input.ReadLine();
            string u = input.ReadLine();

            //  compute
            int score = 0;
            int[,,] backtrack = MultipleLCSBacktrack(v, w, u, ref score);

            List<char> ov = new List<char>();
            List<char> ow = new List<char>();
            List<char> ou = new List<char>();
            OutputMultipleLCS(backtrack, v, w, u, v.Length, w.Length, u.Length, ref ov, ref ow, ref ou);

            //  write result
            StreamWriter result = new StreamWriter("result.txt");
            result.WriteLine(score.ToString());
            foreach (char c in ov) result.Write(c);
            result.WriteLine();
            foreach (char c in ow) result.Write(c);
            result.WriteLine();
            foreach (char c in ou) result.Write(c);
            result.Flush();
        }
        //***********************************************************************************************************
        //  Small parsimony
        /*
            SmallParsimony(T, Character)
             for each node v in tree T
                    Tag(v) ← 0
              if v is a leaf
               Tag(v) ← 1
               for each symbol k in the alphabet
                if Character(v) = k
                 sk(v) ← 0
                else
                 sk(v) ← ∞
             while there exist ripe nodes in T
              v ← a ripe node in T
              Tag(v) ← 1
              for each symbol k in the alphabet
                  sk(v) ← minimum over all symbols i {si(Daughter(v))+δi,k} + minimum over all symbols j {sj(Son(v))+δj,k}
               return minimum over all symbols k {sk(v)}
        */
        static void SmallParsimony(Graph<string, ParsimonyLoad> tree)
        {
            foreach (Node<string, ParsimonyLoad> node in tree.nodes)
            {
                if (node.load != null)    //  leaf
                {
                    for (int i = 0; i < node.load.letter.Length; ++i)
                    {
                        if (node.load.label[0] == node.load.letter[i].nucleo) node.load.letter[i].score = 0;
                        else node.load.letter[i].score = Int32.MaxValue / 2;
                    }
                }
            }

            Node<string, ParsimonyLoad> v = null;
            while (true)
            {
                int v_index = 0;
                for (v_index = 0; v_index < tree.nodes.Count; ++v_index)
                    if (tree.nodes[v_index].load == null)   //  not a leaf
                        if (tree.nodes[v_index].outgoing[0].target.load != null &&
                            tree.nodes[v_index].outgoing[1].target.load != null)
                            break;
                if (v_index == tree.nodes.Count) break;
                v = tree.nodes[v_index];

                v.load = new ParsimonyLoad();

                // for each symbol k in the alphabet
                //  sk(v) ← minimum over all symbols i { si(Daughter(v)) + δi,k}
                //  +minimum over all symbols j { sj(Son(v)) + δj,k}
                for (int k = 0; k < v.load.letter.Length; ++k)
                {
                    int min_left = Int32.MaxValue;
                    foreach (ParsimonyNucleotideInfo left in v.outgoing[0].target.load.letter)
                    {
                        int s = left.score + (left.nucleo == v.load.letter[k].nucleo ? 0 : 1);
                        if (s < min_left)
                        {
                            min_left = s;
                            v.load.letter[k].left = left.nucleo;
                        }
                    }
                    int min_right = Int32.MaxValue;
                    foreach (ParsimonyNucleotideInfo right in v.outgoing[1].target.load.letter)
                    {
                        int s = right.score + (right.nucleo == v.load.letter[k].nucleo ? 0 : 1);
                        if (s < min_right)
                        {
                            min_right = s;
                            v.load.letter[k].right = right.nucleo;
                        }
                    }
                    v.load.letter[k].score = min_left + min_right;
                }
            }

            //  construct the best tree

            int min_score = Int32.MaxValue;
            Node<string, ParsimonyLoad> root = null;
            for (int i = 0; i < tree.nodes.Count; ++i) if (tree.nodes[i].incoming.Count == 0) { root = tree.nodes[i]; break; }
            foreach (ParsimonyNucleotideInfo s in root.load.letter)
            {
                if (s.score < min_score)
                {
                    min_score = s.score;
                    root.load.label = s.nucleo.ToString();
                }
            }

            AssignLetters(root);
        }
        static void AssignLetters(Node<string, ParsimonyLoad> node)
        {
            if (node.outgoing.Count == 0) return;

            foreach (ParsimonyNucleotideInfo info in node.load.letter)
            {
                if (info.nucleo == node.load.label[0])
                {
                    node.outgoing[0].target.load.label = info.left.ToString();
                    node.outgoing[1].target.load.label = info.right.ToString();
                    AssignLetters(node.outgoing[0].target);
                    AssignLetters(node.outgoing[1].target);
                    return;
                }
            }
        }
        struct ParsimonyNucleotideInfo
        {
            public char nucleo;
            public char left;
            public char right;
            public int score;
            public ParsimonyNucleotideInfo(char _nucleo, char _left, char _right, int _score)
            {
                nucleo = _nucleo;
                left = _left;
                right = _right;
                score = _score;
            }
        }
        class ParsimonyLoad
        {
            public ParsimonyNucleotideInfo[] letter;
            public string label;
            public bool tag;
            public ParsimonyLoad()
            {
                letter = new ParsimonyNucleotideInfo[4] {   new ParsimonyNucleotideInfo('A', ' ', ' ', Int32.MaxValue ),
                                                            new ParsimonyNucleotideInfo('C', ' ', ' ', Int32.MaxValue ),
                                                            new ParsimonyNucleotideInfo('G', ' ', ' ', Int32.MaxValue ),
                                                            new ParsimonyNucleotideInfo('T', ' ', ' ', Int32.MaxValue )};
                label = "";
                tag = false;
            }
        }
        static void SmallParsimony(Graph<string, string> tree)
        {
            int leafLength = 0;
            foreach (Node<string, string> node in tree.nodes)
            {
                if (node.outgoing.Count == 0)
                    if (node.label.Length > leafLength) leafLength = node.label.Length;
            }

            List<Graph<string, ParsimonyLoad>> subTrees = new List<Graph<string, ParsimonyLoad>>();
            for (int i = 0; i < leafLength; ++i) subTrees.Add(new Graph<string, ParsimonyLoad>());
            foreach (Edge<string, string> edge in tree.Edges())
            {
                for (int i = 0; i < leafLength; ++i)
                {
                    subTrees[i].AddEdge(edge.source.label, edge.target.label, 0);
                    if (edge.target.outgoing.Count == 0)
                    {
                        subTrees[i].Node(edge.target.label).load = new ParsimonyLoad();
                        subTrees[i].Node(edge.target.label).load.label = edge.target.label[i].ToString();
                    }
                }
            }

            for (int i = 0; i < tree.nodes.Count; ++i)
            {
                tree.nodes[i].load = "";
            }

            foreach (Graph<string, ParsimonyLoad> subTree in subTrees)
            {
                SmallParsimony(subTree);
                for (int i = 0; i < tree.nodes.Count; ++i)
                {
                    tree.nodes[i].load += subTree.Node(tree.nodes[i].label).load.label;
                }
            }

            foreach (Edge<string, string> edge in tree.Edges())
                edge.weight = HammingDistance(edge.source.load, edge.target.load);
        }
        static int HammingDistance(string a, string b)
        {
            int dist = 0;
            for (int i = 0; i < a.Length; ++i)
            {
                if (a[i] != b[i]) ++dist;
            }
            return dist;
        }
        static void SolveSmallParsimony()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            Graph<string, string> graph = new Graph<string, string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    string source = words[0];
                    string target = words[1];

                    graph.AddEdge(source,
                                    target,
                                    0);
                }
            }

            //  compute

            SmallParsimony(graph);

            List<Edge<string, string>> edges = graph.Edges();
            edges.Sort((a, b) => (a.source.load.CompareTo(b.source.load)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(graph.Weight().ToString());
            foreach (Edge<string, string> edge in edges)
                output.WriteLine("{0}->{1}:{2}", edge.source.load, edge.target.load, edge.weight);
            foreach (Edge<string, string> edge in edges)
                output.WriteLine("{1}->{0}:{2}", edge.source.load, edge.target.load, edge.weight);
            output.Flush();
        }
        static Graph<string, string> UnrootedSmallParsimony(Graph<string, string> tree)
        {
            Graph<string, string> bestTree = null;
            double bestScore = Int32.MaxValue;
            foreach (Edge<string, string> edge in tree.Edges())
            {
                string label = edge.source.label + edge.target.label;
                Graph<string, string> rooted = new Graph<string, string>(tree);
                MakeRooted(rooted, edge, label);
                SmallParsimony(rooted);

                Node<string, string> root = rooted.Node(label);
                new Edge<string, string>(root.outgoing[0].target, root.outgoing[1].target,
                                            HammingDistance(root.outgoing[0].target.load, root.outgoing[1].target.load));

                root.outgoing[1].Remove();
                root.outgoing[0].Remove();

                rooted.nodes.Remove(root);

                double score = rooted.Weight();

                if (score < bestScore)
                {
                    bestScore = score;
                    bestTree = rooted;
                }
            }

            return bestTree;
        }
        static void MakeRooted(Graph<string, string> rooted, Edge<string, string> edge, string label)
        {
            Node<string, string> left = rooted.Node(edge.source.label);
            Node<string, string> right = rooted.Node(edge.target.label);

            while (left.incoming.Count > 0) left.incoming[0].Remove();
            while (right.incoming.Count > 0) right.incoming[0].Remove();

            foreach (Edge<string, string> dir in left.outgoing) BuildForward(dir.target);
            foreach (Edge<string, string> dir in right.outgoing) BuildForward(dir.target);

            rooted.AddEdge(label, edge.source.label, 0);
            rooted.AddEdge(label, edge.target.label, 0);
        }
        static void BuildForward(Node<string, string> dir)
        {
            foreach (Edge<string, string> edge in dir.outgoing)
            {
                for (int i = 0; i < edge.target.incoming.Count; ++i)
                {
                    if (edge.target.incoming[i].source != dir) edge.target.incoming[i--].Remove();
                }
                for (int i = 0; i < edge.target.outgoing.Count; ++i)
                {
                    if (edge.target.outgoing[i].target == dir) edge.target.outgoing[i--].Remove();
                }
            }
            foreach (Edge<string, string> edge in dir.outgoing) BuildForward(edge.target);
        }
        static void SolveUnrootedSmallParsimony()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            Graph<string, string> graph = new Graph<string, string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    string source = words[0];
                    string target = words[1];

                    graph.AddEdge(source,
                                    target,
                                    0);
                }
            }

            //  compute

            graph = UnrootedSmallParsimony(graph);

            List<Edge<string, string>> edges = graph.Edges();
            edges.Sort((a, b) => (a.source.load.CompareTo(b.source.load)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(graph.Weight().ToString());
            foreach (Edge<string, string> edge in edges)
                output.WriteLine("{0}->{1}:{2}", edge.source.load, edge.target.load, edge.weight);
            foreach (Edge<string, string> edge in edges)
                output.WriteLine("{1}->{0}:{2}", edge.source.load, edge.target.load, edge.weight);
            output.Flush();
        }
        static void SolveNearestNeighbor()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] nodes = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            string node1 = nodes[0];
            string node2 = nodes[1];

            Graph<string, string> graph = new Graph<string, string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");

                if (line.Length > 0)
                {
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    string source = words[0];
                    string target = words[1];

                    graph.AddEdge(source,
                                    target,
                                    0);
                }
            }

            //  compute
            Graph<string, string> neighbor1, neighbor2;
            NearestNeighbors(graph, node1, node2, out neighbor1, out neighbor2);

            List<Edge<string, string>> edges1 = neighbor2.Edges();
            edges1.Sort((a, b) => (Int32.Parse(a.source.label).CompareTo(Int32.Parse(b.source.label))));

            List<Edge<string, string>> edges2 = neighbor1.Edges();
            edges2.Sort((a, b) => (Int32.Parse(a.source.label).CompareTo(Int32.Parse(b.source.label))));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            foreach (Edge<string, string> edge in edges1)
                output.WriteLine("{0}->{1}", edge.source.label, edge.target.label);
            output.WriteLine();
            foreach (Edge<string, string> edge in edges2)
                output.WriteLine("{0}->{1}", edge.source.label, edge.target.label);
            output.Flush();
        }
        static void NearestNeighbors(Graph<string, string> tree, string node1, string node2,
                                        out Graph<string, string> neighbor1, out Graph<string, string> neighbor2)
        {
            neighbor1 = new Graph<string, string>(tree);
            neighbor2 = new Graph<string, string>(tree);

            int left1 = 0;
            int right1 = 0;
            if (neighbor1.Node(node1).outgoing[0].target.label == node2)
            {
                left1 = 1;
                right1 = 2;
            }
            else if (neighbor1.Node(node1).outgoing[1].target.label == node2)
            {
                left1 = 0;
                right1 = 2;
            }
            else
            {
                left1 = 0;
                right1 = 1;
            }

            int left2 = 0;
            int right2 = 0;
            if (neighbor1.Node(node2).outgoing[0].target.label == node1)
            {
                left2 = 1;
                right2 = 2;
            }
            else if (neighbor1.Node(node2).outgoing[1].target.label == node1)
            {
                left2 = 0;
                right2 = 2;
            }
            else
            {
                left2 = 0;
                right2 = 1;
            }

            Swap(neighbor1.Node(node1), neighbor1.Node(node1).outgoing[right1].target,
                    neighbor1.Node(node2), neighbor1.Node(node2).outgoing[right2].target);


            Swap(neighbor2.Node(node1), neighbor2.Node(node1).outgoing[right1].target,
                    neighbor2.Node(node2), neighbor2.Node(node2).outgoing[left2].target);
        }
        static void Swap(Node<string, string> ra, Node<string, string> a,
                            Node<string, string> rb, Node<string, string> b)
        {
            /*            Node<string, string> tmp = new Node<string, string>();
                        tmp.ShallowCopy(a);
                        a.ShallowCopy(b);
                        b.ShallowCopy(tmp);
              */
            for (int i = 0; i < ra.incoming.Count; ++i)
                if (ra.incoming[i].source == a)
                    ra.incoming[i].Remove();
            for (int i = 0; i < ra.outgoing.Count; ++i)
                if (ra.outgoing[i].target == a)
                    ra.outgoing[i].Remove();

            for (int i = 0; i < rb.incoming.Count; ++i)
                if (rb.incoming[i].source == b)
                    rb.incoming[i].Remove();
            for (int i = 0; i < rb.outgoing.Count; ++i)
                if (rb.outgoing[i].target == b)
                    rb.outgoing[i].Remove();

            new Edge<string, string>(ra, b, 0);
            new Edge<string, string>(b, ra, 0);

            new Edge<string, string>(rb, a, 0);
            new Edge<string, string>(a, rb, 0);
        }
        /*
        NearestNeighborInterchange(Strings)
         score ← ∞
         generate an arbitrary unrooted binary tree Tree with |Strings| leaves
         label the leaves of Tree by arbitrary strings from Strings
         solve the Small Parsimony in an Unrooted Tree Problem for Tree
         label the internal nodes of Tree according to a most parsimonious labeling
         newScore ← the parsimony score of Tree
         newTree ← Tree
         while newScore < score
              score ← newScore
              Tree ← newTree
              for each internal edge e in Tree
                   for each nearest neighbor NeighborTree of Tree with respect to the edge e
                        solve the Small Parsimony in an Unrooted Tree Problem for NeighborTree
                        neighborScore ← the minimum parsimony score of NeighborTree
                        if neighborScore < newScore
                             newScore ← neighborScore
                             newTree ← NeighborTree
         return newTree
        */
        static List<Graph<string, string>> NearestNeighborInterchange(Graph<string, string> tree)
        {
            List<Graph<string, string>> solutions = new List<Graph<string, string>>();

            double score = Double.MaxValue;
            Graph<string, string> ndirected = UnrootedSmallParsimony(tree);
            double newScore = ndirected.Weight();
            Graph<string, string> newTree = tree;

            /*         for (int i=0;i<tree.nodes.Count;++i)
                     {
                         tree.nodes[i].label = tree.nodes[i].load;
                     }
         */
            while (newScore < score)
            {
                score = newScore;
                tree = newTree;
                solutions.Add(ndirected);
                foreach (Edge<string, string> edge in tree.Edges())
                {
                    if (edge.source.outgoing.Count != 3 || edge.target.outgoing.Count != 3) continue;
                    Graph<string, string> neighbor1, neighbor2, ndirected1, ndirected2;
                    NearestNeighbors(tree, edge.source.label, edge.target.label, out neighbor1, out neighbor2);
                    ndirected1 = UnrootedSmallParsimony(neighbor1);
                    ndirected2 = UnrootedSmallParsimony(neighbor2);
                    double neighborScore = ndirected1.Weight();
                    if (neighborScore < newScore)
                    {
                        newScore = neighborScore;
                        ndirected = ndirected1;
                        newTree = neighbor1;
                    }
                    neighborScore = ndirected2.Weight();
                    if (neighborScore < newScore)
                    {
                        newScore = neighborScore;
                        ndirected = ndirected2;
                        newTree = neighbor2;
                    }
                }
            }

            return solutions;
        }
        static void SolveNearestNeighborInterchange()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            int n = Int32.Parse(input.ReadLine());
            Graph<string, string> graph = new Graph<string, string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                line = line.Replace("->", " ");

                if (line.Length > 0)
                {
                    char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                    string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                    string source = words[0];
                    string target = words[1];

                    graph.AddEdge(source,
                                    target,
                                    0);
                }
            }

            //  compute
            List<Graph<string, string>> solutions = NearestNeighborInterchange(graph);

            StreamWriter output = new StreamWriter("result.txt");

            for (int i = 1; i < solutions.Count; ++i)
            {
                List<Edge<string, string>> edges = solutions[i].Edges();
                edges.Sort((a, b) => (a.source.load.CompareTo(b.source.load)));

                //  write result

                output.WriteLine(solutions[i].Weight().ToString());
                foreach (Edge<string, string> edge in edges)
                    output.WriteLine("{0}->{1}:{2}", edge.source.load, edge.target.load, edge.weight);
                foreach (Edge<string, string> edge in edges)
                    output.WriteLine("{1}->{0}:{2}", edge.source.load, edge.target.load, edge.weight);

                if (i < solutions.Count - 1) output.WriteLine();
            }
            output.Flush();
        }
        /*
        GREEDYSORTING(P)
        approxReversalDistance ← 0
        for k = 1 to |P|
            if element k is not sorted
                apply the k-sorting reversal to P
                approxReversalDistance ← approxReversalDistance + 1
            if k-th element of P is −k
                apply the k-sorting reversal to P
                approxReversalDistance ← approxReversalDistance + 1
        return approxReversalDistance
        */
        static List<List<int>> GreedySorting(List<int> P)
        {
            List<List<int>> permutationList = new List<List<int>>();

            for (int k = 0; k < P.Count; ++k)
            {
                bool diff = false;
                for (int st = 0; st <= k; ++st) if (P[st] != (st + 1)) { diff = true; break; }
                if (diff)
                {
                    if (P[k] != -(k + 1))
                    {
                        int j = k;
                        while (P[j] != (k + 1) && P[j] != -(k + 1)) ++j;
                        List<int> PNew = new List<int>();
                        PNew.AddRange(P.GetRange(0, k));
                        for (int i = 0; i <= j - k; ++i) PNew.Add(-P[j - i]);
                        PNew.AddRange(P.GetRange(j + 1, P.Count - j - 1));
                        P = PNew;
                        permutationList.Add(P);
                    }
                }
                if (P[k] == -(k + 1))
                {
                    P = new List<int>(P);
                    P[k] *= -1;
                    permutationList.Add(P);
                }
            }

            return permutationList;
        }
        static void SolveGreedySorting()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();
            line = line.Replace("+", "");
            line = line.Replace("(", "");
            line = line.Replace(")", "");

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] permutation_string = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> permutation = new List<int>();
            foreach (string item in permutation_string)
                permutation.Add(Int32.Parse(item));

            //  compute

            List<List<int>> permutation_list = GreedySorting(permutation);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            for (int j = 0; j < permutation_list.Count; ++j)
            {
                output.Write("(");
                for (int i = 0; i < permutation_list[j].Count; ++i)
                {
                    if (permutation_list[j][i] > 0) output.Write("+");
                    output.Write("{0}", permutation_list[j][i]);
                    if (i != permutation_list[j].Count - 1) output.Write(" ");
                }
                output.Write(")");
                if (j != (permutation_list.Count - 1)) output.Write("\n");
            }

            output.Flush();
        }
        static int PermutationBreakpoints(List<int> permutation)
        {
            permutation.Insert(0, 0);
            permutation.Add(permutation.Count);

            int bp = 0;

            for (int i = 0; i < permutation.Count - 1; ++i)
            {
                if ((permutation[i + 1] - permutation[i]) != 1) ++bp;
            }
            return bp;
        }
        static void SolvePermutationBreakpoints()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();
            line = line.Replace("+", "");
            line = line.Replace("(", "");
            line = line.Replace(")", "");

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] permutation_string = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> permutation = new List<int>();
            foreach (string item in permutation_string)
                permutation.Add(Int32.Parse(item));

            //  compute

            int bp = PermutationBreakpoints(permutation);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            output.WriteLine(bp.ToString());

            output.Flush();
        }

        static char MassToAminoAcid(int mass)
        {
            for (int i = 0; i < AminoAcidMass.Length; ++i) if (AminoAcidMass[i] == mass) return AminoAcid[i];
            return '#';
        }

        static int AminoAcidToMass(char acid)
        {
            for (int i = 0; i < AminoAcid.Length; ++i) if (AminoAcid[i] == acid) return AminoAcidMass[i];
            return 0;
        }

        static Graph<int, int> BuildSpectrumGraph(List<int> Spectrum)
        {
            if (Spectrum[0] != 0) Spectrum.Insert(0, 0);
            Graph<int, int> graph = new Graph<int, int>();

            for (int i = 0; i < Spectrum.Count - 1; ++i)
                for (int j = i; j < Spectrum.Count; ++j)
                {
                    char aminoAcid = MassToAminoAcid(Spectrum[j] - Spectrum[i]);
                    if (aminoAcid != '#')
                    {
                        graph.AddEdge(Spectrum[i], Spectrum[j], 0, aminoAcid.ToString());
                    }
                }
            return graph;
        }

        static void SolveBuildSpectrumGraph()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> Spectrum = new List<int>();

            foreach (string word in words)
            {
                Spectrum.Add(Int32.Parse(word));
            }

            //  compute

            Graph<int, int> graph = BuildSpectrumGraph(Spectrum);

            List<Edge<int, int>> edges = graph.Edges();
            edges.Sort((a, b) => (a.source.label.CompareTo(b.source.label)));

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            foreach (Edge<int, int> edge in edges)
            {
                output.WriteLine("{0}->{1}:{2}", edge.source.label, edge.target.label, edge.label);
            }

            output.Flush();
        }
        /*
         LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)
            PrefixMass(0) ← 0
            for i ← 1 to |Peptide|
                for j ← 1 to 20
                    if AminoAcid(j) =  i-th amino acid in Peptide
                        PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass(j)
            LinearSpectrum ← a list consisting of the single integer 0
            for i ← 0 to |Peptide| − 1
                for j ← i + 1 to |Peptide|
                    add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
            return sorted list LinearSpectrum
        */
        static List<int> LinearSpectrum(string Peptide)
        {
            List<int> spectrum = new List<int>();
            List<int> PrefixMass = new List<int>();
            PrefixMass.Add(0);
            for (int i = 0; i < Peptide.Length; ++i)
                for (int j = 0; j < 20; ++j)
                    if (AminoAcid[j] == Peptide[i])
                        PrefixMass.Add(PrefixMass[i] + AminoAcidMass[j]);

            spectrum.Add(0);
            for (int i = 0; i < Peptide.Length; ++i)
                for (int j = i + 1; j <= Peptide.Length; ++j)
                    spectrum.Add(PrefixMass[j] - PrefixMass[i]);

            spectrum.Sort();

            return spectrum;
        }
        /*
         * CyclicSpectrum(Peptide, AminoAcid, AminoAcidMass)
            PrefixMass(0) ← 0
            for i ← 1 to |Peptide|
                for j ← 1 to 20
                    if AminoAcid(j) =  i-th amino acid in Peptide
                        PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass(j)
            peptideMass ← PrefixMass(|Peptide|)
            CyclicSpectrum ← a list consisting of the single integer 0
            for i ← 0 to |Peptide| − 1
                for j ← i + 1 to |Peptide|
                    add PrefixMass(j) − PrefixMass(i) to CyclicSpectrum
                    if i > 0 and j < |Peptide|
                        add peptideMass - (PrefixMass(j) − PrefixMass(i)) to CyclicSpectrum
            return sorted list CyclicSpectrum
        */
        static List<int> CyclicSpectrum(string Peptide)
        {
            List<int> spectrum = new List<int>();
            List<int> PrefixMass = new List<int>();
            PrefixMass.Add(0);
            for (int i = 0; i < Peptide.Length; ++i)
                for (int j = 0; j < 20; ++j)
                    if (AminoAcid[j] == Peptide[i])
                        PrefixMass.Add(PrefixMass[i] + AminoAcidMass[j]);

            int peptideMass = PrefixMass[Peptide.Length];

            spectrum.Add(0);
            for (int i = 0; i < Peptide.Length; ++i)
                for (int j = i + 1; j <= Peptide.Length; ++j)
                {
                    spectrum.Add(PrefixMass[j] - PrefixMass[i]);
                    if (i > 0 && j < Peptide.Length)
                        spectrum.Add(peptideMass - (PrefixMass[j] - PrefixMass[i]));
                }

            spectrum.Sort();

            return spectrum;
        }
        static string Spell(List<Edge<int, int>> path)
        {
            StringBuilder sb = new StringBuilder();
            foreach (Edge<int, int> edge in path)
            {
                sb.Append(edge.label);
            }
            return sb.ToString();
        }
        /*
            DecodingIdealSpectrum(Spectrum)
             construct Graph(Spectrum)
             for each path Path from source to sink in Graph(Spectrum)
                  Peptide ← the amino acid string spelled by the edge labels of Path
                  if IdealSpectrum(Peptide) = Spectrum
                        return Peptide
        */
        static string DecodingIdealSpectrum(List<int> Spectrum)
        {
            Graph<int, int> graph = BuildSpectrumGraph(Spectrum);

            Node<int, int> sink = null;
            Node<int, int> source = null;
            for (int i = 0; i < graph.nodes.Count; ++i)
            {
                if (graph.nodes[i].outgoing.Count == 0)
                {
                    sink = graph.nodes[i];
                }
                if (graph.nodes[i].incoming.Count == 0)
                {
                    source = graph.nodes[i];
                }
            }

            List<List<Edge<int, int>>> allPaths = graph.AllPaths(source.label, sink.label);

            foreach (List<Edge<int, int>> path in allPaths)
            {
                string peptide = Spell(path);
                List<int> idealSpectrum = LinearSpectrum(peptide);
                int i = 0;
                for (; i < Spectrum.Count; ++i)
                {
                    if (!idealSpectrum.Contains(Spectrum[i]))
                        break;
                }
                if (i == Spectrum.Count) return peptide;
            }

            return "Not found!";
        }
        static void SolveIdealSpectrumDecoding()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> Spectrum = new List<int>();

            foreach (string word in words)
            {
                Spectrum.Add(Int32.Parse(word));
            }

            //  compute

            string peptide = DecodingIdealSpectrum(Spectrum);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(peptide);
            output.Flush();
        }
        static List<int> PeptideVector(string peptide)
        {
            List<int> vector = new List<int>();
            foreach (char acid in peptide)
            {
                int mass = AminoAcidToMass(acid);
                while (--mass > 0) vector.Add(0);
                vector.Add(1);
            }
            return vector;
        }
        static void SolvePeptideVector()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            //  compute

            List<int> peptideVector = PeptideVector(line);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            for (int i = 0; i < peptideVector.Count; ++i)
            {
                output.Write(peptideVector[i]);
                if (i != (peptideVector.Count - 1)) output.Write(" ");
            }
            output.Flush();
        }
        static string Peptide(List<int> peptideVector)
        {
            StringBuilder peptide = new StringBuilder();

            int last = -1;
            for (int i = 0; i < peptideVector.Count; ++i)
            {
                if (peptideVector[i] == 1)
                {
                    peptide.Append(MassToAminoAcid(i - last));
                    last = i;
                }
            }
            return peptide.ToString();
        }
        static void SolvePeptide()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> peptideVector = new List<int>();

            foreach (string word in words)
            {
                peptideVector.Add(Int32.Parse(word));
            }

            //  compute

            string peptide = Peptide(peptideVector);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(peptide);
            output.Flush();
        }

        //********************************************************
        //  Peptide sequencing based on spectral vector
        static string PeptideSequencing(List<int> spectralVector, out int score)
        {
            spectralVector.Insert(0, 0);
            Graph<int, int> graph = new Graph<int, int>();
            for (int i = 0; i < spectralVector.Count; ++i)
            {
                graph.nodes.Add(new Node<int, int>(i, spectralVector[i], spectralVector[i]));
            }
            for (int i = 0; i < spectralVector.Count - 1; ++i)
                for (int j = i; j < spectralVector.Count; ++j)
                {
                    if (AminoAcidMass.Contains(j - i))
                    {
                        graph.AddEdge(i, j, /*j-i*/ spectralVector[j], MassToAminoAcid(j - i).ToString());
                    }
                }

            List<Edge<int, int>> path = graph.LongestPath(0, spectralVector.Count - 1);
            StringBuilder sb = new StringBuilder();
            score = 0;
            for (int i = 0; i < path.Count; ++i)
            {
                sb.Append(path[i].label);
                score += (int)path[i].weight;
            }

            return sb.ToString();
        }
        static string PeptideSequencing(List<int> spectralVector)
        {
            int score;
            return PeptideSequencing(spectralVector, out score);
        }
        static void SolvePeptideSequencing()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> spectralVector = new List<int>();

            foreach (string word in words)
            {
                spectralVector.Add(Int32.Parse(word));
            }

            //  compute

            string peptide = PeptideSequencing(spectralVector);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(peptide);
            output.Flush();
        }
        //******************************************************************************
        //  Peptide identification
        static int Score(string peptide, List<int> spectralVector)
        {
            int score = 0;
            int pos = 0;
            foreach (char acid in peptide)
            {
                pos += AminoAcidToMass(acid);
                if (pos > spectralVector.Count) return Int32.MinValue;
                score += spectralVector[pos - 1];
            }
            if (pos != spectralVector.Count) score = Int32.MinValue;
            return score;
        }
        static string PeptideIdentification(List<int> spectralVector, string Proteome)
        {
            string maxPeptide = "";
            int maxScore = Int32.MinValue;
            int score;
            HashSet<string> visited = new HashSet<string>();
            for (int length = spectralVector.Count / AminoAcidMass.Last(); length < spectralVector.Count / AminoAcidMass.First(); ++length)
            {
                for (int pos = 0; pos < Proteome.Length - length; ++pos)
                {
                    string peptide = Proteome.Substring(pos, length);
                    if (visited.Contains(peptide)) continue;
                    score = Score(peptide, spectralVector);
                    if (score > maxScore)
                    {
                        maxScore = score;
                        maxPeptide = peptide;
                    }

                    visited.Add(peptide);
                }
                visited.Clear();
            }
            return maxPeptide;
        }
        static void SolvePeptideIdentification()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> spectralVector = new List<int>();

            foreach (string word in words)
            {
                spectralVector.Add(Int32.Parse(word));
            }

            string proteome = input.ReadLine();

            //  compute

            string peptide = PeptideIdentification(spectralVector, proteome);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(peptide);
            output.Flush();
        }
        //**************************************************************************
        //  PSMSearch
        /*
            PSMSearch(SpectralVectors, Proteome, threshold).
                PSMSet ← an empty set
                for each vector Spectrum' in SpectralVectors
                      Peptide ← PeptideIdentification(Spectrum', Proteome)
                      if Score(Peptide, Spectrum) ≥ threshold
                          add the PSM (Peptide, Spectrum') to PSMSet
                return PSMSet 
         */
        static List<string> PSMSearch(List<List<int>> SpectralVectors, string Proteome, int threshold)
        {
            List<string> PSMSet = new List<string>();

            foreach (List<int> spectralVector in SpectralVectors)
            {
                string peptide = PeptideIdentification(spectralVector, Proteome);
                if (Score(peptide, spectralVector) >= threshold)
                {
                    PSMSet.Add(peptide);
                }
            }

            return PSMSet;
        }
        static void SolvePSMSearch()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = "";

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };

            List<List<int>> spectralVectors = new List<List<int>>();

            while (true)
            {
                line = input.ReadLine();

                if (!line.Contains('-')) break;

                List<int> spectralVector = new List<int>();

                string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

                foreach (string word in words)
                {
                    spectralVector.Add(Int32.Parse(word));
                }

                spectralVectors.Add(spectralVector);
            }

            string proteome = line;
            int threshold = Int32.Parse(input.ReadLine());

            //  compute

            List<string> peptides = PSMSearch(spectralVectors, proteome, threshold);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            foreach (string peptide in peptides)
                output.WriteLine(peptide);
            output.Flush();
        }
        //***************************************************************************
        //  Spectral Dictionary 
        static int SpectralDictionary(List<int> spectralVector, int threshold, int max_score)
        {
            spectralVector.Insert(0, 0);

            int size_i = spectralVector.Count;
            int size_t = max_score;

            int[,] m = new int[size_i, size_t];

            m[0, 0] = 1;
            for (int t = 1; t < size_t; ++t) m[0, t] = 0;

            for (int i = 1; i < size_i; ++i)
            {
                for (int t = 0; t < size_t; ++t)
                {
                    m[i, t] = 0;
                    if (spectralVector[i] > t || (t - spectralVector[i]) >= size_t) continue;
                    foreach (int a in AminoAcidMass)
                    {
                        if (a > i) continue;
                        m[i, t] += m[i - a, t - spectralVector[i]];
                    }
                }
            }

            int score = 0;
            for (int t = threshold; t < max_score; ++t) score += m[size_i - 1, t];

            ///////////// debug
            /*          for (int j = 0; j < size_t; ++j)
                      {
                          for (int i = 0; i < size_i; ++i)    
                          {
                              Console.Write(m[i, j].ToString());
                          }
                          Console.WriteLine();
                      }
                      Console.WriteLine();
              */
            Console.WriteLine(score.ToString());
            ////////////

            return score;
        }
        static void SolveSpectralDictionary()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> spectralVector = new List<int>();

            foreach (string word in words)
            {
                spectralVector.Add(Int32.Parse(word));
            }

            int threshold = Int32.Parse(input.ReadLine());
            int max_score = Int32.Parse(input.ReadLine());

            //  compute

            int score = SpectralDictionary(spectralVector, threshold, max_score);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(score.ToString());
            output.Flush();
        }


        //********************************************************************
        //  SpectralDictionaryProbability
        static double SpectralDictionaryProbability(List<int> spectralVector, int threshold, int max_score)
        {
            spectralVector.Insert(0, 0);

            int size_i = spectralVector.Count;
            int size_t = max_score;

            double[,] m = new double[size_i, size_t];

            m[0, 0] = 1;
            for (int t = 1; t < size_t; ++t) m[0, t] = 0;

            for (int i = 1; i < size_i; ++i)
            {
                for (int t = 0; t < size_t; ++t)
                {
                    m[i, t] = 0;
                    if (spectralVector[i] > t || (t - spectralVector[i]) >= size_t) continue;
                    foreach (int a in AminoAcidMass)
                    {
                        if (a > i) continue;
                        m[i, t] += m[i - a, t - spectralVector[i]] / 20;
                    }
                }
            }

            double probability = 0;
            for (int t = threshold; t < max_score; ++t) probability += m[size_i - 1, t];

            ///////////// debug
            for (int j = 0; j < size_t; ++j)
            {
                for (int i = 0; i < size_i; ++i)
                {
                    Console.Write(m[i, j].ToString());
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine(probability.ToString());
            ////////////

            return probability;
        }
        static void SolveSpectralDictionaryProbability()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> spectralVector = new List<int>();

            foreach (string word in words)
            {
                spectralVector.Add(Int32.Parse(word));
            }

            int threshold = Int32.Parse(input.ReadLine());
            int max_score = Int32.Parse(input.ReadLine());

            //  compute

            double probability = SpectralDictionaryProbability(spectralVector, threshold, max_score);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(probability.ToString());
            output.Flush();
        }
        //**************************************************************
        //  Spectral Alignment

        static List<int> SpectralAlignment(string peptide, List<int> spectralVector, int k)
        {
            List<int> modifications = new List<int>();

            List<int[,]> score = new List<int[,]>();

            int delta = 0;

            //         int[,]

            return modifications;
        }
        static void SolveSpectralAlignment()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string peptide = input.ReadLine();
            string line = input.ReadLine();

            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            List<int> spectralVector = new List<int>();

            foreach (string word in words)
            {
                spectralVector.Add(Int32.Parse(word));
            }

            int k = Int32.Parse(input.ReadLine());

            //  compute

            List<int> modifications = SpectralAlignment(peptide, spectralVector, k);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            for (int i = 0; i < peptide.Length; ++i)
            {
                output.Write(peptide[i]);
                if (modifications[i] != 0)
                {
                    output.Write(modifications[i] > 0 ? "({+0})" : "({0})", modifications[i]);
                }

            }
            output.Flush();
        }
        //********************************************************* 
        //  Charging station: From genome to breakpoint graph

        static void SolveChromosomeToCycle()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Chromosome chromosome = new Chromosome();
            chromosome.Parse(input.ReadLine());
            //  compute
            Nodes nodes = chromosome.ToCycle();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(nodes.ToString());
            output.Flush();
        }
        
        static void SolveCycleToChromosome()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Nodes nodes = new Nodes();
            nodes.Parse(input.ReadLine());
            //  compute
            Chromosome chromosome = nodes.ToChromosome();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(chromosome);
            output.Flush();
        }
        static void SolveColoredEdges()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Genome genome = new Genome();
            genome.Parse(input.ReadLine());          
            //  compute
            Edges edges = genome.ColoredEdges();
            edges.Sort();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(edges.ToString());
            output.Flush();
        }
        
        static void SolveGraphToGenome()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Edges coloredEdges = new Edges();
            coloredEdges.Parse(input.ReadLine());
            //  compute
            Genome genome = coloredEdges.ToGenome();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(genome.ToString());
            output.Flush();
        }
        static void SolveTwoBreakDistance()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Genome P = new Genome();
            P.Parse(input.ReadLine());
            int blocks = P.BlockNumber();
            Genome Q = new Genome();
            Q.Parse(input.ReadLine());
            //  compute
            Edges edges = P.ColoredEdges();
            edges.Add(Q.ColoredEdges());
            int distance = blocks - edges.Cycles().Count;
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(distance.ToString());
            output.Flush();
        }
        //***************************************************************
        //  Charging station: Solving the 2-Break Sorting Problem
        
        static void SolveTwoBreakOnGenomeGraph()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Edges coloredEdges = new Edges();
            coloredEdges.Parse(input.ReadLine());
            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int i = Int32.Parse(parameter_string[0]);
            int im = Int32.Parse(parameter_string[1]);
            int j = Int32.Parse(parameter_string[2]);
            int jm = Int32.Parse(parameter_string[3]);
            //  compute
            coloredEdges.TwoBreak( i, im, j, jm);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(coloredEdges.ToString());
            output.Flush();
        }
        static void SolveTwoBreakOnGenome()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Genome genome = new Genome();
            genome.Parse(input.ReadLine());
            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int i = Int32.Parse(parameter_string[0]);
            int im = Int32.Parse(parameter_string[1]);
            int j = Int32.Parse(parameter_string[2]);
            int jm = Int32.Parse(parameter_string[3]);
            //  compute
            genome.TwoBreak(i, im, j, jm);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(genome.ToString());
            output.Flush();
        }
        //*********************************************************
        //  Shared kmers

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
        static List<int> SharedKMers(string a, string b, int k)
        {
            List<int> shared = new List<int>();
            string brev = ReverseComplement(b);

            for(int i=0;i <= a.Length - k;++i)
            {
                string kmer = a.Substring(i, k);

                int pos = -1;
                while( (pos = b.IndexOf(kmer, pos+1)) != -1)
                {
                    shared.Add(i);
                    shared.Add(pos);
                }
                pos = -1;
                while ((pos = brev.IndexOf(kmer, pos+1)) != -1)
                {
                    shared.Add(i);
                    shared.Add(brev.Length - pos - k);
                }
            }

            return shared;
        }
        static void SolveSharedKMers()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            int k = Int32.Parse(input.ReadLine());
            string a = input.ReadLine();
            string b = input.ReadLine();

            List<int> shared = SharedKMers(a, b, k);

            //  write result
            StreamWriter output = new StreamWriter("result.txt");

            for (int i = 0; i < shared.Count / 2; ++i)
                output.WriteLine("({0}, {1})", shared[2 * i], shared[2 * i+1]);

            output.Flush();
        }
        static void Main(string[] args)
        {
            SolveUPGMA2();
        }

        static char[] AminoAcid = { 'X', 'Z', 'G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W' };
        static int[] AminoAcidMass = { 4, 5, 57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186 };

    }
}
