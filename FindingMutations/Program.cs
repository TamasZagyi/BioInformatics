using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FindingMutations
{
    class Program
    {
        static void SolveTrieConstruction()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph = new Graph<int>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    graph.AddToTrie(line.Trim());
                }
            }

            //  compute
            
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(graph.ToString());
            output.Flush();
        }
        static void SolvePrefixTrieMatching()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph = new Graph<int>();
            string Text = input.ReadLine().Trim();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    graph.AddToTrie(line.Trim());
                }
            }

            //  compute

            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            for(int i=0;i<Text.Length;++i)
            {
                if(graph.PrefixTrieMatching(Text.Substring(i)))
                {
                    output.Write("{0} ", i);
                }
            }
            
            output.Flush();
        }
        static void SolveSuffixTrie()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph = new Graph<int>();
            string Text = input.ReadLine().Trim();
            for (int i = 0; i < Text.Length; ++i)
            {
                graph.AddToTrie(Text.Substring(i));
                Console.WriteLine("{0}", i);
            }

            //  compute
            graph.Collapse();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            List<Edge<int>> edges = graph.Edges();
            foreach (Edge<int> edge in edges)
                output.WriteLine(edge.label);

            output.WriteLine("leaves: {0}", graph.Leaves());
            output.Flush();
        }
        static void SolveLongestRepeat()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph = new Graph<int>();
            string Text = input.ReadLine().Trim();
            for (int i = 0; i < Text.Length; ++i)
            {
                graph.AddToTrie(Text.Substring(i));
                Console.WriteLine("{0}", i);
            }

            //  compute
            graph.Collapse();
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(graph.LongestRepeat());
        //    output.WriteLine("");
        //    output.Write(graph.ToString());
            output.Flush();
        }
        static void SolveSharedSubstring()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph1 = new Graph<int>();
            Graph<int> graph2 = new Graph<int>();
            string Text1 = input.ReadLine().Trim();
            string Text2 = input.ReadLine().Trim();
            for (int i = 0; i < Text1.Length; ++i)
            {
                graph1.AddToTrie(Text1.Substring(i));
                Console.WriteLine("{0}", i);
            }
            for (int i = 0; i < Text2.Length; ++i)
            {
                graph2.AddToTrie(Text2.Substring(i));
                Console.WriteLine("{0}", i);
            }
            //  compute
            string winner = "";
            for(int i = 0; i < Math.Min( Text1.Length, Text2.Length); ++i)
            {
                List<string> subStrings = graph1.Substring(i);
                foreach(string s in subStrings)
                {
                    if(graph2.TrieMatching(s))
                    {
                        winner = s;
                        break;
                    }
                }
                if (winner.Length < i) break;
            }
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(winner);
            //    output.WriteLine("");
            //    output.Write(graph.ToString());
            output.Flush();
        }
        static void SolveNonSharedSubstring()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            Graph<int> graph1 = new Graph<int>();
            Graph<int> graph2 = new Graph<int>();
            string Text1 = input.ReadLine().Trim();
            string Text2 = input.ReadLine().Trim();
            for (int i = 0; i < Text1.Length; ++i)
            {
                graph1.AddToTrie(Text1.Substring(i));
                Console.WriteLine("{0}", i);
            }
            for (int i = 0; i < Text2.Length; ++i)
            {
                graph2.AddToTrie(Text2.Substring(i));
                Console.WriteLine("{0}", i);
            }
            //  compute
            string winner = "";
            for (int i = 1; i < Math.Min(Text1.Length, Text2.Length); ++i)
            {
                List<string> subStrings = graph1.Substring(i);

                foreach (string s in subStrings)
                {
                    if (!graph2.TrieMatching(s))
                    {
                        winner = s;
                        break;
                    }
                }
                if (winner.Length > 0) break;
            }
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.WriteLine(winner);

            output.Flush();
        }
        ////////////////////////////////////////////////////////////////////
        //  Suffix arrays
        static List<int> ConstructSuffixArray(string Text)
        {
            List<string> Suffixes = new List<string>();
            List<int> SuffixArray = new List<int>();

            for(int i=0;i<Text.Length;++i)
            {
                Suffixes.Add(Text.Substring(i));
            }
            List<string> Sorted = new List<string>(Suffixes);
            Sorted.Sort();

            for (int i = 0; i < Suffixes.Count; ++i)
            {
                SuffixArray.Add( Suffixes.IndexOf( Sorted[i] ) );
            }

            return SuffixArray;
        }
        static void SolveSuffixArrayConstruction()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string Text = input.ReadLine().Trim();

            //  compute
            List<int> SuffixArray = ConstructSuffixArray(Text);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            for (int i = 0; i < SuffixArray.Count; ++i)
            {
                output.Write( i == SuffixArray.Count - 1 ? "{0} " : "{0}, ",SuffixArray[i] );
            }          

            output.Flush();
        }
        static string BurrowsWheelerTransform(string Text)
        {
            List<string> M = new List<string>();
            for (int i = 0; i < Text.Length; ++i)
            {
                M.Add(Text.Substring(i) + Text.Substring(0, i));
            }
            M.Sort();
            StringBuilder sb = new StringBuilder();
            foreach (string s in M) sb.Append(s.Last());
            return sb.ToString();
        }
        static int OccurenceToPos(List<char> text, char c, int index)
        {
            int pos = text.IndexOf(c, 0);
            for (int p = 1; p < index; ++p)
            {
                pos = text.IndexOf(c, pos + 1);
            }
            return pos;
        }
        static void PosToOccurence(List<char> text, int pos, ref char c, ref int index)
        {
            index = 1;
            int lookUp = text.IndexOf(c, 0);
            while (lookUp != pos)
            {
                lookUp = text.IndexOf(c, lookUp + 1);
                ++index;
            }
        }
        static string InverseBurrowsWheelerTransform(string Last)
        {
            List<char> LastAsList = Last.ToList();
            List<char> FirstAsList = new List<char>(LastAsList);
            FirstAsList.Sort();
            
            StringBuilder sb = new StringBuilder();
            int index = 1;
            char c = FirstAsList[0];
            for (int i=1; i < FirstAsList.Count;++i)
            {
                int pos = OccurenceToPos(LastAsList, c, index);
                c = FirstAsList[pos];
                PosToOccurence(FirstAsList, pos, ref c, ref index);
                sb.Append(c);
            }
            sb.Append(FirstAsList[0]);
            return sb.ToString();
        }
        static void SolveBurrowsWheelerTransform()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string Text = input.ReadLine().Trim();
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(BurrowsWheelerTransform(Text));

            output.Flush();
        }
        static void SolveInverseBurrowsWheelerTransform()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string Text = input.ReadLine().Trim();
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(InverseBurrowsWheelerTransform(Text));

            output.Flush();
        }
        static int CountLetters(List<char> where, int until, char what)
        {
            int count = 0;
            for(int i = 0; i < until; ++i)
            {
                if (where[i] == what)
                    ++count;
            }
            return count;
        }
        static List<int> ConstructLastToFirst(string Last)
        {
            List<int> FirstIndices = new List<int>();
            List<char> FirstLetters = Last.ToList();
            FirstLetters.Sort();
            List<int> LastIndices = new List<int>();
            List<char> LastLetters = Last.ToList();
            
            for(int i = 0; i < FirstLetters.Count;++i)
            {
                FirstIndices.Add(CountLetters(FirstLetters, i, FirstLetters[i]) + 1);
            }

            for (int i = 0; i < LastLetters.Count; ++i)
            {
                LastIndices.Add(CountLetters(LastLetters, i, LastLetters[i]) + 1);
            }

            List<int> LastToFirst = new List<int>();
            for (int i = 0; i < LastLetters.Count; ++i)
            {
                int j = 0;
                for(; j < FirstLetters.Count; ++j)
                {
                    if (FirstLetters[j] == LastLetters[i] && FirstIndices[j] == LastIndices[i]) break;
                }
                LastToFirst.Add(j);
            }
            return LastToFirst;
        }
        /*
        BWMATCHING(FirstColumn, LastColumn, Pattern, LastToFirst)
        top ← 0
        bottom ← |LastColumn| − 1
        while top ≤ bottom
            if Pattern is nonempty
                symbol ← last letter in Pattern
                remove last letter from Pattern
                if positions from top to bottom in LastColumn contain an occurrence of symbol
                    topIndex ← first position of symbol among positions from top to bottom in LastColumn
                    bottomIndex ← last position of symbol among positions from top to bottom in LastColumn
                    top ← LastToFirst(topIndex)
                    bottom ← LastToFirst(bottomIndex)
                else
                    return 0
            else
                return bottom − top + 1
        */
        static int BWMATCHING(List<char> FirstColumn, List<char> LastColumn, List<char> Pattern,List<int> LastToFirst)
        {
            int top = 0;
            int bottom = LastColumn.Count - 1;
            while(top<=bottom)
            {
                if (Pattern.Count > 0)
                {
                    char symbol = Pattern.Last();
                    Pattern.RemoveAt(Pattern.Count - 1);
                    int topIndex = top;
                    while (topIndex <= bottom && LastColumn[topIndex] != symbol) ++topIndex;
                    int bottomIndex = bottom;
                    while (bottomIndex >= top && LastColumn[bottomIndex] != symbol) --bottomIndex;
                    if (topIndex <= bottom && bottomIndex >= top)
                    {
                        top = LastToFirst[topIndex];
                        bottom = LastToFirst[bottomIndex];
                    }
                    else return 0;
                }
                else return bottom - top + 1;
            }
            return 0;
        }
        static void SolveBWMatching()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string Text = input.ReadLine().Trim();

            List<char> LastColumn = Text.ToList();
            List<char> FirstColumn = Text.ToList();
            FirstColumn.Sort();
            List<int> LastToFirst = ConstructLastToFirst(Text);
            string Patterns = input.ReadLine();
            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = Patterns.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            
            StreamWriter output = new StreamWriter("result.txt");
            foreach (string word in words)
                output.Write("{0} ", BWMATCHING(FirstColumn,LastColumn,word.ToList(),LastToFirst));

            output.Flush();
        }
        static List<int> BWMATCHING2(List<char> FirstColumn, List<char> LastColumn, List<char> Pattern, List<int> LastToFirst, List<int> SuffixArray)
        {
            List<int> MatchPositions = new List<int>();
            int top = 0;
            int bottom = LastColumn.Count - 1;
            while (top <= bottom)
            {
                if (Pattern.Count > 0)
                {
                    char symbol = Pattern.Last();
                    Pattern.RemoveAt(Pattern.Count - 1);
                    int topIndex = top;
                    while (topIndex <= bottom && LastColumn[topIndex] != symbol) ++topIndex;
                    int bottomIndex = bottom;
                    while (bottomIndex >= top && LastColumn[bottomIndex] != symbol) --bottomIndex;
                    if (topIndex <= bottom && bottomIndex >= top)
                    {
                        top = LastToFirst[topIndex];
                        bottom = LastToFirst[bottomIndex];
                    }
                    else return MatchPositions;
                }
                else
                {
                    //return bottom - top + 1;
              //      top = LastToFirst[top];
              //      bottom = LastToFirst[bottom];
                    for (int pos = top; pos < bottom + 1; ++pos)
                        MatchPositions.Add(SuffixArray[  pos ] );
                    MatchPositions.Sort();
                    return MatchPositions;
                }
                    
            }
            return MatchPositions;
        }
        static void SolveBWMatching2()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string Text = input.ReadLine().Trim();
            List<int> SuffixArray = ConstructSuffixArray(Text);
            Text = BurrowsWheelerTransform(Text);

            List<char> LastColumn = Text.ToList();
            List<char> FirstColumn = Text.ToList();
            FirstColumn.Sort();
            List<int> LastToFirst = ConstructLastToFirst(Text);


            List<int> AllMatches = new List<int>();
            while (!input.EndOfStream)
            {
                string word = input.ReadLine().Trim();
                if (word.Length == 0) break;
                List<int> MatchPositions = BWMATCHING2(FirstColumn, LastColumn, word.ToList(), LastToFirst, SuffixArray);
                AllMatches.AddRange(MatchPositions);
            }
            //       AllMatches.Sort();

            StreamWriter output = new StreamWriter("result.txt");
            List<string> Matches = new List<string>();
            foreach (int pos in AllMatches)
            {
                if (!Matches.Contains(pos.ToString())) Matches.Add(pos.ToString());
            }
            Matches.Sort();

            foreach (string pos in Matches)
                output.Write("{0} ", pos);

            output.Flush();
        }
        static List<int> MismatchTolerantBW(List<char> FirstColumn, List<char> LastColumn, List<char> Pattern, List<int> LastToFirst, List<int> SuffixArray, int d, string Text)
        {
            List<int> MatchPositions = new List<int>();

            int partL = Pattern.Count / (d + 1);
            int partN = d + 1;

            List<int> PartPos = new List<int>();
            for(int i = 0; i < partN; ++i)
            {
                PartPos.Add( i * partL );
            }

            for(int i = 0; i < partN; ++i)
            {
                int partEnd = Math.Min(PartPos[i] + partL, Pattern.Count);
                List<int> PartialMatch = BWMATCHING2(FirstColumn, LastColumn, 
                    Pattern.GetRange(PartPos[i], partEnd - PartPos[i]) , 
                    LastToFirst, SuffixArray);
                foreach(int pos in PartialMatch)
                {
                    int diff = 0;
                    for(int c = 0; c < Pattern.Count; ++c)
                        if(Text[pos - PartPos[i] + c] != Pattern[c])
                        {
                            ++diff;
                            if (diff > d) break;
                        }
                    if (diff <= d)
                    {
                        if(!MatchPositions.Contains(pos - PartPos[i]))
                            MatchPositions.Add(pos - PartPos[i]);
                    }
                }
            }

            MatchPositions.Sort();

            return MatchPositions;
        }
        static void SolveMismatchTolerant()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string TextO = input.ReadLine().Trim() + "$";
            List<int> SuffixArray = ConstructSuffixArray(TextO);
            string Text = BurrowsWheelerTransform(TextO);

            List<char> LastColumn = Text.ToList();
            List<char> FirstColumn = Text.ToList();
            FirstColumn.Sort();
            List<int> LastToFirst = ConstructLastToFirst(Text);
           
            List<int> AllMatches = new List<int>();
            string Patterns = input.ReadLine();
            char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
            string[] words = Patterns.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);

            int d = Int32.Parse(input.ReadLine());

            foreach(string pattern in words)
            {
                List<int> MatchPositions = MismatchTolerantBW(FirstColumn, LastColumn, pattern.ToList(), LastToFirst, SuffixArray,d, TextO);
                AllMatches.AddRange(MatchPositions);
            }

            AllMatches.Sort();
            StreamWriter output = new StreamWriter("result.txt");
            foreach (int pos in AllMatches)
                output.Write("{0} ", pos);

            output.Flush();
        }
        
        static void Main(string[] args)
        {
            //SolveMismatchTolerant();

            //          foreach (int i in ConstructSuffixArray("cocoon$"))
            //            Console.Write("{0} ", i);
            Console.WriteLine(InverseBurrowsWheelerTransform("TTCCATTGGA$"));
            Console.ReadKey();
        }
    }
}
