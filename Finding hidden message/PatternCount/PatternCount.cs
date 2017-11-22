using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace PatternCount
{
    class Program
    {
        string[,] RNA_codon_table = new string[,]{
{"AAA", "K"},{"AAC", "N"},{"AAG", "K"},{"AAU", "N"},{"ACA", "T"},{"ACC", "T"},{"ACG", "T"},{"ACU", "T"},{"AGA", "R"},{"AGC", "S"},{"AGG", "R"},{"AGU", "S"},{"AUA", "I"},
{"AUC", "I"},{"AUG", "M"},{"AUU", "I"},{"CAA", "Q"},{"CAC", "H"},{"CAG", "Q"},{"CAU", "H"},{"CCA", "P"},{"CCC", "P"},{"CCG", "P"},{"CCU", "P"},{"CGA", "R"},{"CGC", "R"},
{"CGG", "R"},{"CGU", "R"},{"CUA", "L"},{"CUC", "L"},{"CUG", "L"},{"CUU", "L"},{"GAA", "E"},{"GAC", "D"},{"GAG", "E"},{"GAU", "D"},{"GCA", "A"},{"GCC", "A"},{"GCG", "A"},
{"GCU", "A"},{"GGA", "G"},{"GGC", "G"},{"GGG", "G"},{"GGU", "G"},{"GUA", "V"},{"GUC", "V"},{"GUG", "V"},{"GUU", "V"},{"UAA", ""},{"UAC", "Y"},{"UAG", ""},{"UAU", "Y"},
{"UCA", "S"},{"UCC", "S"},{"UCG", "S"},{"UCU", "S"},{"UGA", ""},{"UGC", "C"},{"UGG", "W"},{"UGU", "C"},{"UUA", "L"},{"UUC", "F"},{"UUG", "L"},{"UUU", "F"}};

        static int PatternCount(string Text, string Pattern)
        {
            int count = 0;
            for(int i=0; i < Text.Length - Pattern.Length;++i)
            {
                if( Text.Substring(i,Pattern.Length).Equals(Pattern) )
                {
                    ++count;
                }
            }
            return count;
        }
        static List<int> PatternPositions(string Pattern, string Genome)
        {
            List<int> positions = new List<int>();

            int index = Genome.IndexOf(Pattern);
            while(index != -1)
            {
                positions.Add(index);
                index = Genome.IndexOf( Pattern, index + 1 );
            }

            return positions;
        }
        static List<string> FrequentWords(string Text, int k)
        {
            List<string> FrequentPatterns = new List<string>();
            List<int> Count = new List<int>();
            int maxCount = 0;
            for (int i = 0; i < Text.Length - k; ++i )
            {
                string Pattern = Text.Substring(i, k);
                int c = PatternCount(Text, Pattern);
                Count.Add(c);
                if( c > maxCount )
                {
                    maxCount = c;
                }
            }
            for (int i = 0; i < Text.Length - k; ++i )
            {
                if( Count[i] == maxCount )
                {
                    string Pattern = Text.Substring(i, k);
                    if( !FrequentPatterns.Contains(Pattern) )
                        FrequentPatterns.Add( Pattern );
                }
            }
            return FrequentPatterns;
        }
        //  A C G T
        static long PatternToNumber(string Pattern)
        {
            int k = Pattern.Length;
            long pow = 1;
            long val = 0;
            for (int i = k - 1; i >= 0; --i)
            {
                long digit = 0;
                switch(Pattern[i])
                {
                    case 'A': digit = 0; break;
                    case 'C': digit = 1; break;
                    case 'G': digit = 2; break;
                    case 'T': digit = 3; break;
                }
                val += digit * pow;
                pow *= 4;
            }
            return val;
        }
        static string NumberToPattern(long val, long k)
        {
            string Pattern = "";

            for (long i = 0; i < k; ++i)
            {
                switch(val % 4)
                {
                    case 0: Pattern = 'A' + Pattern; break;
                    case 1: Pattern = 'C' + Pattern; break;
                    case 2: Pattern = 'G' + Pattern; break;
                    case 3: Pattern = 'T' + Pattern; break;
                }
                val /= 4;
            }

            return Pattern;
        }
        static List<long> ComputingFrequencies(string Text, int k)
        {
            List<long> FrequencyArray = new List<long>();

            for (long i = 0; i < Math.Pow(4,k); ++i )
            {
                FrequencyArray.Add(0);
            }

            for (int i = 0; i < Text.Length - k+1; ++i )
            {
                FrequencyArray[(int)PatternToNumber(Text.Substring(i, k))]++;
            }
            
            return FrequencyArray;
        }
        //  A - T G - C
        static string ReverseComplement(string Text)
        {
            StringBuilder builder = new StringBuilder();
            for(int i=0; i< Text.Length; ++i)
            {
                char complement = Text[Text.Length-(i+1)];
                switch(complement)
                {
                    case 'A' : complement = 'T';break;
                    case 'T' : complement = 'A';break;
                    case 'G' : complement = 'C';break;
                    case 'C' : complement = 'G';break;
                }
                builder.Append(complement);
            }

            return builder.ToString();
        }
        //  k: k-mer
        //  t: frequency
        //  L: length of window
        static List<string> BetterClumpFinding(string Genome, int k, int t, int L)
        {
            List<int> Clump = new List<int>();
            for (int i = 0; i < Math.Pow(4, k); ++i)
            {
                Clump.Add(0);
            }
            
            //  initial frequency array
            string Text = Genome.Substring(0,L);
            List<long> FrequencyArray = ComputingFrequencies(Text, k);
            for (int i = 0; i < Math.Pow(4, k); ++i)
            {
                if(FrequencyArray[i] >= t)
                {
                    Clump[i] = 1;
                }
            }

            //  shift the window through the genom
            for( int i=1; i <= Genome.Length - L; ++i  )
            {
                string FirstPattern = Genome.Substring(i - 1, k);
                int index = (int)PatternToNumber(FirstPattern);
                FrequencyArray[index] = FrequencyArray[index] - 1;
                string LastPattern = Genome.Substring(i + L - k, k);
                index = (int)PatternToNumber(LastPattern);
                FrequencyArray[index] = FrequencyArray[index] + 1;

                //  conditional add to the clump
                if (FrequencyArray[index] >= t)
                {
                    Clump[index] = 1;
                }
            }

            //  query frequent patterns
            List<string> FrequentPatterns = new List<string>();
            for (int i = 0; i < Math.Pow(4, k); ++i)
            {
                if (Clump[i] == 1)
                {
                    FrequentPatterns.Add(NumberToPattern(i, k));
                }
            }

            return FrequentPatterns;
        }
        static List<int> MinSkew(string Genom)
        {
            List<int> MinSkew = new List<int>();

            int min = 0;
            int act = 0;

            for(int i=0;i<Genom.Length;++i)
            {
                if(Genom[i] == 'C') act--;
                if(Genom[i] == 'G') act++;
                if(act == min)
                {
                    MinSkew.Add(i+1);
                }
                if(act < min)
                {
                    MinSkew.Clear();
                    MinSkew.Add(i+1);
                    min = act;
                }
            }

            return MinSkew;
        }
        static int HammingDistance(string a, string b)
        {
            int dist = 0;
            for(int i = 0; i < a.Length;++i)
            {
                if (a[i] != b[i]) ++dist;
            }
            return dist;
        }
        static List<int> ApproximatePatternMatching(string Pattern,string Text, int max_dist)
        {
            List<int> Matches = new List<int>();

            for(int i=0;i<=Text.Length-Pattern.Length;++i)
            {
                if(HammingDistance(Pattern,Text.Substring(i, Pattern.Length)) <= max_dist)
                {
                    Matches.Add(i);
                }
            }

            return Matches;
        }
        static int ApproximatePatternCount(string Pattern, string Text, int max_dist)
        {
            int Matches = 0;

            for (int i = 0; i <= Text.Length - Pattern.Length; ++i)
            {
                if (HammingDistance(Pattern, Text.Substring(i, Pattern.Length)) <= max_dist)
                {
                    ++Matches;
                }
            }

            return Matches;
        }

        /*
        Neighbors(Pattern, d)
        if d = 0
            return {Pattern}
        if |Pattern| = 1 
            return {A, C, G, T}
        Neighborhood ← an empty set
        SuffixNeighbors ← Neighbors(Suffix(Pattern), d)
        for each string Text from SuffixNeighbors
            if HammingDistance(Suffix(Pattern), Text) < d
                for each nucleotide x
                    add x • Text to Neighborhood
            else
                add FirstSymbol(Pattern) • Text to Neighborhood
        return Neighborhood
        */
        static List<string> Neighbors(string Pattern,int d)
        {
            List<string> Neighborhood = new List<string>();
            char [] nucleotides = { 'A', 'C', 'G', 'T' };

            if( d == 0) { /*do nothing*/ }
            else if( Pattern.Length == 1 )
            {
                foreach (char c in nucleotides)
                    Neighborhood.Add(c.ToString());
            }
            else
            {
                List<string> SuffixNeighbors = Neighbors(Pattern.Substring(1),d);
                foreach(string Text in SuffixNeighbors)
                {
                    if( HammingDistance(Pattern.Substring(1), Text) < d )
                    {
                        foreach(char c in nucleotides)
                        {
                            Neighborhood.Add(c + Text);
                        }
                    }
                    else
                    {
                        Neighborhood.Add(Pattern[0] + Text);
                    }
                }
            }

            return Neighborhood;
        }
        /*
        FrequentWordsWithMismatches(Text, k, d)
        FrequentPatterns ← an empty set
        for i ← 0 to 4k − 1
            Close(i) ← 0
            FrequencyArray ← 0
        for i ← 0 to |Text| − k
            Neighborhood ← Neighbors(Text(i, k), d)
            for each Pattern from Neighborhood
                index ← PatternToNumber(Pattern)
                Close(index) ← 1
        for i ←0 to 4k − 1
            if Close(i) = 1
                Pattern ← NumberToPattern(i, k)
                FrequencyArray(i) ← ApproximatePatternCount(Text, Pattern, d)
        maxCount ← maximal value in FrequencyArray
        for i ←0 to 4k − 1
            if FrequencyArray(i) = maxCount
                Pattern ← NumberToPattern(i, k)
                add Pattern to the set FrequentPatterns
       return FrequentPatterns 
        */
        static List<string> FrequentWordsWithMismatches(string Text, int k,int d)
        {
            List<string> FrequentPatterns = new List<string>();

            List<int> FrequencyArray = new List<int>();
            List<int> Close = new List<int>();

            for (int i = 0; i <= Math.Pow(4, k) - 1; ++i)
            {
                FrequencyArray.Add(0);
                Close.Add(0);
            }

            for(int i = 0; i <= Text.Length - k; ++i)
            {
                foreach(string Pattern in Neighbors( Text.Substring(i,k), d ) )
                {
                    Close[(int)PatternToNumber(Pattern)] = 1;
                }
            }

            long maxCount = 0;
            for (int i = 0; i <= Math.Pow(4, k) - 1; ++i)
            {
                if( Close[i] == 1 )
                {
                    FrequencyArray[i] = ApproximatePatternCount( NumberToPattern(i, k), Text, d)
                                        + ApproximatePatternCount( ReverseComplement( NumberToPattern(i, k) ), Text, d);
                    if (FrequencyArray[i] > maxCount) maxCount = FrequencyArray[i];
                }
            }

            for (int i = 0; i <= Math.Pow(4, k) - 1; ++i)
            {
                if (FrequencyArray[i] == maxCount)
                {
                    FrequentPatterns.Add(NumberToPattern(i, k));
                }
            }

            return FrequentPatterns;
        }
        /*
        MOTIFENUMERATION(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in Dna
            for each k-mer Pattern’ differing from Pattern by at most d
              mismatches
                if Pattern' appears in each string from Dna with at most d
                mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
        */
        static List<string> MotifEnumeration(List<string> Dna, int k, int d)
        {
            List<string> Patterns = new List<string>();

            foreach(string source in Dna)
            {
                for (int i = 0; i <= source.Length - k; ++i)
                {
                    foreach (string Pattern in Neighbors(source.Substring(i, k), d))
                    {
                        bool bingo = true;
                        foreach (string check in Dna)
                        {
                            if (ApproximatePatternCount(Pattern, check,d) == 0)
                            {
                                bingo = false;
                                break;
                            }
                        }
                        if(bingo)
                        {
                            if(!Patterns.Contains(Pattern))
                                Patterns.Add(Pattern);
                        }
                    }
                }
            }

            return Patterns;
        }
        /*
         * DistanceBetweenPatternAndStrings(Pattern, Dna)
            k ← |Pattern|
            distance ← 0
            for each string Text in Dna
                HammingDistance ← ∞
                for each k-mer Pattern’ in Text
                    if HammingDistance > HammingDistance(Pattern, Pattern’)
                        HammingDistance ← HammingDistance(Pattern, Pattern’)
                distance ← distance + HammingDistance
            return distance
         */
        static int DistanceBetweenPatternAndStrings(string Pattern,List<string> Dna)
        {
            int k = Pattern.Length;
            int distance = 0;
            foreach(string Text in Dna)
            {
                int MinDistance = Int32.MaxValue;
                for(int i = 0; i <= Text.Length - k; ++i)
                {
                    int Distance = HammingDistance(Pattern, Text.Substring(i,k));
                    if(MinDistance > Distance)
                    {
                        MinDistance = Distance;
                    }
                }
                distance += MinDistance;
            }
            return distance;
        }
        /*
         MedianString(Dna, k)
    distance ← ∞
    for i ←0 to 4k −1
        Pattern ← NumberToPattern(i, k)
        if distance > DistanceBetweenPatternAndStrings(Pattern, Dna)
            distance ← DistanceBetweenPatternAndStrings(Pattern, Dna)
            Median ← Pattern
    return Median
         */
        static string MedianString(List<string> Dna,int k)
        {
            int distance = Int32.MaxValue;
            string Median = "";
            for(int i = 0; i < Math.Pow(4,k); ++i )
            {
                string Pattern = NumberToPattern(i, k);
                int d = DistanceBetweenPatternAndStrings(Pattern, Dna);
                if(distance > d)
                {
                    distance = d;
                    Median = Pattern;
                }
            }
            return Median;
        }
        static int NucleotideToIndex(char c)
        {
            switch(c)
            {
                case 'a':
                case 'A': return 0;
                case 'c':
                case 'C': return 1;
                case 'g':
                case 'G': return 2;
                case 't':
                case 'T': return 3;
            }
            return -1;
        }
        static char IndexToNucleotide(int i)
        {
            switch(i)
            {
                case 0: return 'A';
                case 1: return 'C';
                case 2: return 'G';
                case 3: return 'T';
            }
            return ' ';
        }
        static string ProfileMostProbableKMer(string Text, int k, List<List<double>> Profile)
        {
            double MaxProbability = 0;
            string kmer = Text.Substring(0, k);
            for(int i=0;i<=Text.Length-k;++i)
            {
                //  A,C,G,T
                string pattern = Text.Substring(i, k);
                double probability = 1;
                for(int j = 0; j < k; ++j)
                {
                    probability *= Profile[NucleotideToIndex(pattern[j])][j];
                }
                if(probability > MaxProbability)
                {
                    MaxProbability = probability;
                    kmer = pattern;
                }
            }
            return kmer;
        }
        /*
         * GREEDYMOTIFSEARCH(Dna, k, t)
        BestMotifs ← motif matrix formed by first k-mers in each string
                      from Dna
        for each k-mer Motif in the first string from Dna
            Motif1 ← Motif
            for i = 2 to t
                form Profile from motifs Motif1, …, Motifi - 1
                Motifi ← Profile-most probable k-mer in the i-th string
                          in Dna
            Motifs ← (Motif1, …, Motift)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
         */
        static List<string> GreedyMotifSearch(List<string> Dna, int t, int k)
        {
            List<string> BestMotifs = new List<string>();
            int BestScore = Int32.MaxValue;
            foreach(string Text in Dna)
            {
                BestMotifs.Add(Text.Substring(0, k));
            }

            if(Dna.Count == 0) return BestMotifs;

            for (int j = 0; j <= Dna[0].Length - k; ++j )
            {
                List<string> Motifs = new List<string>();
                Motifs.Add(Dna[0].Substring(j, k));
                for(int i=1; i < t;++i)
                {
                    List<List<double>> Profile = RawProfile(Motifs);
                    for (int x = 0; x < Profile.Count; ++x)
                        for (int y = 0; y < Profile[x].Count; ++y)
                            Profile[x][y] /= i*2;
                    Motifs.Add(ProfileMostProbableKMer(Dna[i],k,Profile));
                }
                int score = Score(Motifs);
                if (score < BestScore)
                {
                    BestMotifs = Motifs;
                    BestScore = score;
                }
            }

            return BestMotifs;
        }
        static List<List<double>> RawProfile(List<string> Motifs)
        {
            List<List<double>> Profile = new List<List<double>>();
            for (int x = 0; x < 4; ++x)
            {
                Profile.Add(new List<double>());
                for (int y = 0; y < Motifs[0].Length; ++y)
                {
                    Profile[x].Add(1);
                }
            }
            foreach (string Motif in Motifs)
            {
                for (int x = 0; x < Motif.Length; ++x)
                    Profile[NucleotideToIndex(Motif[x])][x]++;
            }
            return Profile;
        }
        static int Score(List<string> Motifs)
        {
            List<List<double>> Profile = RawProfile(Motifs);
            string Consensus = "";
            for (int i = 0; i < Motifs[0].Length;++i )
            {
                double max = 0;
                int nucleotide = 0;
                for (int j = 0; j < 4; ++j)
                {
                    if (Profile[j][i] > max)
                    {
                        max = Profile[j][i];
                        nucleotide = j;
                    }
                }
                Consensus += IndexToNucleotide(nucleotide);
            }

            int distance = 0;
            foreach(string Motif in Motifs)
            {
                distance += HammingDistance(Motif, Consensus);
            }

            return distance;
        }
        static void Main(string[] args)
        {
            Console.WriteLine("Where to read from? Console of File? (C/F)");
            string input = Console.ReadLine().ToUpper().Trim();

            List<string> Dna = new List<string>();
            int k = 0;
            int t = 0;

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");

                k = Int32.Parse(sr.ReadLine());
                t = Int32.Parse(sr.ReadLine());

                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();

                    if (line.Length > 0)
                    {
                        Dna.Add(line);
                    }
                }
            }
            else
            {
                return;
            }

            List<string> Motifs = GreedyMotifSearch(Dna, t, k);

            StreamWriter sw = new StreamWriter("result.txt");

            foreach (string motif in Motifs)
                sw.WriteLine(motif);

            sw.Flush();


   /*         string Text = "";
            int k = 0;
            List<List<double>> Profile = new List<List<double>>();
            Profile.Add(new List<double>());
            Profile.Add(new List<double>());
            Profile.Add(new List<double>());
            Profile.Add(new List<double>());

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");

                Text = sr.ReadLine();
                k = Int32.Parse(sr.ReadLine());

                for (int i = 0; i < 4;++i )
                {
                    string line = sr.ReadLine();

                    if (line.Length > 0)
                    {
                        char[] delimiterChars = { ' ', '\t' };
                        string[] words = line.Split(delimiterChars,StringSplitOptions.RemoveEmptyEntries);
                        foreach (string word in words)
                            Profile[i].Add(Double.Parse(word));
                    }
                }
            }
            else
            {
                return;
            }

            Console.WriteLine(ProfileMostProbableKMer(Text, k, Profile));
            */
 /*           List<string> Dna = new List<string>();
            int k=0;

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");

                k = Int32.Parse( sr.ReadLine() );
                
                while (!sr.EndOfStream)
                {
                    string line = sr.ReadLine();
                    
                    if (line.Length > 0)
                    {
                 //       char[] delimiterChars = { ' ', ',', '.', ':', '\t' };
                   //     string[] words = line.Split(delimiterChars,StringSplitOptions.RemoveEmptyEntries);
                     //   foreach(string word in words)
                            Dna.Add(line);
                    } 
                }
            }
            else
            {
                return;
            }

            Console.WriteLine(MedianString(Dna,k));
*/
      /*      StreamWriter sw = new StreamWriter("result.txt");

            foreach (string motif in MotifEnumeration(Patterns, k, d))
                sw.Write(motif + " ");

            sw.Flush();
            */
    /*        int k = 0;
            int d = 0;
            List<string> Patterns = new List<string>();

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");

                k = Int32.Parse(sr.ReadLine());
                d = Int32.Parse(sr.ReadLine());

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

            foreach (string motif in MotifEnumeration(Patterns, k, d))
                sw.Write(motif + " ");

            sw.Flush();
            */
   /*         string Pattern = "";
            int k = 0;
            int d = 0;

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                Pattern = sr.ReadLine();
                k = Int32.Parse(sr.ReadLine());
                d = Int32.Parse(sr.ReadLine());
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");

            foreach (string neighbor in FrequentWordsWithMismatches(Pattern, k, d))
                sw.Write(neighbor + " ");

            sw.Flush();
            */
            /*          string Pattern = "";
                      int k = 0;

                      if (input == "C")
                      {

                      }
                      else if (input == "F")
                      {
                          StreamReader sr = new StreamReader("dataset.txt");
                          Pattern = sr.ReadLine();
                          k = Int32.Parse(sr.ReadLine());
                      }
                      else
                      {
                          return;
                      }

                      StreamWriter sw = new StreamWriter("result.txt");

                      foreach(string neighbor in Neighbors(Pattern,k))
                          sw.WriteLine(neighbor);

                      sw.Flush();
          */
            /*
            string Pattern = "";
            string Text = "";
            int k = 0;

            if (input == "C")
            {

            }
            else if (input == "F")
            {
                StreamReader sr = new StreamReader("dataset.txt");
                Pattern = sr.ReadLine();
                Text = sr.ReadLine();
                k = Int32.Parse(sr.ReadLine());
            }
            else
            {
                return;
            }

            StreamWriter sw = new StreamWriter("result.txt");

           
                sw.Write("{0}", ApproximatePatternCount(Pattern,Text,k));
            

            sw.Flush();
            */
            /*        string Pattern = "";
                    string Text = "";
                    int k = 0;

                    if (input == "C")
                    {

                    }
                    else if (input == "F")
                    {
                        StreamReader sr = new StreamReader("dataset.txt");
                        Pattern = sr.ReadLine();
                        Text = sr.ReadLine();
                        k = Int32.Parse(sr.ReadLine());
                    }
                    else
                    {
                        return;
                    }

                    StreamWriter sw = new StreamWriter("result.txt");

                    foreach (int i in ApproximatePatternMatching(Pattern,Text, k))
                    {
                        sw.Write("{0} ", i);
                    }

                    sw.Flush();
                    */
            /*           string TextA = "";
                       string TextB = "";

                       if (input == "C")
                       {

                       }
                       else if(input == "F")
                       {
                           StreamReader sr = new StreamReader("dataset.txt");
                           TextA = sr.ReadLine();
                           TextB = sr.ReadLine();
                       }
                       else
                       {
                           return;
                       }

                       StreamWriter sw = new StreamWriter("result.txt");

                           sw.WriteLine("{0}", HammingDistance(TextA, TextB));

                       sw.Flush();
           */
            /*           string Text;

                       if (input == "C")
                       {
                           Console.WriteLine("Enter Text: ");
                           Text = Console.ReadLine();
                       }
                       else if (input == "F")
                       {
                           StreamReader sr = new StreamReader("dataset.txt");
                           Text = sr.ReadLine();
                       }
                       else
                       {
                           return;
                       }

                       StreamWriter sw = new StreamWriter("result.txt");
                       foreach (int i in MinSkew(Text))
                       {
                           sw.Write("{0} ", i);
                       }
                       sw.Flush();
           */
            /*           Console.WriteLine("Enter the genome: ");
                       string Text = Console.ReadLine();
                       Console.WriteLine("Enter the pattern length: ");
                       int k = Int32.Parse(Console.ReadLine());
                       Console.WriteLine("Enter window length: ");
                       int L = Int32.Parse(Console.ReadLine());
                       Console.WriteLine("Enter minimum frequency: ");
                       int t = Int32.Parse(Console.ReadLine());

                       try
                       {
                           StreamReader sr = new StreamReader("dataset.txt");
                           string Text = sr.ReadLine();
                           int k = Int32.Parse(sr.ReadLine());
                           int L = Int32.Parse(sr.ReadLine());
                           int t = Int32.Parse(sr.ReadLine());

                           List<string> FrequentPatterns = BetterClumpFinding(Text, k, t, L);
                           Console.WriteLine("The following k-mers have been found\n");
                           foreach (string kmer in FrequentPatterns)
                               Console.Write("{0} ", kmer);

                           Console.WriteLine("Number of elements in the list {0}", FrequentPatterns.Count);
            */
            /*               List<int> Frequencies = PatternPositions(Pattern, Text);
            //               Console.WriteLine("The following frequencies have been found\n");
                           foreach (int f in Frequencies)
                               Console.Write("{0} ", f);

             //              Console.WriteLine(ReverseComplement(Text));
                       }
                       catch(Exception e)
                       {
                           Console.WriteLine(e.StackTrace);
                       }
           */

            //            Console.WriteLine("{0}", NumberToPattern(7020, 9));


            /*            List<string> FrequentPatterns = FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA",3);
                        Console.WriteLine("The following k-mers have been found\n");
                        foreach (string kmer in FrequentPatterns)
                            Console.Write("{0} ", kmer);
            */


            Console.ReadKey();
        }
    }
}
