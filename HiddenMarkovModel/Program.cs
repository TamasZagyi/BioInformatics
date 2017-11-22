using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HiddenMarkovModel
{
    class States : List<string>
    {
        public States(string line)
        {
            Parse(line);
        }
        public States()
        {

        }
        public States(States other)
        {
            AddRange(other);
        }
        public void Parse(string line)
        {
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            foreach (string item in line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries))
                Add(item);
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < Count; ++i)
            {
                sb.AppendFormat("{0}", this[i]);
                if (i != Count - 1) sb.Append(" ");
            }

            return sb.ToString();
        }
      
    }

    class Transition : List<List<double>>
    {
        public Transition(States states)
        {
            _states = new States(states);
        }
        public double Pr(string A, string B)
        {
            return this[_states.IndexOf(A)][_states.IndexOf(B)];
        }
        public void Parse(List<string> matrix)
        {
            if (matrix.Count != _states.Count + 1) throw new Exception("Dimension missmatch!");
            Clear();
            for (int i = 1; i < matrix.Count; ++i)
            {
                List<double> row = new List<double>();
                char[] delimiterChars = { ' ', ',', ':', '\t' };
                string[] items = matrix[i].Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                for (int j = 1; j < items.Length; ++j)
                {
                    row.Add(Double.Parse(items[j]));
                }
                Add(row);
            }
        }
        public void Estimate(string path)
        {
            for(int i=0;i<_states.Count;++i)
            {
                List<double> row = new List<double>();
                for (int j = 0; j < _states.Count; ++j)
                {
                    row.Add(0);
                }
                Add(row);
            }

            for(int i=1;i<path.Length;++i)
            {
                this[_states.IndexOf(path[i - 1].ToString())][_states.IndexOf(path[i].ToString())] += 1;
            }

            for (int i = 0; i < _states.Count; ++i)
            {
                double sum = 0;
                for (int j = 0; j < _states.Count; ++j)
                {
                    sum += this[i][j];
                }
                
                for (int j = 0; j < _states.Count; ++j)
                {
                    if (sum != 0) this[i][j] /= sum;
                    else this[i][j] = 1.0 / (double)_states.Count;
                }
            }
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < Count; ++i)
            {
                sb.AppendFormat("\t{0}",_states[i]);
            }
            sb.AppendLine();
            for (int i = 0; i < Count; ++i)
            {
                sb.AppendFormat("{0}", _states[i]);
                for (int j = 0; j < this[i].Count; ++j)
                {
                    sb.AppendFormat("\t{0:0.000}", this[i][j]);
                } 
                sb.AppendLine();
            }

            return sb.ToString();
        }
        private States _states;
    }

    class Alphabet : List<char>
    {
        public Alphabet(string line)
        {
            Parse(line);
        }
        public Alphabet()
        {

        }
        public Alphabet(Alphabet other)
        {
            AddRange(other);
        }
        public void Parse(string line)
        {
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            foreach (string item in line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries))
                Add(item[0]);
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < Count; ++i)
            {
                sb.AppendFormat("{0}", this[i]);
                if (i != Count - 1) sb.Append(" ");
            }

            return sb.ToString();
        }

    }

    class Emission : List<List<double>>
    {
        public Emission(States states, Alphabet alphabet)
        {
            _states = new States(states);
            _alphabet = new Alphabet(alphabet);
        }
        public double Pr(string state, char letter)
        {
            return this[_states.IndexOf(state)][_alphabet.IndexOf(letter)];
        }
        public void Parse(List<string> matrix)
        {
            if (matrix.Count != _states.Count + 1) throw new Exception("Dimension missmatch!");
            Clear();
            for (int i = 1; i < matrix.Count; ++i)
            {
                List<double> row = new List<double>();
                char[] delimiterChars = { ' ', ',', ':', '\t' };
                string[] items = matrix[i].Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
                if(items.Length != _alphabet.Count + 1) throw new Exception("Dimension missmatch!");
                for (int j = 1; j < items.Length; ++j)
                {
                    row.Add(Double.Parse(items[j]));
                }
                Add(row);
            }
        }
        public void Estimate(string outcome, string path)
        {
            for (int i = 0; i < _states.Count; ++i)
            {
                List<double> row = new List<double>();
                for (int j = 0; j < _alphabet.Count; ++j)
                {
                    row.Add(0);
                }
                Add(row);
            }

            for (int i = 0; i < outcome.Length; ++i)
            {
                this[_states.IndexOf(path[i].ToString())]
                    [_alphabet.IndexOf(outcome[i])] += 1;
            }

            for (int i = 0; i < _states.Count; ++i)
            {
                double sum = 0;
                for (int j = 0; j < _alphabet.Count; ++j)
                {
                    sum += this[i][j];
                }
                for (int j = 0; j < _alphabet.Count; ++j)
                {
                    if (sum != 0) this[i][j] /= sum;
                    else this[i][j] = 1.0 / (double)_alphabet.Count;
                }
            }
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < _alphabet.Count; ++i)
            {
                sb.AppendFormat("\t{0}", _alphabet[i]);
            }
            sb.AppendLine();
            for (int i = 0; i < _states.Count; ++i)
            {
                sb.AppendFormat("{0}", _states[i]);
                for (int j = 0; j < _alphabet.Count; ++j)
                {
                    sb.AppendFormat("\t{0:0.000}", this[i][j]);
                }
                sb.AppendLine();
            }

            return sb.ToString();
        }
        private States _states;
        private Alphabet _alphabet;
    }

    class HiddenMarkovModel
    {
        public void ParseStates(string line)
        {
            _states = new States(line);
        }
        public void ParseAlphabet(string line)
        {
            _alphabet = new Alphabet(line);
        }
        public void ParseTransition(List<string> matrix)
        {
            _transition = new Transition(_states);
            _transition.Parse(matrix);
        }
        public void ParseEmission(List<string> matrix)
        {
            _emission = new Emission(_states, _alphabet);
            _emission.Parse(matrix);
        }
        public double PrPath(string path)
        {
            double pr = 1.0 / _states.Count;
            for(int i = 1; i < path.Length; ++i)
            {
                pr *= _transition.Pr(path[i - 1].ToString(), path[i].ToString());
            }
            return pr;
        }
        public double PrOutcomeGivenPath(string outcome, string path)
        {
            double pr = 1.0;
            for (int i = 0; i < path.Length; ++i)
            {
                pr *= _emission.Pr(path[i].ToString(), outcome[i]);
            }
            return pr;
        }
        public void EstimateParameters(string outcome, string path)
        {
            _transition = new Transition(_states);
            _transition.Estimate(path);

            _emission = new Emission(_states, _alphabet);
            _emission.Estimate(outcome,path);
        }

        public string DecodePath(string outcome, ref double prob)
        {
            double[,] s = new double[_states.Count, outcome.Length];
            double[,] forward = new double[_states.Count, outcome.Length];
            double[,] backward = new double[_states.Count, outcome.Length];
            double[,] pr = new double[_states.Count, outcome.Length];
            int[,] backtrack = new int[_states.Count, outcome.Length];

            for (int i = 0; i < _states.Count; ++i)
            {
                s[i, 0] = 1.0 / _states.Count * _emission.Pr(_states[i], outcome[0]);
                forward[i, 0] = s[i, 0];
                backward[i, outcome.Length-1] = 1.0 / _states.Count * _emission.Pr(_states[i], outcome[outcome.Length - 1]);
            }

            for (int j = 1; j < outcome.Length; ++j)
                for (int i = 0; i < _states.Count; ++i)
                {
                    int i_max = -1;
                    double d_max = 0.0;

                    for(int k = 0; k < _states.Count; ++k)
                    {
                        double weight = _transition.Pr(_states[k], _states[i]) * _emission.Pr(_states[i], outcome[j]);
                        double d = s[k, j - 1] * weight;
                        forward[i, j] += forward[k, j - 1] * weight;
                        if( d > d_max )
                        {
                            d_max = d;
                            i_max = k;
                        }
                    }
                    s[i, j] = d_max;
                    backtrack[i, j - 1] = i_max;
                }

            for (int j = outcome.Length - 2; j >= 0; --j)
                for (int i = 0; i < _states.Count; ++i)
                {
                    for (int k = 0; k < _states.Count; ++k)
                    {
                        double weight = _transition.Pr(_states[i], _states[k]) * _emission.Pr(_states[k], outcome[j+1]);
                        backward[i, j] += backward[k, j + 1] * weight;
                    }
                }

            double d_last = 0;
            int i_last = -1;
            prob = 0;
            for(int i = 0; i < _states.Count; ++i)
            {
                if(s[i,outcome.Length-1] > d_last)
                {
                    d_last = s[i, outcome.Length - 1];
                    i_last = i;
                }
                prob += forward[i, outcome.Length - 1];
            }

            for (int j = 0; j < outcome.Length; ++j)
                for (int i = 0; i < _states.Count; ++i)
                {
                    pr[i, j] = forward[i, j] * backward[i, j] / prob * 4.3469;
                }

            StringBuilder sb = new StringBuilder();
            
            for(int i = outcome.Length - 2; i >= -1; --i)
            {
                sb.Append(_states[i_last]);
                if (i == -1) break;
                i_last = backtrack[i_last, i];
            }

            string result = sb.ToString();
            char[] charArray = result.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }
        class ProfileState : List<char>
        {
            public int CountOf(char c)
            {
                int count = 0;
                for(int i=0;i<Count;++i)
                {
                    if (this[i] == c) ++count;
                }
                return count;
            }
            public double Pr(char c)
            {
                return (double)CountOf(c) / (double)Count;
            }
            public bool gray;
        }
        public void BuildProfile(List<string> alignment, double theta)
        {
            //  build profile matrix
            List<ProfileState> profile = new List<ProfileState>();
            for(int i=0;i < alignment[0].Length;++i)
            {
                ProfileState ps = new ProfileState();
                for(int j = 0; j < alignment.Count;++j)
                {
                    ps.Add(alignment[j][i]);
                }
                ps.gray = ps.Pr('-') > theta;
                profile.Add(ps);
            }

            //  build states
            _states.Add("S");
            _states.Add("I0");
            int m = 1;
            for(int i = 0; i < profile.Count; ++i)
            {
                if(!profile[i].gray)
                {     
                    _states.Add("M" + m.ToString());
                    _states.Add("D" + m.ToString());
                    _states.Add("I" + m.ToString());

                    ++m;
                }
            }
            _states.Add("E");

            //  build transition matrix
            _transition = new Transition(_states);
            
            //  build emission matrix
            _emission = new Emission(_states, _alphabet);
        } 
        private States _states;
        private Alphabet _alphabet;
        public Transition _transition;
        public Emission _emission;
    }
    class Program
    {
        static void ReadStates(StreamReader input, HiddenMarkovModel HMM)
        {
            string states = input.ReadLine().Trim();
            HMM.ParseStates(states);
            input.ReadLine();
        }
        static void ReadAlphabet(StreamReader input, HiddenMarkovModel HMM)
        {
            string alphabet = input.ReadLine().Trim();
            HMM.ParseAlphabet(alphabet);
            input.ReadLine();
        }
        static void ReadTransition(StreamReader input, HiddenMarkovModel HMM)
        {
            List<string> transition = new List<string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine().Trim();
                if (line.Contains('-')) break;
                transition.Add(line);
            }
            HMM.ParseTransition(transition);
        }
        static void ReadEmission(StreamReader input, HiddenMarkovModel HMM)
        {
            List<string> emission = new List<string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine().Trim();
                if (line.Contains('-')) break;
                emission.Add(line);
            }
            HMM.ParseEmission(emission);
        }
        static void SolveProbabilityOfPath()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");

            string path = input.ReadLine().Trim();
            input.ReadLine();

            HiddenMarkovModel HMM = new HiddenMarkovModel();

            ReadStates(input, HMM);
            ReadTransition(input, HMM);

            StreamWriter output = new StreamWriter("result.txt");
            output.Write("{0} ", HMM.PrPath(path));
            output.Flush();
        }
        static void SolveProbabilityOfOutcomeGivenPath()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();
            ReadAlphabet(input, HMM);

            string path = input.ReadLine().Trim();
            input.ReadLine();
            ReadStates(input, HMM);

            ReadEmission(input, HMM);

            StreamWriter output = new StreamWriter("result.txt");
            output.Write("{0} ", HMM.PrOutcomeGivenPath(outcome, path));
            output.Flush();
        }
        static void SolveDecodingProblem()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();

            ReadAlphabet(input, HMM);
            ReadStates(input, HMM);
            ReadTransition(input, HMM);
            ReadEmission(input, HMM);

            StreamWriter output = new StreamWriter("result.txt");
            double prob = 0;
            output.Write("{0} ", HMM.DecodePath(outcome, ref prob));
            output.Flush();
        }
        static void SolveDecodingProbability()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();

            ReadAlphabet(input, HMM);
            ReadStates(input, HMM);
            ReadTransition(input, HMM);
            ReadEmission(input, HMM);

            StreamWriter output = new StreamWriter("result.txt");
            double prob = 0;
            HMM.DecodePath(outcome, ref prob);
            output.Write("{0}", prob);
            output.Flush();
        }
        static void SolveProfileHMM()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();
            double theta = Double.Parse(input.ReadLine().Trim());
            input.ReadLine();
            ReadAlphabet(input, HMM);
            List<string> alignment = new List<string>();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine().Trim();
                if (line.Length == 0) break;
                alignment.Add(line);
            }
            HMM.BuildProfile(alignment, theta);
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(HMM._transition.ToString());
            output.Write("--------");
            output.Write(HMM._emission.ToString());
            output.Flush();
        }
        static void SolveParameterEstimation()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();
            ReadAlphabet(input, HMM);

            string path = input.ReadLine().Trim();
            input.ReadLine();
            ReadStates(input, HMM);

            HMM.EstimateParameters(outcome, path);
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(HMM._transition.ToString());
            output.WriteLine("--------");
            output.Write(HMM._emission.ToString());
            output.Flush();
        }
        static void SolveViterbiLearning()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            int iterations = Int32.Parse(input.ReadLine().Trim());
            input.ReadLine();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();

            ReadAlphabet(input, HMM);
            ReadStates(input, HMM);
            ReadTransition(input, HMM);
            ReadEmission(input, HMM);

            for(int i = 0; i < iterations; ++i)
            {
                double prob = 0.0;
                string path = HMM.DecodePath(outcome, ref prob);
                HMM.EstimateParameters(outcome, path);
            }

            StreamWriter output = new StreamWriter("result.txt");
            output.Write(HMM._transition.ToString());
            output.WriteLine("--------");
            output.Write(HMM._emission.ToString());
            output.Flush();
        }
        static void SolveSoftDecoding()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            HiddenMarkovModel HMM = new HiddenMarkovModel();

            string outcome = input.ReadLine().Trim();
            input.ReadLine();

            ReadAlphabet(input, HMM);
            ReadStates(input, HMM);
            ReadTransition(input, HMM);
            ReadEmission(input, HMM);

            StreamWriter output = new StreamWriter("result.txt");
            double prob = 0;
            HMM.DecodePath(outcome, ref prob);
            output.Write("{0}", prob);
            output.Flush();
        }
        static void Main(string[] args)
        {
            SolveSoftDecoding();
        }
    }
}
