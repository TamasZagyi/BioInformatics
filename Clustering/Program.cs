using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Clustering
{
    class Program
    {
        class Point : List<double>
        {
            public Point(string line)
            {
                Parse(line);
            }
            public Point()
            {
                
            }
            public void Add(Point o)
            {
                while (Count < o.Count) Add(0.0);
                for(int i=0;i<Count;++i)
                {
                    this[i] += o[i];
                }
            }
            public void Divide(double c)
            {
                for (int i = 0; i < Count; ++i)
                {
                    this[i] /= c;
                }
            }
            public void Parse(string line)
            {
                line = line.Replace("+", "");
                line = line.Replace("(", " ");
                line = line.Replace(")", " ");

                char[] delimiterChars = { ' ', ',', ':', '\t' };

                foreach (string item in line.Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries))
                    Add(Double.Parse(item));
            }
            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();

                for (int i = 0; i < Count; ++i)
                {
                    sb.AppendFormat("{0}", this[i].ToString("F3"));
                    if (i != Count - 1) sb.Append(" ");
                }

                return sb.ToString();
            }
            public double d(Point o)
            {
                if (o.Count != Count) throw new Exception("Cannot handle different dimensions");

                double dist = 0.0;
                for (int i = 0; i < Count; ++i) dist += (this[i] - o[i]) * (this[i] - o[i]);

                return Math.Sqrt(dist);
            }
            public double d(Cluster cluster)
            {
                double min_dist = Double.MaxValue;

                foreach(Point o in cluster)
                {
                    if (d(o) < min_dist) min_dist = d(o);
                }

                return min_dist;
            }
            public double d(Cluster cluster, ref int index)
            {
                double min_dist = Double.MaxValue;

                for (int i = 0; i < cluster.Count;++i)
                {
                    if (d(cluster[i]) < min_dist)
                    {
                        min_dist = d(cluster[i]);
                        index = i;
                    }
                }

                return min_dist;
            }
        }
        class Cluster : List<Point>
        {
            public double MaxDistance(Cluster cluster, ref int max_index)
            {
                max_index = -1;
                double max_dist = Double.MinValue;

                for (int i=0;i<Count;++i)
                {
                    double dist = this[i].d(cluster);
                    if (dist > max_dist)
                    {
                        max_dist = dist;
                        max_index = i;
                    }
                }

                return max_dist;
            }
            public double Distortion(Cluster centers)
            {
                double distortion = 0.0;

                for (int i = 0; i < Count; ++i)
                {
                    double d = this[i].d(centers);
                    distortion += d * d;
                }

                return distortion / Count;
            }
            public List<Cluster> CentersToClusters(Cluster centers)
            {
                List<Cluster> clusters = new List<Cluster>();
                foreach (Point c in centers) clusters.Add(new Cluster());
                for (int i = 0; i < Count; ++i)
                {
                    int min_index = -1;
                    this[i].d(centers,ref min_index);
                    clusters[min_index].Add(this[i]);
                }

                return clusters;
            }
            public Cluster ClustersToCenters(List<Cluster> clusters)
            {
                Cluster centers = new Cluster();
                
                for (int i = 0; i < clusters.Count; ++i)
                {
                    centers.Add(clusters[i].CenterOfGravity());
                }

                return centers;
            }
            public double[,] CentersToClustersSoft(Cluster centers, double stiffness)
            {
                int k = centers.Count;
                int n = Count;

                double [,] HiddenMatrix = new double[k, n];

                for(int i = 0; i < k;++i)
                {
                    for(int j = 0; j < n;++j)
                    {
                        double div = 0.0;
                        for(int ii = 0; ii < centers.Count; ++ii)
                        {
                            div += Math.Pow(Math.E, -stiffness * this[j].d(centers[ii]));
                        }

                        HiddenMatrix[i,j] = Math.Pow(Math.E, -stiffness * this[j].d(centers[i])) / div;
                    }
                }

                return HiddenMatrix;
            }
            public Cluster ClustersToCentersSoft(double[,] HiddenMatrix)
            {
                Cluster centers = new Cluster();

                for (int i = 0; i < HiddenMatrix.GetLength(0); ++i)
                {
                    Point center = new Point();
                    for(int j=0;j<this[0].Count;++j)
                    {
                        double counter = 0.0;
                        double divisor = 0.0;
                        for(int ii=0;ii<Count;++ii)
                        {
                            counter += HiddenMatrix[i, ii] * this[ii][j];
                            divisor += HiddenMatrix[i, ii];
                        }
                        center.Add(counter / divisor);
                    }
                    centers.Add(center);
                }

                return centers;
            }
            public void Parse(string line)
            {
                Add(new Point(line));
            }
            public override string ToString()
            {
                StringBuilder sb = new StringBuilder();

                for (int i = 0; i < Count; ++i)
                {
                    sb.AppendFormat("{0}", this[i].ToString());
                    if (i != Count - 1) sb.AppendLine();
                }

                return sb.ToString();
            }
            /*
            FarthestFirstTraversal(Data, k) 
             Centers ← the set consisting of a single randomly chosen point from Data
              while |Centers| < k 
               DataPoint ← the point in Data maximizing d(DataPoint, Centers) 
               add DataPoint to Centers 
             return Centers 
            */
            public Cluster FarthestFirstTraversal(int k)
            {
                Cluster Centers = new Cluster();
                Centers.Add(this[0]);
                while(Centers.Count < k)
                {
                    int max_index = -1;
                    MaxDistance(Centers, ref max_index);
                    Centers.Add(this[max_index]);
                }

                return Centers;
            }
            public Cluster Lloyd(int k)
            {
                Cluster Centers = new Cluster();
                for(int i=0;i<k;++i) Centers.Add(this[i]);

                for(int i=0;i<20*k;++i)
                {
                    List<Cluster> clusters = CentersToClusters(Centers);
                    Centers = ClustersToCenters(clusters);
                }

                return Centers;
            }
            public Cluster SoftLloyd(int k, double stiffness)
            {
                Cluster Centers = new Cluster();
                for (int i = 0; i < k; ++i) Centers.Add(this[i]);

                for (int i = 0; i < 100; ++i)
                {
                    double[,] HiddenMatrix = CentersToClustersSoft(Centers, stiffness);
                    Centers = ClustersToCentersSoft(HiddenMatrix);
                }

                return Centers;
            }
            public Point CenterOfGravity()
            {
                Point c = new Point();
                foreach(Point p in this)
                {
                    c.Add(p);
                }
                c.Divide(Count);
                return c;
            }
        } 
        static void SolveFarthestFirstTraversal()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int k = Int32.Parse(parameter_string[0]);
            int m = Int32.Parse(parameter_string[1]);
            Cluster Data = new Cluster();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    Data.Parse(line);
                }
            }

            //  compute
            Cluster Centers = Data.FarthestFirstTraversal(k);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(Centers.ToString());
            output.Flush();
        }
        static void SolveDistortion()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int k = Int32.Parse(parameter_string[0]);
            int m = Int32.Parse(parameter_string[1]);
            Cluster Centers = new Cluster();
            for (int i = 0; i < k; ++i) Centers.Parse(input.ReadLine());
            input.ReadLine();   //  eat separator
            Cluster Data = new Cluster();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    Data.Parse(line);
                }
            }

            //  compute
            double distortion = Data.Distortion(Centers);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(distortion.ToString());
            output.Flush();
        }
        static void SolveLloyd()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int k = Int32.Parse(parameter_string[0]);
            int m = Int32.Parse(parameter_string[1]);
            Cluster Data = new Cluster();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    Data.Parse(line);
                }
            }

            //  compute
            Cluster Centers = Data.Lloyd(k);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(Centers.ToString());
            output.Flush();
        }
        static void SolveSoftLloyd()
        {
            //  read data
            StreamReader input = new StreamReader("data.txt");
            char[] delimiterChars = { ' ', ',', ':', '\t' };
            string[] parameter_string = input.ReadLine().Split(delimiterChars, StringSplitOptions.RemoveEmptyEntries);
            int k = Int32.Parse(parameter_string[0]);
            int m = Int32.Parse(parameter_string[1]);
            double stiffness = Double.Parse(input.ReadLine());
            Cluster Data = new Cluster();
            while (!input.EndOfStream)
            {
                string line = input.ReadLine();
                if (line.Length > 0)
                {
                    Data.Parse(line);
                }
            }

            //  compute
            Cluster Centers = Data.SoftLloyd(k, stiffness);
            //  write result
            StreamWriter output = new StreamWriter("result.txt");
            output.Write(Centers.ToString());
            output.Flush();
        }
        static void Main(string[] args)
        {
            SolveSoftLloyd();
     /*       Cluster data = new Cluster();
            data.Add(new Point("17 0 -4"));
            data.Add(new Point("3 14 23"));
            data.Add(new Point("9 7 16"));
            data.Add(new Point("7 3 5"));

            Cluster centers = new Cluster();
            centers.Add(new Point("4 5"));
            centers.Add(new Point("7 4"));

            Point c = data.CenterOfGravity();
            Console.WriteLine( c.ToString()  );

            Console.ReadKey();*/
        }
    }
}
