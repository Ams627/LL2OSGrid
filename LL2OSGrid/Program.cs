using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LL2OSGrid
{
    internal class Program
    {
        private static void Main(string[] args)
        {
            try
            {
                var loc = new CvLocation();
                loc.SetLatLongDegrees(50.715812, -2.011931);
                var (eastings, northings) = loc.OsGrid;
                Console.WriteLine($"{eastings} {northings}");
            }
            catch (Exception ex)
            {
                var codeBase = System.Reflection.Assembly.GetExecutingAssembly().CodeBase;
                var progname = Path.GetFileNameWithoutExtension(codeBase);
                Console.Error.WriteLine(progname + ": Error: " + ex.Message);
            }

        }
    }
}
