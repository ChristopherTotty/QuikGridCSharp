//http://www.galiander.ca/quikgrid/index.html
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuickGrid
{
    public class Xpand
    {
        // Shared Variables. 
        static int NumXcoords,
               NumYcoords,
               xIndex,
               yIndex,
               PercentEdgeFactor = 100,
               PercentDensityRatio = 150;

        static float[] PointDistSquared = new float[8];   // Distance**2 from point.
        static float ClosestSquared, EdgeSenseDistSq, Radius, Diameter, DiameterSq,
                    xyStopDistSq, xGridMax, yGridMax, xGridMin, yGridMin;
        static float ScanStopRatio = 16.0f;
        static float DensityStopRatio = 1.5f;
        static float EdgeSenseFactor = 1.0f;
        static float UndefinedZ = -99999;
        
        static long[] PointInOctant = new long[8];  // Number of point closest in Octant
        static long TotalRange,
                 NumberDone,
                 TotalShells,
                 NumDataPoints,
                 NoFound,
                 TotalGridPoints,
                 nSample = 1;          // Sample every n'th point.

        static GridXType Locator = new GridXType(0, 0, 0);


        //**************************************************************
        //        X p a n d
        //                     Generate the entire grid in one call. 
        //*************************************************************
        public Xpand(ref SurfaceGrid Zgrid, ref ScatteredData RandomData)
        {
            XpandInit(ref Zgrid, ref RandomData);
            while (XpandPoint(ref Zgrid, ref RandomData)) ;
        }

        //******************************************************************
        //                    X p a n d   P o i n t
        //
        //  Evaluates one grid intersection only and advances indexes so
        //       the next call will evaluate the next one. 
        //******************************************************************
        bool XpandPoint(ref SurfaceGrid Zgrid, ref ScatteredData RandomData)
        {
            float Zvalue;
            //assert( xIndex < NumXcoords );
            //assert( yIndex < NumYcoords );
            if (NumDataPoints < 3) return false;

            // Select the closest point in each octant surround the grid point.
            SelectPoints(ref RandomData, ref Zgrid);

            // Check if point can be included and if so.. calculate it. 
            if (IncludeGridPoint() > 0) Zvalue = WeightedAverage(ref RandomData);
            else Zvalue = UndefinedZ;

            Zgrid.zset(xIndex, yIndex, Zvalue);

            // Move to next grid intersection....
            ++xIndex;
            if (xIndex >= NumXcoords) { ++yIndex; xIndex = 0; }

            if (yIndex >= NumYcoords) { Locator.New(0, NumXcoords, NumYcoords); return false; }
            ++NumberDone;
            //assert( NumberDone <= TotalGridPoints);

            return true;
        }

        //***************************************************************
        //     Set and Query  Xpand  Status  and parameters
        //**************************************************************

        // The following three functions are only of interest if you spin out
        // a separate thread to generate the grid.

        int XpandPercentDone() { return Convert.ToInt32((NumberDone * 100) / TotalGridPoints) ; }

        int XpandBandWidth()
        {
            if (NumberDone == 0) return 0;
            return Convert.ToInt32(TotalRange / NumberDone);
        }

        int XpandPercentShell()
        {
            if (NumberDone == 0) return 0;
            return Convert.ToInt32(TotalShells / NumberDone);
        }

        // The following allow you to control the grid generation parameters
        // See the QuikGrid help file for an explanation of what they do.

        int XpandScanRatio() { return (int)ScanStopRatio; }
        int XpandScanRatio( int NewRatio )
        {
            float OldRatio = ScanStopRatio;
            ScanStopRatio=NewRatio;
            if( ScanStopRatio< 1.0 ) ScanStopRatio = 1.0f;
            if( ScanStopRatio > 100.0 ) ScanStopRatio = 100.0f;
            return (int) OldRatio;
        }

        int XpandDensityRatio() { return PercentDensityRatio; }
        int XpandDensityRatio(  int NewRatio)
        {
            int OldRatio = PercentDensityRatio;
            PercentDensityRatio = NewRatio;
            if (PercentDensityRatio < 1) PercentDensityRatio = 1;
            if (PercentDensityRatio > 10000) PercentDensityRatio = 10000;
            DensityStopRatio = (float)PercentDensityRatio * 0.01f;
            return OldRatio;
        }

        int XpandEdgeFactor() { return PercentEdgeFactor; }
        int XpandEdgeFactor(int NewFactor)
        {
            int OldFactor = PercentEdgeFactor;
            PercentEdgeFactor = NewFactor;
            if (PercentEdgeFactor < 1) PercentEdgeFactor = 1;
            if (PercentEdgeFactor > 10000) PercentEdgeFactor = 10000;
            EdgeSenseFactor = (float)PercentEdgeFactor * 0.01f;
            return OldFactor;
        }
        float XpandUndefinedZ() { return UndefinedZ; }
        float XpandUndefinedZ(float z)
        {
            float OldZ = UndefinedZ;
            UndefinedZ = z;
            return OldZ;
        }
        long XpandSample() { return nSample; }
        long XpandSample( long i)
        {
            long oldSample = nSample;
            nSample = i;
            if (nSample < 1) nSample = 1;
            return oldSample;
        }

        //***************************************************************
        //            X p a n d   I n i t
        //***************************************************************
        void XpandInit(ref SurfaceGrid Zgrid, ref ScatteredData Data)
        {
            NumXcoords = Convert.ToInt32(Zgrid.xsize);
            NumYcoords = Convert.ToInt32(Zgrid.ysize);

            TotalGridPoints = (long)NumXcoords * (long)NumYcoords;
            xIndex = yIndex = 0;
            NumberDone = TotalRange = TotalShells = 0;
            NumDataPoints = Data.Size;
            if (NumDataPoints < 3) return;

            xGridMin = Zgrid.x(0);
            xGridMax = Zgrid.x(NumXcoords - 1);
            yGridMin = Zgrid.y(0);
            yGridMax = Zgrid.y(NumYcoords - 1);
            float xRange = Data.xMax - Data.xMin;
            float yRange = Data.yMax - Data.yMin;
            float Volume = xRange * yRange;
            float VolPerPoint = Volume / ((float)NumDataPoints / (float)nSample);
            Radius = (float)Math.Sqrt(VolPerPoint / Math.PI);
            Diameter = Radius * 2;
            DiameterSq = Diameter * Diameter;
            xyStopDistSq = DiameterSq * DensityStopRatio * DensityStopRatio;
            EdgeSenseDistSq = DiameterSq * EdgeSenseFactor * EdgeSenseFactor;

            LocateInit(ref Data, ref Zgrid);

        }
        //**********************************************************************
        //             S c a n    O n e   G r i d
        //**********************************************************************
        static void ScanOneGrid(int i, int j, ref ScatteredData Data, ref SurfaceGrid Grid)
        {
            long GotOne;
            for (int n = 0; true; n++)
            {
                GotOne = Locator.Search(i, j, n);
                if (GotOne < 0) return;
                PutInOctant(Grid.x(xIndex), Grid.y(yIndex), ref Data, GotOne);
                TotalRange++;
            }
        }

        //****************************************************************
        //                 P u t   I n    O c t a n t
        //****************************************************************
        static void PutInOctant( float x, float y, ref ScatteredData Data, long DataIndex)
        {  
           float Xdistance = Data.X[DataIndex] - x;
                float Ydistance = Data.Y[DataIndex] - y;

                // Select the octant the data point is in. They are labeled:
                //        31
                //       2  0
                //       6  4
                //        75

                int Octant = 0;
           if( Ydistance< 0.0 ) Octant = 4;
           if (Xdistance< 0.0 ) Octant += 2;
           if (Math.Abs(Xdistance) < Math.Abs(Ydistance)) Octant += 1;

           // and get the distance squared and save that information
           // for that Octant iff this is the first point found or closer than
           // the point previously saved for that octant.

           float DistSquared = Xdistance * Xdistance + Ydistance * Ydistance;
           if( NoFound == 0 ) ClosestSquared = DistSquared;
           if( DistSquared<ClosestSquared) ClosestSquared = DistSquared;
           NoFound++;

           if ( PointInOctant[Octant] == -1 ||
	        DistSquared<PointDistSquared[Octant] )
	        {
	          PointInOctant[Octant] = DataIndex ;
	          PointDistSquared[Octant] = DistSquared;
	        }
        }


        // ***************************************************************
        //                    S e l e c t    P o i n t s
        //*****************************************************************
        static void SelectPoints(ref ScatteredData Data, ref SurfaceGrid Grid)

        // This routine will search the array of Data points looking for
        // the closest point in each of the 8 octants surrounding the grid
        // coordinate currently being evaluated (at location x, y).

        {
            int Start,
                End,
                Row,
                Column;

            float TestDist;

            bool TopDone = false;      // These are logical flags controlling the 
            bool BottomDone = false;   // Shelling process. 
            bool LeftDone = false;
            bool RightDone = false;
            int i;

            // Zero out the arrays which keep track of closest point and its distance.
            for (i = 0; i < 8; i++)
            {
                PointInOctant[i] = -1;
                PointDistSquared[i] = 0.0f;
            }
            NoFound = 0;

            ScanOneGrid(xIndex, yIndex, ref Data, ref Grid); // Do home grid first.

            for (int shell = 1; true; shell++) // Now shell outwards from home. 
            {
                Start = xIndex - shell; if (Start < 0) Start = 0; // Thanks to Matt Gessner
                End = xIndex + shell + 1; if (End > NumXcoords) End = NumXcoords;
                // Do top row.
                if (!TopDone)
                {
                    Row = yIndex + shell;
                    if (Row >= NumYcoords) TopDone = true;
                    else
                    {
                        TestDist = Grid.y(Row) - Grid.y(yIndex);
                        TestDist = TestDist * TestDist;
                        if (((NoFound > 0) && ((TestDist > ClosestSquared * ScanStopRatio)) ||
                          (TestDist > xyStopDistSq))) TopDone = true;
                        else
                            for (i = Start; i < End; i++) ScanOneGrid(i, Row, ref Data, ref Grid);
                    }
                }

                // Do bottom row.
                if (!BottomDone)
                {
                    Row = yIndex - shell;
                    if (Row < 0) BottomDone = true;
                    else
                    {
                        TestDist = Grid.y(yIndex) - Grid.y(Row);
                        TestDist = TestDist * TestDist;
                        if (((NoFound > 0) && ((TestDist > ClosestSquared * ScanStopRatio)) ||
                            (TestDist > xyStopDistSq))) BottomDone = true;
                        else
                            for (i = Start; i < End; i++) ScanOneGrid(i, Row, ref Data, ref Grid);
                    }
                }

                Start = yIndex - shell + 1; if (Start < 0) Start = 0; // Thanks to Matt Gessner
                End = yIndex + shell; if (End > NumYcoords) End = NumYcoords;
                // Do left column.
                if (!LeftDone)
                {
                    Column = xIndex - shell;
                    if (Column < 0) LeftDone = true;
                    else
                    {
                        TestDist = Grid.x(xIndex) - Grid.x(Column);
                        TestDist = TestDist * TestDist;
                        if (((NoFound > 0) && ((TestDist > ClosestSquared * ScanStopRatio)) ||
                            (TestDist > xyStopDistSq))) LeftDone = true;
                        else
                            for (i = Start; i < End; i++) ScanOneGrid(Column, i, ref Data, ref Grid);
                    }
                }
                // Do right column.
                if (!RightDone)
                {
                    Column = xIndex + shell;
                    if (Column >= NumXcoords) RightDone = true;
                    else
                    {
                        TestDist = Grid.x(Column) - Grid.x(xIndex);
                        TestDist = TestDist * TestDist;
                        if (((NoFound > 0) && ((TestDist > ClosestSquared * ScanStopRatio)) ||
                            (TestDist > xyStopDistSq))) RightDone = true;
                        else
                            for (i = Start; i < End; i++) ScanOneGrid(Column, i, ref Data, ref Grid);
                    }
                }
                if (TopDone && BottomDone && LeftDone && RightDone)
                { TotalShells += shell; break; }
            }

        }

        //************************************************************
        //   I n c l u d e   G r i d   P o i n t
        //************************************************************
        static int IncludeGridPoint()
        {
            // Routine to decide whether to evaluate a Grid point.
            // Return 0 = reject; 1 = accept.

            int[] LookUpTable = new int[11] { 0, 1, 3, 2, 6, 7, 5, 4, 0, 1, 3 };

            if (NoFound == 0) return 0; // No points found.
            if (ClosestSquared > xyStopDistSq) return 0; // Closest too far away
            if (ClosestSquared < EdgeSenseDistSq) return 1; // 

            int NumConsecutiveOctants = 0; // code looks for 4 consecutive
            for (int i = 0; i < 11; i++)    // empty octants.
            {
                if (PointInOctant[LookUpTable[i]] > -1) NumConsecutiveOctants = 0;
                else
                {
                    NumConsecutiveOctants += 1;
                    if (NumConsecutiveOctants == 4) return 0;
                }
            }
            return 1;
        }

        //************************************************************
        //       W e i g h t e d    A v e r a g e
        //************************************************************
        static float WeightedAverage(ref ScatteredData Data)

        // function calculates the weighted average value for the
        // specified grid point. The needed data is all in PointInOctant
        // and PointDistSquared.


        {
            float SumSquared = 0.0f;
            float SumWeights = 0.0f;
            for (int i = 0; i < 8; i++)
            {
                long Point = PointInOctant[i];
                if (Point > -1)
                {
                    float DistanceSquared = PointDistSquared[i];
                    // Check for point right on a grid intersection.
                    if (DistanceSquared == 0.0) return Data.Z[Point];
                    SumSquared += 1.0f / DistanceSquared;
                    SumWeights += (1.0f / DistanceSquared) * Data.Z[Point];
                }
            }
            return SumWeights / SumSquared;
        }

        //**************************************************************
        //             L o c a t e    I n i t
        //**************************************************************
        static void LocateInit(ref ScatteredData Data, ref SurfaceGrid Grid)
        {
            int ix = 0, iy = 0;
            Locator.New(NumDataPoints, NumXcoords, NumYcoords);

            for (long i = 0; i < NumDataPoints; i += nSample)
            {
                if (LocateGridX(ref ix, ref iy, Data.X[i], Data.Y[i], ref Grid))
                    Locator.setnext(i, ix, iy);
            }
            Locator.Sort();
        }

        //**************************************************************
        //             L o c a t e   G r i d   X
        //**************************************************************
        static bool LocateGridX(ref int ix, ref int iy, float xLocn, float yLocn, ref SurfaceGrid Grid )
        {
            // Finds the closest grid intersection to data point i.
            // Return's 1 if successful,
            //          0 if data point is too far away from the intersection. 
            int test = 0,
                    nx, ny,
                    top,
                    bot;
            float distance;

            // check to see if point too far outside the grid perimeter.
            distance = xGridMin - xLocn;
            if (xLocn < xGridMin && distance * distance > xyStopDistSq) return false;
            distance = xGridMax - xLocn;
            if (xLocn > xGridMax && distance * distance > xyStopDistSq) return false;
            distance = yGridMin - yLocn;
            if (yLocn < yGridMin && distance * distance > xyStopDistSq) return false;
            distance = yGridMax - yLocn;
            if (yLocn > yGridMax && distance * distance > xyStopDistSq) return false;


            // Binary search for closest in x. 
            nx = NumXcoords;
            top = nx - 1;
            bot = 0;

            while ((top - bot) > 1)
            {
                test = (top - bot) / 2 + bot;
                if (xLocn > Grid.x(test)) bot = test;
                else top = test;
            }

            float dist = Math.Abs(xLocn - Grid.x(test));
            if (test < (nx - 1))
                if (Math.Abs(xLocn - Grid.x(test + 1)) < dist) test++;
            if (test > 0)
                if (Math.Abs(xLocn - Grid.x(test - 1)) < dist) test--;
            //assert( test >= 0 && test < nx );
            ix = test;

            // Do the same for closest in y. 
            ny = NumYcoords;
            top = ny - 1;
            bot = 0;
            while ((top - bot) > 1)
            {
                test = (top - bot) / 2 + bot;
                if (yLocn > Grid.y(test)) bot = test;
                else top = test;
            }
            dist = Math.Abs(yLocn - Grid.y(test));
            if (test < (ny - 1))
                if (Math.Abs(yLocn - Grid.y(test + 1)) < dist) test++;
            if (test > 0)
                if (Math.Abs(yLocn - Grid.y(test - 1)) < dist) test--;
            //assert( test >= 0 && test < ny);
            iy = test;
            return true;
        }
    }

    public class ScatteredData
    {
        public float xMax, xMin, yMax, yMin, zMax, zMin;
        public float[] X;
        public float[] Y;
        public float[] Z;
        public int Size;
        public int MaxPoints;

        public ScatteredData(int MaxSize)
        {
            MaxPoints = MaxSize;
            X = new float[MaxPoints];
            Y = new float[MaxPoints];
            Z = new float[MaxPoints];
        }

        public void SetNext (float x, float y, float z)
        {
            if (Size >= MaxPoints) { return; }

            X[Size] = x;
            Y[Size] = y;
            Z[Size] = z;

            if (Size == 0)
            {
                xMax = xMin = x;
                yMax = yMin = y;
                zMax = zMin = z;
            }
            else
            {
                if (x > xMax) xMax = x;
                if (x < xMin) xMin = x;
                if (y > yMax) yMax = y;
                if (y < yMin) yMin = y;
                if (z > zMax) zMax = z;
                if (z < zMin) zMin = z;
            }

            ++Size;
        }
    }

    public class SurfaceGrid
    {

        float[] zgrid;       // z[i][j]
        float[] xvect;       // x[i]
		float[] yvect;       // y[j]

        public long xsize, ysize;     // size of the array.

        public SurfaceGrid(int i, int j)
        {
            xsize = i;
            ysize = j;
            zgrid = new float[i * j];
            xvect = new float[xsize];
            yvect = new float[ysize];
            for (var k = 0; k < xsize; k++) { xvect[k] = k; }
            for (var k = 0; k < ysize; k++) { yvect[k] = k; }
        }

        public void zset (int i, int j, float a)
        {
            long offset = (long)j * xsize + (long)i;
            zgrid[offset] = a;
        }

        public void xset (int i, float a)
        {
            xvect[i] = a;
        }

        public void yset (int i, float a)
        {
            yvect[i] = a;
        }

        public float x(int i)
        {
            return xvect[i];
        }

        public float y(int i)
        {
            return yvect[i];
        }

        public float z(int i, int j)
        {
            long offset = (long)j * xsize + (long)i;
            return zgrid[offset];
        }
    }

    public class GridXType
    {
        private int Shift = 8192;
        private struct LocateStructure
        {
            public long intersection;
            public long DataLocation;
        }

        private LocateStructure[] FindPoints;
        long Size, nx, ny, LookupSize, PreviousSearch, np;
        long[] Lookup;

        public long x(long i)
        {
            return FindPoints[i].intersection / Shift;
        }

        public long y(long i)
        {
            return FindPoints[i].intersection % Shift;
        }

        public long location(long i)
        {
            return FindPoints[i].DataLocation;
        }

        public void setnext(long i, int ix, int iy)
        {
            FindPoints[np].DataLocation = i;
            FindPoints[np].intersection = (long)ix * (long)Shift + (long)iy;
            np++;
        }
        

        public GridXType(long iSize, long ax, long ay)
        {
            New(iSize, ax, ay);
        }

        public void New(long iSize, long ax, long ay)
        {
            if (iSize == Size && ax == nx && ay == ny) { return; }

            nx = ax; ny = ay;
            Size = iSize;
            np = 0;
            PreviousSearch = -1;
            if (Size <= 0) { return; }
            FindPoints = new LocateStructure[Size];
            LookupSize = nx * ny;
            Lookup = new long[LookupSize];
        }
        
        private static void qsort(LocateStructure[] elements, int left, int right)
        {
            int i = left, j = right;
            LocateStructure pivot = elements[(left + right) / 2];

            while (i <= j)
            {
                while (i < elements.Length - 1 && elements[i].intersection < pivot.intersection)
                {
                    i++;
                }

                while (j > 0 && elements[j].intersection > pivot.intersection)
                {
                    j--;
                }

                if (i <= j)
                {
                    LocateStructure tmp = elements[i];
                    elements[i] = elements[j];
                    elements[j] = tmp;

                    i++;
                    j--;
                }
            }
            
            if (left < j)
            {
                qsort(elements, left, j);
            }

            if (i < right)
            {
                qsort(elements, i, right);
            }
        }
        public void Sort()
        {
            long i, ix, iy, ixiy, ixiyold, ixiyLook;

            qsort(FindPoints, 0, FindPoints.Length - 1);

            // Build lookup table. -1 means no data points attached to that node.
            for (i = 0; i < LookupSize; i++) { Lookup[i] = -1; }

            ixiyold = -1;
            for (i = 0; i < np; i++)
            {
                ixiy = FindPoints[i].intersection;
                if (ixiy == ixiyold) continue;
                ixiyold = ixiy;
                ix = ixiy / Shift;
                iy = ixiy % Shift;
                ixiyLook = iy * nx + ix;
                Lookup[ixiyLook] = i;
            }
        }

        public long Search(int i, int j, int n)
        {
            long ixiy, newij, result;
            ixiy = (long)j * nx + (long)i;
            result = Lookup[ixiy];

            if (result == -1) return -1;
            newij = (long)i * (long)Shift + (long)j;

            result += n;
            if (result >= np) return -1;

            if (newij == FindPoints[result].intersection)
                return FindPoints[result].DataLocation;
            return -1;
        }
    }
}
