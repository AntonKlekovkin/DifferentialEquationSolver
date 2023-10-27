using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace DifferentialEquationSolver
{
    public partial class Form1 : Form
    {
        const int tMax = 10;
        const double stepDE = 0.01;
        const int numberPoints = (int)(tMax / stepDE) + 1;
        

        const double g = 9.8;

        double[] dx = new double[numberPoints];
        double[] dy = new double[numberPoints];
        double[] x = new double[numberPoints];
        double[] y = new double[numberPoints];
        double[] t = new double[numberPoints];

        int[] xRange = new int[2] { 0, tMax };
        int[] yRange = new int[2] { -7, 7 };

        int[] xRangeTr = new int[2] { -4, 4 };
        int[] yRangeTr = new int[2] { -4, 4 };

        int markerSizeForPlot = 2;
        int markerSizeForAnimation = 10;

        int stepPointForTimer;

        ChartArea c;
        //Series series;
        //Series seriesForTimer;

        Random r = new Random(1);

        double[][] randomNoise = new double[10][];
        double[][] xNoise = new double[5][];
        double[][] yNoise = new double[5][];
        double[] xMean = new double[numberPoints];
        double[] yMean = new double[numberPoints];

        double[] bigDeltaX = new double[numberPoints];
        double[] bigDeltaY = new double[numberPoints];

        int numberNameSeries = 0;

        SeriesChartType standartPlotType = SeriesChartType.Line;

        public Form1()
        {
            InitializeComponent();

            lbNumPoints.Text = "Points: " + numberPoints.ToString();
            lbTime.Text = "Time: " + tMax.ToString();
            lbStep.Text = "Step: " + stepDE.ToString();

            stepPointForTimer = (int)(numberPoints * (timer1.Interval / 1000.0) / tMax);

            progressBar1.Maximum = numberPoints / stepPointForTimer;

            for (int i = 0; i < randomNoise.Length; i++)
            {
                randomNoise[i] = new double[numberPoints];                
            }

            for (int i = 0; i < xNoise.Length; i++)
            {
                xNoise[i] = new double[numberPoints];
                yNoise[i] = new double[numberPoints];
            }
        }
        private double ConvertToDoubleFromNumeric(NumericUpDown num, double defaultValue)
        {
            double ret;

            try
            {
                ret = Convert.ToDouble(num.Value);
            }
            catch
            {
                ret = defaultValue;
            }

            return ret;
        }

        public void InitPlot(int[] xRange, int[] yRange)
        {
            chart1.Series.Clear();
            c = chart1.ChartAreas[0];
            c.AxisX.Minimum = xRange[0];
            c.AxisX.Maximum = xRange[1];

            c.AxisY.Minimum = yRange[0];
            c.AxisY.Maximum = yRange[1];

            Grid gridX = new Grid();
            gridX.Enabled = true;
            //gridX.Interval = 10;
            //gridX.IntervalOffsetType = DateTimeIntervalType.Seconds;
            c.AxisX.MajorGrid = gridX;

            Grid gridY = new Grid();
            gridY.Enabled = true;
            c.AxisY.MajorGrid = gridY;
            
            c.AxisX.Interval = (xRange[1] - xRange[0]) / 8;
            c.AxisY.Interval = (yRange[1] - yRange[0]) / 8;
            //c.AxisY.Crossing = 0;

            //series = chart1.Series.Add("Point");
            //series.ChartType = SeriesChartType.Line;
            //series.IsVisibleInLegend = false;
            //series.Color = Color.Red;
            //series.MarkerSize = markerSize;
            //series.MarkerStyle = MarkerStyle.Circle;
        }

        public void Plot(double[] x, double[] y, int markerSize, Color col, SeriesChartType type)
        {
            Series series;
            series = chart1.Series.Add("Point" + numberNameSeries++.ToString());
            series.ChartType = type; // SeriesChartType.Point;
            series.IsVisibleInLegend = false;
            series.Color = col;
            series.MarkerSize = markerSize;
            series.MarkerStyle = MarkerStyle.Circle;

            if (x.Length != y.Length)
            {
                return;
            }            
            
            for (int i = 0; i < x.Length; i++)
            {
                series.Points.AddXY(x[i], y[i]);
            }
        }

        void SolveEquations()
        {

            if(rbEuler.Checked)
            {
                SolveEquatinEuler();
            }
            else if(rbHeun.Checked)
            {
                SolveEquatinHeun();
            }
            else if(rbRK4.Checked)
            {
                SolveEquationRK4();
            }
            else
            {
                return;
            }
            
            
        }

        void SolveEquatinEuler()
        {
            Vector4d zeroBias = Vector4d.zero;

            x[0] = 2;
            y[0] = 0;
            dx[0] = 0;
            dy[0] = 0;
            t[0] = 0;

            for (int i = 0; i < numberPoints - 1; i++)
            {
                t[i+1] = t[i] + stepDE;

                dx[i+1] = StepEuler(dx[i], fdX, i);
                dy[i+1] = StepEuler(dy[i], fdY, i);
                x[i+1] = StepEuler(x[i], fX, i);
                y[i+1] = StepEuler(y[i], fY, i);
            }
        }
        private double StepEuler(double oldValue, FuncDelegate func, int iter)
        {
            double newValue = 0;
            Vector4d zeroBias = Vector4d.zero;

            newValue = oldValue + stepDE * func(dx[iter], dy[iter], x[iter], y[iter], zeroBias);

            return newValue;
        }
        void SolveEquatinHeun()
        {
            Vector4d zeroBias = Vector4d.zero;

            x[0] = 2;
            y[0] = 0;
            dx[0] = 0;
            dy[0] = 0;
            t[0] = 0;

            for (int i = 0; i < numberPoints - 1; i++)
            {
                t[i + 1] = t[i] + stepDE;

                //var dxTemp = dx[i] + stepDE * fdX(dx[i], dy[i], x[i], y[i], zeroBias);
                //var dyTemp = dy[i] + stepDE * fdY(dx[i], dy[i], x[i], y[i], zeroBias);
                //var xTemp = x[i] + stepDE * dx[i];
                //var yTemp = y[i] + stepDE * dy[i];

                //dx[i + 1] = dx[i] + 0.5 * stepDE * fdX(dx[i], dy[i], x[i], y[i], zeroBias) + 0.5 * stepDE * fdX(dxTemp, dyTemp, xTemp, yTemp, zeroBias);
                //dy[i + 1] = dy[i] + 0.5 * stepDE * fdY(dx[i], dy[i], x[i], y[i], zeroBias) + 0.5 * stepDE * fdY(dxTemp, dyTemp, xTemp, yTemp, zeroBias);
                //x[i + 1] = x[i] + 0.5 * stepDE * dx[i] + 0.5 * stepDE * dxTemp;
                //y[i + 1] = y[i] + 0.5 * stepDE * dy[i] + 0.5 * stepDE * dyTemp;

                //FuncDelegate[] functionsPtr = new FuncDelegate[4];  //dx, dy, x ,y
                //functionsPtr[0] = fdX;
                //functionsPtr[1] = fdY;
                //functionsPtr[2] = fX;
                //functionsPtr[3] = fX;

                Vector4d bias = new Vector4d();

                bias[0] = stepDE * fdX(dx[i], dy[i], x[i], y[i], zeroBias);
                bias[1] = stepDE * fdY(dx[i], dy[i], x[i], y[i], zeroBias);
                bias[2] = stepDE * fX(dx[i], dy[i], x[i], y[i], zeroBias);
                bias[3] = stepDE * fY(dx[i], dy[i], x[i], y[i], zeroBias);


                dx[i + 1] = dx[i] + 0.5 * stepDE * (fdX(dx[i], dy[i], x[i], y[i], zeroBias) + fdX(dx[i], dy[i], x[i], y[i], bias));
                dy[i + 1] = dy[i] + 0.5 * stepDE * (fdY(dx[i], dy[i], x[i], y[i], zeroBias) + fdY(dx[i], dy[i], x[i], y[i], bias));
                x[i + 1] = x[i] + 0.5 * stepDE * (fX(dx[i], dy[i], x[i], y[i], zeroBias) + fX(dx[i], dy[i], x[i], y[i], bias));
                y[i + 1] = y[i] + 0.5 * stepDE * (fY(dx[i], dy[i], x[i], y[i], zeroBias) + fY(dx[i], dy[i], x[i], y[i], bias));
            }
        }

        void SolveEquationRK4()
        {
            Vector4d zeroBias = Vector4d.zero;

            x[0] = 2;
            y[0] = 0;
            dx[0] = 0;
            dy[0] = 0;
            t[0] = 0;

            for (int i = 0; i < numberPoints - 1; i++)
            {
                t[i + 1] = t[i] + stepDE;

                Vector4d k1 = new Vector4d();
                Vector4d k2 = new Vector4d();
                Vector4d k3 = new Vector4d();
                Vector4d k4 = new Vector4d();

                k1[0] = fdX(dx[i], dy[i], x[i], y[i], zeroBias);
                k1[1] = fdY(dx[i], dy[i], x[i], y[i], zeroBias);
                k1[2] = fX(dx[i], dy[i], x[i], y[i], zeroBias);
                k1[3] = fY(dx[i], dy[i], x[i], y[i], zeroBias);


                k2[0] = fdX(dx[i], dy[i], x[i], y[i], k1 * (stepDE * 0.5) );
                k2[1] = fdY(dx[i], dy[i], x[i], y[i], k1 * (stepDE * 0.5) );
                k2[2] = fX(dx[i], dy[i], x[i], y[i], k1 * (stepDE * 0.5) );
                k2[3] = fY(dx[i], dy[i], x[i], y[i], k1 * (stepDE * 0.5) );

                k3[0] = fdX(dx[i], dy[i], x[i], y[i], k2 * (stepDE * 0.5));
                k3[1] = fdY(dx[i], dy[i], x[i], y[i], k2 * (stepDE * 0.5));
                k3[2] = fX(dx[i], dy[i], x[i], y[i], k2 * (stepDE * 0.5));
                k3[3] = fY(dx[i], dy[i], x[i], y[i], k2 * (stepDE * 0.5));

                k4[0] = fdX(dx[i], dy[i], x[i], y[i], k3 * stepDE);
                k4[1] = fdY(dx[i], dy[i], x[i], y[i], k3 * stepDE);
                k4[2] = fX(dx[i], dy[i], x[i], y[i], k3 * stepDE);
                k4[3] = fY(dx[i], dy[i], x[i], y[i], k3 * stepDE);

                dx[i + 1] = dx[i] + stepDE * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
                dy[i + 1] = dy[i] + stepDE * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
                x[i + 1] = x[i] + stepDE * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0;
                y[i + 1] = y[i] + stepDE * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]) / 6.0;
            }
        }

        private delegate double FuncDelegate(double dx, double dy, double x, double y, Vector4d bias );

        double fdX(double dx_, double dy_, double x_, double y_, Vector4d bias)
        {
            double ret = 0;

            dx_ += bias[0];
            dy_ += bias[1];
            x_ += bias[2];
            y_ += bias[3];

            ret = -x_ * (dx_ * dx_ + dy_ * dy_ - g * y_) / (x_ * x_ + y_ * y_);

            return ret;
        }

        double fdY(double dx_, double dy_, double x_, double y_, Vector4d bias)
        {
            double ret = 0;

            dx_ += bias[0];
            dy_ += bias[1];
            x_ += bias[2];
            y_ += bias[3];

            ret = -g - y_ * (dx_ * dx_ + dy_ * dy_ - g * y_) / (x_ * x_ + y_ * y_);

            return ret;
        }
        double fX(double dx_, double dy_, double x_, double y_, Vector4d bias)
        {
            double ret = 0;

            dx_ += bias[0];
            dy_ += bias[1];
            x_ += bias[2];
            y_ += bias[3];
            
            ret = dx_;

            return ret;
        }

        double fY(double dx_, double dy_, double x_, double y_, Vector4d bias)
        {
            double ret = 0;

            dx_ += bias[0];
            dy_ += bias[1];
            x_ += bias[2];
            y_ += bias[3];

            ret = dy_;

            return ret;
        }


        private void btnDx_Click(object sender, EventArgs e)
        {
            SolveEquations();
            InitPlot(xRange, yRange);
            Plot(t, dx, markerSizeForPlot, Color.Red, standartPlotType);
        }
        private void btnDy_Click(object sender, EventArgs e)
        {
            SolveEquations();
            InitPlot(xRange, yRange);
            Plot(t, dy, markerSizeForPlot, Color.Red, standartPlotType);
        }

        private void btnX_Click(object sender, EventArgs e)
        {
            SolveEquations();
            InitPlot(xRange, yRange);
            Plot(t, x, markerSizeForPlot, Color.Red, standartPlotType);
        }

        private void btnY_Click(object sender, EventArgs e)
        {
            SolveEquations();
            InitPlot(xRange, yRange);
            Plot(t, y, markerSizeForPlot, Color.Red, standartPlotType);
        }

        

        private void btnTrajectory_Click(object sender, EventArgs e)
        {
            if (!timer1.Enabled)
            {
                numPointFotTimer = 0;
                SolveEquations();
                timer1.Enabled = true;
                progressBar1.Value = 0;

                InitPlot(xRangeTr, yRangeTr);
            }
            else
            {
                numPointFotTimer = 0;
                timer1.Enabled = false;
                progressBar1.Value = 0;
            }

            
        }

        int numPointFotTimer = 0;
        private void timer1_Tick(object sender, EventArgs e)
        {
            InitPlot(xRangeTr, yRangeTr);
            Plot(new double[] { x[numPointFotTimer] }, new double[] { y[numPointFotTimer] }, markerSizeForAnimation, Color.Red, standartPlotType);
            progressBar1.Increment(1);

            numPointFotTimer+= stepPointForTimer;
            if(numPointFotTimer >= numberPoints)
            {
                numPointFotTimer = 0;
                timer1.Enabled = false;
            }
        }

        private void btnRandom_Click(object sender, EventArgs e)
        {
            double scale = ConvertToDoubleFromNumeric(numScaleRandom, 0.1);
            

            for (int i = 0; i < randomNoise.Length; i++)
            {
                randomNoise[i][0] = 0;
                for (int j = 1; j < numberPoints; j++)
                {
                    randomNoise[i][j] = randomNoise[i][j-1] + (r.NextDouble() * scale) - scale/2;
                    //randomNoise[i][j] = (r.NextDouble() * scale) - scale / 2;   // it is bad
                }
            }

            for (int i = 0; i < xNoise.Length; i++)
            {
                for (int j = 0; j < numberPoints; j++)
                {
                    xNoise[i][j] = x[j] + randomNoise[i * 2][j];
                    yNoise[i][j] = y[j] + randomNoise[i * 2 + 1][j];
                }
            }

            
            

        }

        private void btnXNoise_Click(object sender, EventArgs e)
        {
            //InitPlot(xRange, yRange);

            Color[] colArr = {Color.Red, Color.Blue, Color.Green, Color.Aqua, Color.Gold};

            for (int i = 0; i < xNoise.Length; i++)
            {
                Plot(t, xNoise[i], markerSizeForPlot, colArr[i], standartPlotType);
            }
        }

        private void btnYNoise_Click(object sender, EventArgs e)
        {
            //InitPlot(xRange, yRange);

            Color[] colArr = { Color.Red, Color.Blue, Color.Green, Color.Aqua, Color.Gold };

            for (int i = 0; i < xNoise.Length; i++)
            {
                Plot(t, yNoise[i], markerSizeForPlot, colArr[i], standartPlotType);
            }
        }

        private void btnMeanX_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < numberPoints; i++)
            {
                double sum = 0;
                for (int j = 0; j < xNoise.Length; j++)
                {
                    sum += xNoise[j][i];
                }

                xMean[i] = sum / xNoise.Length;
            }

            Plot(t, xMean, markerSizeForPlot, Color.Black, standartPlotType);
        }

        private void btnMeanY_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < numberPoints; i++)
            {
                double sum = 0;
                for (int j = 0; j < yNoise.Length; j++)
                {
                    sum += yNoise[j][i];
                }

                yMean[i] = sum / yNoise.Length;
            }

            Plot(t, yMean, markerSizeForPlot, Color.Black, standartPlotType);
        }

        private void btnConfIntX_Click(object sender, EventArgs e)
        {
            double squareSigma, sum, bigDelta;
            const double z99 = 5.598;  // 99.5%
            const double z95 = 2.776;  // 95%

            double z;

            if(rb95.Checked)
            {
                z = z95;
            }
            else
            {
                z = z99;
            }

            for (int i = 0; i < xMean.Length; i++) // numberPoints - all points
            {
                sum = 0;
                for (int j = 0; j < xNoise.Length; j++) // 5
                {
                    sum += Math.Pow((xNoise[j][i] - xMean[i]), 2);
                }

                squareSigma = sum / xNoise.Length;
                bigDelta = z * Math.Sqrt(squareSigma / xNoise.Length);

                bigDeltaX[i] = bigDelta;
            }

            for (int i = 0; i < xMean.Length; i++)
            {
                Plot(new double[] { t[i], t[i] }, new double[] {xMean[i] - bigDeltaX[i], xMean[i] + bigDeltaX[i] }, 1, Color.LightGreen, standartPlotType);
            }

            btnXNoise.PerformClick();
            btnMeanX.PerformClick();
        }

        private void btnConfIntY_Click(object sender, EventArgs e)
        {
            double squareSigma, sum, bigDelta;
            //const double z = 5.89;  // 99.8%
            const double z = 2.571;  // 95%

            for (int i = 0; i < yMean.Length; i++) // numberPoints - all points
            {
                sum = 0;
                for (int j = 0; j < yNoise.Length; j++) // 5
                {
                    sum += Math.Pow((yNoise[j][i] - yMean[i]), 2);
                }

                squareSigma = sum / yNoise.Length;
                bigDelta = z * Math.Sqrt(squareSigma / 5);

                bigDeltaY[i] = bigDelta;
            }

            for (int i = 0; i < yMean.Length; i++)
            {
                Plot(new double[] { t[i], t[i] }, new double[] { yMean[i] - bigDeltaY[i], yMean[i] + bigDeltaY[i] }, 1, Color.LightGreen, standartPlotType);
            }


        }

        private void btnClearPlot_Click(object sender, EventArgs e)
        {
            InitPlot(xRange, yRange);
        }
    }
}
