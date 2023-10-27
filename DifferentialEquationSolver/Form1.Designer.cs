namespace DifferentialEquationSolver
{
    partial class Form1
    {
        /// <summary>
        /// Обязательная переменная конструктора.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Освободить все используемые ресурсы.
        /// </summary>
        /// <param name="disposing">истинно, если управляемый ресурс должен быть удален; иначе ложно.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Код, автоматически созданный конструктором форм Windows

        /// <summary>
        /// Требуемый метод для поддержки конструктора — не изменяйте 
        /// содержимое этого метода с помощью редактора кода.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            this.chart1 = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.btnDx = new System.Windows.Forms.Button();
            this.btnDy = new System.Windows.Forms.Button();
            this.btnX = new System.Windows.Forms.Button();
            this.btnY = new System.Windows.Forms.Button();
            this.lbNumPoints = new System.Windows.Forms.Label();
            this.btnTrajectory = new System.Windows.Forms.Button();
            this.lbTime = new System.Windows.Forms.Label();
            this.timer1 = new System.Windows.Forms.Timer(this.components);
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            this.btnRandom = new System.Windows.Forms.Button();
            this.btnXNoise = new System.Windows.Forms.Button();
            this.numScaleRandom = new System.Windows.Forms.NumericUpDown();
            this.btnYNoise = new System.Windows.Forms.Button();
            this.btnMeanX = new System.Windows.Forms.Button();
            this.btnMeanY = new System.Windows.Forms.Button();
            this.btnConfIntX = new System.Windows.Forms.Button();
            this.btnConfIntY = new System.Windows.Forms.Button();
            this.rb99 = new System.Windows.Forms.RadioButton();
            this.rb95 = new System.Windows.Forms.RadioButton();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.rbHeun = new System.Windows.Forms.RadioButton();
            this.rbRK4 = new System.Windows.Forms.RadioButton();
            this.rbEuler = new System.Windows.Forms.RadioButton();
            this.lbStep = new System.Windows.Forms.Label();
            this.btnClearPlot = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.numScaleRandom)).BeginInit();
            this.groupBox1.SuspendLayout();
            this.SuspendLayout();
            // 
            // chart1
            // 
            chartArea2.Name = "ChartArea1";
            this.chart1.ChartAreas.Add(chartArea2);
            this.chart1.Location = new System.Drawing.Point(12, 12);
            this.chart1.Name = "chart1";
            this.chart1.Size = new System.Drawing.Size(550, 550);
            this.chart1.TabIndex = 0;
            this.chart1.Text = "chart1";
            // 
            // btnDx
            // 
            this.btnDx.Location = new System.Drawing.Point(568, 12);
            this.btnDx.Name = "btnDx";
            this.btnDx.Size = new System.Drawing.Size(75, 23);
            this.btnDx.TabIndex = 1;
            this.btnDx.Text = "dX";
            this.btnDx.UseVisualStyleBackColor = true;
            this.btnDx.Click += new System.EventHandler(this.btnDx_Click);
            // 
            // btnDy
            // 
            this.btnDy.Location = new System.Drawing.Point(568, 41);
            this.btnDy.Name = "btnDy";
            this.btnDy.Size = new System.Drawing.Size(75, 23);
            this.btnDy.TabIndex = 1;
            this.btnDy.Text = "dY";
            this.btnDy.UseVisualStyleBackColor = true;
            this.btnDy.Click += new System.EventHandler(this.btnDy_Click);
            // 
            // btnX
            // 
            this.btnX.Location = new System.Drawing.Point(568, 70);
            this.btnX.Name = "btnX";
            this.btnX.Size = new System.Drawing.Size(75, 23);
            this.btnX.TabIndex = 1;
            this.btnX.Text = "X";
            this.btnX.UseVisualStyleBackColor = true;
            this.btnX.Click += new System.EventHandler(this.btnX_Click);
            // 
            // btnY
            // 
            this.btnY.Location = new System.Drawing.Point(568, 99);
            this.btnY.Name = "btnY";
            this.btnY.Size = new System.Drawing.Size(75, 23);
            this.btnY.TabIndex = 1;
            this.btnY.Text = "Y";
            this.btnY.UseVisualStyleBackColor = true;
            this.btnY.Click += new System.EventHandler(this.btnY_Click);
            // 
            // lbNumPoints
            // 
            this.lbNumPoints.AutoSize = true;
            this.lbNumPoints.Location = new System.Drawing.Point(648, 437);
            this.lbNumPoints.Name = "lbNumPoints";
            this.lbNumPoints.Size = new System.Drawing.Size(35, 13);
            this.lbNumPoints.TabIndex = 2;
            this.lbNumPoints.Text = "label1";
            // 
            // btnTrajectory
            // 
            this.btnTrajectory.Location = new System.Drawing.Point(568, 128);
            this.btnTrajectory.Name = "btnTrajectory";
            this.btnTrajectory.Size = new System.Drawing.Size(75, 23);
            this.btnTrajectory.TabIndex = 1;
            this.btnTrajectory.Text = "Trajectory";
            this.btnTrajectory.UseVisualStyleBackColor = true;
            this.btnTrajectory.Click += new System.EventHandler(this.btnTrajectory_Click);
            // 
            // lbTime
            // 
            this.lbTime.AutoSize = true;
            this.lbTime.Location = new System.Drawing.Point(648, 414);
            this.lbTime.Name = "lbTime";
            this.lbTime.Size = new System.Drawing.Size(35, 13);
            this.lbTime.TabIndex = 3;
            this.lbTime.Text = "label1";
            // 
            // timer1
            // 
            this.timer1.Tick += new System.EventHandler(this.timer1_Tick);
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(568, 158);
            this.progressBar1.Name = "progressBar1";
            this.progressBar1.Size = new System.Drawing.Size(75, 23);
            this.progressBar1.TabIndex = 4;
            // 
            // btnRandom
            // 
            this.btnRandom.Location = new System.Drawing.Point(569, 188);
            this.btnRandom.Name = "btnRandom";
            this.btnRandom.Size = new System.Drawing.Size(75, 23);
            this.btnRandom.TabIndex = 5;
            this.btnRandom.Text = "Random";
            this.btnRandom.UseVisualStyleBackColor = true;
            this.btnRandom.Click += new System.EventHandler(this.btnRandom_Click);
            // 
            // btnXNoise
            // 
            this.btnXNoise.Location = new System.Drawing.Point(569, 217);
            this.btnXNoise.Name = "btnXNoise";
            this.btnXNoise.Size = new System.Drawing.Size(75, 23);
            this.btnXNoise.TabIndex = 7;
            this.btnXNoise.Text = "xNoise";
            this.btnXNoise.UseVisualStyleBackColor = true;
            this.btnXNoise.Click += new System.EventHandler(this.btnXNoise_Click);
            // 
            // numScaleRandom
            // 
            this.numScaleRandom.DecimalPlaces = 2;
            this.numScaleRandom.Increment = new decimal(new int[] {
            1,
            0,
            0,
            131072});
            this.numScaleRandom.Location = new System.Drawing.Point(651, 190);
            this.numScaleRandom.Maximum = new decimal(new int[] {
            101,
            0,
            0,
            131072});
            this.numScaleRandom.Name = "numScaleRandom";
            this.numScaleRandom.Size = new System.Drawing.Size(73, 20);
            this.numScaleRandom.TabIndex = 8;
            this.numScaleRandom.Value = new decimal(new int[] {
            20,
            0,
            0,
            131072});
            // 
            // btnYNoise
            // 
            this.btnYNoise.Location = new System.Drawing.Point(569, 246);
            this.btnYNoise.Name = "btnYNoise";
            this.btnYNoise.Size = new System.Drawing.Size(75, 23);
            this.btnYNoise.TabIndex = 9;
            this.btnYNoise.Text = "yNoise";
            this.btnYNoise.UseVisualStyleBackColor = true;
            this.btnYNoise.Click += new System.EventHandler(this.btnYNoise_Click);
            // 
            // btnMeanX
            // 
            this.btnMeanX.Location = new System.Drawing.Point(569, 275);
            this.btnMeanX.Name = "btnMeanX";
            this.btnMeanX.Size = new System.Drawing.Size(75, 23);
            this.btnMeanX.TabIndex = 10;
            this.btnMeanX.Text = "Mean X";
            this.btnMeanX.UseVisualStyleBackColor = true;
            this.btnMeanX.Click += new System.EventHandler(this.btnMeanX_Click);
            // 
            // btnMeanY
            // 
            this.btnMeanY.Location = new System.Drawing.Point(569, 304);
            this.btnMeanY.Name = "btnMeanY";
            this.btnMeanY.Size = new System.Drawing.Size(75, 23);
            this.btnMeanY.TabIndex = 11;
            this.btnMeanY.Text = "Mean Y";
            this.btnMeanY.UseVisualStyleBackColor = true;
            this.btnMeanY.Click += new System.EventHandler(this.btnMeanY_Click);
            // 
            // btnConfIntX
            // 
            this.btnConfIntX.Location = new System.Drawing.Point(569, 333);
            this.btnConfIntX.Name = "btnConfIntX";
            this.btnConfIntX.Size = new System.Drawing.Size(75, 23);
            this.btnConfIntX.TabIndex = 12;
            this.btnConfIntX.Text = "ConfIntX";
            this.btnConfIntX.UseVisualStyleBackColor = true;
            this.btnConfIntX.Click += new System.EventHandler(this.btnConfIntX_Click);
            // 
            // btnConfIntY
            // 
            this.btnConfIntY.Location = new System.Drawing.Point(569, 362);
            this.btnConfIntY.Name = "btnConfIntY";
            this.btnConfIntY.Size = new System.Drawing.Size(75, 23);
            this.btnConfIntY.TabIndex = 13;
            this.btnConfIntY.Text = "ConfIntY";
            this.btnConfIntY.UseVisualStyleBackColor = true;
            this.btnConfIntY.Click += new System.EventHandler(this.btnConfIntY_Click);
            // 
            // rb99
            // 
            this.rb99.AutoSize = true;
            this.rb99.Checked = true;
            this.rb99.Location = new System.Drawing.Point(651, 336);
            this.rb99.Name = "rb99";
            this.rb99.Size = new System.Drawing.Size(54, 17);
            this.rb99.TabIndex = 14;
            this.rb99.TabStop = true;
            this.rb99.Text = "99.5%";
            this.rb99.UseVisualStyleBackColor = true;
            // 
            // rb95
            // 
            this.rb95.AutoSize = true;
            this.rb95.Location = new System.Drawing.Point(650, 365);
            this.rb95.Name = "rb95";
            this.rb95.Size = new System.Drawing.Size(45, 17);
            this.rb95.TabIndex = 14;
            this.rb95.Text = "95%";
            this.rb95.UseVisualStyleBackColor = true;
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.rbHeun);
            this.groupBox1.Controls.Add(this.rbRK4);
            this.groupBox1.Controls.Add(this.rbEuler);
            this.groupBox1.Location = new System.Drawing.Point(649, 12);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(75, 100);
            this.groupBox1.TabIndex = 15;
            this.groupBox1.TabStop = false;
            // 
            // rbHeun
            // 
            this.rbHeun.AutoSize = true;
            this.rbHeun.Location = new System.Drawing.Point(6, 65);
            this.rbHeun.Name = "rbHeun";
            this.rbHeun.Size = new System.Drawing.Size(51, 17);
            this.rbHeun.TabIndex = 1;
            this.rbHeun.TabStop = true;
            this.rbHeun.Text = "Heun";
            this.rbHeun.UseVisualStyleBackColor = true;
            // 
            // rbRK4
            // 
            this.rbRK4.AutoSize = true;
            this.rbRK4.Location = new System.Drawing.Point(6, 42);
            this.rbRK4.Name = "rbRK4";
            this.rbRK4.Size = new System.Drawing.Size(46, 17);
            this.rbRK4.TabIndex = 0;
            this.rbRK4.Text = "RK4";
            this.rbRK4.UseVisualStyleBackColor = true;
            // 
            // rbEuler
            // 
            this.rbEuler.AutoSize = true;
            this.rbEuler.Checked = true;
            this.rbEuler.Location = new System.Drawing.Point(6, 19);
            this.rbEuler.Name = "rbEuler";
            this.rbEuler.Size = new System.Drawing.Size(49, 17);
            this.rbEuler.TabIndex = 0;
            this.rbEuler.TabStop = true;
            this.rbEuler.Text = "Euler";
            this.rbEuler.UseVisualStyleBackColor = true;
            // 
            // lbStep
            // 
            this.lbStep.AutoSize = true;
            this.lbStep.Location = new System.Drawing.Point(648, 460);
            this.lbStep.Name = "lbStep";
            this.lbStep.Size = new System.Drawing.Size(35, 13);
            this.lbStep.TabIndex = 16;
            this.lbStep.Text = "label1";
            // 
            // btnClearPlot
            // 
            this.btnClearPlot.Location = new System.Drawing.Point(569, 537);
            this.btnClearPlot.Name = "btnClearPlot";
            this.btnClearPlot.Size = new System.Drawing.Size(75, 23);
            this.btnClearPlot.TabIndex = 17;
            this.btnClearPlot.Text = "Clear Plot";
            this.btnClearPlot.UseVisualStyleBackColor = true;
            this.btnClearPlot.Click += new System.EventHandler(this.btnClearPlot_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(736, 572);
            this.Controls.Add(this.btnClearPlot);
            this.Controls.Add(this.lbStep);
            this.Controls.Add(this.groupBox1);
            this.Controls.Add(this.rb95);
            this.Controls.Add(this.rb99);
            this.Controls.Add(this.btnConfIntY);
            this.Controls.Add(this.btnConfIntX);
            this.Controls.Add(this.btnMeanY);
            this.Controls.Add(this.btnMeanX);
            this.Controls.Add(this.btnYNoise);
            this.Controls.Add(this.numScaleRandom);
            this.Controls.Add(this.btnXNoise);
            this.Controls.Add(this.btnRandom);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.lbTime);
            this.Controls.Add(this.lbNumPoints);
            this.Controls.Add(this.btnTrajectory);
            this.Controls.Add(this.btnY);
            this.Controls.Add(this.btnX);
            this.Controls.Add(this.btnDy);
            this.Controls.Add(this.btnDx);
            this.Controls.Add(this.chart1);
            this.Name = "Form1";
            this.Text = "Form1";
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.numScaleRandom)).EndInit();
            this.groupBox1.ResumeLayout(false);
            this.groupBox1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart chart1;
        private System.Windows.Forms.Button btnDx;
        private System.Windows.Forms.Button btnDy;
        private System.Windows.Forms.Button btnX;
        private System.Windows.Forms.Button btnY;
        private System.Windows.Forms.Label lbNumPoints;
        private System.Windows.Forms.Button btnTrajectory;
        private System.Windows.Forms.Label lbTime;
        private System.Windows.Forms.Timer timer1;
        private System.Windows.Forms.ProgressBar progressBar1;
        private System.Windows.Forms.Button btnRandom;
        private System.Windows.Forms.Button btnXNoise;
        private System.Windows.Forms.NumericUpDown numScaleRandom;
        private System.Windows.Forms.Button btnYNoise;
        private System.Windows.Forms.Button btnMeanX;
        private System.Windows.Forms.Button btnMeanY;
        private System.Windows.Forms.Button btnConfIntX;
        private System.Windows.Forms.Button btnConfIntY;
        private System.Windows.Forms.RadioButton rb99;
        private System.Windows.Forms.RadioButton rb95;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.RadioButton rbRK4;
        private System.Windows.Forms.RadioButton rbEuler;
        private System.Windows.Forms.RadioButton rbHeun;
        private System.Windows.Forms.Label lbStep;
        private System.Windows.Forms.Button btnClearPlot;
    }
}

