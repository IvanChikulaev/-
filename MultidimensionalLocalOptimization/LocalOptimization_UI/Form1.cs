using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using MultidimensionalLocalOptimization;


namespace LocalOptimization_UI
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            int n = 0;
            n = Convert.ToInt32(tb_n.Text);
            double Xmin1=0.0;
            double Xmin2=0.0;
            int fdigits = -1;
            double[] X0 = new double[n];
            double[] Xmin = new double[n];
            X0[0] = Convert.ToDouble(tb_x0.Text);
            X0[1] = Convert.ToDouble(tb_x1.Text);
            n = Convert.ToInt32(tb_n.Text);
            double[] Xp = new double[n];
            Xp[0] = Convert.ToDouble(tb_Xp.Text);
            Xp[1] = Convert.ToDouble(tb_Xp.Text);

            double[] Xc = new double[n];
            Xc[0] = Convert.ToDouble(tb_Xc.Text);
            Xc[1] = Convert.ToDouble(tb_Xc.Text);

            //UMDRIVER(n,Xc,4,fdigits,false);
            UMDRIVER(n, X0, Xmin, 1, fdigits, false);
            Xmin1 = Xmin[0];
            Xmin2 = Xmin[1];
            tb_BFGS0.Text = Xmin1.ToString();
            tb_BFGS1.Text = Xmin2.ToString();
        }

        public void FN(int n,double[] x, out double fx)
        {
            fx = ((1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]));
            //fx = (x[0]*x[0])+(x[1]*x[1]);
        }

        public void FN1(int n, double[] x,double ei,double stepsize, out double fneighbor)
        {
            for (int i = 0; i < n; i++) x[i] += stepsize * ei;
            fneighbor = (x[0] * x[0]) + (x[1] * x[1]);
        }

        public void FN2(int n, double[] x, double ei, double stepsize, out double Fii)
        {
            for (int i = 0; i < n; i++) x[i] += 2*stepsize * ei;
            Fii = ((1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]));
        }

        public void FN(int n, double[] x, double ei, double ej, double stepsizei, double stepsizej, out double Fij)
        {
            for (int i = 0; i < n; i++)
                x[i] += stepsizei*ei+stepsizej*ej;
            Fij = ((1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]));
        }

        public double[] UMDRIVER(int n, double[] X0, double[] Xmin, int global,int fdigits, bool analgrad = true)
        {
            bool maxtaken=false; // посмотреть пправильное присвоение
            double Fcurrent = 0.0, Fplus = 0.0,macheps = 0.0, eta = 0.0, gradtol = 0.0, typf=0.0; 
            int retcode = 0, itnlimit = 0, termcode = 0, consecmax = 0, itncount = 0;
            double[] Sx = new double[n];
            double[] Xf = new double[n];
            double[] Gcurrent = new double[n];
            double[] Gplus = new double[n];
            double[] g = new double[n];
            double[] Xcurrent = new double[n];
            double[] p = new double[n];
            double[,] Dx = new double[n,n];
            double[,] H = new double[n, n];
            double[] Xplus = new double[n];
            double maxstep = 0.0;
            // 1
            macheps = MACHINEPS();
            double steptol = 0.0;
            steptol = Math.Pow(macheps, 2.0 / 3.0);
            for (int i = 0; i < n; i++)
            {
                Sx[i] = 0.0;
                for (int j=0; j<n;j++)
                    Dx[i,j] = 0.0;
            }
            // 2
            UMINCK(X0,n,fdigits, Sx, Dx,macheps,
                   out eta,out termcode,out typf,out gradtol, 
                   out steptol, out itnlimit, out maxstep);
            // 3
            if (termcode < 0)
                for (int i = 0; i < n; i++) Xmin[i] = X0[i];
            // 5
            FN(n, X0,out Fcurrent); //Fcurrent=4.5
            /*6*/
            if (analgrad == false)
                FDGRAD(n, X0, Sx, eta, g); // g = 3 & 3
            else
                 ANAL_GRAD(n, X0,g);
            //7 
            UMSTOP0(n, X0, Fcurrent, g, Sx, typf, gradtol,
                out termcode, out consecmax); //termcode=0, consecmax=0 
            //8
            if (termcode > 0) for (int i = 0; i < n; i++) Xmin[i] = X0[i];
            else
            {
                //FDHESSF(n, Xcurrent, Fcurrent, Sx, eta, H);
                HESS(n, Xcurrent, H);
            }
            /*9*/
            for (int i = 0; i < n; i++)
                Xcurrent[i] = X0[i];
            /*10*/
            while (termcode == 0)
            {
                double Grad_Norma2 = 0.0;
                for (int i = 0; i < n; i++)
                {
                    Grad_Norma2 += Math.Abs(g[i] * g[i]);
                    if (i == n - 1)
                        Math.Pow(Grad_Norma2, 1.0 / 2.0);
                }
                for (int i = 0; i < n; i++) p[i] = -(g[i] / Grad_Norma2);
                // 10.1
                itncount += 1;
                if (global == 1)
                {
                    LINESEARCH(n, Xcurrent, Fcurrent, g, p, Sx,Dx, maxstep, steptol, Xplus,
                        out retcode,out Fplus,out maxtaken);
                    FDGRAD(n,Xplus,Sx, eta, Gplus);
                    FDGRAD(n, Xcurrent, Sx, eta, Gcurrent);
                    //if (analgrad == true) //call GRlf AD
                    //else //call FDGRAD
                }
                // 10.6 
                UMSTOP(n,Xcurrent, Xplus, Fplus, g, Sx, typf, retcode, gradtol, steptol, itncount, itnlimit,maxtaken,
                       out termcode, out consecmax);
                
                if (termcode > 0)
                { 
                    for(int i=0; i<n; i++)  Xf[i] = Xplus[i];
                }
                else
                {
                    BFGS(n,Xcurrent,Xplus,Gcurrent,Gplus,macheps, eta, true,H);
                }
                for (int i=0; i<n; i++)
                {
                    Xcurrent[i] = Xplus[i];
                    Xmin[i] = Xcurrent[i];
                    Fcurrent = Fplus;
                    Gcurrent[i] = Gplus[i];
                }
            }

            double[] driver = new double[n];
            for (int i = 0; i < n; i++)
                driver[i] = Xcurrent[i];
            //driver[1, 0] = Fcurrent;
            return driver;
        }

        public void LINESEARCH(int n, double[] Xcurrent, double Fcurrent, double[] Gcurrent, double[] p, double[] Sx, double[,] Dx, double maxstep, double steptol, double[] Xplus,
            out int retcode, out double FNplus, out bool maxtaken)
        {
            double alpha = 10e-4, Newtlen = 0.0, initslope = 0.0, rellenght = 0.0,
                   minlambda = 0.0, lambda = 1.0, lambdatemp = 0.0, lambdaPrev = 0.0, 
                   FplusPrev = 0.0, disc = 0.0;
            double[] array_max_el = new double[n];
            FNplus = 0.0;
            // 1
            maxtaken = false;
            // 2
            retcode = 2;

            for (int i = 0; i < n; i++)
            {
                Newtlen += Math.Pow(Dx[i, i], 2) * Math.Pow(p[i], 2);
                if (i == n - 1)
                    Newtlen = Math.Pow(Newtlen,1.0/2.0);
            }
            //5
            if (Newtlen > maxstep)
            {
                for (int i = 0; i < n; i++) p[i] = p[i] * (maxstep / Newtlen);
                Newtlen = maxstep;
            }
            //6
            for (int i = 0; i < n; i++)
                initslope += Gcurrent[i] * p[i];
            // 7
            for (int i = 0; i < n; i++)
            {
                double max_el = 0.0;
                array_max_el[i] = (Math.Abs(p[i]) / (max(Xcurrent[i], 1 / Sx[i])));
                max_el = array_max_el[0];
                rellenght = max_el;
                if (max_el > array_max_el[i])
                {
                    max_el = array_max_el[i];
                    rellenght = max_el;
                }
            }
            minlambda = steptol / rellenght;
            /*10*/
            while (retcode >= 2) // посмотреть условие 
            {
                for (int i = 0; i < n; i++)
                    Xplus[i] = Xcurrent[i] + lambda * p[i];
               
                FN(n, Xplus,out FNplus);
                // 10.3a
                if (FNplus <= Fcurrent + (alpha * lambda * initslope))
                {
                    retcode = 0;
                    if (lambda == 1.0 & (Newtlen > 0.99 * maxstep)) 
                        maxtaken = true;
                }
                // 10.3b
                else if (lambda < minlambda) 
                { 
                    retcode = 1; 
                    Xplus = Xcurrent; 
                }
                // 10.3c
                else
                {
                    double[] rezult = new double[2];
                    if (lambda == 1.0)
                    {
                        lambdatemp = -initslope / (2 * (FNplus - Fcurrent - initslope));
                    }
                    else
                    {
                        double[] matrix1_2 = new double[2] { FNplus - Fcurrent - lambda * initslope, FplusPrev - Fcurrent - lambdaPrev * initslope };
                        double[,] matrix2_2 = new double[2, 2] { { 1 / lambda * lambda, -(1 / lambdaPrev * lambdaPrev) }, { -(lambdaPrev / lambda * lambda), lambda / lambdaPrev } };
                        for (int i = 0; i < 2; i++)
                        {
                            for (int j = 0; j < 2; j++)
                            {
                                rezult[j] = 1.0 / (lambda - lambdaPrev) * matrix2_2[i, j] * matrix1_2[j];
                            }
                        }
                        double a = rezult[0], b = rezult[1];
                        disc = b * b - 3 * a * initslope;
                        if (a == 0)
                        {
                            lambdatemp = -initslope / (2 * b);
                        }
                        else
                        {
                            lambdatemp = (-b + Math.Pow(disc, 1.0 / 2.0) / 3 * a);
                            if (lambdatemp > 0.5 * lambda) 
                                lambdatemp = 0.5 * lambda;
                        }

                    }
                    lambdaPrev = lambda;
                    FplusPrev = FNplus;
                    if (lambdatemp <= 0.1 * lambda) 
                        lambda = 0.1 * lambda;
                    else lambda = lambdatemp;
                }

            }
        }

        public void UMSTOP(int n, double[] Xcurrent, double[] Xplus, double f,
            double[] g, double[] Sx, double typf, int retcode, 
            double gradtol, double steptol, int itcount, int itnlimit, bool maxtaken,
            out int termcode, out int consecmax)
        {
            termcode = 0;
            consecmax = 0;
            if (retcode == 1) termcode = 3;
            else if (MaxElementUMSTOP2b(g, Xplus, Sx, f,typf, n, typf) <= gradtol) termcode = 1;
            else if (MaxElementUMSTOP2c(n, Xplus, Xcurrent, Sx) <= steptol) termcode = 2;
            else if (itcount >= itnlimit) termcode = 4;
            else if (maxtaken == true)
            {
                consecmax += 1;
                if (consecmax == 5) termcode = 5;
            }
            else consecmax = 0;
        }

        public double MaxElementUMSTOP2b(double[] g, double[] Xplus, double[] Sx, double f,double typf, int n, double typef)
        {
            double max_element = 0.0;
            double[] array = new double[n];
            for (int i = 0; i < n; i++)
            {
                double el1 = Math.Abs(max(Math.Abs(Xplus[i]), 1.0 / Sx[i]));
                double el2 = max(Math.Abs(f), typf);
                array[i] = g[i] * (el1 / el2);
                if (i == 0) max_element = array[0];
                else
                {
                    if (array[i] > max_element)
                    {
                        max_element = array[i];
                    }
                }
            }
            return max_element;
        }

        public double max_el_from_array(int n,double[] array)
        {
            double max_element = 0.0;
            max_element = array[0];
            for (int i = 0; i < n; i++)
            {
                if (array[i] > max_element)
                {
                    max_element = array[i];
                }
            }
            return max_element;
        }

        public double MaxElementUMSTOP2c(int n, double[] Xplus, double[] Xcurrent, double[] Sx)
        {
            double max_element = 0.0;
            double[] array = new double[n];
            for (int i = 0; i < n; i++)
            {
                double el1 = Math.Abs(Xplus[i] - Xcurrent[i]);
                double el2 = max(Math.Abs(Xplus[i]), 1 / Sx[i]);
                array[i] = el1/el2;
                if (i == 0) max_element = array[0];
                else
                {
                    if (array[i] > max_element)
                    {
                        max_element = array[i];
                    }
                }
            }
            return max_element;
        }

        public double MACHINEPS()
        {
            double macheps = 1.0;
            while (1 + macheps > 1)
                macheps /= 2.0;
            return macheps * 2;
        }

        public void UMSTOP0(int n, double[] X0, double f, double[] g, double[] Sx, double typf, double gradtol, 
            out int termcode, out int consecmax)
        {
            consecmax = 0;
            if (MaxElementUMSTOP02(n, g, X0, Sx, f, typf) <= (10e-3 * gradtol))
                termcode = 1;
            else termcode = 0;
        }

        public double MaxElementUMSTOP02(int n, double[] g, double[] X0, double[] Sx, double fun, double typf)
        {
            double max_element = 0.0;
            double[] array = new double[n];
            for (int i = 0; i < n; i++)
            {
                double el1 = max(Math.Abs(X0[i]), 1.0/Sx[i]);
                double el2 = max(Math.Abs(fun),typf);
                array[i] =Math.Abs(g[i])* (el1 / el2);
                if(i==0) max_element = array[0];
                else
                {
                    if (array[i] > max_element)
                    {
                        max_element = array[i];
                    }
                }
            }
            return max_element;
        }

        public void UMINCK(double[] X0,int n, double fdigits, double[] Sx, double[,] Dx,double macheps,
            out double eta, out int termcode,out double typf,
            out double gradtol,out double steptol,out int itnlimit, out double maxstep)
        {
            typf = 1;
            gradtol = Math.Pow(macheps, 1.0 / 3.0);
            steptol = Math.Pow(macheps, 2.0 / 3.0);
            maxstep = 0.0;
            double el1 = 0.0, el2 = 0.0;
            itnlimit = 100;
            int global = 1, delta = 0;
            if (n < 1)
                termcode = -1;
            // если typf введен пользователем, то Sx[i] = 1/typf[i]
            // иначе 
            for (int i = 0; i < n; i++) Sx[i] = 1.0;
            for (int i = 0; i < n; i++) Dx[i,i] = Sx[i];

            for (int i = 0; i < n; i++)
            {
                el1 += (Math.Abs(Dx[i, i]) * Math.Abs(Dx[i, i])) * (Math.Abs(X0[i]) * Math.Abs(X0[i]));
                if (i == n - 1)
                    el1 = Math.Pow(el1, (1.0 / 2.0));
            }
            for (int i = 0; i < n; i++)
            {
                el2 += Math.Abs(Dx[i, i]) * Math.Abs(Dx[i, i]);
                if (i == n - 1)
                    el2 = Math.Pow(el2, (1.0 / 2.0));
            }
            maxstep = Math.Pow(10,3)*max(el1, el2);
            // если задано пользователем fdigits то eta = max(machieps, 10pow(-fdigits))
            // если fdigits =-1 то eta = machieps
            eta = macheps;
            if (global == 2 | global == 3) delta = -1;    // и значение дельты не введено
            termcode = 1;        
        }

        static public double sign(double Xc)
        {
            if (Xc > 0) return 1.0;
            else return -1.0;
        }
        public void ANAL_GRAD(int n, double[] Xc, double[] g)
        {
            for (int i = 0; i < n; i++) 
                g[i] = 2 * Xc[i] + 2*Xc[i];
        }

        public void FDGRAD(int n, double[] Xc, double[] Sx, double eta, double[] g)
        {
            double StepSizeJ=0.0, tempJ=0.0, fx=0.0,fp=0.0,fm=0.0;
            double cuberteta = Math.Pow(eta, 1.0 / 3.0);
            for (int j = 0; j < n; j++)
            {
                StepSizeJ = max(Xc[j], 1.0 / Sx[j]) * sign(Xc[j]);
                tempJ = Xc[j];
                Xc[j] = Xc[j] + StepSizeJ;
                StepSizeJ = Xc[j] - tempJ;
                FN(n, Xc,out fp);
                Xc[j] = tempJ - StepSizeJ;
                FN(n, Xc,out fm);
                g[j] = (fp - fm) / (2 * StepSizeJ);
                Xc[j] = tempJ;
            }
        }

        public void FDHESSF(int n, double[] Xcurrent, double Fcurrent, double[] Sx, double eta, double[,] H)
        {
            double cuberteta = 0.0;
            double[] stepsize = new double[n];
            double[] fneighbor = new double[n];
            double tempi = 0.0,tempj=0.0, ei = 1.0,ej=0.0, Fii=0.0,Fij=0.0,hi=0.0;
            //hi = Math.Pow(eta,1.0/3.0)*max(Math.Abs())
            //1
            cuberteta = Math.Pow(eta, 1.0 / 3.0);
            //2
            for (int i = 0; i < n; i++)
            {
                stepsize[i] = cuberteta * max(Math.Abs(Xcurrent[i]), 1.0 / Sx[i]) * sign(Xcurrent[i]);
                tempi = Xcurrent[i];
                Xcurrent[i] = Xcurrent[i] + stepsize[i];
                stepsize[i] = Xcurrent[i] - tempi;
                FN1(n, Xcurrent,ei,stepsize[i],out fneighbor[i]); // подумать над реализацией
                Xcurrent[i] = tempi;
            }
            // 3
            for(int i=0; i<n; i++)
            {
                tempi = Xcurrent[i];
                Xcurrent[i] = Xcurrent[i] + 2 * stepsize[i];
                FN2(n, Xcurrent,ei,stepsize[i], out Fii);
                H[i, i] = ((Fcurrent - fneighbor[i]) + (Fii - fneighbor[i]));
                Xcurrent[i] = tempi + stepsize[i];
                // 3.6
                for (int j=i+1;j<n;j++)
                {
                    tempj = Xcurrent[j];
                    Xcurrent[j] = Xcurrent[j] + stepsize[j]; //stepsize ij
                    FN(n, Xcurrent,ei,ej,stepsize[i], stepsize[i],out Fij);
                    H[i, j] = ((Fcurrent - fneighbor[i]) + (Fij - fneighbor[j])) / (stepsize[i]*stepsize[j]);
                    Xcurrent[j] = tempj;
                }
                Xcurrent[i] = tempi;
            }
                
        }

        public void HESS(int n, double[] Xcurrent, double[,] H)
        {
            for(int i=0; i<n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j) 
                        (-400.0) * Xcurrent[0];
                    else H[i, j] = 0.0; 
                }
            }
        }

        public void BFGS(int n, double[] Xc, double[] Xp, double[] Gc, double[] Gp, double macheps, double eta, bool grad, double[,] H)
        {
            double[] s = new double[n];
            double[] y = new double[n];
            double temp1 = 0.0;
            double tol = 0.0, s_norma2 = 0.0, y_norma2 = 0.0;
            bool skipupdate = true;
            for (int i = 0; i < n; i++)
                /*1*/
                s[i] = Xp[i] - Xc[i];
            for (int i = 0; i < n; i++)
                /*2*/
                y[i] = Gp[i] - Gc[i];
            //double temp = 1.0;
            for (int i = 0; i < n; i++)
                /*3*/
                temp1 += y[i] * s[i];
            for (int i = 0; i < n; i++)
            {
                s_norma2 += Math.Abs(s[i])* Math.Abs(s[i]);
                y_norma2 += Math.Abs(y[i])* Math.Abs(y[i]);
                if(i==n-1)
                {
                    s_norma2 = Math.Pow(s_norma2, 1.0 / 2.0);
                    y_norma2 = Math.Pow(y_norma2, 1.0 / 2.0);
                }
            }
            /*4*/
            if (temp1 >= (Math.Pow(macheps, 1.0/2.0)) * s_norma2 * y_norma2)
            {
                /*4.1*/
                if (grad == true) tol = eta;
                else tol = Math.Pow(eta, 1.0 / 2.0);
                /*4.2*/
                skipupdate = true;
                /*4.3*/
                double[] t = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0.0, sum2 = 0.0;
                    for (int j = 0; j < i; j++)
                    {
                        sum1 += H[j, i] * s[j];
                    }
                    for (int j = i + 1; j < n; j++)
                    {
                        sum2 += H[i, j] * s[j];
                    }
                    /*4.3.1*/
                    t[i] = sum1 + sum2;
                    /*4.3.2*/
                    if (Math.Abs(y[i] - t[i]) >= (max(Math.Abs(Gc[i]), Math.Abs(Gp[i]))*tol))
                        skipupdate = false;
                }
                /*4.4*/
                if (skipupdate == false)
                {
                    double temp2 = 0.0;
                    double[] prod = new double[n];
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                            prod[i] = s[i] * H[j, i];

                    for (int i = 0; i < n; i++)
                        temp2 += prod[i] * s[i];

                    for (int i = 0; i < n; i++)
                        for (int j = i; j < n; j++)
                            H[i, j] = H[i, j] + ((y[i] * y[j]) / temp1) - ((t[i] * t[j]) / temp2);
                    
                }
            }
        }

        public double max(double a, double b)
        {
            if (a > b)
                return (a);
            else
                return (b);
        }
        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void tb_x0_TextChanged(object sender, EventArgs e)
        {

        }
    }
}

