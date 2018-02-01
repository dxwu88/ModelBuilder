using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Microsoft.SolverFoundation.Common;
using Microsoft.SolverFoundation.Solvers;
using Microsoft.SolverFoundation.Services;

namespace ModelTools
{
    public class Optimizer
    {
        public int _numberOfVariables;

        public int _termOrder;
        public int _modelItemCount;
        public int[,]_variableTerms = null;

        public int _dataRecordCount;
        public double[,] _independentVariables;
        public double[] _dependentVariable;
        public double[] _rawDependentVariable;
        public double[] _parameters = new double[50];

        // skewclear
        public int[] _modelItemCounts = null;
        public int[,,] _itemVariableTerms = null;
        public double[,] _itemParameters = null;

        // PP
        public double[] _parametersPP = new double[25];

        public Optimizer(int dataRecordCount, double[,] independentVariables, double[] dependentVariable)
        {
            _dataRecordCount = dataRecordCount;
            int row = independentVariables.GetLength(0);
            int col = independentVariables.GetLength(1);
            _independentVariables = new double[row, col];
            _independentVariables = independentVariables;
            _dependentVariable = new double[row];
            _dependentVariable = dependentVariable;

            _variableTerms = new int[50, 3];

            _modelItemCounts = new int[8];
            _itemVariableTerms = new int[8, 50, 3];
            _itemParameters = new double[8, 50];
        }

        public void FitSkewModel()
        {
            var solverParams = new NelderMeadSolverParams();
            NelderMeadSolver solver = new NelderMeadSolver();
            int vidRow;
            int[] vidVaribales = new int[8];

            //add variables
            for (int i = 0; i < 8; ++i)
            {
                solver.AddVariable(null, out vidVaribales[i]);
            }

            //set bounds
            solver.SetBounds(vidVaribales[0], 0.0, 1.0E5);
            for (int i = 1; i < 8; ++i)
            {
                solver.SetBounds(vidVaribales[i], -1.0E5, 1.0E5);
            }

            //set initial values
            for (int i = 0; i < 8; ++i)
            {
                solver.SetProperty(SolverProperties.VariableStartValue, vidVaribales[i], 1.0);
            }

            //add a row and set it as the goal
            solver.AddRow(null, out vidRow);
            solver.AddGoal(vidRow, 0, true);

            solver.FunctionEvaluator = EvalFunction;
            //solver.GradientEvaluator = EvalGradient;
            var nms = solver.Solve(solverParams);
            Console.WriteLine(solver.ToString());

            // Delete the file if it exists.
            string strFile = @"C:\data.dat";
            if (File.Exists(strFile))
            {
                File.Delete(strFile);
            }

            FileStream fs = new FileStream(strFile, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);

            double x = 0.0;
            double[] res = new double[1354];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double T = _independentVariables[i, 1]; // - 185.0) / 1200.0;
                double yp = _parameters[0] * Math.Exp(_parameters[1] * (T + _parameters[2]));
                yp /= Math.Exp(-_parameters[3] * (T - _parameters[4])) + _parameters[5] * Math.Exp(_parameters[6] * (T - _parameters[7]));

                res[i] = _dependentVariable[i] - yp;
                x += Math.Abs(res[i]);

                sw.WriteLine("{0}, {1}, {2}, {3}", _independentVariables[i,1], _dependentVariable[i], yp, res[i]);
            }
            x /= _dataRecordCount;

            sw.Close();
        }

        private double EvalFunction(INonlinearModel model, int rowVid, ValuesByIndex values, bool newValues)
        {
            double obj = 0.0;

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                // T
                double T = _independentVariables[i, 1];

                double yp = values[1] * Math.Exp(values[2] * (T + values[3]));
                yp /= Math.Exp(-values[4] * (T - values[5])) + values[6] * Math.Exp(values[7] * (T - values[8]));

                double y = _dependentVariable[i];
                obj += (y - yp) * (y - yp);
            }

            for (int i = 0; i < 8; ++i)
            {
                _parameters[i] = values[i + 1];
            }

            return obj;
        }

        private double EvalReFitFunction(INonlinearModel model, int rowVid, ValuesByIndex values, bool newValues)
        {
            double obj = 0.0;

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double z = 0.0;
                for (int j = 0; j < _modelItemCount; ++j)
                {
                    int k1 = _variableTerms[j, 0];
                    int k2 = _variableTerms[j, 1];
                    int k3 = _variableTerms[j, 2];

                    double D1 = _independentVariables[i, k1];
                    double D2 = _independentVariables[i, k2];
                    double D3 = _independentVariables[i, k3];
                    z += _parameters[j] * D1 * D2 * D3;
                }
                
                double res = z - _dependentVariable[i];
                if (_independentVariables[i, 1] == 0.0 ||_independentVariables[i, 1] == (350.0 - 185.0) / 1200.0)
                    obj += 20.0*Math.Abs(res);
                else
                    obj += Math.Abs(res);
            }

            for (int i = 0; i < _modelItemCount; ++i)
            {
                _parameters[i] = values[i + 1];
            }

            return obj;
        }
        public void FitExpModel()
        {
            var solverParams = new NelderMeadSolverParams();
            NelderMeadSolver solver = new NelderMeadSolver();
            int vidRow;
            int[] vidVaribales = new int[3];

            //add variables
            for (int i = 0; i < 3; ++i)
            {
                solver.AddVariable(null, out vidVaribales[i]);
            }

            //set bounds
            for (int i = 0; i < 3; ++i)
            {
                solver.SetBounds(vidVaribales[i], -1.0E5, 1.0E5);
            }

            //set initial values
            for (int i = 0; i < 3; ++i)
            {
                solver.SetProperty(SolverProperties.VariableStartValue, vidVaribales[i], 1.0);
            }

            //add a row and set it as the goal
            solver.AddRow(null, out vidRow);
            solver.AddGoal(vidRow, 0, true);

            solver.FunctionEvaluator = EvalExpFunction;
            //solver.GradientEvaluator = EvalGradient;
            var nms = solver.Solve(solverParams);
            Console.WriteLine(solver.ToString());

            // Delete the file if it exists.
            string strFile = @"C:\data.dat";
            if (File.Exists(strFile))
            {
                File.Delete(strFile);
            }

            FileStream fs = new FileStream(strFile, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);

            double x = 0.0;
            double[] res = new double[1354];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double T = _independentVariables[i, 1]; // - 185.0) / 1200.0;
                double yp = 1.0 / (_parameters[0] + Math.Exp(_parameters[1] * (T - _parameters[2])));

                res[i] = _dependentVariable[i] - yp;
                x += Math.Abs(res[i]);

                sw.WriteLine("{0}, {1}, {2}, {3}", _independentVariables[i, 1], _dependentVariable[i], yp, res[i]);
            }
            x /= _dataRecordCount;

            sw.Close();
        }

        private double EvalExpFunction(INonlinearModel model, int rowVid, ValuesByIndex values, bool newValues)
        {
            double obj = 0.0;

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                // T
                double T = _independentVariables[i, 1];

                double yp = 1.0 / (values[1] + Math.Exp(values[2] * (T - values[3])));

                double y = _dependentVariable[i];
                obj += (y - yp) * (y - yp);
            }

            for (int i = 0; i < 3; ++i)
            {
                _parameters[i] = values[i + 1];
            }

            return obj;
        }

        public double CalFunction(int parameterIndex, double T, double Y, double[] parameters)
        {
            double D, D1, D2, N, yprime, p0, p1, p2, p3, p4, p5, p6, p7;

            N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
            D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));

            switch (parameterIndex)
            {
                case 0:
                    N = Math.Exp(parameters[1] * (T + parameters[2]));
                    D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    p0 = Y * D / N;
                    return p0;

                case 1:
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    p1 = Y * D / parameters[0];
                    p1 = Math.Log(p1) / (T + parameters[2]);
                    return p1;

                case 2:
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    p2 = Y * D / parameters[0];
                    p2 = Math.Log(p2) / parameters[1] - T;
                    return p2;

                case 3:
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                    D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    p3 = N / Y - D2;
                    p3 = - Math.Log(p3) / (T - parameters[4]);
                    return p3;

                case 4:
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                    D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    p4 = N / Y - D2;
                    p4 = (Math.Log(p4) - D2) / parameters[3] + T;
                    return p4;

                case 5:
                    D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                    D2 = Math.Exp(parameters[6] * (T - parameters[7]));
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    yprime = N / Y;
                    yprime = (Math.Log(yprime) - D1) / D2;
                    return yprime;

                case 6:
                    D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                    D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    yprime = (N / Y - D1) / parameters[5];
                    yprime = Math.Log(yprime) / (T - parameters[7]);
                    return yprime;

                case 7:
                    D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                    D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                    N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                    yprime = (N / Y - D1) / parameters[5];
                    yprime = Math.Log(yprime) / parameters[6] - T;
                    return yprime;

            }

            return 0.0;
        }

        public double[] CalculateNewDependentVariable(int item)
        {
            double[] parameters = new double[8];
            double[] y = new double[_dataRecordCount];
            double N, D, D1, D2;

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double T = _independentVariables[i, 1]; 
                double Y = _dependentVariable[i];

                switch (item)
                {
                    case 0:
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);

                        N = Math.Exp(parameters[1] * (T + parameters[2]));
                        D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p0 = Y * D / N;
                        y[i] = p0;

                        break;

                    case 1:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);
                        if (parameters[0] <= 0.0)
                        {
                            int nnnnnn = 0;
                        }

                        D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p1 = Y * D / parameters[0];
                        p1 = Math.Log(p1) / (T + parameters[2]);
                        y[i] = p1;

                        break;

                    case 2:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);

                        D = Math.Exp(-parameters[3] * (T - parameters[4])) + parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p2 = Y * D / parameters[0];
                        p2 = Math.Log(p2) / parameters[1] - T;
                        y[i] = p2;

                        break;

                    case 3:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);

                        N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                        //D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                        D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p3 = N / Y - D2;
                        p3 = - Math.Log(p3) / (T - parameters[4]);
                        y[i] = p3;

                        break;

                    case 4:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);

                        N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                        //D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                        D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p4 = N / Y - D2;
                        p4 = T + Math.Log(p4) / parameters[3];
                        y[i] = p4;

                        break;

                    case 5:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[6] = CalculateParameter(6, i);
                        parameters[7] = CalculateParameter(7, i);

                        N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                        D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                        //D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p5 = N / Y - D1;
                        p5 = p5 / Math.Exp(parameters[6] * (T - parameters[7]));
                        y[i] = p5;

                        break;

                    case 6:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[7] = CalculateParameter(7, i);

                        N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                        D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                        //D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p6 = N / Y - D1;
                        p6 = Math.Log(p6 / parameters[5]) / (T - parameters[7]);
                        y[i] = p6;

                        break;

                    case 7:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);
                        parameters[3] = CalculateParameter(3, i);
                        parameters[4] = CalculateParameter(4, i);
                        parameters[5] = CalculateParameter(5, i);
                        parameters[6] = CalculateParameter(6, i);

                        N = parameters[0] * Math.Exp(parameters[1] * (T + parameters[2]));
                        D1 = Math.Exp(-parameters[3] * (T - parameters[4]));
                        //D2 = parameters[5] * Math.Exp(parameters[6] * (T - parameters[7]));
                        double p7 = N / Y - D1;
                        p7 = T - Math.Log(p7 / parameters[5]) / parameters[6];
                        y[i] = p7;

                        break;

                    default:
                        break;

                }
            }

            return y;
        }

        public double[] CalculateNewExpDependentVariable(int item)
        {
            double[] parameters = new double[3];
            double[] y = new double[_dataRecordCount];

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double T = _independentVariables[i, 1];
                double Y = _dependentVariable[i];

                switch (item)
                {
                    case 0:
                        parameters[1] = CalculateParameter(1, i);
                        parameters[2] = CalculateParameter(2, i);

                        double p0 = 1.0 / Y - Math.Exp(T - parameters[1]);
                        y[i] = p0;

                        break;

                    case 1:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[2] = CalculateParameter(2, i);

                        double p1 = Math.Log(1.0 / Y  - parameters[0]) / (T - parameters[1]);
                        y[i] = p1;

                        break;

                    case 2:
                        parameters[0] = CalculateParameter(0, i);
                        parameters[1] = CalculateParameter(1, i);

                        double p2 = T - Math.Log(1.0 / Y - parameters[0]) / parameters[1];
                        y[i] = p2;

                        break;

                    default:
                        break;

                }
            }

            return y;
        }

        public double CalculateParameter(int item, int dataRecord)
        {
            double x = 0.0;
            for (int k = 0; k <= _modelItemCounts[item]; ++k)
            {
                int k1 = _itemVariableTerms[item, k, 0];
                int k2 = _itemVariableTerms[item, k, 1];
                int k3 = _itemVariableTerms[item, k, 2];

                double d1 = _independentVariables[dataRecord, k1];
                double d2 = _independentVariables[dataRecord, k2];
                double d3 = _independentVariables[dataRecord, k3];

                x += _itemParameters[item, k] * d1 * d2 * d3;
            }


            return x <= 0.0 ? 0.000001 : x;
        }

        public void FitNewModel()
        {
            _itemParameters = new double[8, 50];
            for (int item = 0; item < 8; ++item)
            {
                _itemParameters[item, 0] = _parameters[item];
            }

            int count = 0;
            while (count < 5)
            {
                ++count;

                for (int item = 0; item <= 8 && item != 3 && item != 5 && item != 5 && item != 6; ++item)
                {
                    double x;
                    int modelItemCount = 0;
                    double AIC = double.MaxValue;

                    double aicCreterion0 = double.MaxValue;
                    double aicCreterion = double.MaxValue / 2.0;
                    int[] oneTerm = { 0, 0, 0 };

                    double[] initialGuess = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

                    double[] parameters = new double[50];
                    double residual = double.MaxValue;

                    double[] pY = CalculateNewDependentVariable(item);

                    while (aicCreterion < aicCreterion0 && modelItemCount < 49)
                    {
                        aicCreterion0 = aicCreterion;

                        bool findMinAIC = false;

                        // first order terms
                        for (int i = 0; i < _numberOfVariables && i != 1; ++i)
                        {
                            int[] oneTerm0 = { i, 0, 0 };
                            if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                            {
                                SetModelTermCount(item, modelItemCount);
                                SetVariableTerms(item, oneTerm0);

                                double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];

                                SetupModelData(item, pData);

                                double[] param = new double[_modelItemCounts[item] + 1];
                                x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                x = Math.Log(x) + 0.0001 * (modelItemCount + 1.0);
                                if (x < AIC)
                                {
                                    for (int k = 0; k <= modelItemCount; ++k)
                                    {
                                        initialGuess[k] = param[k];
                                    }

                                    AIC = x;
                                    Array.Copy(oneTerm0, oneTerm, 3);
                                    findMinAIC = true;
                                }
                            }
                        }

                        // second order terms
                        for (int i = 2; i < _numberOfVariables; ++i)
                        {
                            for (int j = i; j < _numberOfVariables; ++j)
                            {
                                int[] oneTerm0 = { i, j, 0 };
                                if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                {
                                    SetModelTermCount(item, modelItemCount);
                                    SetVariableTerms(item, oneTerm0);

                                    double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                    SetupModelData(item, pData);

                                    double[] param = new double[_modelItemCounts[item] + 1];
                                    x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                    x = Math.Log(x) + 0.0001 * (modelItemCount + 1.0);
                                    if (x < AIC)
                                    {
                                        for (int k = 0; k <= modelItemCount; ++k)
                                        {
                                            initialGuess[k] = param[k];
                                        }

                                        AIC = x;
                                        Array.Copy(oneTerm0, oneTerm, 3);
                                        findMinAIC = true;
                                    }
                                }
                            }
                        }

                        // third order terms
                        //if (modelItemCount > 5)
                        {
                            for (int i = 2; i < _numberOfVariables; ++i)
                            {
                                for (int j = i; j < _numberOfVariables; ++j)
                                {
                                    for (int k = j; k < _numberOfVariables; ++k)
                                    {
                                        int[] oneTerm0 = { i, j, k };
                                        if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                        {
                                            SetModelTermCount(item, modelItemCount);
                                            SetVariableTerms(item, oneTerm0);

                                            double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                            SetupModelData(item, pData);

                                            double[] param = new double[_modelItemCounts[item] + 1];
                                            x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                            x = Math.Log(x) + 0.0001 * (modelItemCount + 1.0);
                                            if (x < AIC)
                                            {
                                                for (int n = 0; n <= modelItemCount; ++n)
                                                {
                                                    initialGuess[n] = param[n];
                                                }

                                                AIC = x;
                                                Array.Copy(oneTerm0, oneTerm, 3);
                                                findMinAIC = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (findMinAIC)
                        {
                            Array.Copy(initialGuess, parameters, modelItemCount + 1);
                            aicCreterion = AIC;
                            residual = Math.Sqrt(Math.Exp(AIC - 0.0001 * (modelItemCount + 1.0))) / _dataRecordCount;

                            _itemVariableTerms[item, modelItemCount, 0] = oneTerm[0];
                            _itemVariableTerms[item, modelItemCount, 1] = oneTerm[1];
                            _itemVariableTerms[item, modelItemCount, 2] = oneTerm[2];

                            for (int n = 0; n < modelItemCount + 1; ++n)
                            {
                                _itemParameters[item, n] = parameters[n];
                            }

                            modelItemCount++;
                        }
                    }

                } // for 8
            } // while

            double x00 = 0.0;
            double[] res = new double[_dataRecordCount];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                    res[i] = CalculateParameter(5, i);
            }

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                //double parameter0 = CalculateParameter(0, i);
                for (int j = 0; j < 8; ++j)
                {
                    _parameters[j] = CalculateParameter(j, i);
                }

                double T = _independentVariables[i, 1]; 
                double N = Math.Exp(_parameters[1] * (T + _parameters[2]));
                double D = Math.Exp(-_parameters[3] * (T - _parameters[4])) + _parameters[5] * Math.Exp(_parameters[6] * (T - _parameters[7]));
                res[i] = _dependentVariable[i] - _parameters[0] * N / D;
                x00 += Math.Abs(res[i]);
            }
            x00 /= _dataRecordCount;

        }

        public void FitExtSkewModel()
        {
            _itemParameters = new double[8, 50];
            for (int item = 0; item < 8; ++item)
            {
                _itemParameters[item, 0] = _parameters[item];
            }

            int count = 0;
            //while (count < 5)
            {
                ++count;

                for (int item = 0; item < 1; ++item)
                {
                    double x;
                    int modelItemCount = 0;
                    double AIC = double.MaxValue;

                    double aicCreterion0 = double.MaxValue;
                    double aicCreterion = double.MaxValue / 2.0;
                    int[] oneTerm = { 0, 0, 0 };

                    double[] initialGuess = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

                    double[] parameters = new double[50];
                    double residual = double.MaxValue;

                    double[] pY = CalculateNewDependentVariable(item);

                    while (aicCreterion < aicCreterion0)
                    {
                        aicCreterion0 = aicCreterion;

                        bool findMinAIC = false;

                        // first order terms
                        for (int i = 0; i < _numberOfVariables; ++i)
                        {
                            int[] oneTerm0 = { i, 0, 0 };
                            if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                            {
                                SetModelTermCount(item, modelItemCount);
                                SetVariableTerms(item, oneTerm0);

                                double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];

                                SetupModelData(item, pData);

                                double[] param = new double[_modelItemCounts[item] + 1];
                                x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                if (x < AIC)
                                {
                                    for (int k = 0; k <= modelItemCount; ++k)
                                    {
                                        initialGuess[k] = param[k];
                                    }

                                    AIC = x;
                                    Array.Copy(oneTerm0, oneTerm, 3);
                                    findMinAIC = true;
                                }
                            }
                        }

                        // second order terms
                        for (int i = 1; i < _numberOfVariables; ++i)
                        {
                            for (int j = i; j < _numberOfVariables; ++j)
                            {
                                int[] oneTerm0 = { i, j, 0 };
                                if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                {
                                    SetModelTermCount(item, modelItemCount);
                                    SetVariableTerms(item, oneTerm0);

                                    double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                    SetupModelData(item, pData);

                                    double[] param = new double[_modelItemCounts[item] + 1];
                                    x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                    x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                    if (x < AIC)
                                    {
                                        for (int k = 0; k <= modelItemCount; ++k)
                                        {
                                            initialGuess[k] = param[k];
                                        }

                                        AIC = x;
                                        Array.Copy(oneTerm0, oneTerm, 3);
                                        findMinAIC = true;
                                    }
                                }
                            }
                        }

                        // third order terms
                        //if (modelItemCount > 5)
                        {
                            for (int i = 1; i < _numberOfVariables; ++i)
                            {
                                for (int j = i; j < _numberOfVariables; ++j)
                                {
                                    for (int k = j; k < _numberOfVariables; ++k)
                                    {
                                        int[] oneTerm0 = { i, j, k };
                                        if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                        {
                                            SetModelTermCount(item, modelItemCount);
                                            SetVariableTerms(item, oneTerm0);

                                            double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                            SetupModelData(item, pData);

                                            double[] param = new double[_modelItemCount + 1];
                                            x = mcgs(_dataRecordCount, _modelItemCount + 1, pY, pData, param);

                                            x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                            if (x < AIC)
                                            {
                                                for (int n = 0; n <= modelItemCount; ++n)
                                                {
                                                    initialGuess[n] = param[n];
                                                }

                                                AIC = x;
                                                Array.Copy(oneTerm0, oneTerm, 3);
                                                findMinAIC = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (findMinAIC)
                        {
                            Array.Copy(initialGuess, parameters, modelItemCount + 1);
                            aicCreterion = AIC;
                            residual = Math.Sqrt(Math.Exp(AIC - 0.001 * (modelItemCount + 1.0))) / _dataRecordCount;

                            _itemVariableTerms[item, modelItemCount, 0] = oneTerm[0];
                            _itemVariableTerms[item, modelItemCount, 1] = oneTerm[1];
                            _itemVariableTerms[item, modelItemCount, 2] = oneTerm[2];

                            for (int n = 0; n < modelItemCount + 1; ++n)
                            {
                                _itemParameters[item, n] = parameters[n];
                            }

                            modelItemCount++;
                        }
                    }

                } // for 8
            } // while

            double[] res = new double[_dataRecordCount];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                //double parameter0 = CalculateParameter(0, i);
                for (int j = 0; j < 8; ++j)
                {
                    _parameters[j] = CalculateParameter(j, i);
                }

                double T = _independentVariables[i, 1]; // - 185.0) / 1200.0;
                double N = Math.Exp(_parameters[1] * (T + _parameters[2]));
                double D = Math.Exp(-_parameters[3] * (T - _parameters[4])) + _parameters[5] * Math.Exp(_parameters[6] * (T - _parameters[7]));
                res[i] = _dependentVariable[i] - _parameters[0] * N / D;
            }
        }

        public void FitModel()
        {
            double paramWeight = 0.11;

            double x;
            int modelItemCount = 0;
            double AIC = double.MaxValue;

            double aicCreterion0 = double.MaxValue;
            double aicCreterion = double.MaxValue / 2.0;
            int[] oneTerm = { 0, 0, 0 };

            double[] initialGuess = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

            double[] parameters = new double[50];
            double residual = double.MaxValue;

            // build model
            while (aicCreterion < aicCreterion0 && modelItemCount < 50)
            {
                aicCreterion0 = aicCreterion;

                bool findMinAIC = false;

                // first order terms
                for (int i = 0; i < _numberOfVariables; ++i)
                {
                    int[] oneTerm0 = { i, 0, 0 };
                    if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                    {
                        SetModelTermCount(modelItemCount);
                        SetVariableTerms(oneTerm0);

                        double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];

                        SetupModelData(pData);

                        double[] param = new double[_modelItemCount + 1];
                        x = mcgs(_dataRecordCount, _modelItemCount + 1, _dependentVariable, pData, param);

                        x = _dataRecordCount * Math.Log(x / _dataRecordCount) + paramWeight * (modelItemCount + 1.0);
                        if (x < AIC)
                        {
                            for (int k = 0; k <= modelItemCount; ++k)
                            {
                                initialGuess[k] = param[k];
                            }

                            AIC = x;
                            Array.Copy(oneTerm0, oneTerm, 3);
                            findMinAIC = true;
                        }
                    }
                }

                // second order terms
                for (int i = 1; i < _numberOfVariables; ++i)
                {
                    for (int j = i; j < _numberOfVariables; ++j)
                    {
                        int[] oneTerm0 = { i, j, 0 };
                        if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                        {
                            SetModelTermCount(modelItemCount);
                            SetVariableTerms(oneTerm0);

                            double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];
                            SetupModelData(pData);

                            double[] param = new double[_modelItemCount + 1];
                            x = mcgs(_dataRecordCount, _modelItemCount + 1, _dependentVariable, pData, param);

                            x = _dataRecordCount * Math.Log(x / _dataRecordCount) + paramWeight * (modelItemCount + 1.0);
                            if (x < AIC)
                            {
                                for (int k = 0; k <= modelItemCount; ++k)
                                {
                                    initialGuess[k] = param[k];
                                }

                                AIC = x;
                                Array.Copy(oneTerm0, oneTerm, 3);
                                findMinAIC = true;
                            }
                        }
                    }
                }

                // third order terms
                //if (modelItemCount > 5)
                {
                    for (int i = 1; i < _numberOfVariables; ++i)
                    {
                        for (int j = i; j < _numberOfVariables; ++j)
                        {
                            for (int k = j; k < _numberOfVariables; ++k)
                            {
                                int[] oneTerm0 = { i, j, k };
                                if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                                {
                                    SetModelTermCount(modelItemCount);
                                    SetVariableTerms(oneTerm0);

                                    double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];
                                    SetupModelData(pData);

                                    double[] param = new double[_modelItemCount + 1];
                                    x = mcgs(_dataRecordCount, _modelItemCount + 1, _dependentVariable, pData, param);

                                    x = _dataRecordCount * Math.Log(x / _dataRecordCount) + paramWeight * (modelItemCount + 1.0);
                                    if (x < AIC)
                                    {
                                        for (int n = 0; n <= modelItemCount; ++n)
                                        {
                                            initialGuess[n] = param[n];
                                        }

                                        AIC = x;
                                        Array.Copy(oneTerm0, oneTerm, 3);
                                        findMinAIC = true;
                                    }
                                }
                            }
                        }
                    }

                    if (findMinAIC)
                    {
                        Array.Copy(initialGuess, parameters, modelItemCount + 1);
                        aicCreterion = AIC;
                        residual = Math.Sqrt(Math.Exp((AIC / _dataRecordCount - paramWeight * (modelItemCount + 1.0)) / _dataRecordCount) * _dataRecordCount) / _dataRecordCount;

                        _variableTerms[modelItemCount, 0] = oneTerm[0];
                        _variableTerms[modelItemCount, 1] = oneTerm[1];
                        _variableTerms[modelItemCount, 2] = oneTerm[2];

                        for (int n = 0; n < modelItemCount + 1; ++n)
                        {
                            _parameters[n] = parameters[n];
                        }

                        modelItemCount++;
                    }
                }
            }

            // remove items
            int rc = 0;
            while (rc < 15)
            {
                double tempAIC = AIC;
                int savedModelItemCount = _modelItemCount;
                int[,] savedVariableTerms = new int[50, 3];
                for (int i = 0; i < 50; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        savedVariableTerms[i, j] = _variableTerms[i, j];
                    }
                }

                int ITEM = -1;
                //AIC = double.MaxValue;
                for (int n = 0; n < savedModelItemCount; ++n)
                {
                    SetModelTermCount(_modelItemCount - 1);
                    RemoveVariableTerm(n);

                    double[,] pData = new double[_dataRecordCount, _modelItemCount];

                    int tempModelItemCount = _modelItemCount;
                    --_modelItemCount;
                    SetupModelData(pData);
                    _modelItemCount = tempModelItemCount;

                    double[] param = new double[_modelItemCount];
                    x = mcgs(_dataRecordCount, _modelItemCount, _dependentVariable, pData, param);

                    x = _dataRecordCount * Math.Log(x / _dataRecordCount) + paramWeight * _modelItemCount;
                    if (x < AIC)
                    {
                        AIC = x;
                        ITEM = n;
                        for (int n1 = 0; n1 < _modelItemCount; ++n1)
                        {
                            _parameters[n1] = param[n1];
                        }
                    }

                    _modelItemCount = savedModelItemCount;
                    for (int i = 0; i < 50; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            _variableTerms[i, j] = savedVariableTerms[i, j];
                        }
                    }
                }

                if (ITEM >= 0)
                {
                    --savedModelItemCount;
                    --_modelItemCount;

                    RemoveVariableTerm(ITEM);
                }
                else
                    break;

                ++rc;
            }

            //ReFitMe();
        }

        public void ReFitMe()
        {
            var solverParams = new NelderMeadSolverParams();
            NelderMeadSolver solver = new NelderMeadSolver();
            int vidRow;
            int[] vidVaribales = new int[_modelItemCount];

            //add variables
            for (int i = 0; i < _modelItemCount; ++i)
            {
                solver.AddVariable(null, out vidVaribales[i]);
            }

            //set bounds
            for (int i = 0; i < _modelItemCount; ++i)
            {
                solver.SetBounds(vidVaribales[i], -1.0E5, 1.0E5);
            }

            //set initial values
            for (int i = 0; i < _modelItemCount; ++i)
            {
                solver.SetProperty(SolverProperties.VariableStartValue, vidVaribales[i], _parameters[i]);
            }

            //add a row and set it as the goal
            solver.AddRow(null, out vidRow);
            solver.AddGoal(vidRow, 0, true);

            solver.FunctionEvaluator = EvalReFitFunction;
            //solver.GradientEvaluator = EvalGradient;
            var nms = solver.Solve(solverParams);
            Console.WriteLine(solver.ToString());
        }

        public void FitNewExpModel()
        {
            _itemParameters = new double[3, 50];
            for (int item = 0; item < 3; ++item)
            {
                _itemParameters[item, 0] = _parameters[item];
            }

            int count = 0;
            //while (count < 5)
            {
                ++count;

                for (int item = 0; item < 3; ++item)
                {
                    double x;
                    int modelItemCount = 0;
                    double AIC = double.MaxValue;

                    double aicCreterion0 = double.MaxValue;
                    double aicCreterion = double.MaxValue / 2.0;
                    int[] oneTerm = { 0, 0, 0 };

                    double[] initialGuess = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

                    double[] parameters = new double[50];
                    double residual = double.MaxValue;

                    double[] pY = CalculateNewExpDependentVariable(item);

                    while (aicCreterion < aicCreterion0)
                    {
                        aicCreterion0 = aicCreterion;

                        bool findMinAIC = false;

                        // first order terms
                        for (int i = 0; i < _numberOfVariables; ++i)
                        {
                            int[] oneTerm0 = { i, 0, 0 };
                            if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                            {
                                SetModelTermCount(item, modelItemCount);
                                SetVariableTerms(item, oneTerm0);

                                double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];

                                SetupModelData(item, pData);

                                double[] param = new double[_modelItemCounts[item] + 1];
                                x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                if (x < AIC)
                                {
                                    for (int k = 0; k <= modelItemCount; ++k)
                                    {
                                        initialGuess[k] = param[k];
                                    }

                                    AIC = x;
                                    Array.Copy(oneTerm0, oneTerm, 3);
                                    findMinAIC = true;
                                }
                            }
                        }

                        // second order terms
                        for (int i = 1; i < _numberOfVariables; ++i)
                        {
                            for (int j = i; j < _numberOfVariables; ++j)
                            {
                                int[] oneTerm0 = { i, j, 0 };
                                if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                {
                                    SetModelTermCount(item, modelItemCount);
                                    SetVariableTerms(item, oneTerm0);

                                    double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                    SetupModelData(item, pData);

                                    double[] param = new double[_modelItemCounts[item] + 1];
                                    x = mcgs(_dataRecordCount, _modelItemCounts[item] + 1, pY, pData, param);

                                    x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                    if (x < AIC)
                                    {
                                        for (int k = 0; k <= modelItemCount; ++k)
                                        {
                                            initialGuess[k] = param[k];
                                        }

                                        AIC = x;
                                        Array.Copy(oneTerm0, oneTerm, 3);
                                        findMinAIC = true;
                                    }
                                }
                            }
                        }

                        // third order terms
                        //if (modelItemCount > 5)
                        {
                            for (int i = 1; i < _numberOfVariables; ++i)
                            {
                                for (int j = i; j < _numberOfVariables; ++j)
                                {
                                    for (int k = j; k < _numberOfVariables; ++k)
                                    {
                                        int[] oneTerm0 = { i, j, k };
                                        if (!IsTermIncluded(item, oneTerm0, _itemVariableTerms, modelItemCount))
                                        {
                                            SetModelTermCount(item, modelItemCount);
                                            SetVariableTerms(item, oneTerm0);

                                            double[,] pData = new double[_dataRecordCount, _modelItemCounts[item] + 1];
                                            SetupModelData(item, pData);

                                            double[] param = new double[_modelItemCount + 1];
                                            x = mcgs(_dataRecordCount, _modelItemCount + 1, pY, pData, param);

                                            x = Math.Log(x) + 0.001 * (modelItemCount + 1.0);
                                            if (x < AIC)
                                            {
                                                for (int n = 0; n <= modelItemCount; ++n)
                                                {
                                                    initialGuess[n] = param[n];
                                                }

                                                AIC = x;
                                                Array.Copy(oneTerm0, oneTerm, 3);
                                                findMinAIC = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (findMinAIC)
                        {
                            Array.Copy(initialGuess, parameters, modelItemCount + 1);
                            aicCreterion = AIC;
                            residual = Math.Sqrt(Math.Exp(AIC - 0.001 * (modelItemCount + 1.0))) / _dataRecordCount;

                            _itemVariableTerms[item, modelItemCount, 0] = oneTerm[0];
                            _itemVariableTerms[item, modelItemCount, 1] = oneTerm[1];
                            _itemVariableTerms[item, modelItemCount, 2] = oneTerm[2];

                            for (int n = 0; n < modelItemCount + 1; ++n)
                            {
                                _itemParameters[item, n] = parameters[n];
                            }

                            modelItemCount++;
                        }
                    }

                } // for 8
            } // while

            double[] res = new double[_dataRecordCount];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                //double parameter0 = CalculateParameter(0, i);
                for (int j = 0; j < 3; ++j)
                {
                    _parameters[j] = CalculateParameter(j, i);
                }

                double T = _independentVariables[i, 1]; // - 185.0) / 1200.0;
                res[i] = _dependentVariable[i] - 1.0 / (_parameters[0] + Math.Exp(_parameters[1] * (T - _parameters[2])));
            }

        }

        private double EvalExtFunction(INonlinearModel model, int rowVid, ValuesByIndex values, bool newValues)
        {
            double obj = 0.0;

            for (int i = 0; i < _dataRecordCount; ++i)
            {
                // T
                double T = (_independentVariables[0, i] - 185.0) / 1200.0;

                double yp = values[1] * Math.Exp(values[2] * (T + values[3]));
                yp /= Math.Exp(-values[4] * (T - values[5])) + values[6] * Math.Exp(values[7] * (T - values[8]));

                double y = 100.0 * _dependentVariable[i];
                obj += (y - yp) * (y - yp);

                //return x[0] * exp(x[1] * (T + x[2]) + x[8] * T * T + x[9] * T*T*T) / (exp(-x[3] * (T + x[4])) + x[7] * exp(x[5] * (T + x[6])));
                //5.9823920361110536 - 0.021416168386291994	0.4236341547994758	0.067851588693342887 - 0.10541848089095807	0.052406523716965019	0.3474224248234688	0.536865899999497 - 0.16578961041102355	0.039836182836672077 
            }
            for (int i = 0; i < 8; ++i)
            {
                _parameters[i] = values[i + 1];
            }

            return obj;
        }

        public double mcgs(int dataRow, int dataCol, double[] pY, double[,] pData, double[] x)
        {
            int i, j, k, k1;
            double temp = 0.0;

            double[] g = new double[dataCol];

            double[,] alpha = new double[dataRow, dataCol];

            double[,] w = new double[dataRow, dataCol];

            double[] res = new double[dataRow];
            double[] theta = new double[dataCol];

            double[] pTempY = new double[dataRow];
            for (i = 0; i < dataRow; i++)
            {
                pTempY[i] = pY[i];
            }

            double[,] pTempData = new double[dataRow, dataCol];

            for (i = 0; i < dataRow; i++)
            {
                for (j = 0; j < dataCol; j++)
                {
                    pTempData[i,j] = pData[i,j];
                }
            }

            for (k = 0; k < dataCol; k++)
            {
                if (k < dataCol - 1)
                {
                    for (j = 0; j < dataRow; j++)
                    {
                        w[j,k] = pTempData[j,k];
                    }

                    for (i = k + 1; i < dataCol; i++)
                    {
                        temp = 0.0;
                        alpha[k,i] = 0.0;

                        for (k1 = 0; k1 < dataRow; k1++)
                        {
                            temp = temp + w[k1,k] * w[k1,k];
                            alpha[k,i] = alpha[k,i] + w[k1,k] * pTempData[k1,i];
                        }
                        if (temp < 1.0e-20)
                        {
                            // printf("MCGS numerical error !\n");
                            return (-1.0); /* protecting overflow of dividing */
                        }

                        alpha[k,i] = alpha[k,i] / temp;

                        for (k1 = 0; k1 < dataRow; k1++)
                        {
                            pTempData[k1,i] = pTempData[k1,i] - alpha[k,i] * w[k1,k];
                        }
                    }
                }

                if (k == dataCol - 1)
                {
                    temp = 0.0;
                    for (k1 = 0; k1 < dataRow; k1++)
                    {
                        w[k1,dataCol - 1] = pTempData[k1,dataCol - 1];
                        temp = temp + w[k1,k] * w[k1,k];
                    }
                }

                g[k] = 0.0;
                for (k1 = 0; k1 < dataRow; k1++)
                {
                    g[k] = g[k] + w[k1,k] * pTempY[k1];
                 }
                if (temp < 1.0e-20)
                    return (-1.0); /* protecting overflow of dividing */
                g[k] = g[k] / temp;

                for (k1 = 0; k1 < dataRow; k1++)
                {
                    pTempY[k1] = pTempY[k1] - g[k] * w[k1,k];
                }
            }

            for (i = dataCol - 1; i > -1; i--)
            {
                theta[i] = g[i];
                for (j = i; j < dataCol - 1; j++)
                {
                    theta[i] = theta[i] - theta[j + 1] * alpha[i,j + 1];
                }
                x[i] = theta[i];
            }

            temp = 0.0;
            for (i = 0; i < dataRow; i++)
            {
                res[i] = pY[i];
                for (j = 0; j < dataCol; j++)
                {
                    res[i] = res[i] - theta[j] * pData[i,j];
                }
                temp = temp + res[i] * res[i];
            }

            return temp;
        }

        void SetModelingData(int row, int col, double[] pData1, double[] pData2)
        {
            Array.Copy(pData2, _dependentVariable, row);
            Array.Copy(pData1, _independentVariables, row*col);
        }
        void SetRawData(double[] data)
        {
            _rawDependentVariable = data;
        }
        void SetTermOrder(int termOrder)
        {
            _termOrder = termOrder;
        }

        void SetModelTermCount(int modelItemCount)
        {
            _modelItemCount = modelItemCount;
        }
        void SetModelTermCount(int item, int modelItemCount)
        {
            _modelItemCounts[item] = modelItemCount;
        }

        void SetNumberOfVariables(int numberOfVariables)
        {
            _numberOfVariables = numberOfVariables;
        }

        void SetDataRecordCount(int count)
        {
            _dataRecordCount = count;
        }

        void SetVariableTerms(int[] pOneTerm)
        {
            _variableTerms[_modelItemCount,0] = pOneTerm[0];
            _variableTerms[_modelItemCount,1] = pOneTerm[1];
            _variableTerms[_modelItemCount,2] = pOneTerm[2];
        }
        void SetVariableTerms(int item, int[] pOneTerm)
        {
            _itemVariableTerms[item, _modelItemCounts[item], 0] = pOneTerm[0];
            _itemVariableTerms[item, _modelItemCounts[item], 1] = pOneTerm[1];
            _itemVariableTerms[item, _modelItemCounts[item], 2] = pOneTerm[2];
        }

        void RemoveVariableTerm(int item)
        {
            for (int i = item+1; i < 50; ++i)
            {
                _variableTerms[i-1, 0] = _variableTerms[i, 0];
                _variableTerms[i-1, 1] = _variableTerms[i, 1];
                _variableTerms[i-1, 2] = _variableTerms[i, 2];
            }
        }

        void SetupModelData(double[,] pData)
        {
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                for (int k = 0; k <= _modelItemCount; ++k)
                {
                    int k1 = _variableTerms[k, 0];
                    int k2 = _variableTerms[k, 1];
                    int k3 = _variableTerms[k, 2];

                    double d1 = _independentVariables[i, k1];
                    double d2 = _independentVariables[i, k2];
                    double d3 = _independentVariables[i, k3];

                    pData[i, k] = d1 * d2 * d3;
                }
            }
        }

        void SetupModelData(int item, double[,] pData)
        {
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                for (int k = 0; k <= _modelItemCounts[item]; ++k)
                {
                    int k1 = _itemVariableTerms[item, k, 0];
                    int k2 = _itemVariableTerms[item, k, 1];
                    int k3 = _itemVariableTerms[item, k, 2];

                    double d1 = _independentVariables[i, k1];
                    double d2 = _independentVariables[i, k2];
                    double d3 = _independentVariables[i, k3];

                    pData[i, k] = d1 * d2 * d3;
                }
            }
        }

        public bool IsTermIncluded(int[] oneTerm, int[,] variableTerms, int modelItemCount)
        {
            if (modelItemCount == 0)
                return false;

            for (int i = 0; i < modelItemCount; ++i)
            {
                if (oneTerm[0] == variableTerms[i,0] && oneTerm[1] == variableTerms[i,1] && oneTerm[2] == variableTerms[i,2])
                    return true;
            }
            return false;
        }
        public bool IsTermIncluded(int item, int[] oneTerm, int[,,] variableTerms, int modelItemCount)
        {
            if (modelItemCount == 0)
                return false;

            for (int i = 0; i < modelItemCount; ++i)
            {
                if (oneTerm[0] == variableTerms[item,i, 0] && oneTerm[1] == variableTerms[item,i, 1] && oneTerm[2] == variableTerms[item,i, 2])
                    return true;
            }
            return false;
        }


        public void Start()
        {
            double x;
            int modelItemCount = 0;
            double AIC = double.MaxValue;

            double aicCreterion0 = double.MaxValue;
            double aicCreterion = double.MaxValue / 2.0;
            int[] oneTerm = { 0, 0, 0 };

            double[] initialGuess = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            
            double[] parameters = new double[50];
            double residual = double.MaxValue;

            while (aicCreterion < aicCreterion0)
            {
                aicCreterion0 = aicCreterion;

                bool findMinAIC = false;

                // first order terms
                for (int i = 0; i < _numberOfVariables; ++i)
                {
                    int[] oneTerm0 = { i, 0, 0 };
                    if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                    {
                        SetModelTermCount(modelItemCount);
                        SetVariableTerms(oneTerm0);

                        double[] pY = new double[_dataRecordCount];
                        Array.Copy(_dependentVariable, pY, _dataRecordCount);
                       
                        double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];

                        SetupModelData(pData);

                        double[] param = new double[_modelItemCount + 1];
                        x = mcgs(_dataRecordCount, _modelItemCount + 1, pY, pData, param);

                        x = Math.Log(x) + 2.0 * (modelItemCount + 1.0);
                        if (x < AIC)
                        {
                            for (int k = 0; k <= modelItemCount; ++k)
                            {
                                initialGuess[k] = param[k];
                            }

                            AIC = x;
                            Array.Copy(oneTerm0, oneTerm, 3);
                            findMinAIC = true;
                        }
                    }
                }

                // second order terms
                for (int i = 1; i < _numberOfVariables; ++i)
                {
                    for (int j = i; j < _numberOfVariables; ++j)
                    {
                        int[] oneTerm0 = { i, j, 0 };
                        if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                        {
                            SetModelTermCount(modelItemCount);
                            SetVariableTerms(oneTerm0);

                            double[] pY = new double[_dataRecordCount];
                            Array.Copy(_dependentVariable, pY, _dataRecordCount);

                            double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];

                            SetupModelData(pData);

                            double[] param = new double[_modelItemCount + 1];
                            x = mcgs(_dataRecordCount, _modelItemCount + 1, pY, pData, param);


                            x = Math.Log(x) + 2.0 * (modelItemCount + 1.0);
                            if (x < AIC)
                            {
                                for (int k = 0; k <= modelItemCount; ++k)
                                {
                                    initialGuess[k] = param[k];
                                }

                                AIC = x;
                                Array.Copy(oneTerm0, oneTerm, 3);
                                findMinAIC = true;
                            }
                        }
                    }
                }

                // third order terms
                if (modelItemCount > 5)
                {
                    for (int i = 1; i < _numberOfVariables; ++i)
                    {
                        for (int j = i; j < _numberOfVariables; ++j)
                        {
                            for (int k = j; k < _numberOfVariables; ++k)
                            {
                                int[] oneTerm0 = { i, j, k };
                                if (!IsTermIncluded(oneTerm0, _variableTerms, modelItemCount))
                                {
                                    SetModelTermCount(modelItemCount);
                                    SetVariableTerms(oneTerm0);

                                    double[] pY = new double[_dataRecordCount];
                                    Array.Copy(_dependentVariable, pY, _dataRecordCount);

                                    double[,] pData = new double[_dataRecordCount, _modelItemCount + 1];

                                    SetupModelData(pData);

                                    double[] param = new double[_modelItemCount + 1];
                                    x = mcgs(_dataRecordCount, _modelItemCount + 1, pY, pData, param);

                                    x = Math.Log(x) + 2.0 * (modelItemCount + 1.0);
                                    if (x < AIC)
                                    {
                                        for (int n = 0; n <= modelItemCount; ++n)
                                        {
                                            initialGuess[n] = param[n];
                                        }

                                        AIC = x;
                                        Array.Copy(oneTerm0, oneTerm, 3);
                                        findMinAIC = true;
                                    }
                                }
                            }
                        }
                    }
                }

                if (findMinAIC)
                {
                    Array.Copy(initialGuess, parameters, modelItemCount + 1);
                    aicCreterion = AIC;
                    residual = Math.Sqrt(Math.Exp(AIC - 2.0 * (modelItemCount + 1.0))) / _dataRecordCount;

                    _variableTerms[modelItemCount, 0] = oneTerm[0];
                    _variableTerms[modelItemCount, 1] = oneTerm[1];
                    _variableTerms[modelItemCount, 2] = oneTerm[2];

                    modelItemCount++;
                }
            }
        }

        public void FitPPModel(int numberOfCuts)
        {
            var solverParams = new NelderMeadSolverParams();
            NelderMeadSolver solver = new NelderMeadSolver();
            int vidRow;
            int[] vidVaribales = new int[numberOfCuts];

            //add variables
            for (int i = 0; i < numberOfCuts; ++i)
            {
                solver.AddVariable(null, out vidVaribales[i]);
            }

            //set bounds
            for (int i = 0; i < numberOfCuts; ++i)
            {
                solver.SetBounds(vidVaribales[i], 0.8, 1.2);
            }

            //set initial values
            for (int i = 0; i < numberOfCuts; ++i)
            {
                solver.SetProperty(SolverProperties.VariableStartValue, vidVaribales[i], 1.0);
            }

            //add a row and set it as the goal
            solver.AddRow(null, out vidRow);
            solver.AddGoal(vidRow, 0, true);

            solver.FunctionEvaluator = EvalExpFunction;
            //solver.GradientEvaluator = EvalGradient;
            var nms = solver.Solve(solverParams);
            Console.WriteLine(solver.ToString());

            // Delete the file if it exists.
            string strFile = @"C:\data.dat";
            if (File.Exists(strFile))
            {
                File.Delete(strFile);
            }

            FileStream fs = new FileStream(strFile, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);

            double x = 0.0;
            double[] res = new double[1354];
            for (int i = 0; i < _dataRecordCount; ++i)
            {
                double T = _independentVariables[i, 1]; // - 185.0) / 1200.0;
                double yp = 1.0 / (_parameters[0] + Math.Exp(_parameters[1] * (T - _parameters[2])));

                res[i] = _dependentVariable[i] - yp;
                x += Math.Abs(res[i]);

                sw.WriteLine("{0}, {1}, {2}, {3}", _independentVariables[i, 1], _dependentVariable[i], yp, res[i]);
            }
            x /= _dataRecordCount;

            sw.Close();
        }

        private double EvalPPFunction(INonlinearModel model, int rowVid, ValuesByIndex values, bool newValues)
        {
            double obj = 0.0;
            int numberOfCuts = 5;
            int dataRecordCount = 5;

            double m1 = 5.21847;
            double m2 = 4.75195;
            double[] r2Yields = {12.67957933, 13.06306916, 15.68171093, 3.568119376, 56.22002506};
            double[] r3Yields = {18.90663233, 22.51575479, 18.34551976, 2.830675389, 37.19562219};
            double[] ovYields = {28.15965709, 31.11094062, 25.2072688, 3.105008641, 6.363347871};

            for (int i = 0; i < dataRecordCount; ++i)
            {

                //obj += (values[i+1] * r2Yields[i] + values[i + 1] * (m2/m1)*r3Yields[i] - ovYields[i]) * (y - yp);
            }

            for (int i = 0; i < numberOfCuts; ++i)
            {
                _parametersPP[i] = values[i + 1];
            }

            return obj;
        }


    }
}
