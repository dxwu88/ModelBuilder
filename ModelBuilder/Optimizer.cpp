// Optimizer.cpp: implementation of the COptimizer class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Optimizer.h"
//#include "frontkey.h"

#ifdef _DEBUG
//#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

//#pragma warning(disable : 4244 4800 4305)
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

COptimizer *pOptimizerPtr=NULL;

COptimizer::COptimizer()
{
	m_DataRecordCount = 10;
	pOptimizerPtr = this;
}

COptimizer::~COptimizer()
{
	delete pOptimizerPtr;
}

COptimizer::COptimizer(int termOrder, int modelItemCount, int numberOfVariables, int* pVariableTerms, double* pData1, double* pData2)
{
	pOptimizerPtr = this;

	m_NumberOfVariables = numberOfVariables;

	m_TermOrder = termOrder;
	m_ModelItemCount = modelItemCount;
	memcpy(&m_VariableTerms[0][0], pVariableTerms, sizeof(m_VariableTerms));

	memcpy(&m_IndependentVariables[0][0], pData1, sizeof(m_IndependentVariables));
	memcpy(&m_DependentVariable[0], pData2, sizeof(m_DependentVariable));
}

void COptimizer::SetModelingData(int row, int col, double** pData1, double* pData2)
{
	memcpy(&m_DependentVariable[0], pData2, sizeof(m_DependentVariable));	
	//memcpy(&m_IndependentVariables[0][0], pData1, sizeof(m_IndependentVariables));
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			m_IndependentVariables[i][j] = pData1[i][j];
		}
	}
}

void COptimizer::SetTermOrder(int termOrder)
{
	m_TermOrder = termOrder;
}

void COptimizer::SetModelTermCount(int modelItemCount)
{
	m_ModelItemCount = modelItemCount;
}

void COptimizer::SetNumberOfVariables(int numberOfVariables)
{
	m_NumberOfVariables = numberOfVariables;
}

void COptimizer::SetDataRecordCount(int count)
{
	m_DataRecordCount = count;
}

void COptimizer::SetVariableTerms(int* pOneTerm)
{
	m_VariableTerms[m_ModelItemCount][0] = pOneTerm[0];
	m_VariableTerms[m_ModelItemCount][1] = pOneTerm[1];
	m_VariableTerms[m_ModelItemCount][2] = pOneTerm[2];
}

bool COptimizer::IsTermIncluded(int oneTerm[3], int variableTerms[50][3], int modelItemCount)
{
	if (modelItemCount == 0)
		return false;

	for (int i = 0; i < modelItemCount; ++i)
	{
		if (oneTerm[0] == variableTerms[i][0] && oneTerm[1] == variableTerms[i][1] && oneTerm[2] == variableTerms[i][2])
			return true;
	}
	return false;
}

double COptimizer::FindNextModelFit(int oneTerm[3], double* pInitialGuess)
{
	int modelItemCount = m_ModelItemCount + 1;

	double* objDer = new double[modelItemCount];
	double* lbDer = new double[modelItemCount];
	double* ubDer = new double[modelItemCount];
	double* x = new double[modelItemCount];

	for (int i = 0; i < modelItemCount; ++i)
	{
		lbDer[i] = -INFBOUND;
		ubDer[i] = INFBOUND;

		x[i] = pInitialGuess[i];
	}

	double rhs[1] = { 0.0 };
	char sense[1];
	strncpy(sense, "E", 1);
	double matval[225];

	long statDer;
	double objvalDer;
	double* djDer = new double[modelItemCount];
	double* pioutDer = new double[modelItemCount];
	double* slackDer = new double[modelItemCount];

	HPROBLEM lpDer = NULL;
	setintparam(lpDer, PARAM_ARGCK, 1);
	int maxloopDer = 0;

	lpDer = loadnlp(PROBNAME, 2, 0, 1, objDer, NULL, (unsigned char*)sense,
		NULL, NULL, NULL, NULL, x, lbDer, ubDer, NULL, NULL,
		EvalFunc, NULL);

	if (!lpDer) return false;

	setintparam(lpDer, PARAM_DERIV, 1);
	setdblparam(lpDer, PARAM_EPCONV, 1.0E-8);
	if (optimize(lpDer)) { // Did not find the solution
		unloadprob(&lpDer);
	}
	else {
		solution(lpDer, &statDer, &objvalDer, x, pioutDer, slackDer, djDer);
		unloadprob(&lpDer);

		for (int i = 0; i < modelItemCount; i++)
		{
			double z = x[i];
			pInitialGuess[i] = x[i];
		}
	}

	memcpy(pInitialGuess, &x[0], modelItemCount);
	return objvalDer;
}

void COptimizer::StartNLP()
{
	double x;
	int modelItemCount = 0;
	double AIC = DBL_MAX;

	double aicCreterion0 = DBL_MAX;
	double aicCreterion = DBL_MAX / 2.0;
	int oneTerm[] = { 0, 0, 0 };

	double InitialGuess[10] = { 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
	double parameters[10];

	while (aicCreterion < aicCreterion0)
	{
		aicCreterion0 = aicCreterion;

		bool findMinAIC = false;

		// first order terms
		for (int i = 0; i < m_NumberOfVariables; ++i)
		{
			int oneTerm0[] = { i, 0, 0 };
			if ( !IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount) )
			{ 
				//SetTermOrder(1);
				SetModelTermCount(modelItemCount);
				SetVariableTerms(&oneTerm0[0]);

				x = FindNextModelFit(oneTerm0, &InitialGuess[0]);
				x = log(x) + 0.0 * (modelItemCount + 1);
				if (x < AIC)
				{
					AIC = x;
					memcpy(oneTerm, oneTerm0, 3);
					findMinAIC = true;
				}
			}
		}

		// second order terms
		for (int i = 1; i < m_NumberOfVariables; ++i)
		{
			for (int j = i; j < m_NumberOfVariables; ++j)
			{
				int oneTerm0[] = { i, j, 0 };
				if (!IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount))
				{
					//SetTermOrder(2);
					SetModelTermCount(modelItemCount);
					SetVariableTerms(&oneTerm0[0]);

					x = FindNextModelFit(oneTerm0, &InitialGuess[0]);
					x = log(x) + 0.0 * (modelItemCount + 1);
					if (x < AIC)
					{
						AIC = x;
						memcpy(oneTerm, oneTerm0, 3);
						findMinAIC = true;
					}
				}
			}
		}
/***
		// third order terms
		for (int i = 1; i < m_NumberOfVariables; ++i)
		{
			for (int j = i; j < m_NumberOfVariables; ++j)
			{
				for (int k = j; k < m_NumberOfVariables; ++k)
				{
					int oneTerm0[] = { i, j, k };
					if (!IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount))
					{
						//SetTermOrder(3);
						SetModelTermCount(modelItemCount);
						SetVariableTerms(&oneTerm0[0]);

						x = FindNextModelFit(oneTerm0, &InitialGuess[0]);
						x = log(x) + 2.0 * (modelItemCount + 1);
						if (x < AIC)
						{
							AIC = x;
							memcpy(oneTerm, oneTerm0, 3);
							findMinAIC = true;
						}
					}
				}
			}
		}
***/
		if (findMinAIC)
		{
			memcpy(&parameters[0], InitialGuess, modelItemCount);

			aicCreterion = AIC;

			memcpy(&m_VariableTerms[modelItemCount][0], oneTerm, 3);

			modelItemCount++;
		}
	}

}


double COptimizer::FitSkewModel()
{
	int modelItemCount = 8;

	double* objDer = new double[modelItemCount];
	double* lbDer = new double[modelItemCount];
	double* ubDer = new double[modelItemCount];
	double* x = new double[modelItemCount] { 
		0.623238235,
		1.937953622,
		1.919573719,
		2.538293111,
		0.41273197,
		3.392042252,
		6.110383618,
		1.008861977
	};

	for (int i = 0; i < modelItemCount; ++i)
	{
		lbDer[i] = 0.0;
		ubDer[i] = 7.0;

		//x[i] = 1.0;
	}
	

	double rhs[1] = { 0.0 };
	char sense[1];
	strncpy(sense, "E", 1);
	double matval[225];

	long statDer;
	double objvalDer;
	double* djDer = new double[modelItemCount];
	double* pioutDer = new double[modelItemCount];
	double* slackDer = new double[modelItemCount];

	HPROBLEM lpDer = NULL;
	setintparam(lpDer, PARAM_ARGCK, 1);
	int maxloopDer = 0;

	lpDer = loadnlp(PROBNAME, 8, 0, 1, objDer, NULL, NULL,
		NULL, NULL, NULL, NULL, x, lbDer, ubDer, NULL, NULL,
		EvalSkewFunc, NULL);

	if (!lpDer) return false;

	setintparam(lpDer, PARAM_DERIV, 1);
	setdblparam(lpDer, PARAM_EPCONV, 1.0E-8);
	if (optimize(lpDer)) { // Did not find the solution
		unloadprob(&lpDer);
	}
	else {
		solution(lpDer, &statDer, &objvalDer, x, pioutDer, slackDer, djDer);
		unloadprob(&lpDer);

		for (int i = 0; i < modelItemCount; i++)
		{
			double z = x[i];
		}
	}

	return objvalDer;
}

long COptimizer::EvalSkewFunc(HPROBLEM lp, INTARG numcols,
   INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG x,
   INTARG varone, INTARG vartwo)
{
	objval[0] = 0.0;

	for (int i = 0; i < pOptimizerPtr->m_DataRecordCount; ++i)
	{
		// T
		double T = pOptimizerPtr->m_IndependentVariables[i][1]; // (pOptimizerPtr->m_IndependentVariables[i][0] - 185.0) / 1200.0;

		double yp = x[0] *exp(x[1]*(T + x[2]));
		yp /= exp(x[3] * (T - x[4])) + x[5] * exp(x[6] * (T - x[7]));

		double y = pOptimizerPtr->m_DependentVariable[i];
		objval[0] += (y - yp) * (y - yp);
	
		//return x[0] * exp(x[1] * (T + x[2]) + x[8] * T * T + x[9] * T*T*T) / (exp(-x[3] * (T + x[4])) + x[7] * exp(x[5] * (T + x[6])));
		//5.9823920361110536 - 0.021416168386291994	0.4236341547994758	0.067851588693342887 - 0.10541848089095807	0.052406523716965019	0.3474224248234688	0.536865899999497 - 0.16578961041102355	0.039836182836672077 
	}

	return 0;
}

long COptimizer::EvalFunc(HPROBLEM lp, INTARG numcols,
	INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG x,
	INTARG varone, INTARG vartwo)
{
	objval[0] = 0.0;

	for (int i = 0; i < pOptimizerPtr->m_DataRecordCount; ++i)
	{
		double yp = 0.0;
		for (int k = 0; k <= pOptimizerPtr->m_ModelItemCount; ++k)
		{
			int k1 = pOptimizerPtr->m_VariableTerms[k][0];
			int k2 = pOptimizerPtr->m_VariableTerms[k][1];
			int k3 = pOptimizerPtr->m_VariableTerms[k][2];

			double d1 = pOptimizerPtr->m_IndependentVariables[i][k1];
			double d2 = pOptimizerPtr->m_IndependentVariables[i][k2];
			double d3 = pOptimizerPtr->m_IndependentVariables[i][k3];

			yp += x[k] * d1 * d2 * d3;
		}

		double y = pOptimizerPtr->m_DependentVariable[i];
		objval[0] += (y - yp) * (y - yp);
	}

	return 0;
}


double COptimizer::mcgs(int dataRow, int dataCol, double* pY, double** pData, double x[])
{
	int i, j, k, k1;
	double temp = 0.0;

	double* g = new double[dataCol];
		
	double** alpha = new double*[dataRow];
	for (int i = 0; i < dataRow; ++i)
	{
		alpha[i] = new double[dataCol];
	}

	double** w = new double*[dataRow];
	for (int i = 0; i < dataRow; ++i)
	{
		w[i] = new double[dataCol];
	}

	double* res = new double[dataRow];
	double*	theta = new double[dataCol];

	double* pTempY = new double[dataRow];
	memcpy(pTempY, pY, dataRow * sizeof(double));

	double **pTempData = new double*[dataRow];
	for (int i = 0; i < dataRow; ++i)
	{
		pTempData[i] = new double[dataCol];
	}

	for (i = 0; i < dataRow; i++)
	{
		for (j = 0; j < dataCol; j++)
		{
			pTempData[i][j] = pData[i][j];
		}
	}

	for ( k = 0; k < dataCol; k++)
	{
		if (k < dataCol - 1)
		{
			for (j = 0; j < dataRow; j++)
			{
				w[j][k] = pTempData[j][k];
			}

			for (i = k + 1; i < dataCol; i++) {
				temp = 0.0; 
				alpha[k][i] = 0.0;

				for (k1 = 0; k1 < dataRow; k1++)
				{
					temp = temp + w[k1][k] * w[k1][k];
					alpha[k][i] = alpha[k][i] + w[k1][k] * pTempData[k1][i];
				}
				if (temp < 1.0e-20) 
				{
					// printf("MCGS numerical error !\n");
					return (DBL_MAX); /* protecting overflow of dividing */
				}

				alpha[k][i] = alpha[k][i] / temp;

				for (k1 = 0; k1 < dataRow; k1++)
				{
					pTempData[k1][i] = pTempData[k1][i] - alpha[k][i] * w[k1][k];
				}
			}
		}

		if (k == dataCol - 1)
		{
			temp = 0.0;
			for (k1 = 0; k1 < dataRow; k1++)
			{
				w[k1][dataCol - 1] = pTempData[k1][dataCol - 1];
				temp = temp + w[k1][k] * w[k1][k];
			}
		}

		g[k] = 0.0;
		for (k1 = 0; k1 < dataRow; k1++)
		{
			g[k] = g[k] + w[k1][k] * pTempY[k1];
		}
		if (temp < 1.0e-20) 
			return (DBL_MAX); /* protecting overflow of dividing */
		g[k] = g[k] / temp;

		for (k1 = 0; k1 < dataRow; k1++)
		{
			pTempY[k1] = pTempY[k1] - g[k] * w[k1][k];
		}
	}

	for (i = dataCol - 1; i > -1; i--)
	{
		theta[i] = g[i];
		for (j = i; j < dataCol - 1; j++)
		{
			theta[i] = theta[i] - theta[j + 1] * alpha[i][j + 1];
		}
		x[i] = theta[i];
	}

	temp = 0.0;
	for (i = 0; i < dataRow; i++)
	{
		res[i] = pY[i];
		for (j = 0; j < dataCol; j++)
		{
			res[i] = res[i] - theta[j] * pData[i][j];
		}
		temp = temp + res[i] * res[i];
	}

	delete [] g;

	for (int i = 0; i < dataRow; ++i)
	{
		delete[] alpha[i];
	}
	delete [] alpha;

	for (int i = 0; i < dataRow; ++i)
	{
		delete[] w[i];
	}
	delete[] w;

	delete [] res;
	delete [] theta;

	return (temp);
}

void COptimizer::SetupModelData(double** pData)
{
	for (int i = 0; i < pOptimizerPtr->m_DataRecordCount; ++i)
	{
		for (int k = 0; k <= pOptimizerPtr->m_ModelItemCount; ++k)
		{
			int k1 = pOptimizerPtr->m_VariableTerms[k][0];
			int k2 = pOptimizerPtr->m_VariableTerms[k][1];
			int k3 = pOptimizerPtr->m_VariableTerms[k][2];

			double d1 = pOptimizerPtr->m_IndependentVariables[i][k1];
			double d2 = pOptimizerPtr->m_IndependentVariables[i][k2];
			double d3 = pOptimizerPtr->m_IndependentVariables[i][k3];

			pData[i][k] = d1 * d2 * d3;
		}
	}
}

void COptimizer::Start()
{
	double x;
	int modelItemCount = 0;
	double AIC = DBL_MAX;

	double aicCreterion0 = DBL_MAX;
	double aicCreterion = DBL_MAX / 2.0;
	int oneTerm[] = { 0, 0, 0 };

	double InitialGuess[50] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
	double parameters[50];
	double residual = DBL_MAX;

	while (aicCreterion < aicCreterion0)
	{
		aicCreterion0 = aicCreterion;

		bool findMinAIC = false;

		// first order terms
		for (int i = 0; i < m_NumberOfVariables; ++i)
		{
			int oneTerm0[] = { i, 0, 0 };
			if (!IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount))
			{
				SetModelTermCount(modelItemCount);
				SetVariableTerms(&oneTerm0[0]);

				double* pY = new double[m_DataRecordCount];
				memcpy(pY, m_DependentVariable, m_DataRecordCount * sizeof(double));

				double **pData = AllocateDynamicArray<double>(m_DataRecordCount, m_ModelItemCount + 1);

				SetupModelData(pData);

				double* param = new double[m_ModelItemCount + 1];
				x = mcgs(m_DataRecordCount, m_ModelItemCount + 1, pY, pData, param);


				x = log(x) + 2.0 * (modelItemCount + 1);
				if (x < AIC)
				{
					for (int n = 0; n <= modelItemCount; ++n)
					{
						InitialGuess[n] = param[n];
					}

					AIC = x;
					memcpy(oneTerm, oneTerm0, 3);
					findMinAIC = true;
				}

				delete[] pY;
				delete[] param;
				FreeDynamicArray<double>(pData);
			}
		}

		// second order terms
		for (int i = 1; i < m_NumberOfVariables; ++i)
		{
			for (int j = i; j < m_NumberOfVariables; ++j)
			{
				int oneTerm0[] = { i, j, 0 };
				if (!IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount))
				{
					SetModelTermCount(modelItemCount);
					SetVariableTerms(&oneTerm0[0]);

					double* pY = new double[m_DataRecordCount];
					memcpy(pY, m_DependentVariable, m_DataRecordCount * sizeof(double));

					double **pData = AllocateDynamicArray<double>(m_DataRecordCount, m_ModelItemCount + 1);

					SetupModelData(pData);

					double* param = new double[m_ModelItemCount + 1];
					x = mcgs(m_DataRecordCount, m_ModelItemCount + 1, pY, pData, param);


					x = log(x) + 2.0 * (modelItemCount + 1);
					if (x < AIC)
					{
						for (int n = 0; n <= modelItemCount; ++n)
						{
							InitialGuess[n] = param[n];
						}

						AIC = x;
						memcpy(oneTerm, oneTerm0, 3);
						findMinAIC = true;
					}

					delete[] pY;
					delete[] param;
					FreeDynamicArray<double>(pData);
				}
			}
		}

		// third order terms
		if (modelItemCount > 5)
		{
			for (int i = 1; i < m_NumberOfVariables; ++i)
			{
				for (int j = i; j < m_NumberOfVariables; ++j)
				{
					for (int k = j; k < m_NumberOfVariables; ++k)
					{
						int oneTerm0[] = { i, j, k };
						if (!IsTermIncluded(oneTerm0, m_VariableTerms, modelItemCount))
						{
							SetModelTermCount(modelItemCount);
							SetVariableTerms(&oneTerm0[0]);

							double* pY = new double[m_DataRecordCount];
							memcpy(pY, m_DependentVariable, m_DataRecordCount * sizeof(double));

							double **pData = AllocateDynamicArray<double>(m_DataRecordCount, m_ModelItemCount + 1);

							SetupModelData(pData);

							double* param = new double[m_ModelItemCount + 1];
							x = mcgs(m_DataRecordCount, m_ModelItemCount + 1, pY, pData, param);

							x = log(x) + 2.0 * (modelItemCount + 1);
							if (x < AIC)
							{
								for (int n = 0; n <= modelItemCount; ++n)
								{
									InitialGuess[n] = param[n];
								}

								AIC = x;
								memcpy(oneTerm, oneTerm0, 3);
								findMinAIC = true;
							}

							delete[] pY;
							delete[] param;
							FreeDynamicArray<double>(pData);
						}
					}
				}
			}
		}

		if (findMinAIC)
		{
			memcpy(&parameters[0], InitialGuess, (modelItemCount + 1) * sizeof(double));

			aicCreterion = AIC;
			residual = sqrt(exp(AIC - 2.0 * (modelItemCount + 1.0))) / m_DataRecordCount;

			memcpy(&m_VariableTerms[modelItemCount][0], oneTerm, 3);

			modelItemCount++;
		}
	}
}

// x0 * exp(x1*(T-x2))
double COptimizer::CalculateExp1(double x0, double x1, double x2, double T)
{
	return x0 * exp(x1*(T - x2));
}

// exp(-x3*(T-x4))
double COptimizer::CalculateExp2(double x3, double x4, double T)
{
	return exp(-x3*(T - x4));
}

// x5 * exp(x6*(T-x7))
double COptimizer::CalculateExp3(double x5, double x6, double x7, double T)
{
	return x5 * exp(x6*(T - x7));
}

// x0 * exp(x1*(T-x2))
double COptimizer::CalculateX1(double T, double x0, double x2, double reminder)
{
	return log(reminder / x0) / (T - x2);
}

// x0 * exp(x1*(T-x2))
double COptimizer::CalculateX2(double T, double x0, double x1, double reminder)
{
	return T - log(reminder / x0) / x1;
}

// exp(-x3*(T-x4))
double COptimizer::CalculateX3(double T, double x4, double reminder)
{
	return - log(reminder) / (T + x4);
}

// exp(-x3*(T-x4))
double COptimizer::CalculateX4(double T, double x3, double reminder)
{
	return -log(reminder) / x3 - T;
}

// x5 * exp(x6*(T-x7))
double COptimizer::CalculateX6(double T, double x5, double x7, double reminder)
{
	return log(reminder / x5) / (T + x7);
}

// x5 * exp(x6*(T-x7))
double COptimizer::CalculateX7(double T, double x5, double x6, double reminder)
{
	return log(reminder / x5) / x6 - T;

}


