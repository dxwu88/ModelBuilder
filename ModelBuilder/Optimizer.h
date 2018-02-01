// Optimizer.h: interface for the COptimizer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_Optimizer_H__7F089FD1_9DDA_4986_B240_E2031EDE6084__INCLUDED_)
#define AFX_Optimizer_H__7F089FD1_9DDA_4986_B240_E2031EDE6084__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

template <typename T>
T **AllocateDynamicArray(int nRows, int nCols)
{
	T **dynamicArray;

	dynamicArray = new T*[nRows];
	for (int i = 0; i < nRows; i++)
		dynamicArray[i] = new T[nCols];

	return dynamicArray;
}

template <typename T>
void FreeDynamicArray(T** dArray)
{
	delete[] * dArray;
	delete[] dArray;
}

class COptimizer  
{
public:
	COptimizer();
	COptimizer(int termOrder, int modelItemCount, int numberOfVariables, int* pVariableTerms, double* pData1, double* pData2);

	virtual ~COptimizer();

// Attributes
	int m_NumberOfVariables;

	int m_TermOrder;
	int m_ModelItemCount;
	int m_VariableTerms[50][3];

	int m_DataRecordCount;
	double m_IndependentVariables[1353][8];
	double m_DependentVariable[1353];

//Operations
	void SetModelingData(int row, int col, double** pData1, double* pData2);
	void COptimizer::SetTermOrder(int termOrder);
	void SetModelTermCount(int modelItemCount);
	void SetNumberOfVariables(int numberOfVariables);
	void SetVariableTerms(int* pOneTerm);
	void SetDataRecordCount(int count);
	double mcgs(int dataRow, int dataCol, double* pY, double** pData, double x[]);
	void SetupModelData(double** pData);

	double FindNextModelFit(int oneTerm[3], double* pInitialGuess);
	bool IsTermIncluded(int oneTerm[3], int variableTerms[50][3], int modelItemCount);
	void StartNLP();
	void Start();

	double CalculateExp1(double x0, double x1, double x2, double T);
	double CalculateExp2(double x3, double x4, double T);
	double CalculateExp3(double x5, double x6, double x7, double T);
	double CalculateX1(double T, double x0, double x2, double reminder);
	double CalculateX2(double T, double x0, double x2, double reminder);
	double CalculateX3(double T, double x4, double reminder);
	double CalculateX4(double T, double x3, double reminder);
	double CalculateX6(double T, double x5, double x7, double reminder);
	double CalculateX7(double T, double x5, double x6, double reminder);

	double FitSkewModel();

	static long __stdcall EvalFunc(HPROBLEM lp, INTARG numcols,
	   INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG var,
	   INTARG varone, INTARG vartwo);
	static long __stdcall EvalSkewFunc(HPROBLEM lp, INTARG numcols,
		INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG var,
		INTARG varone, INTARG vartwo);


	double OptimizeModelParameter(int k, double fYieldData[300][31], double fParamsForParam[34]);
	static long __stdcall ParameterEvalFunc(HPROBLEM lp, INTARG numcols,
	   INTARG numrows, LPREALARG objval, LPREALARG lhs, LPREALARG var,
	   INTARG varone, INTARG vartwo);

	double fp(int iParam, int iData, int iCut, double fParams[34]);
	double* MCGS(int row, int col, double* pfY, double* pfX);

	double getSkew15F1(double std, double fParams[]);
	double getSkew15F2(double std, double fParams[]);
	double getSkew15F3(double std, double fParams[]);
	double getSkew15F1D(double std, double fParams[]);
	double getSkew15F2D(double std, double fParams[]);
	double getSkew15F3D(double std, double fParams[]);
	double getSkew15F1DD(double std, double fParams[]);
	double getSkew15F2DD(double std, double fParams[]);
	double getSkew15F3DD(double std, double fParams[]);
	double getSkew15(double std, double fParams[]);
	double getSkew15D(double std, double fParams[]);
	double getSkew15DD(double std, double fParams[]);
	double* ReCalculateCoeffs15(double* fParams, double* coeffs);  
	double CalculateVolatility15(double std, double* fParams, double* fCoeffs, bool bForceRecalcCoeffs);
};

#endif // !defined(AFX_Optimizer_H__7F089FD1_9DDA_4986_B240_E2031EDE6084__INCLUDED_)
