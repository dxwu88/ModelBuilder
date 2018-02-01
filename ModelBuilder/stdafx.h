// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include "frontmip.h"
#include "frontkey.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <utility>
#include <type_traits>
#include <typeinfo>
//#include <cxxabi.h>
using namespace std;

#define maxvar 2
#define maxcol 2
#define maxrow 100

//template <typename T>
//T **AllocateDynamicArray(int nRows, int nCols)
//{
//	T **dynamicArray;
//
//	dynamicArray = new T*[nRows];
//	for (int i = 0; i < nRows; i++)
//		dynamicArray[i] = new T[nCols];
//
//	return dynamicArray;
//}
//
//template <typename T>
//void FreeDynamicArray(T** dArray)
//{
//	delete[] * dArray;
//	delete[] dArray;
//}

// TODO: reference additional headers your program requires here
