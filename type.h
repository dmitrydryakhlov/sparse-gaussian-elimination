#ifndef __TYPE__
#define __TYPE__


//#define _DEBUG_MEMORY //define if want to use malloc wrapper with memory table
#ifdef _DEBUG_MEMORY
#include <string>
#include <stdarg.h>
#include <stdio.h>
#endif


#define SINGLE_PRECISION 0
#define DOUBLE_PRECISION 1

#define TYPE_OF_FLOAT  DOUBLE_PRECISION

#ifndef TYPE_OF_FLOAT 
#define TYPE_OF_FLOAT DOUBLE_PRECISION
#endif

#if TYPE_OF_FLOAT == DOUBLE_PRECISION
#define EPSILON 1E-15
#define FLOAT_TYPE double
#elif TYPE_OF_FLOAT == SINGLE_PRECISION
#define EPSILON 1E-8
#define FLOAT_TYPE float
#endif


#define ERROR_TYPE int

#define Int32 int
#define uInt32 unsigned int
#define Int64 long long
#define uInt64 unsigned long long

#define indexType Int32
#define CountType uInt64


#define PRINT_TYPE double

enum PermutationType
{
	Metis,
	UserSupplied
};

struct CRSMatrix
{
	indexType n;
	CountType nz;
	CountType* row;
	indexType* column;
	FLOAT_TYPE* val;
};

struct RHS
{
	indexType n;
	indexType m;
	FLOAT_TYPE* b;
};

struct SolverParameters
{
	PermutationType pType;//тип перестановки
	int* iperm;
	int threshold;
	uInt32 Relaxation;
};

#ifdef _DEBUG_MEMORY
struct memoryRecord
{
	uInt32 fileLine;
	std::string fileName;
	Int64 size;
};
#endif



#endif
