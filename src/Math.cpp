#include "Math.hpp"
#include <limits.h>

namespace Math
{
#define INTEGER_ABS(T)\
	template<>\
	T Abs(const T &pVal)\
	{\
		return pVal & ~(1 << (sizeof(T) - 1));\
	}

	INTEGER_ABS(long long)
	INTEGER_ABS(long)
	INTEGER_ABS(int)
	INTEGER_ABS(short)
	INTEGER_ABS(signed char)
#undef INTEGER_ABS

	template<>
	int EuclidianRemainder(const int &pLhs, const int &pRhs)
	{
#if (-1 >> 1) == -1
		#define OFFSET ((pLhs % pRhs >> (sizeof(int) * CHAR_BIT - 1)) & pRhs)
#else
		#define OFFSET ((a % b < 0) * b)
#endif //ARITHMETIC SHIFT CHECK

		return pLhs % pRhs + OFFSET;
#undef OFFSET
	}
}

