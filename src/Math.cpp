#include "Math.hpp"
#include <limits.h>

namespace Math
{
#define INTEGER_ABS(T)\
	template<>\
	T Abs(const T &v)\
	{\
		return v & ~(1 << (sizeof(T) - 1));\
	}

	INTEGER_ABS(long long)
	INTEGER_ABS(long)
	INTEGER_ABS(int)
	INTEGER_ABS(short)
	INTEGER_ABS(signed char)
#undef INTEGER_ABS

	template<>
	int EuclidianRemainder(const int &a, const int &b)
	{
#if (-1 >> 1) == -1
		#define OFFSET ((a % b >> (sizeof(int) * CHAR_BIT - 1)) & b)
#else
		#define OFFSET ((a % b < 0) * b)
#endif //ARITHMETIC SHIFT CHECK

		return a % b + OFFSET;
#undef OFFSET
	}
}

