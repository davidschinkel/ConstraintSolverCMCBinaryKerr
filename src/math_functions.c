#include "main.h"

ftype csch(ftype x)
{
	return 1/sinh(x);
}

ftype csc(ftype x)
{
	return 1/sin(x);
}

ftype sech(ftype x)
{
	return 1/cosh(x);
}

ftype cot(ftype x)
{
	return 1/tan(x);
}

ftype ArcTan(ftype e1, ftype e2)
{//ArcTan-Funktion von Mathematica, dort sind die Argumente vertauscht
	return atan2(e2,e1);
}

ftype maximum(ftype e1, ftype e2)
{
	return (e1>e2)? e1 : e2;
}

ftype minimum(ftype e1, ftype e2)
{
	return (e1<e2)? e1 : e2;
}

int sgn(ftype x)
{
	{
	    if (x==0)
	        return 0;
	    else
	        return (x>0) ? 1 : -1;
	}
}
