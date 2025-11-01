
#define abs128(a)		(a>0? a : -1*a)
__int128 __umodti3(__int128 a, __int128 b)
{
    __int128 sign = b > 0 ? 1 : -1;
    __int128 absa = abs128(a);
    __int128 absb = abs128(b);
    __int128 mod = 0;
    while (absa > mod*absb) {
        mod++;
    }
    return mod * sign;
}

double __floattidf(__int128 a)
{
    int sign = a > 0 ? 1 : -1;
    a = abs128(a);
    long long low = (long long)(a & 0xFFFFFFFFFFFFFFFF);
    long long high = (long long)((a >> 64) & 0xFFFFFFFFFFFFFFFF);
    double result = (double)high * 18446744073709551616.0 + (double)low;
    return result * sign;
}