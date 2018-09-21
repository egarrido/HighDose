typedef uint64_t u64;
 
#define QSIZE 0xfff
#define CNG (cng = 6906969069ULL * cng + 13579)
#define XS (xs ^= (xs << 13), xs ^= (xs >> 17), xs ^= (xs << 43))
#define KISS (B64MWC() + CNG + XS)
 
static u64 QARY[QSIZE];
static int j;
static u64 carry;
static u64 xs;
static u64 cng;
 
u64 B64MWC(void)
{
	u64 t, x;
	j = (j + 1) & (QSIZE - 1);
	x = QARY[j];
	t = (x << 28) + carry;
	carry = (x >> 36) - (t < x);
	return (QARY[j] = t - x);
}
 
/* Initialize PRNG with default seed */
void randk_seed(u64 time)
{
	u64 i;
	j = QSIZE - 1;
	carry = 0;
	xs = 362436069362436069ULL;
	cng = 123456789987654321ULL;
	cng = time;
	/* Seed QARY[] with CNG+XS: */
	for (i = 0; i < QSIZE; i++)
		QARY[i] = CNG + XS;
}
 
/* Generate a pseudorandom 64-bit unsigned integer. */
u64 randk(void)
{
	return KISS;
}
