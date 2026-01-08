#include "sample.h"
#include "fips202.h"
#include "aes.h"
#include <string.h>
#include "config.h"
// This code is based on the implementation of FrodoKEM
void scloudplus_mul_add_as_e(const uint8_t *seedA, const uint16_t *S,
							 const uint16_t *E, uint16_t *B)
{
	memcpy(B, E, 2 * scloudplus_m * scloudplus_nbar);
	ALIGN_HEADER(32)
	uint16_t AROWOUT[4 * scloudplus_n] ALIGN_FOOTER(32) = {0};
	ALIGN_HEADER(32)
	uint32_t AROWIN[4 * scloudplus_block_rowlen] ALIGN_FOOTER(32) = {0};
	uint8_t aes_key_schedule[16 * 11];
	AES128_load_schedule(seedA, aes_key_schedule);
	for (int i = 0; i < scloudplus_m; i += 4)
	{

		for (int j = 0; j < scloudplus_block_number; j += 1)
		{
			AROWIN[scloudplus_block_size * j + 0 * scloudplus_block_rowlen] =
				i * scloudplus_block_number + j;
			AROWIN[scloudplus_block_size * j + 1 * scloudplus_block_rowlen] =
				(i + 1) * scloudplus_block_number + j;
			AROWIN[scloudplus_block_size * j + 2 * scloudplus_block_rowlen] =
				(i + 2) * scloudplus_block_number + j;
			AROWIN[scloudplus_block_size * j + 3 * scloudplus_block_rowlen] =
				(i + 3) * scloudplus_block_number + j;
		}
		AES128_CTR_enc_sch((uint8_t *)AROWIN, 4 * scloudplus_n * sizeof(uint16_t),
						   aes_key_schedule, (uint8_t *)AROWOUT);

		for (int k = 0; k < scloudplus_nbar; k++)
		{
			uint16_t sum[4] = {0};
			for (int j = 0; j < scloudplus_n; j++)
			{
				uint16_t sp = S[k * scloudplus_n + j];
				sum[0] += AROWOUT[0 * scloudplus_n + j] * sp;
				sum[1] += AROWOUT[1 * scloudplus_n + j] * sp;
				sum[2] += AROWOUT[2 * scloudplus_n + j] * sp;
				sum[3] += AROWOUT[3 * scloudplus_n + j] * sp;
			}
			B[(i + 0) * scloudplus_nbar + k] += sum[0];
			B[(i + 1) * scloudplus_nbar + k] += sum[1];
			B[(i + 2) * scloudplus_nbar + k] += sum[2];
			B[(i + 3) * scloudplus_nbar + k] += sum[3];
		}
	}
	AES128_free_schedule(aes_key_schedule);
}
// This code is based on the implementation of FrodoKEM
void scloudplus_mul_add_sa_e(const uint8_t *seedA, const uint16_t *S,
							 uint16_t *E, uint16_t *C)
{

	ALIGN_HEADER(32)
	uint16_t AROWOUT[8 * scloudplus_n] ALIGN_FOOTER(32) = {0};

	uint8_t aes_key_schedule[16 * 11];
	AES128_load_schedule(seedA, aes_key_schedule);

	ALIGN_HEADER(32)
	uint32_t AROWIN[8 * scloudplus_block_rowlen] ALIGN_FOOTER(32) = {0};

	for (int i = 0; i < scloudplus_m; i += 8)
	{
		for (int q = 0; q < 8; q++)
		{
			for (int p = 0; p < scloudplus_block_number; p += 1)
			{
				AROWIN[q * scloudplus_block_rowlen + scloudplus_block_size * p] =
					(i + q) * scloudplus_block_number + p;
			}
		}
		AES128_CTR_enc_sch((uint8_t *)AROWIN, 8 * scloudplus_n * sizeof(uint16_t),
						   aes_key_schedule, (uint8_t *)AROWOUT);

		for (int j = 0; j < scloudplus_mbar; j++)
		{
			uint16_t sum = 0;
			uint16_t sp[8];
			for (int p = 0; p < 8; p++)
			{
				sp[p] = S[j * scloudplus_m + i + p];
			}
			for (int q = 0; q < scloudplus_n; q++)
			{
				sum = E[j * scloudplus_n + q];
				for (int p = 0; p < 8; p++)
				{
					sum += sp[p] * AROWOUT[p * scloudplus_n + q];
				}
				E[j * scloudplus_n + q] = sum;
			}
		}
	}
	memcpy((unsigned char *)C, (unsigned char *)E,
		   2 * scloudplus_mbar * scloudplus_n);
	AES128_free_schedule(aes_key_schedule);
}

static inline uint32_t read3bytestou32(const uint8_t *ptr)
{
	return ((uint32_t)ptr[0] << 0) | ((uint32_t)ptr[1] << 8) |
		   ((uint32_t)ptr[2] << 16);
}
static inline uint32_t read4bytestou32(const uint8_t *ptr)
{
	return ((uint32_t)ptr[0]) | ((uint32_t)ptr[1] << 8) |
		   ((uint32_t)ptr[2] << 16) | ((uint32_t)ptr[3] << 24);
}

static inline uint64_t read7bytestou64(const uint8_t *ptr)
{
	return ((uint64_t)ptr[0] << 0) | ((uint64_t)ptr[1] << 8) |
		   ((uint64_t)ptr[2] << 16) | ((uint64_t)ptr[3] << 24) |
		   ((uint64_t)ptr[4] << 32) | ((uint64_t)ptr[5] << 40) |
		   ((uint64_t)ptr[6] << 48);
}

static inline void cbd1(uint8_t in, uint16_t *out)
{
	uint8_t b, b0, b1;
	b = in;
	for (size_t j = 0; j < 4; j++)
	{
		b0 = b & 1;
		b1 = (b >> 1) & 1;
		*(out + j) = (uint16_t)(b0 - b1);
		b = b >> 2;
	}
}

static inline void cbd2(uint8_t in, uint16_t *out)
{
	uint8_t b = 0;
	b += in & 0x55;
	b += (in >> 1) & 0x55;
	*out = (uint16_t)((b & 0x03) - ((b >> 2) & 0x03));
	*(out + 1) = (uint16_t)(((b >> 4) & 0x03) - ((b >> 6) & 0x03));
}

static inline void cbd3(uint32_t in, uint16_t *out)
{
	uint32_t b = 0;
	b += in & 0x00249249;
	b += (in >> 1) & 0x00249249;
	b += (in >> 2) & 0x00249249;
	for (int i = 0; i < 4; i++)
	{
		out[i] = ((b >> (6 * i)) & 0x07) - ((b >> (6 * i + 3)) & 0x07);
	}
}

static inline void cbd7(uint64_t in, uint16_t *out)
{
	uint64_t b0 = 0;
	b0 += in & 0x2040810204081;
	b0 += (in >> 1) & 0x2040810204081;
	b0 += (in >> 2) & 0x2040810204081;
	b0 += (in >> 3) & 0x2040810204081;
	b0 += (in >> 4) & 0x2040810204081;
	b0 += (in >> 5) & 0x2040810204081;
	b0 += (in >> 6) & 0x2040810204081;
	for (int i = 0; i < 4; i++)
	{
		out[i] = ((b0 >> (14 * i)) & 0x7F) - ((b0 >> (14 * i + 7)) & 0x7F);
	}
}

void scloudplus_sampleeta1(uint8_t *seed, uint16_t *matrixe)
{
	size_t hashlen =
		((scloudplus_m * scloudplus_nbar) * (2 * scloudplus_eta1)) >> 3;
	uint8_t *tmp = (uint8_t *)malloc(hashlen * sizeof(uint8_t));
	uint8_t *ptrtmp = tmp;
	uint16_t *ptrmatrix = matrixe;
	memset(matrixe, 0, scloudplus_m * scloudplus_nbar * sizeof(uint16_t));
	shake256(tmp, hashlen, seed, 32);
#if (scloudplus_eta1 == 2)
	for (size_t i = 0; i < scloudplus_m * scloudplus_nbar; i = i + 2)
	{
		cbd2(*ptrtmp, ptrmatrix);
		ptrtmp = ptrtmp + 1;
		ptrmatrix = ptrmatrix + 2;
	}
#elif (scloudplus_eta1 == 3)
	for (size_t i = 0; i < scloudplus_m * scloudplus_nbar; i = i + 4)
	{
		cbd3(read3bytestou32(ptrtmp) & 0xFFFFFF, ptrmatrix);
		ptrtmp = ptrtmp + 3;
		ptrmatrix = ptrmatrix + 4;
	}
#elif (scloudplus_eta1 == 7)
	for (size_t i = 0; i < scloudplus_m * scloudplus_nbar; i = i + 4)
	{
		cbd7(read7bytestou64(ptrtmp) & 0xFFFFFFFFFFFFFF, ptrmatrix);
		ptrtmp = ptrtmp + 7;
		ptrmatrix = ptrmatrix + 4;
	}
#endif
	free(tmp);
}

void scloudplus_sampleeta2(uint8_t *seed, uint16_t *matrixe1,
						   uint16_t *matrixe2)
{
	size_t hash1len =
		((scloudplus_mbar * scloudplus_n) * (2 * scloudplus_eta2)) >> 3;
	size_t hash2len =
		((scloudplus_mbar * scloudplus_nbar) * (2 * scloudplus_eta2) + 7) >> 3;
	uint8_t *tmp = (uint8_t *)malloc((hash1len + hash2len) * sizeof(uint8_t));
	uint8_t *ptrtmp1 = tmp;
	uint8_t *ptrtmp2 = tmp + hash1len;
	uint16_t *ptrmatrix1 = matrixe1;
	uint16_t *ptrmatrix2 = matrixe2;
	memset(matrixe1, 0, scloudplus_mbar * scloudplus_n * 2);
	memset(matrixe2, 0, scloudplus_mbar * scloudplus_nbar * 2);
	shake256(tmp, hash1len + hash2len, seed, 32);
#if (scloudplus_eta2 == 1)
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_n; i = i + 4)
	{
		cbd1(*ptrtmp1, ptrmatrix1);
		ptrtmp1 = ptrtmp1 + 1;
		ptrmatrix1 = ptrmatrix1 + 4;
	}
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_nbar; i = i + 4)
	{
		cbd1(*ptrtmp2, ptrmatrix2);
		ptrtmp2 = ptrtmp2 + 1;
		ptrmatrix2 = ptrmatrix2 + 4;
	}
#elif (scloudplus_eta2 == 2)
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_n; i = i + 2)
	{
		cbd2(*ptrtmp1, ptrmatrix1);
		ptrtmp1 = ptrtmp1 + 1;
		ptrmatrix1 = ptrmatrix1 + 2;
	}
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_nbar; i = i + 2)
	{
		cbd2(*ptrtmp2, ptrmatrix2);
		ptrtmp2 = ptrtmp2 + 1;
		ptrmatrix2 = ptrmatrix2 + 2;
	}
#elif (scloudplus_eta2 == 7)
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_n; i = i + 4)
	{
		cbd7(read7bytestou64(ptrtmp1) & 0xFFFFFFFFFFFFFF, ptrmatrix1);
		ptrtmp1 = ptrtmp1 + 7;
		ptrmatrix1 = ptrmatrix1 + 4;
	}
	for (size_t i = 0; i < scloudplus_mbar * scloudplus_nbar; i = i + 4)
	{
		cbd7(read7bytestou64(ptrtmp2) & 0xFFFFFFFFFFFFFF, ptrmatrix2);
		ptrtmp2 = ptrtmp2 + 7;
		ptrmatrix2 = ptrmatrix2 + 4;
	}
#endif
	free(tmp);
}


static inline uint32_t div600_u32(uint32_t x) {
  return (uint64_t)x * 0x1b4e81b5 >> 38;
}

static inline uint64_t div1120_u64(uint64_t x)
{
  uint64_t q =
    ( (__uint128_t)x * 0xea0ea0ea0ea0ea0full ) >> 74;
  return q;
}

static inline uint64_t div1136_u64(uint64_t x)
{
  uint64_t q = ( (__uint128_t)(x>>4) * 0xe6c2b4481cd8569ull ) >> 66;
  return q;
}

void readu8ton(uint8_t *in, int n, uint16_t *out, int *outlen)
{
	uint8_t *ptrin = in;
	uint16_t *ptrout = out;
	*outlen = 0;
#if (scloudplus_l == 128)
	uint32_t tmp, q, r;
	for (int i = 0; i < n; i = i + 7)
	{
		tmp = read4bytestou32(ptrin) & 0xFFFFFFF;
		if (tmp < scloudplus_n3)
		{
			// scloudplus_n==600
			//*ptrout = tmp % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*ptrout = r;
			tmp = q;
			//*(ptrout + 1) = tmp / scloudplus_n % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*(ptrout + 1) = r;
			//*(ptrout + 2) = tmp / scloudplus_n2 % scloudplus_n;
			*(ptrout + 2) = q;
			ptrout = ptrout + 3;
			*outlen += 3;
		}
		tmp = (read4bytestou32(ptrin + 3) >> 4) & 0xFFFFFFF;
		if (tmp < scloudplus_n3)
		{
			// scloudplus_n==600
			//*ptrout = tmp % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*ptrout = r;
			tmp = q;
			//*(ptrout + 1) = tmp / scloudplus_n % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*(ptrout + 1) = r;
			//*(ptrout + 2) = tmp / scloudplus_n2 % scloudplus_n;
			*(ptrout + 2) = q;
			ptrout = ptrout + 3;
			*outlen += 3;
		}
		ptrin = ptrin + 7;
	}
#elif (scloudplus_l == 192)
	uint64_t t64, x;
	for (int i = 0; i < n/15; i++)
	{
		t64 = *(uint64_t *)ptrin & 0x0FFFFFFFFFFFFFFFull;
		for (int j = 0; j < 6; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_n) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
		t64 = *(uint64_t *)(ptrin + 7) >> 4;
		for (int j = 0; j < 6; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_n) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
		ptrin += 15;
	}
	if (n%15 != 0) {
		t64 = *(uint64_t *)(ptrin - (8-n%15)) >> 8*(8-n%15);
		for (int j = 0; j < 4; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_n) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
	}
#elif (scloudplus_l == 256)
	uint64_t A[8], q, r;
	for (int i = 0; i < scloudplus_rejblocks*136/51; i++)
	{
		A[0] = *(uint64_t *)ptrin & 0x7FFFFFFFFFFFF;
		A[1] = (*(uint64_t *)(ptrin + 6) >> 3) & 0x7FFFFFFFFFFFF;
		A[2] = (*(uint64_t *)(ptrin + 12) >> 6) & 0x7FFFFFFFFFFFF;
		A[3] = (*(uint64_t *)(ptrin + 19) >> 1) & 0x7FFFFFFFFFFFF;
		A[4] = (*(uint64_t *)(ptrin + 25) >> 4) & 0x7FFFFFFFFFFFF;
		A[5] = (*(uint64_t *)(ptrin + 31) >> 7) & 0x7FFFFFFFFFFFF;
		A[6] = (*(uint64_t *)(ptrin + 38) >> 2) & 0x7FFFFFFFFFFFF;
		A[7] = (*(uint64_t *)(ptrin + 44) >> 5) & 0x7FFFFFFFFFFFF;
		for (int j = 0; j < 8; j++)
		{
			if (A[j] < scloudplus_n5)
			{
				// scloudplus_n==1120
				//*ptrout = A[j] % scloudplus_n;
				q = div1120_u64(A[j]);
				r = A[j] - q*1120;
				*ptrout = r;
				A[j] = q;
				//*(ptrout + 1) = A[j] / scloudplus_n % scloudplus_n;
				q = div1120_u64(A[j]);
				r = A[j] - q*1120;
				*(ptrout + 1) = r;
				A[j] = q;
				//*(ptrout + 2) = A[j] / scloudplus_n2 % scloudplus_n;
				q = div1120_u64(A[j]);
				r = A[j] - q*1120;
				*(ptrout + 2) = r;
				A[j] = q;
				//*(ptrout + 3) = A[j] / scloudplus_n3 % scloudplus_n;
				q = div1120_u64(A[j]);
				r = A[j] - q*1120;
				*(ptrout + 3) = r;
				//*(ptrout + 4) = A[j] / scloudplus_n4 % scloudplus_n;
				*(ptrout + 4) = q;
				ptrout = ptrout + 5;
				*outlen += 5;
			}
		}
		ptrin = ptrin + 51;
	}
#endif
}

void readu8tom(uint8_t *in, int n, uint16_t *out, int *outlen)
{
	uint8_t *ptrin = in;
	uint16_t *ptrout = out;
	*outlen = 0;
#if (scloudplus_l == 128)
	uint32_t tmp, q, r;
	for (int i = 0; i < n; i = i + 7)
	{
		tmp = read4bytestou32(ptrin) & 0xFFFFFFF;
		if (tmp < scloudplus_m3)
		{
			// scloudplus_m==600
			//*ptrout = tmp % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*ptrout = r;
			tmp = q;
			//*(ptrout + 1) = tmp / scloudplus_n % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*(ptrout + 1) = r;
			//*(ptrout + 2) = tmp / scloudplus_n2 % scloudplus_n;
			*(ptrout + 2) = q;
			ptrout = ptrout + 3;
			*outlen += 3;
		}
		tmp = (read4bytestou32(ptrin + 3) >> 4) & 0xFFFFFFF;
		if (tmp < scloudplus_m3)
		{
			// scloudplus_m==600
			//*ptrout = tmp % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*ptrout = r;
			tmp = q;
			//*(ptrout + 1) = tmp / scloudplus_n % scloudplus_n;
			q = div600_u32(tmp);
			r = tmp - q * 600;
			*(ptrout + 1) = r;
			//*(ptrout + 2) = tmp / scloudplus_n2 % scloudplus_n;
			*(ptrout + 2) = q;
			ptrout = ptrout + 3;
			*outlen += 3;
		}
		ptrin = ptrin + 7;
	}
#elif (scloudplus_l == 192)
	uint64_t t64, x;
	for (int i = 0; i < n/15; i++)
	{
		t64 = *(uint64_t *)ptrin & 0x0FFFFFFFFFFFFFFF;
		for (int j = 0; j < 6; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_m) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
		t64 = *(uint64_t *)(ptrin + 7) >> 4;
		for (int j = 0; j < 6; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_m) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
		ptrin += 15;
	}
	if (n%15 != 0) {
		t64 = *(uint64_t *)(ptrin - (8-n%15)) >> 8*(8-n%15);
		for (int j = 0; j < 4; j++) {
			x = t64 & 0x3FF;
			t64 >>= 10;
			if (x < scloudplus_m) {
				*ptrout = (uint16_t) x;
				ptrout += 1;
				*outlen += 1;
			}
		}
	}
#elif (scloudplus_l == 256)
	uint64_t A[8], q, r;
	for (int i = 0; i < scloudplus_rejblocks*136/51; i++)
	{
		A[0] = *(uint64_t *)ptrin & 0x7FFFFFFFFFFFF;
		A[1] = (*(uint64_t *)(ptrin + 6) >> 3) & 0x7FFFFFFFFFFFF;
		A[2] = (*(uint64_t *)(ptrin + 12) >> 6) & 0x7FFFFFFFFFFFF;
		A[3] = (*(uint64_t *)(ptrin + 19) >> 1) & 0x7FFFFFFFFFFFF;
		A[4] = (*(uint64_t *)(ptrin + 25) >> 4) & 0x7FFFFFFFFFFFF;
		A[5] = (*(uint64_t *)(ptrin + 31) >> 7) & 0x7FFFFFFFFFFFF;
		A[6] = (*(uint64_t *)(ptrin + 38) >> 2) & 0x7FFFFFFFFFFFF;
		A[7] = (*(uint64_t *)(ptrin + 44) >> 5) & 0x7FFFFFFFFFFFF;
		for (int j = 0; j < 8; j++)
		{
			if (A[j] < scloudplus_m5)
			{
				// scloudplus_m==1136
				//*ptrout = A[j] % scloudplus_m;
				q = div1136_u64(A[j]);
				r = A[j] - q*1136;
				*ptrout = r;
				A[j] = q;
				//*(ptrout + 1) = A[j] / scloudplus_n % scloudplus_m;
				q = div1136_u64(A[j]);
				r = A[j] - q*1136;
				*(ptrout + 1) = r;
				A[j] = q;
				//*(ptrout + 2) = A[j] / scloudplus_n2 % scloudplus_m;
				q = div1136_u64(A[j]);
				r = A[j] - q*1136;
				*(ptrout + 2) = r;
				A[j] = q;
				//*(ptrout + 3) = A[j] / scloudplus_n3 % scloudplus_m;
				q = div1136_u64(A[j]);
				r = A[j] - q*1136;
				*(ptrout + 3) = r;
				//*(ptrout + 4) = A[j] / scloudplus_n4 % scloudplus_m;
				*(ptrout + 4) = q;
				ptrout = ptrout + 5;
				*outlen += 5;
			}
		}
		ptrin = ptrin + 51;
	}
#endif
}

void scloudplus_samplepsi(uint8_t *seed, uint16_t *matrixs)
{
	memset(matrixs, 0, scloudplus_n * scloudplus_nbar * 2);
	uint8_t hash[scloudplus_rejblocks*136];
	uint16_t tmp[scloudplus_mnout];
	keccak_state state;
	int outlen, condition, k = 0;
	uint16_t location, mask;
	shake256_absorb_once(&state, seed, 32);
	shake256_squeezeblocks(hash, scloudplus_rejblocks, &state);
	readu8ton(hash, 136*scloudplus_rejblocks, tmp, &outlen);
	for (int i = 0; i < scloudplus_nbar; i++)
	{
		int j = 0;
		while (j < scloudplus_h1 * 2)
		{
			if (k == outlen)
			{
				shake256_squeezeblocks(hash, scloudplus_rejblocks, &state);
				readu8ton(hash, 136*scloudplus_rejblocks, tmp, &outlen);
				k = 0;
			}
			location = tmp[k];
			condition = (matrixs[i * scloudplus_n + location] == 0);
			mask = -condition;
			matrixs[i * scloudplus_n + location] =
				(matrixs[i * scloudplus_n + location] & ~mask) |
				((1 - 2 * (j & 1)) & mask);
			j += condition;
			k++;
		}
	}
}

void scloudplus_samplephi(uint8_t *seed, uint16_t *matrixs)
{
	memset(matrixs, 0, scloudplus_mbar * scloudplus_m * 2);
	uint8_t hash[scloudplus_rejblocks*136];
	uint16_t tmp[scloudplus_mnout];
	keccak_state state;
	int outlen, condition, k = 0;
	uint16_t location, mask;
	shake256_absorb_once(&state, seed, 32);
	shake256_squeezeblocks(hash, scloudplus_rejblocks, &state);
	readu8tom(hash, 136*scloudplus_rejblocks, tmp, &outlen);
	for (int i = 0; i < scloudplus_mbar; i++)
	{
		int j = 0;
		while (j < scloudplus_h2 * 2)
		{
			if (k == outlen)
			{
				shake256_squeezeblocks(hash, scloudplus_rejblocks, &state);
				readu8tom(hash, 136*scloudplus_rejblocks, tmp, &outlen);
				k = 0;
			}
			location = tmp[k];
			condition = (matrixs[i * scloudplus_m + location] == 0);
			mask = -condition;
			matrixs[i * scloudplus_m + location] =
				(matrixs[i * scloudplus_m + location] & ~mask) |
				((1 - 2 * (j & 1)) & mask);
			j += condition;
			k++;
		}
	}
}
