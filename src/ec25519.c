/*
  Copyright (c) 2012-2015, Matthias Schiffer <mschiffer@universe-factory.net>
  Partly based on public domain code by Matthew Dempsky and D. J. Bernstein.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/** \file
 * EC group operations for Twisted Edwards Curve \f$ ax^2 + y^2 = 1 + dx^2y^2 \f$ with
 *    \f$ a = 486664 \f$ and
 *    \f$ d = 486660 \f$
 * on prime field \f$ p = 2^{255} - 19 \f$.
 *
 * The curve is equivalent to the Montgomery Curve used in D. J. Bernstein's
 * Curve25519 Diffie-Hellman algorithm.
 *
 * See http://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html for add and
 * double operations.
 *
 * Doxygen comments for public APIs can be found in the public header file.
 *
 * Invariant that must be held by all public API: the components of an
 * \ref ecc_25519_work_t are always in the range \f$ [0, 2p) \f$.
 * Integers in this range will be called \em squeezed in the following.
 */

#include <libuecc/ecc.h>


const ecc_25519_work_t ecc_25519_work_identity = {{0}, {1}, {1}, {0}};

const ecc_25519_work_t ecc_25519_work_default_base = {
	{0xd4, 0x6b, 0xfe, 0x7f, 0x39, 0xfa, 0x8c, 0x22,
	 0xe1, 0x96, 0x23, 0xeb, 0x26, 0xb7, 0x8e, 0x6a,
	 0x34, 0x74, 0x8b, 0x66, 0xd6, 0xa3, 0x26, 0xdd,
	 0x19, 0x5e, 0x9f, 0x21, 0x50, 0x43, 0x7c, 0x54},
	{0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
	 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
	 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
	 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66},
	{1},
	{0x47, 0x56, 0x98, 0x99, 0xc7, 0x61, 0x0a, 0x82,
	 0x1a, 0xdf, 0x82, 0x22, 0x1f, 0x2c, 0x72, 0x88,
	 0xc3, 0x29, 0x09, 0x52, 0x78, 0xe9, 0x1e, 0xe4,
	 0x47, 0x4b, 0x4c, 0x81, 0xa6, 0x02, 0xfd, 0x29}
};


static const uint32_t zero[32] = {0};
static const uint32_t one[32] = {1};


/** Adds two unpacked integers (modulo p) */
static void add(uint32_t out[32], const uint32_t a[32], const uint32_t b[32]) {
	unsigned int j;
	uint32_t u;

	u = 0;

	for (j = 0; j < 31; j++) {
		u += a[j] + b[j];
		out[j] = u & 255;
		u >>= 8;
	}

	u += a[31] + b[31];
	out[31] = u;
}

/**
 * Subtracts two unpacked integers (modulo p)
 *
 * b must be \em squeezed.
 */
static void sub(uint32_t out[32], const uint32_t a[32], const uint32_t b[32]) {
	unsigned int j;
	uint32_t u;

	u = 218;

	for (j = 0;j < 31;++j) {
		u += a[j] + UINT32_C(65280) - b[j];
		out[j] = u & 255;
		u >>= 8;
	}

	u += a[31] - b[31];
	out[31] = u;
}

/**
 * Performs carry and reduce on an unpacked integer
 *
 * The result is not always fully reduced, but it will be significantly smaller than \f$ 2p \f$.
 */
static void squeeze(uint32_t a[32]) {
	unsigned int j;
	uint32_t u;

	u = 0;

	for (j = 0;j < 31;++j) {
		u += a[j];
		a[j] = u & 255;
		u >>= 8;
	}

	u += a[31];
	a[31] = u & 127;
	u = 19 * (u >> 7);

	for (j = 0;j < 31;++j) {
		u += a[j];
		a[j] = u & 255;
		u >>= 8;
	}

	u += a[31];
	a[31] = u;
}

/**
 * Ensures that the output of a previous \ref squeeze is fully reduced
 *
 * After a \ref freeze, only the lower byte of each integer part holds a meaningful value.
 */
static void freeze(uint32_t a[32]) {
	static const uint32_t minusp[32] = {
		19, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 128
	};

	uint32_t aorig[32];
	unsigned int j;
	uint32_t negative;

	for (j = 0; j < 32; j++)
		aorig[j] = a[j];
	add(a, a, minusp);
	negative = -((a[31] >> 7) & 1);

	for (j = 0; j < 32; j++)
		a[j] ^= negative & (aorig[j] ^ a[j]);
}

/**
 * Multiplies two unpacked integers (modulo p)
 *
 * The result will be \em squeezed.
 */
static void mult(uint32_t out[32], const uint32_t a[32], const uint32_t b[32]) {
	unsigned int i, j;
	uint32_t u;

	for (i = 0; i < 32; ++i) {
		u = 0;

		for (j = 0; j <= i; j++)
			u += a[j] * b[i - j];

		for (j = i + 1; j < 32; j++)
			u += 38 * a[j] * b[i + 32 - j];

		out[i] = u;
	}

	squeeze(out);
}

/**
 * Multiplies an unpacked integer with a small integer (modulo p)
 *
 * The result will be \em squeezed.
 */
static void mult_int(uint32_t out[32], uint32_t n, const uint32_t a[32]) {
	unsigned int j;
	uint32_t u;

	u = 0;

	for (j = 0; j < 31; j++) {
		u += n * a[j];
		out[j] = u & 255;
		u >>= 8;
	}

	u += n * a[31]; out[31] = u & 127;
	u = 19 * (u >> 7);

	for (j = 0; j < 31; j++) {
		u += out[j];
		out[j] = u & 255;
		u >>= 8;
	}

	u += out[j];
	out[j] = u;
}

/**
 * Squares an unpacked integer
 *
 * The result will be sqeezed.
 */
static void square(uint32_t out[32], const uint32_t a[32]) {
	unsigned int i, j;
	uint32_t u;

	for (i = 0; i < 32; i++) {
		u = 0;

		for (j = 0; j < i - j; j++)
			u += a[j] * a[i - j];

		for (j = i + 1; j < i + 32 - j; j++)
			u += 38 * a[j] * a[i + 32 - j];

		u *= 2;

		if ((i & 1) == 0) {
			u += a[i / 2] * a[i / 2];
			u += 38 * a[i / 2 + 16] * a[i / 2 + 16];
		}

		out[i] = u;
	}

	squeeze(out);
}

/** Checks for the equality of two unpacked integers */
static int check_equal(const uint32_t x[32], const uint32_t y[32]) {
	uint32_t differentbits = 0;
	int i;

	for (i = 0; i < 32; i++) {
		differentbits |= ((x[i] ^ y[i]) & 0xffff);
		differentbits |= ((x[i] ^ y[i]) >> 16);
	}

	return (1 & ((differentbits - 1) >> 16));
}

/**
 * Checks if an unpacked integer equals zero (modulo p)
 *
 * The interger must be squeezed before.
 */
static int check_zero(const uint32_t x[32]) {
	static const uint32_t p[32] = {
		0xed, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
	};

	return (check_equal(x, zero) | check_equal(x, p));
}

/** Copies r to out when b == 0, s when b == 1 */
static void selectw(ecc_25519_work_t *out, const ecc_25519_work_t *r, const ecc_25519_work_t *s, uint32_t b) {
	unsigned int j;
	uint32_t t;
	uint32_t bminus1;

	bminus1 = b - 1;
	for (j = 0; j < 32; ++j) {
		t = bminus1 & (r->X[j] ^ s->X[j]);
		out->X[j] = s->X[j] ^ t;

		t = bminus1 & (r->Y[j] ^ s->Y[j]);
		out->Y[j] = s->Y[j] ^ t;

		t = bminus1 & (r->Z[j] ^ s->Z[j]);
		out->Z[j] = s->Z[j] ^ t;

		t = bminus1 & (r->T[j] ^ s->T[j]);
		out->T[j] = s->T[j] ^ t;
	}
}

/** Copies r to out when b == 0, s when b == 1 */
static void select(uint32_t out[32], const uint32_t r[32], const uint32_t s[32], uint32_t b) {
	unsigned int j;
	uint32_t t;
	uint32_t bminus1;

	bminus1 = b - 1;
	for (j = 0;j < 32;++j) {
		t = bminus1 & (r[j] ^ s[j]);
		out[j] = s[j] ^ t;
	}
}

/**
 * Computes the square root of an unpacked integer (in the prime field modulo p)
 *
 * If the given integer has no square root, 0 is returned, 1 otherwise.
 */
static int square_root(uint32_t out[32], const uint32_t z[32]) {
	static const uint32_t minus1[32] = {
		0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
	};

	static const uint32_t rho_s[32] = {
		0xb0, 0xa0, 0x0e, 0x4a, 0x27, 0x1b, 0xee, 0xc4,
		0x78, 0xe4, 0x2f, 0xad, 0x06, 0x18, 0x43, 0x2f,
		0xa7, 0xd7, 0xfb, 0x3d, 0x99, 0x00, 0x4d, 0x2b,
		0x0b, 0xdf, 0xc1, 0x4f, 0x80, 0x24, 0x83, 0x2b
	};

	/* raise z to power (2^252-2), check if power (2^253-5) equals -1 */

	uint32_t z2[32];
	uint32_t z9[32];
	uint32_t z11[32];
	uint32_t z2_5_0[32];
	uint32_t z2_10_0[32];
	uint32_t z2_20_0[32];
	uint32_t z2_50_0[32];
	uint32_t z2_100_0[32];
	uint32_t t0[32];
	uint32_t t1[32];
	uint32_t z2_252_1[32];
	uint32_t z2_252_1_rho_s[32];
	int i;

	/* 2 */ square(z2, z);
	/* 4 */ square(t1, z2);
	/* 8 */ square(t0, t1);
	/* 9 */ mult(z9, t0, z);
	/* 11 */ mult(z11, z9, z2);
	/* 22 */ square(t0, z11);
	/* 2^5 - 2^0 = 31 */ mult(z2_5_0, t0, z9);

	/* 2^6 - 2^1 */ square(t0, z2_5_0);
	/* 2^7 - 2^2 */ square(t1, t0);
	/* 2^8 - 2^3 */ square(t0, t1);
	/* 2^9 - 2^4 */ square(t1, t0);
	/* 2^10 - 2^5 */ square(t0, t1);
	/* 2^10 - 2^0 */ mult(z2_10_0, t0, z2_5_0);

	/* 2^11 - 2^1 */ square(t0, z2_10_0);
	/* 2^12 - 2^2 */ square(t1, t0);
	/* 2^20 - 2^10 */ for (i = 2; i < 10; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^20 - 2^0 */ mult(z2_20_0, t1, z2_10_0);

	/* 2^21 - 2^1 */ square(t0, z2_20_0);
	/* 2^22 - 2^2 */ square(t1, t0);
	/* 2^40 - 2^20 */ for (i = 2; i < 20; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^40 - 2^0 */ mult(t0, t1, z2_20_0);

	/* 2^41 - 2^1 */ square(t1, t0);
	/* 2^42 - 2^2 */ square(t0, t1);
	/* 2^50 - 2^10 */ for (i = 2; i < 10; i += 2) { square(t1, t0); square(t0, t1); }
	/* 2^50 - 2^0 */ mult(z2_50_0, t0, z2_10_0);

	/* 2^51 - 2^1 */ square(t0, z2_50_0);
	/* 2^52 - 2^2 */ square(t1, t0);
	/* 2^100 - 2^50 */ for (i = 2; i < 50; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^100 - 2^0 */ mult(z2_100_0, t1, z2_50_0);

	/* 2^101 - 2^1 */ square(t1, z2_100_0);
	/* 2^102 - 2^2 */ square(t0, t1);
	/* 2^200 - 2^100 */ for (i = 2; i < 100; i += 2) { square(t1, t0); square(t0, t1); }
	/* 2^200 - 2^0 */ mult(t1, t0, z2_100_0);

	/* 2^201 - 2^1 */ square(t0, t1);
	/* 2^202 - 2^2 */ square(t1, t0);
	/* 2^250 - 2^50 */ for (i = 2; i < 50; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^250 - 2^0 */ mult(t0, t1, z2_50_0);

	/* 2^251 - 2^1 */ square(t1, t0);
	/* 2^252 - 2^2 */ square(t0, t1);
	/* 2^252 - 2^1 */ mult(z2_252_1, t0, z2);

	/* 2^253 - 2^3 */ square(t1, t0);
	/* 2^253 - 6 */ mult(t0, t1, z2);
	/* 2^253 - 5 */ mult(t1, t0, z);

	mult(z2_252_1_rho_s, z2_252_1, rho_s);

	select(out, z2_252_1, z2_252_1_rho_s, check_equal(t1, minus1));

	/* Check the root */
	square(t0, out);
	return check_equal(t0, z);
}

/** Computes the reciprocal of an unpacked integer (in the prime field modulo p) */
static void recip(uint32_t out[32], const uint32_t z[32]) {
	uint32_t z2[32];
	uint32_t z9[32];
	uint32_t z11[32];
	uint32_t z2_5_0[32];
	uint32_t z2_10_0[32];
	uint32_t z2_20_0[32];
	uint32_t z2_50_0[32];
	uint32_t z2_100_0[32];
	uint32_t t0[32];
	uint32_t t1[32];
	int i;

	/* 2 */ square(z2, z);
	/* 4 */ square(t1, z2);
	/* 8 */ square(t0, t1);
	/* 9 */ mult(z9, t0, z);
	/* 11 */ mult(z11, z9, z2);
	/* 22 */ square(t0, z11);
	/* 2^5 - 2^0 = 31 */ mult(z2_5_0, t0, z9);

	/* 2^6 - 2^1 */ square(t0, z2_5_0);
	/* 2^7 - 2^2 */ square(t1, t0);
	/* 2^8 - 2^3 */ square(t0, t1);
	/* 2^9 - 2^4 */ square(t1, t0);
	/* 2^10 - 2^5 */ square(t0, t1);
	/* 2^10 - 2^0 */ mult(z2_10_0, t0, z2_5_0);

	/* 2^11 - 2^1 */ square(t0, z2_10_0);
	/* 2^12 - 2^2 */ square(t1, t0);
	/* 2^20 - 2^10 */ for (i = 2; i < 10; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^20 - 2^0 */ mult(z2_20_0, t1, z2_10_0);

	/* 2^21 - 2^1 */ square(t0, z2_20_0);
	/* 2^22 - 2^2 */ square(t1, t0);
	/* 2^40 - 2^20 */ for (i = 2; i < 20; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^40 - 2^0 */ mult(t0, t1, z2_20_0);

	/* 2^41 - 2^1 */ square(t1, t0);
	/* 2^42 - 2^2 */ square(t0, t1);
	/* 2^50 - 2^10 */ for (i = 2; i < 10; i += 2) { square(t1, t0); square(t0, t1); }
	/* 2^50 - 2^0 */ mult(z2_50_0, t0, z2_10_0);

	/* 2^51 - 2^1 */ square(t0, z2_50_0);
	/* 2^52 - 2^2 */ square(t1, t0);
	/* 2^100 - 2^50 */ for (i = 2; i < 50; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^100 - 2^0 */ mult(z2_100_0, t1, z2_50_0);

	/* 2^101 - 2^1 */ square(t1, z2_100_0);
	/* 2^102 - 2^2 */ square(t0, t1);
	/* 2^200 - 2^100 */ for (i = 2; i < 100; i += 2) { square(t1, t0); square(t0, t1); }
	/* 2^200 - 2^0 */ mult(t1, t0, z2_100_0);

	/* 2^201 - 2^1 */ square(t0, t1);
	/* 2^202 - 2^2 */ square(t1, t0);
	/* 2^250 - 2^50 */ for (i = 2; i < 50; i += 2) { square(t0, t1); square(t1, t0); }
	/* 2^250 - 2^0 */ mult(t0, t1, z2_50_0);

	/* 2^251 - 2^1 */ square(t1, t0);
	/* 2^252 - 2^2 */ square(t0, t1);
	/* 2^253 - 2^3 */ square(t1, t0);
	/* 2^254 - 2^4 */ square(t0, t1);
	/* 2^255 - 2^5 */ square(t1, t0);
	/* 2^255 - 21 */ mult(out, t1, z11);
}

int ecc_25519_load_xy(ecc_25519_work_t *out, const ecc_int256_t *x, const ecc_int256_t *y) {
	int i;
	uint32_t X2[32], Y2[32], aX2[32], dX2[32], dX2Y2[32], aX2_Y2[32], _1_dX2Y2[32], r[32];

	for (i = 0; i < 32; i++) {
		out->X[i] = x->p[i];
		out->Y[i] = y->p[i];
		out->Z[i] = (i == 0);
	}

	/* Check validity */
	square(X2, out->X);
	square(Y2, out->Y);
	mult_int(aX2, UINT32_C(486664), X2);
	mult_int(dX2, UINT32_C(486660), X2);
	mult(dX2Y2, dX2, Y2);
	add(aX2_Y2, aX2, Y2);
	add(_1_dX2Y2, one, dX2Y2);
	sub(r, aX2_Y2, _1_dX2Y2);
	squeeze(r);

	if (!check_zero(r))
	    return 0;

	mult(out->T, out->X, out->Y);

	return 1;
}

void ecc_25519_store_xy(ecc_int256_t *x, ecc_int256_t *y, const ecc_25519_work_t *in) {
	uint32_t X[32], Y[32], Z[32];
	int i;

	recip(Z, in->Z);

	if (x) {
		mult(X, Z, in->X);
		freeze(X);
		for (i = 0; i < 32; i++)
			x->p[i] = X[i];
	}

	if (y) {
		mult(Y, Z, in->Y);
		freeze(Y);
		for (i = 0; i < 32; i++)
			y->p[i] = Y[i];
	}
}

int ecc_25519_load_packed(ecc_25519_work_t *out, const ecc_int256_t *in) {
	int i;
	uint32_t X2[32] /* X^2 */, aX2[32] /* aX^2 */, dX2[32] /* dX^2 */, _1_aX2[32] /* 1-aX^2 */, _1_dX2[32] /* 1-aX^2 */;
	uint32_t _1_1_dX2[32]  /* 1/(1-aX^2) */, Y2[32] /* Y^2 */, Y[32], Yt[32];

	for (i = 0; i < 32; i++) {
		out->X[i] = in->p[i];
		out->Z[i] = (i == 0);
	}

	out->X[31] &= 0x7f;

	square(X2, out->X);
	mult_int(aX2, UINT32_C(486664), X2);
	mult_int(dX2, UINT32_C(486660), X2);
	sub(_1_aX2, one, aX2);
	sub(_1_dX2, one, dX2);
	recip(_1_1_dX2, _1_dX2);
	mult(Y2, _1_aX2, _1_1_dX2);

	if (!square_root(Y, Y2))
		return 0;

	/* No squeeze is necessary after subtractions from zero if the subtrahend is squeezed */
	sub(Yt, zero, Y);

	select(out->Y, Y, Yt, (in->p[31] >> 7) ^ (Y[0] & 1));

	mult(out->T, out->X, out->Y);

	return 1;
}

void ecc_25519_store_packed(ecc_int256_t *out, const ecc_25519_work_t *in) {
	ecc_int256_t y;

	ecc_25519_store_xy(out, &y, in);
	out->p[31] |= (y.p[0] << 7);
}

int ecc_25519_is_identity(const ecc_25519_work_t *in) {
	uint32_t Y_Z[32];

	sub(Y_Z, in->Y, in->Z);
	squeeze(Y_Z);

	return (check_zero(in->X)&check_zero(Y_Z));
}

void ecc_25519_negate(ecc_25519_work_t *out, const ecc_25519_work_t *in) {
	int i;

	for (i = 0; i < 32; i++) {
		out->Y[i] = in->Y[i];
		out->Z[i] = in->Z[i];
	}

	/* No squeeze is necessary after subtractions from zero if the subtrahend is squeezed */
	sub(out->X, zero, in->X);
	sub(out->T, zero, in->T);
}

void ecc_25519_double(ecc_25519_work_t *out, const ecc_25519_work_t *in) {
	uint32_t A[32], B[32], C[32], D[32], E[32], F[32], G[32], H[32], t0[32], t1[32], t2[32], t3[32];

	square(A, in->X);
	square(B, in->Y);
	square(t0, in->Z);
	mult_int(C, 2, t0);
	mult_int(D, UINT32_C(486664), A);
	add(t1, in->X, in->Y);
	square(t2, t1);
	sub(t3, t2, A);
	sub(E, t3, B);
	add(G, D, B);
	sub(F, G, C);
	sub(H, D, B);
	mult(out->X, E, F);
	mult(out->Y, G, H);
	mult(out->T, E, H);
	mult(out->Z, F, G);
}

void ecc_25519_add(ecc_25519_work_t *out, const ecc_25519_work_t *in1, const ecc_25519_work_t *in2) {
	uint32_t A[32], B[32], C[32], D[32], E[32], F[32], G[32], H[32], t0[32], t1[32], t2[32], t3[32], t4[32], t5[32];

	mult(A, in1->X, in2->X);
	mult(B, in1->Y, in2->Y);
	mult_int(t0, UINT32_C(486660), in2->T);
	mult(C, in1->T, t0);
	mult(D, in1->Z, in2->Z);
	add(t1, in1->X, in1->Y);
	add(t2, in2->X, in2->Y);
	mult(t3, t1, t2);
	sub(t4, t3, A);
	sub(E, t4, B);
	sub(F, D, C);
	add(G, D, C);
	mult_int(t5, UINT32_C(486664), A);
	sub(H, B, t5);
	mult(out->X, E, F);
	mult(out->Y, G, H);
	mult(out->T, E, H);
	mult(out->Z, F, G);
}

void ecc_25519_sub(ecc_25519_work_t *out, const ecc_25519_work_t *in1, const ecc_25519_work_t *in2) {
	ecc_25519_work_t in2_neg;

	ecc_25519_negate(&in2_neg, in2);
	ecc_25519_add(out, in1, &in2_neg);
}

void ecc_25519_scalarmult_bits(ecc_25519_work_t *out, const ecc_int256_t *n, const ecc_25519_work_t *base, unsigned bits) {
	ecc_25519_work_t Q2, Q2p;
	ecc_25519_work_t cur = ecc_25519_work_identity;
	int b, pos;

	if (bits > 256)
		bits = 256;

	for (pos = bits - 1; pos >= 0; --pos) {
		b = n->p[pos / 8] >> (pos & 7);
		b &= 1;

		ecc_25519_double(&Q2, &cur);
		ecc_25519_add(&Q2p, &Q2, base);
		selectw(&cur, &Q2, &Q2p, b);
	}

	*out = cur;
}

void ecc_25519_scalarmult(ecc_25519_work_t *out, const ecc_int256_t *n, const ecc_25519_work_t *base) {
	ecc_25519_scalarmult_bits(out, n, base, 256);
}

void ecc_25519_scalarmult_base_bits(ecc_25519_work_t *out, const ecc_int256_t *n, unsigned bits) {
	ecc_25519_scalarmult_bits(out, n, &ecc_25519_work_default_base, bits);
}

void ecc_25519_scalarmult_base(ecc_25519_work_t *out, const ecc_int256_t *n) {
	ecc_25519_scalarmult(out, n, &ecc_25519_work_default_base);
}
