/*
  Copyright (c) 2012, Matthias Schiffer <mschiffer@universe-factory.net>
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

/*
  EC group operations for Twisted Edwards Curve ax^2 + y^2 + 1 + dx^2y^2 with
    a = 486664
    d = 486660
  on prime field p = 2^255 - 19.

  The curve is equivalent to the Montgomery Curve used in D. J. Bernstein's
  Curve25519 Diffie-Hellman algorithm

  See http://hyperelliptic.org/EFD/g1p/auto-twisted-inverted.html for add and
  double operations
*/

#include <libuecc/ecc.h>


static void add(unsigned int out[32], const unsigned int a[32], const unsigned int b[32]) {
	unsigned int j;
	unsigned int u;
	u = 0;
	for (j = 0;j < 31;++j) { u += a[j] + b[j]; out[j] = u & 255; u >>= 8; }
	u += a[31] + b[31]; out[31] = u;
}

static void sub(unsigned int out[32], const unsigned int a[32], const unsigned int b[32]) {
	unsigned int j;
	unsigned int u;
	u = 218;
	for (j = 0;j < 31;++j) {
		u += a[j] + 65280 - b[j];
		out[j] = u & 255;
		u >>= 8;
	}
	u += a[31] - b[31];
	out[31] = u;
}

static void squeeze(unsigned int a[32]) {
	unsigned int j;
	unsigned int u;
	u = 0;
	for (j = 0;j < 31;++j) { u += a[j]; a[j] = u & 255; u >>= 8; }
	u += a[31]; a[31] = u & 127;
	u = 19 * (u >> 7);
	for (j = 0;j < 31;++j) { u += a[j]; a[j] = u & 255; u >>= 8; }
	u += a[31]; a[31] = u;
}

static const unsigned int minusp[32] = {
	19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128
} ;

static void freeze(unsigned int a[32]) {
	unsigned int aorig[32];
	unsigned int j;
	unsigned int negative;

	for (j = 0; j < 32; j++) aorig[j] = a[j];
	add(a, a, minusp);
	negative = -((a[31] >> 7) & 1);
	for (j = 0; j < 32; j++) a[j] ^= negative & (aorig[j] ^ a[j]);
}

static void mult(unsigned int out[32], const unsigned int a[32], const unsigned int b[32]) {
	unsigned int i;
	unsigned int j;
	unsigned int u;

	for (i = 0; i < 32; ++i) {
		u = 0;
		for (j = 0;j <= i;++j) u += a[j] * b[i - j];
		for (j = i + 1;j < 32;++j) u += 38 * a[j] * b[i + 32 - j];
		out[i] = u;
	}
	squeeze(out);
}

static void mult_int(unsigned int out[32], const unsigned int n, const unsigned int a[32]) {
	unsigned int j;
	unsigned int u;

	u = 0;
	for (j = 0;j < 31;++j) { u += n * a[j]; out[j] = u & 255; u >>= 8; }
	u += n * a[31]; out[31] = u & 127;
	u = 19 * (u >> 7);
	for (j = 0;j < 31;++j) { u += out[j]; out[j] = u & 255; u >>= 8; }
	u += out[j]; out[j] = u;
}

static void square(unsigned int out[32], const unsigned int a[32]) {
	unsigned int i;
	unsigned int j;
	unsigned int u;

	for (i = 0; i < 32; ++i) {
		u = 0;
		for (j = 0;j < i - j;++j) u += a[j] * a[i - j];
		for (j = i + 1;j < i + 32 - j;++j) u += 38 * a[j] * a[i + 32 - j];
		u *= 2;
		if ((i & 1) == 0) {
			u += a[i / 2] * a[i / 2];
			u += 38 * a[i / 2 + 16] * a[i / 2 + 16];
		}
		out[i] = u;
	}
	squeeze(out);
}

static int check_zero(const unsigned int x[32]) {
	unsigned int differentbits = 0;
	int i;

	for (i = 0; i < 32; i++) {
		differentbits |= (x[i] & 0xffff);
		differentbits |= (x[i] >> 16);
	}

	return (1-(1 & ((differentbits - 1) >> 16)));
}

static int check_equal(const unsigned int x[32], const unsigned int y[32]) {
	unsigned int differentbits = 0;
	int i;

	for (i = 0; i < 32; i++) {
		differentbits |= ((x[i] ^ y[i]) & 0xffff);
		differentbits |= ((x[i] ^ y[i]) >> 16);
	}

	return (1-(1 & ((differentbits - 1) >> 16)));
}

static void selectw(ec_25519_work *out, const ec_25519_work *r, const ec_25519_work *s, unsigned int b) {
	unsigned int j;
	unsigned int t;
	unsigned int bminus1;

	bminus1 = b - 1;
	for (j = 0;j < 32;++j) {
		t = bminus1 & (r->X[j] ^ s->X[j]);
		out->X[j] = s->X[j] ^ t;

		t = bminus1 & (r->Y[j] ^ s->Y[j]);
		out->Y[j] = s->Y[j] ^ t;

		t = bminus1 & (r->Z[j] ^ s->Z[j]);
		out->Z[j] = s->Z[j] ^ t;
	}
}

static void select(unsigned int out[32], const unsigned int r[32], const unsigned int s[32], unsigned int b) {
	unsigned int j;
	unsigned int t;
	unsigned int bminus1;

	bminus1 = b - 1;
	for (j = 0;j < 32;++j) {
		t = bminus1 & (r[j] ^ s[j]);
		out[j] = s[j] ^ t;
	}
}

static const unsigned int rho_s[32] = {
	0xb0, 0xa0, 0x0e, 0x4a, 0x27, 0x1b, 0xee, 0xc4,
	0x78, 0xe4, 0x2f, 0xad, 0x06, 0x18, 0x43, 0x2f,
	0xa7, 0xd7, 0xfb, 0x3d, 0x99, 0x00, 0x4d, 0x2b,
	0x0b, 0xdf, 0xc1, 0x4f, 0x80, 0x24, 0x83, 0x2b
};

static const unsigned int zero[32] = {0};
static const unsigned int minus1[32] = {
	0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
	0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
};

static void square_root(unsigned int out[32], const unsigned int z[32]) {
	/* raise z to the (2^252-2)th power */

	unsigned int z2[32];
	unsigned int z9[32];
	unsigned int z11[32];
	unsigned int z2_5_0[32];
	unsigned int z2_10_0[32];
	unsigned int z2_20_0[32];
	unsigned int z2_50_0[32];
	unsigned int z2_100_0[32];
	unsigned int t0[32];
	unsigned int t1[32];
	unsigned int rt_sq[32];
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

	/* 2^252 - 2 */ mult(t1, t0, z2);

	mult(t0, t1, rho_s);

	square(rt_sq, t1);

	select(out, t0, t1, check_equal(rt_sq, minus1));
}

static void recip(unsigned int out[32], const unsigned int z[32]) {
	unsigned int z2[32];
	unsigned int z9[32];
	unsigned int z11[32];
	unsigned int z2_5_0[32];
	unsigned int z2_10_0[32];
	unsigned int z2_20_0[32];
	unsigned int z2_50_0[32];
	unsigned int z2_100_0[32];
	unsigned int t0[32];
	unsigned int t1[32];
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

void ec_25519_inflate(ec_25519_work *out, const ec_public_key_256 *in) {
	int i;
	unsigned int X2[32], d_X2[32] = {0x04, 0x6d, 0x07} /* 486660 */, a_X2[32] = {0x08, 0x6d, 0x07} /* 486664 */, _1_a_X2[32], d_X2_a_X2[32], Y[32], Yt[32];

	for (i = 0; i < 32; i++) {
		out->X[i] = in->p[i];
		out->Z[i] = (i == 0);
	}

	out->X[31] &= 0x7f;

	square(X2, out->X);
	sub(d_X2, d_X2, X2);
	sub(a_X2, a_X2, X2);
	recip(_1_a_X2, a_X2);
	mult(d_X2_a_X2, d_X2, _1_a_X2);
	square_root(Y, d_X2_a_X2);
	sub(Yt, zero, Y);

	select(out->Y, Y, Yt, in->p[31] >> 7);
}

void ec_25519_deflate(ec_public_key_256 *out, ec_25519_work *in) {
	unsigned int x[32], y[32], z[32];
	int i;

	recip(z, in->Z);

	mult(x, z, in->X);
	mult(y, z, in->Y);

	freeze(x);
	freeze(y);

	for (i = 0; i < 32; i++)
		out->p[i] = x[i];

	out->p[31] |= (y[0] << 7);
}

static const ec_25519_work infty = {{0}, {0}, {1}};

void ec_25519_double(ec_25519_work *out, const ec_25519_work *in) {
	unsigned int A[32], B[32], C[32], D[32], E[32], U[32], t0[32], t1[32], t2[32], t3[32], t4[32], t5[32];

	square(A, in->X);
	square(B, in->Y);
	mult_int(U, 486664, B);
	add(C, A, U);
	sub(D, A, U);
	add(t0, in->X, in->Y);
	square(t1, t0);
	sub(t2, t1, A);
	sub(E, t2, B);
	mult(out->X, C, D);
	square(t3, in->Z);
	mult_int(t4, 973320, t3);
	sub(t5, C, t4);
	mult(out->Y, E, t5);
	mult(out->Z, D, E);
	selectw(out, &infty, out, check_zero(out->X)*check_zero(out->Y));
}

void ec_25519_add(ec_25519_work *out, const ec_25519_work *in1, const ec_25519_work *in2) {
	unsigned int A[32], B[32], C[32], D[32], E[32], H[32], I[32], t0[32], t1[32], t2[32], t3[32], t4[32], t5[32], t6[32], t7[32], t8[32];

	mult(A, in1->Z, in2->Z);
	square(t0, A);
	mult_int(B, 486660, t0);
	mult(C, in1->X, in2->X);
	mult(D, in1->Y, in2->Y);
	mult(E, C, D);
	mult_int(t1, 486664, D);
	sub(H, C, t1);
	add(t2, in1->X, in1->Y);
	add(t3, in2->X, in2->Y);
	mult(t4, t2, t3);
	sub(t5, t4, C);
	sub(I, t5, D);
	add(t6, E, B);
	mult(out->X, t6, H);
	sub(t7, E, B);
	mult(out->Y, t7, I);
	mult(t8, H, I);
	mult(out->Z, A, t8);
	selectw(out, in1, out, check_zero(t3));
	selectw(out, in2, out, check_zero(t2));
}

void ec_25519_scalarmult(ec_25519_work *out, const ec_secret_key_256 *n, const ec_25519_work *base) {
	ec_25519_work Q2, Q2p, cur;
	int i, b, pos;

	for (i = 0; i < 32; i++) {
		cur.X[i] = 0;
		cur.Y[i] = 0;
		cur.Z[i] = (i == 0);
	}

	for (pos = 254;pos >= 0;--pos) {
		b = n->s[pos / 8] >> (pos & 7);
		b &= 1;

		ec_25519_double(&Q2, &cur);
		ec_25519_add(&Q2p, &Q2, base);
		selectw(&cur, &Q2, &Q2p, b);
	}

	for (i = 0; i < 32; i++) {
		out->X[i] = cur.X[i];
		out->Y[i] = cur.Y[i];
		out->Z[i] = cur.Z[i];
	}
}
