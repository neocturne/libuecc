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
  Simple finite field operations on the prime field F_q for
  q = 2^252 + 27742317777372353535851937790883648493, which
  is the order of the base point used for ec25519
*/

#include <libuecc/ecc.h>


#define IS_NEGATIVE(n) ((int)((((unsigned)n) >> (8*sizeof(n)-1))&1))


static const unsigned char q[32] = {
	0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
	0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
};


static void select(unsigned char out[32], const unsigned char r[32], const unsigned char s[32], unsigned int b) {
	unsigned int j;
	unsigned int t;
	unsigned int bminus1;

	bminus1 = b - 1;
	for (j = 0;j < 32;++j) {
		t = bminus1 & (r[j] ^ s[j]);
		out[j] = s[j] ^ t;
	}
}

int ecc_25519_secret_is_zero(const ecc_secret_key_256 *in) {
	int i;
	ecc_secret_key_256 r;
	unsigned int bits;

	ecc_25519_secret_reduce(&r, in);

	for (i = 0; i < 32; i++)
		bits |= r.s[i];

	return (((bits-1)>>8) & 1);
}

void ecc_25519_secret_add(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2) {
	unsigned int j;
	int u1, u2, u3;
	unsigned char out1[32], out2[32], out3[32];

	u1 = u2 = u3 = 0;
	for (j = 0; j < 31; ++j) {
		u1 += in1->s[j] + in2->s[j];
		u2 += in1->s[j] + in2->s[j] - 8*q[j];
		u3 += in1->s[j] + in2->s[j] - 16*q[j];

		out1[j] = u1; out2[j] = u2; out3[j] = u3;
		u1 = (u1+IS_NEGATIVE(u1))/256 - IS_NEGATIVE(u1);
		u2 = (u2+IS_NEGATIVE(u2))/256 - IS_NEGATIVE(u2);
		u3 = (u3+IS_NEGATIVE(u3))/256 - IS_NEGATIVE(u3);
	}
	u1 += in1->s[31] + in2->s[31];
	u2 += in1->s[31] + in2->s[31] - 8*q[31];
	u3 += in1->s[31] + in2->s[31] - 16*q[31];
	out1[31] = u1; out2[31] = u2; out3[31] = u3;

	select(out->s, out1, out2, (u1 >> 8) & 1);
	select(out->s, out->s, out3, ((u1 >> 8) & (u2 >> 8)) & 1);
}

void ecc_25519_secret_sub(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2) {
	unsigned int j;
	int u1, u2, u3;
	unsigned char out1[32], out2[32], out3[32];

	u1 = u2 = u3 = 0;
	for (j = 0; j < 31; ++j) {
		u1 += in1->s[j] - in2->s[j] + 16*q[j];
		u2 += in1->s[j] - in2->s[j] + 8*q[j];
		u3 += in1->s[j] - in2->s[j];

		out1[j] = u1; out2[j] = u2; out3[j] = u3;
		u1 = (u1+IS_NEGATIVE(u1))/256 - IS_NEGATIVE(u1);
		u2 = (u2+IS_NEGATIVE(u2))/256 - IS_NEGATIVE(u2);
		u3 = (u3+IS_NEGATIVE(u3))/256 - IS_NEGATIVE(u3);
	}
	u1 += in1->s[31] - in2->s[31] + 16*q[31];
	u2 += in1->s[31] - in2->s[31] + 8*q[31];
	u3 += in1->s[31] - in2->s[31];
	out1[31] = u1; out2[31] = u2; out3[31] = u3;

	select(out->s, out1, out2, (u1 >> 8) & 1);
	select(out->s, out->s, out3, ((u1 >> 8) & (u2 >> 8)) & 1);
}

static void reduce(unsigned char a[32]) {
	unsigned int j;
	int nq = a[31] >> 4;
	int u1, u2;
	unsigned char out1[32], out2[32];

	u1 = u2 = 0;
	for (j = 0; j < 31; ++j) {
		u1 += a[j] - nq*q[j];
		u2 += a[j] - (nq-1)*q[j];

		out1[j] = u1; out2[j] = u2;
		u1 = (u1+IS_NEGATIVE(u1))/256 - IS_NEGATIVE(u1);
		u2 = (u2+IS_NEGATIVE(u2))/256 - IS_NEGATIVE(u2);
	}
	u1 += a[31] - nq*q[31];
	u2 += a[31] - (nq-1)*q[31];
	out1[31] = u1; out2[31] = u2;

	select(a, out1, out2, IS_NEGATIVE(u1));
}

void ecc_25519_secret_reduce(ecc_secret_key_256 *out, const ecc_secret_key_256 *in) {
	int i;

	for (i = 0; i < 32; i++)
		out->s[i] = in->s[i];

	reduce(out->s);
}

/* Montgomery modular multiplication algorithm */
static void montgomery(unsigned char out[32], const unsigned char a[32], const unsigned char b[32]) {
	unsigned int i, j;
	unsigned int nq;
	unsigned int u;

	for (i = 0; i < 32; i++)
		out[i] = 0;

	for (i = 0; i < 32; i++) {
		u = out[0] + a[i]*b[0];
		nq = (u*27) & 255;
		u += nq*q[0];

		for (j = 1; j < 32; ++j) {
			u += (out[j] + a[i]*b[j] + nq*q[j]) << 8;
			u >>= 8;
			out[j-1] = u;
		}

		out[31] = u >> 8;
	}
}


void ecc_25519_secret_mult(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2) {
	/* 2^512 mod q */
	static const unsigned char C[32] = {
		0x01, 0x0f, 0x9c, 0x44, 0xe3, 0x11, 0x06, 0xa4,
		0x47, 0x93, 0x85, 0x68, 0xa7, 0x1b, 0x0e, 0xd0,
		0x65, 0xbe, 0xf5, 0x17, 0xd2, 0x73, 0xec, 0xce,
		0x3d, 0x9a, 0x30, 0x7c, 0x1b, 0x41, 0x99, 0x03
	};

	unsigned char B[32];
	unsigned char R[32];
	unsigned int i;

	for (i = 0; i < 32; i++)
		B[i] = in2->s[i];

	reduce(B);

	montgomery(R, in1->s, B);
	montgomery(out->s, R, C);
}

