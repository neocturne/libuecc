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

#ifndef _LIBUECC_ECC_H_
#define _LIBUECC_ECC_H_

typedef union _ecc_int_256 {
	unsigned char p[32];

	/* old name */
	unsigned char s[32];
} ecc_int_256, ecc_secret_key_256, ecc_public_key_256;

/* a point on the curve unpacked for efficient calculation */
typedef struct _ecc_25519_work {
	unsigned int X[32];
	unsigned int Y[32];
	unsigned int Z[32];
	unsigned int T[32];
} ecc_25519_work;


void ecc_25519_load_xy(ecc_25519_work *out, const ecc_int_256 *x, const ecc_int_256 *y);
void ecc_25519_store_xy(ecc_int_256 *x, ecc_int_256 *y, const ecc_25519_work *in);

void ecc_25519_load_packed(ecc_25519_work *out, const ecc_int_256 *in);
void ecc_25519_store_packed(ecc_int_256 *out, const ecc_25519_work *in);

int ecc_25519_is_identity(const ecc_25519_work *in);
void ecc_25519_add(ecc_25519_work *out, const ecc_25519_work *in1, const ecc_25519_work *in2);
void ecc_25519_double(ecc_25519_work *out, const ecc_25519_work *in);
void ecc_25519_scalarmult(ecc_25519_work *out, const ecc_int_256 *n, const ecc_25519_work *base);
void ecc_25519_scalarmult_base(ecc_25519_work *out, const ecc_int_256 *n);

/* operations on elements of the prime field F_q for q = 2^252 + 27742317777372353535851937790883648493 */
int ecc_25519_gf_is_zero(const ecc_int_256 *in);
void ecc_25519_gf_add(ecc_int_256 *out, const ecc_int_256 *in1, const ecc_int_256 *in2);
void ecc_25519_gf_sub(ecc_int_256 *out, const ecc_int_256 *in1, const ecc_int_256 *in2);
void ecc_25519_gf_reduce(ecc_int_256 *out, const ecc_int_256 *in);
void ecc_25519_gf_mult(ecc_int_256 *out, const ecc_int_256 *in1, const ecc_int_256 *in2);
void ecc_25519_gf_recip(ecc_int_256 *out, const ecc_int_256 *in);

void ecc_25519_gf_sanitize_secret(ecc_int_256 *out, const ecc_int_256 *in);

/* defines for the old names */
#define ecc_25519_load            ecc_25519_load_packed
#define ecc_25519_store           ecc_25519_store_packed

#define ecc_25519_secret_is_zero  ecc_25519_gf_is_zero
#define ecc_25519_secret_add      ecc_25519_gf_add
#define ecc_25519_secret_sub      ecc_25519_gf_sub
#define ecc_25519_secret_reduce   ecc_25519_gf_reduce
#define ecc_25519_secret_mult     ecc_25519_gf_mult
#define ecc_25519_secret_sanitize ecc_25519_gf_sanitize_secret

#endif /* _LIBUECC_ECC_H_ */
