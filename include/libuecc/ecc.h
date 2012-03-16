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

typedef struct _ecc_public_key_256 {
	unsigned char p[32];
} ecc_public_key_256;

typedef struct _ecc_secret_key_256 {
	unsigned char s[32];
} ecc_secret_key_256;

typedef struct _ecc_25519_work {
	unsigned int X[32];
	unsigned int Y[32];
	unsigned int Z[32];
} ecc_25519_work;


void ecc_25519_load(ecc_25519_work *out, const ecc_public_key_256 *in);
void ecc_25519_store(ecc_public_key_256 *out, const ecc_25519_work *in);

int ecc_25519_is_identity(const ecc_25519_work *in);
void ecc_25519_add(ecc_25519_work *out, const ecc_25519_work *in1, const ecc_25519_work *in2);
void ecc_25519_double(ecc_25519_work *out, const ecc_25519_work *in);
void ecc_25519_scalarmult(ecc_25519_work *out, const ecc_secret_key_256 *n, const ecc_25519_work *base);
void ecc_25519_scalarmult_base(ecc_25519_work *out, const ecc_secret_key_256 *n);

int ecc_25519_secret_is_zero(const ecc_secret_key_256 *in);
void ecc_25519_secret_add(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2);
void ecc_25519_secret_sub(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2);
void ecc_25519_secret_reduce(ecc_secret_key_256 *out, const ecc_secret_key_256 *in);
void ecc_25519_secret_mult(ecc_secret_key_256 *out, const ecc_secret_key_256 *in1, const ecc_secret_key_256 *in2);

#endif /* _LIBUECC_ECC_H_ */
