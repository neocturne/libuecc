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

#ifndef DEPRECATED
#define DEPRECATED __attribute__((deprecated))
#endif

/**
 * A 256 bit integer
 *
 * All functions of libuecc treat \ref ecc_int256_t as unsigned little-endian.
 */
typedef union _ecc_int256 {
	/** Data bytes */
	unsigned char p[32];

	/**
	 * Old name of p
	 * \deprecated Use \ref ecc_int256_t::p instead.
	 */
	unsigned char s[32] DEPRECATED;
} ecc_int256_t;

/**
 * A point on the curve unpacked for efficient calculation
 *
 * The internal representation of an unpacked point isn't unique, so for serialization
 * it should always be packed.
 */
typedef struct _ecc_25519_work {
	unsigned int X[32];
	unsigned int Y[32];
	unsigned int Z[32];
	unsigned int T[32];
} ecc_25519_work_t;

/**
 * \defgroup curve_ops Operations on points of the Elliptic Curve
 * @{
 */

void ecc_25519_load_xy(ecc_25519_work_t *out, const ecc_int256_t *x, const ecc_int256_t *y);
void ecc_25519_store_xy(ecc_int256_t *x, ecc_int256_t *y, const ecc_25519_work_t *in);

void ecc_25519_load_packed(ecc_25519_work_t *out, const ecc_int256_t *in);
void ecc_25519_store_packed(ecc_int256_t *out, const ecc_25519_work_t *in);

int ecc_25519_is_identity(const ecc_25519_work_t *in);
void ecc_25519_double(ecc_25519_work_t *out, const ecc_25519_work_t *in);
void ecc_25519_add(ecc_25519_work_t *out, const ecc_25519_work_t *in1, const ecc_25519_work_t *in2);
void ecc_25519_scalarmult(ecc_25519_work_t *out, const ecc_int256_t *n, const ecc_25519_work_t *base);
void ecc_25519_scalarmult_base(ecc_25519_work_t *out, const ecc_int256_t *n);

/**@}*/

/**
 * \defgroup gf_ops Prime field operations for the order of the base point of the Elliptic Curve
 * @{
 */

extern const ecc_int256_t ecc_25519_gf_order;

int ecc_25519_gf_is_zero(const ecc_int256_t *in);
void ecc_25519_gf_add(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);
void ecc_25519_gf_sub(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);
void ecc_25519_gf_reduce(ecc_int256_t *out, const ecc_int256_t *in);
void ecc_25519_gf_mult(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);
void ecc_25519_gf_recip(ecc_int256_t *out, const ecc_int256_t *in);
void ecc_25519_gf_sanitize_secret(ecc_int256_t *out, const ecc_int256_t *in);

/**@}*/

/* declarations for the old names */

/**
 * Old name of \ref ecc_int256_t
 * \deprecated Use \ref ecc_int256_t instead.
 */
typedef ecc_int256_t ecc_secret_key_256 DEPRECATED;

/**
 * Old name of \ref ecc_int256_t
 * \deprecated Use \ref ecc_int256_t instead.
 */
typedef ecc_int256_t ecc_public_key_256 DEPRECATED;

/**
 * Old name of \ref ecc_25519_work_t
 * \deprecated Use \ref ecc_25519_work_t instead.
 */
typedef ecc_25519_work_t ecc_25519_work DEPRECATED;


/**
 * Loads a packed point into its unpacked representation
 *
 * \deprecated Use \ref ecc_25519_load_packed instead.
 */
DEPRECATED static inline void ecc_25519_load(ecc_25519_work_t *out, const ecc_int256_t *in) {
	ecc_25519_load_packed(out, in);
}

/**
 * Stores a point into its packed representation
 *
 * \deprecated Use \ref ecc_25519_store_packed instead.
 */
DEPRECATED static inline void ecc_25519_store(ecc_int256_t *out, const ecc_25519_work_t *in) {
	ecc_25519_store_packed(out, in);
}

/**
 * Checks if an integer is equal to zero (after reduction)
 *
 * \deprecated Use \ref ecc_25519_gf_is_zero instead.
 */
DEPRECATED static inline int ecc_25519_secret_is_zero(const ecc_int256_t *in) {
	return ecc_25519_gf_is_zero(in);
}

/**
 * Adds two integers as Galois field elements
 *
 * \deprecated Use \ref ecc_25519_gf_add instead.
 */
DEPRECATED static inline void ecc_25519_secret_add(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2) {
	ecc_25519_gf_add(out, in1, in2);
}

/**
 * Subtracts two integers as Galois field elements
 *
 * \deprecated Use \ref ecc_25519_gf_sub instead.
 */
DEPRECATED static inline void ecc_25519_secret_sub(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2) {
	ecc_25519_gf_sub(out, in1, in2);
}

/**
 * Reduces an integer to a unique representation in the range \f$ [0,q-1] \f$
 *
 * \deprecated Use \ref ecc_25519_gf_reduce instead.
 */
DEPRECATED static inline void ecc_25519_secret_reduce(ecc_int256_t *out, const ecc_int256_t *in) {
	ecc_25519_gf_reduce(out, in);
}

/**
 * Multiplies to integers as Galois field elements
 *
 * \deprecated Use \ref ecc_25519_gf_mult instead.
 */
DEPRECATED static inline void ecc_25519_secret_mult(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2) {
	ecc_25519_gf_mult(out, in1, in2);
}

/**
 * Ensures some properties of a Galois field element to make it fit for use as a secret key
 *
 * \deprecated Use \ref ecc_25519_gf_sanitize_secret instead.
 */
DEPRECATED static inline void ecc_25519_secret_sanitize(ecc_int256_t *out, const ecc_int256_t *in) {
	ecc_25519_gf_sanitize_secret(out, in);
}

#endif /* _LIBUECC_ECC_H_ */
