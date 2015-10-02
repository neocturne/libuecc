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

#ifndef _LIBUECC_ECC_H_
#define _LIBUECC_ECC_H_

/**
 * A 256 bit integer
 *
 * All functions of libuecc treat \ref ecc_int256_t as unsigned little-endian.
 */
typedef union _ecc_int256 {
	/** Data bytes */
	unsigned char p[32];
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

/** The identity element */
extern const ecc_25519_work_t ecc_25519_work_identity;

/** The ec25519 default base */
extern const ecc_25519_work_t ecc_25519_work_default_base;



/** Loads a point with given coordinates into its unpacked representation */
int ecc_25519_load_xy(ecc_25519_work_t *out, const ecc_int256_t *x, const ecc_int256_t *y);

/**
 * Stores a point's x and y coordinates
 *
 * \param x Returns the x coordinate of the point. May be NULL.
 * \param y Returns the y coordinate of the point. May be NULL.
 * \param in The unpacked point to store.
 */
void ecc_25519_store_xy(ecc_int256_t *x, ecc_int256_t *y, const ecc_25519_work_t *in);


/** Loads a packed point into its unpacked representation */
int ecc_25519_load_packed(ecc_25519_work_t *out, const ecc_int256_t *in);

/** Stores a point into its packed representation */
void ecc_25519_store_packed(ecc_int256_t *out, const ecc_25519_work_t *in);


/** Checks if a point is the identity element of the Elliptic Curve group */
int ecc_25519_is_identity(const ecc_25519_work_t *in);

/**
 * Negates a point of the Elliptic Curve
 *
 * The same pointer may be given for input and output
 */
void ecc_25519_negate(ecc_25519_work_t *out, const ecc_25519_work_t *in);

/**
 * Doubles a point of the Elliptic Curve
 *
 * ecc_25519_double(out, in) is equivalent to ecc_25519_add(out, in, in), but faster.
 *
 * The same pointer may be given for input and output.
 */
void ecc_25519_double(ecc_25519_work_t *out, const ecc_25519_work_t *in);

/**
 * Adds two points of the Elliptic Curve
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_add(ecc_25519_work_t *out, const ecc_25519_work_t *in1, const ecc_25519_work_t *in2);

/**
 * Subtracts two points of the Elliptic Curve
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_sub(ecc_25519_work_t *out, const ecc_25519_work_t *in1, const ecc_25519_work_t *in2);

/**
 * Does a scalar multiplication of a point of the Elliptic Curve with an integer of a given bit length
 *
 * To speed up scalar multiplication when it is known that not the whole 256 bits of the scalar
 * are used. The bit length should always be a constant and not computed at runtime to ensure
 * that no timing attacks are possible.
 *
 * The same pointer may be given for input and output.
 **/
void ecc_25519_scalarmult_bits(ecc_25519_work_t *out, const ecc_int256_t *n, const ecc_25519_work_t *base, unsigned bits);

/**
 * Does a scalar multiplication of a point of the Elliptic Curve with an integer
 *
 * The same pointer may be given for input and output.
 **/
void ecc_25519_scalarmult(ecc_25519_work_t *out, const ecc_int256_t *n, const ecc_25519_work_t *base);

/**
 * Does a scalar multiplication of the default base point (generator element) of the Elliptic Curve with an integer of a given bit length
 *
 * The order of the base point is \f$ 2^{252} + 27742317777372353535851937790883648493 \f$.
 *
 * See the notes about \ref ecc_25519_scalarmult_bits before using this function.
 */
void ecc_25519_scalarmult_base_bits(ecc_25519_work_t *out, const ecc_int256_t *n, unsigned bits);

/**
 * Does a scalar multiplication of the default base point (generator element) of the Elliptic Curve with an integer
 *
 * The order of the base point is \f$ 2^{252} + 27742317777372353535851937790883648493 \f$.
 */
void ecc_25519_scalarmult_base(ecc_25519_work_t *out, const ecc_int256_t *n);

/**@}*/

/**
 * \defgroup gf_ops Prime field operations for the order of the base point of the Elliptic Curve
 * @{
 */

/**
 * The order of the prime field
 *
 * The order is \f$ 2^{252} + 27742317777372353535851937790883648493 \f$.
 */
extern const ecc_int256_t ecc_25519_gf_order;


/** Checks if an integer is equal to zero (after reduction) */
int ecc_25519_gf_is_zero(const ecc_int256_t *in);

/**
 * Adds two integers as Galois field elements
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_gf_add(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);

/**
 * Subtracts two integers as Galois field elements
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_gf_sub(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);

/**
 * Reduces an integer to a unique representation in the range \f$ [0,q-1] \f$
 *
 * The same pointer may be given for input and output.
 */
void ecc_25519_gf_reduce(ecc_int256_t *out, const ecc_int256_t *in);

/**
 * Multiplies two integers as Galois field elements
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_gf_mult(ecc_int256_t *out, const ecc_int256_t *in1, const ecc_int256_t *in2);

/**
 * Computes the reciprocal of a Galois field element
 *
 * The same pointers may be given for input and output.
 */
void ecc_25519_gf_recip(ecc_int256_t *out, const ecc_int256_t *in);

/**
 * Ensures some properties of a Galois field element to make it fit for use as a secret key
 *
 * This sets the 255th bit and clears the 256th and the bottom three bits (so the key
 * will be a multiple of 8). See Daniel J. Bernsteins paper "Curve25519: new Diffie-Hellman speed records."
 * for the rationale of this.
 *
 * The same pointer may be given for input and output.
 */
void ecc_25519_gf_sanitize_secret(ecc_int256_t *out, const ecc_int256_t *in);

/**@}*/

#endif /* _LIBUECC_ECC_H_ */
