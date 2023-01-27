#include "blst.h"

typedef blst_p1 g1_t;         /**< Internal G1 group element type */
typedef blst_p2 g2_t;         /**< Internal G2 group element type */
typedef blst_fr fr_t;         /**< Internal Fr field element type */

/**
 * Stores the setup and parameters needed for performing FFTs.
 */
typedef struct {
    uint64_t max_width;            /**< The maximum size of FFT these settings support, a power of 2. */
    fr_t *expanded_roots_of_unity; /**< Ascending powers of the root of unity, size `width + 1`. */
    fr_t *reverse_roots_of_unity;  /**< Descending powers of the root of unity, size `width + 1`. */
    fr_t *roots_of_unity;          /**< Powers of the root of unity in bit-reversal permutation, size `width`. */
} FFTSettings;

/**
 * The common return type for all routines in which something can go wrong.
 */
typedef enum {
    C_KZG_OK = 0,  /**< Success! */
    C_KZG_BADARGS, /**< The supplied data is invalid in some way */
    C_KZG_ERROR,   /**< Internal error - this should never occur and may indicate a bug in the library */
    C_KZG_MALLOC,  /**< Could not allocate memory */
} C_KZG_RET;

/** This is 1 in Blst's `blst_fr` limb representation. Crazy but true. */
static const fr_t fr_one = {0x00000001fffffffeL, 0x5884b7fa00034802L, 0x998c4fefecbc4ff5L, 0x1824b159acc5056fL};

C_KZG_RET new_fr_array(fr_t **x, size_t n);
C_KZG_RET c_kzg_malloc(void **x, size_t n);

int log_2_byte(byte b);

bool fr_is_one(const fr_t *p);

bool fr_equal(const fr_t *aa, const fr_t *bb);

void fr_div(fr_t *out, const fr_t *a, const fr_t *b);

void fr_pow(fr_t *out, const fr_t *a, uint64_t n);

void fr_from_uint64(fr_t *out, uint64_t n);

C_KZG_RET fr_batch_inv(fr_t *out, const fr_t *a, size_t len);

void g1_mul(g1_t *out, const g1_t *a, const fr_t *b);

void g1_sub(g1_t *out, const g1_t *a, const g1_t *b);

void g2_sub(g2_t *out, const g2_t *a, const g2_t *b);

void g2_mul(g2_t *out, const g2_t *a, const fr_t *b);

bool pairings_verify(const g1_t *a1, const g2_t *a2, const g1_t *b1, const g2_t *b2);

int log2_pow2(uint32_t n);
