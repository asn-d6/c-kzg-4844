#define CHECK(cond) \
    if (!(cond)) return C_KZG_BADARGS


/**
 * The first 32 roots of unity in the finite field F_r.
 *
 * For element `{A, B, C, D}`, the field element value is `A + B * 2^64 + C * 2^128 + D * 2^192`. This format may be
 * converted to an `fr_t` type via the #blst_fr_from_uint64 function.
 *
 * The decimal values may be calculated with the following Python code:
 * @code{.py}
 * MODULUS = 52435875175126190479447740508185965837690552500527637822603658699938581184513
 * PRIMITIVE_ROOT = 7
 * [pow(PRIMITIVE_ROOT, (MODULUS - 1) // (2**i), MODULUS) for i in range(32)]
 * @endcode
 *
 * Note: Being a "primitive root" in this context means that r^k != 1 for any k < q-1 where q is the modulus. So
 * powers of r generate the field. This is also known as being a "primitive element".
 *
 * This is easy to check for: we just require that r^((q-1)/2) != 1. Instead of 5, we could use 7, 10, 13, 14, 15, 20...
 * to create the roots of unity below. There are a lot of primitive roots:
 * https://crypto.stanford.edu/pbc/notes/numbertheory/gen.html
 */
static const uint64_t SCALE2_ROOT_OF_UNITY[][4] = {
    {0x0000000000000001L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L},
    {0xffffffff00000000L, 0x53bda402fffe5bfeL, 0x3339d80809a1d805L, 0x73eda753299d7d48L},
    {0x0001000000000000L, 0xec03000276030000L, 0x8d51ccce760304d0L, 0x0000000000000000L},
    {0x7228fd3397743f7aL, 0xb38b21c28713b700L, 0x8c0625cd70d77ce2L, 0x345766f603fa66e7L},
    {0x53ea61d87742bcceL, 0x17beb312f20b6f76L, 0xdd1c0af834cec32cL, 0x20b1ce9140267af9L},
    {0x360c60997369df4eL, 0xbf6e88fb4c38fb8aL, 0xb4bcd40e22f55448L, 0x50e0903a157988baL},
    {0x8140d032f0a9ee53L, 0x2d967f4be2f95155L, 0x14a1e27164d8fdbdL, 0x45af6345ec055e4dL},
    {0x5130c2c1660125beL, 0x98d0caac87f5713cL, 0xb7c68b4d7fdd60d0L, 0x6898111413588742L},
    {0x4935bd2f817f694bL, 0x0a0865a899e8deffL, 0x6b368121ac0cf4adL, 0x4f9b4098e2e9f12eL},
    {0x4541b8ff2ee0434eL, 0xd697168a3a6000feL, 0x39feec240d80689fL, 0x095166525526a654L},
    {0x3c28d666a5c2d854L, 0xea437f9626fc085eL, 0x8f4de02c0f776af3L, 0x325db5c3debf77a1L},
    {0x4a838b5d59cd79e5L, 0x55ea6811be9c622dL, 0x09f1ca610a08f166L, 0x6d031f1b5c49c834L},
    {0xe206da11a5d36306L, 0x0ad1347b378fbf96L, 0xfc3e8acfe0f8245fL, 0x564c0a11a0f704f4L},
    {0x6fdd00bfc78c8967L, 0x146b58bc434906acL, 0x2ccddea2972e89edL, 0x485d512737b1da3dL},
    {0x034d2ff22a5ad9e1L, 0xae4622f6a9152435L, 0xdc86b01c0d477fa6L, 0x56624634b500a166L},
    {0xfbd047e11279bb6eL, 0xc8d5f51db3f32699L, 0x483405417a0cbe39L, 0x3291357ee558b50dL},
    {0xd7118f85cd96b8adL, 0x67a665ae1fcadc91L, 0x88f39a78f1aeb578L, 0x2155379d12180caaL},
    {0x08692405f3b70f10L, 0xcd7f2bd6d0711b7dL, 0x473a2eef772c33d6L, 0x224262332d8acbf4L},
    {0x6f421a7d8ef674fbL, 0xbb97a3bf30ce40fdL, 0x652f717ae1c34bb0L, 0x2d3056a530794f01L},
    {0x194e8c62ecb38d9dL, 0xad8e16e84419c750L, 0xdf625e80d0adef90L, 0x520e587a724a6955L},
    {0xfece7e0e39898d4bL, 0x2f69e02d265e09d9L, 0xa57a6e07cb98de4aL, 0x03e1c54bcb947035L},
    {0xcd3979122d3ea03aL, 0x46b3105f04db5844L, 0xc70d0874b0691d4eL, 0x47c8b5817018af4fL},
    {0xc6e7a6ffb08e3363L, 0xe08fec7c86389beeL, 0xf2d38f10fbb8d1bbL, 0x0abe6a5e5abcaa32L},
    {0x5616c57de0ec9eaeL, 0xc631ffb2585a72dbL, 0x5121af06a3b51e3cL, 0x73560252aa0655b2L},
    {0x92cf4deb77bd779cL, 0x72cf6a8029b7d7bcL, 0x6e0bcd91ee762730L, 0x291cf6d68823e687L},
    {0xce32ef844e11a51eL, 0xc0ba12bb3da64ca5L, 0x0454dc1edc61a1a3L, 0x019fe632fd328739L},
    {0x531a11a0d2d75182L, 0x02c8118402867ddcL, 0x116168bffbedc11dL, 0x0a0a77a3b1980c0dL},
    {0xe2d0a7869f0319edL, 0xb94f1101b1d7a628L, 0xece8ea224f31d25dL, 0x23397a9300f8f98bL},
    {0xd7b688830a4f2089L, 0x6558e9e3f6ac7b41L, 0x99e276b571905a7dL, 0x52dd465e2f094256L},
    {0x474650359d8e211bL, 0x84d37b826214abc6L, 0x8da40c1ef2bb4598L, 0x0c83ea7744bf1beeL},
    {0x694341f608c9dd56L, 0xed3a181fabb30adcL, 0x1339a815da8b398fL, 0x2c6d4e4511657e1eL},
    {0x63e7cb4906ffc93fL, 0xf070bb00e28a193dL, 0xad1715b02e5713b5L, 0x4b5371495990693fL}
};


/**
 * Wrapped `malloc()` that reports failures to allocate.
 *
 * @param[out] x Pointer to the allocated space
 * @param[in]  n The number of bytes to be allocated
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_MALLOC  Memory allocation failed
 */
static C_KZG_RET c_kzg_malloc(void **x, size_t n) {
    if (n > 0) {
        *x = malloc(n);
        return *x != NULL ? C_KZG_OK : C_KZG_MALLOC;
    }
    *x = NULL;
    return C_KZG_OK;
}

/**
 * Allocate memory for an array of G1 group elements.
 *
 * @remark Free the space later using `free()`.
 *
 * @param[out] x Pointer to the allocated space
 * @param[in]  n The number of G1 elements to be allocated
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_MALLOC  Memory allocation failed
 */
static C_KZG_RET new_g1_array(g1_t **x, size_t n) {
    return c_kzg_malloc((void **)x, n * sizeof **x);
}

/**
 * Allocate memory for an array of G2 group elements.
 *
 * @remark Free the space later using `free()`.
 *
 * @param[out] x Pointer to the allocated space
 * @param[in]  n The number of G2 elements to be allocated
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_MALLOC  Memory allocation failed
 */
static C_KZG_RET new_g2_array(g2_t **x, size_t n) {
    return c_kzg_malloc((void **)x, n * sizeof **x);
}

/**
 * Discrete fourier transforms over arrays of G1 group elements.
 *
 * Also known as [number theoretic
 * transforms](https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general)#Number-theoretic_transform).
 *
 * @remark Functions here work only for lengths that are a power of two.
 */

/**
 * Fast Fourier Transform.
 *
 * Recursively divide and conquer.
 *
 * @param[out] out          The results (array of length @p n)
 * @param[in]  in           The input data (array of length @p n * @p stride)
 * @param[in]  stride       The input data stride
 * @param[in]  roots        Roots of unity (array of length @p n * @p roots_stride)
 * @param[in]  roots_stride The stride interval among the roots of unity
 * @param[in]  n            Length of the FFT, must be a power of two
 */
static void fft_g1_fast(g1_t *out, const g1_t *in, uint64_t stride, const fr_t *roots, uint64_t roots_stride,
                        uint64_t n) {
    uint64_t half = n / 2;
    if (half > 0) { // Tunable parameter
        fft_g1_fast(out, in, stride * 2, roots, roots_stride * 2, half);
        fft_g1_fast(out + half, in + stride, stride * 2, roots, roots_stride * 2, half);
        for (uint64_t i = 0; i < half; i++) {
            g1_t y_times_root;
            g1_mul(&y_times_root, &out[i + half], &roots[i * roots_stride]);
            g1_sub(&out[i + half], &out[i], &y_times_root);
            blst_p1_add_or_double(&out[i], &out[i], &y_times_root);
        }
    } else {
        *out = *in;
    }
}

/**
 * The main entry point for forward and reverse FFTs over the finite field.
 *
 * @param[out] out     The results (array of length @p n)
 * @param[in]  in      The input data (array of length @p n)
 * @param[in]  inverse `false` for forward transform, `true` for inverse transform
 * @param[in]  n       Length of the FFT, must be a power of two
 * @param[in]  fs      Pointer to previously initialised FFTSettings structure with `max_width` at least @p n.
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_BADARGS Invalid parameters were supplied
 */
static C_KZG_RET fft_g1(g1_t *out, const g1_t *in, bool inverse, uint64_t n, const FFTSettings *fs) {
    uint64_t stride = fs->max_width / n;
    CHECK(n <= fs->max_width);
    CHECK(is_power_of_two(n));
    if (inverse) {
        fr_t inv_len;
        fr_from_uint64(&inv_len, n);
        blst_fr_eucl_inverse(&inv_len, &inv_len);
        fft_g1_fast(out, in, 1, fs->reverse_roots_of_unity, stride, n);
        for (uint64_t i = 0; i < n; i++) {
            g1_mul(&out[i], &out[i], &inv_len);
        }
    } else {
        fft_g1_fast(out, in, 1, fs->expanded_roots_of_unity, stride, n);
    }
    return C_KZG_OK;
}

/**
 * Generate powers of a root of unity in the field for use in the FFTs.
 *
 * @remark @p root must be such that @p root ^ @p width is equal to one, but no smaller power of @p root is equal to
 * one.
 *
 * @param[out] out   The generated powers of the root of unity (array size @p width + 1)
 * @param[in]  root  A root of unity
 * @param[in]  width One less than the size of @p out
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_BADARGS Invalid parameters were supplied
 */
static C_KZG_RET expand_root_of_unity(fr_t *out, const fr_t *root, uint64_t width) {
    out[0] = fr_one;
    out[1] = *root;

    for (uint64_t i = 2; !fr_is_one(&out[i - 1]); i++) {
        CHECK(i <= width);
        blst_fr_mul(&out[i], &out[i - 1], root);
    }
    CHECK(fr_is_one(&out[width]));

    return C_KZG_OK;
}

/**
 * Initialise an FFTSettings structure.
 *
 * Space is allocated for, and arrays are populated with, powers of the roots of unity. The two arrays contain the same
 * values in reverse order for convenience in inverse FFTs.
 *
 * `max_width` is the maximum size of FFT that can be calculated with these settings, and is a power of two by
 * construction. The same settings may be used to calculated FFTs of smaller power sizes.
 *
 * @remark As with all functions prefixed `new_`, this allocates memory that needs to be reclaimed by calling the
 * corresponding `free_` function. In this case, #free_fft_settings.
 * @remark These settings may be used for FFTs on both field elements and G1 group elements.
 *
 * @param[out] fs        The new settings
 * @param[in]  max_scale Log base 2 of the max FFT size to be used with these settings
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_BADARGS Invalid parameters were supplied
 * @retval C_CZK_ERROR   An internal error occurred
 * @retval C_CZK_MALLOC  Memory allocation failed
 */
static C_KZG_RET new_fft_settings(FFTSettings *fs, unsigned int max_scale) {
    C_KZG_RET ret;
    fr_t root_of_unity;

    fs->max_width = (uint64_t)1 << max_scale;
    fs->expanded_roots_of_unity = NULL;
    fs->reverse_roots_of_unity = NULL;
    fs->roots_of_unity = NULL;

    CHECK((max_scale < sizeof SCALE2_ROOT_OF_UNITY / sizeof SCALE2_ROOT_OF_UNITY[0]));
    blst_fr_from_uint64(&root_of_unity, SCALE2_ROOT_OF_UNITY[max_scale]);

    // Allocate space for the roots of unity
    ret = new_fr_array(&fs->expanded_roots_of_unity, fs->max_width + 1);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_fr_array(&fs->reverse_roots_of_unity, fs->max_width + 1);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_fr_array(&fs->roots_of_unity, fs->max_width);
    if (ret != C_KZG_OK) goto out_error;

    // Populate the roots of unity
    ret = expand_root_of_unity(fs->expanded_roots_of_unity, &root_of_unity, fs->max_width);
    if (ret != C_KZG_OK) goto out_error;

    // Populate reverse roots of unity
    for (uint64_t i = 0; i <= fs->max_width; i++) {
        fs->reverse_roots_of_unity[i] = fs->expanded_roots_of_unity[fs->max_width - i];
    }

    // Permute the roots of unity
    memcpy(fs->roots_of_unity, fs->expanded_roots_of_unity, sizeof(fr_t) * fs->max_width);
    ret = bit_reversal_permutation(fs->roots_of_unity, sizeof(fr_t), fs->max_width);
    if (ret != C_KZG_OK) goto out_error;

    goto out_success;

out_error:
    free(fs->expanded_roots_of_unity);
    free(fs->reverse_roots_of_unity);
    free(fs->roots_of_unity);
out_success:
    return ret;
}

/**
 * Free the memory that was previously allocated by #new_fft_settings.
 *
 * @param fs The settings to be freed
 */
static void free_fft_settings(FFTSettings *fs) {
    free(fs->expanded_roots_of_unity);
    free(fs->reverse_roots_of_unity);
    free(fs->roots_of_unity);
    fs->max_width = 0;
}

/**
 * Free the memory that was previously allocated by #new_kzg_settings.
 *
 * @param ks The settings to be freed
 */
static void free_kzg_settings(KZGSettings *ks) {
    free((FFTSettings*)ks->fs);
    free(ks->g1_values);
    free(ks->g2_values);
}

/**
 * Load trusted setup into a KZGSettings.
 *
 * @remark Free after use with #free_trusted_setup.
 *
 * @param[out] out Pointer to the stored trusted setup data
 * @param g1_bytes Array of G1 elements
 * @param n1       Length of `g1`
 * @param g2_bytes Array of G2 elements
 * @param n2       Length of `g2`
 * @retval C_CZK_OK      All is well
 * @retval C_CZK_BADARGS Invalid parameters were supplied
 * @retval C_CZK_ERROR   An internal error occurred
 * @retval C_CZK_MALLOC  Memory allocation failed
 */
C_KZG_RET load_trusted_setup(KZGSettings *out, const uint8_t *g1_bytes, size_t n1, const uint8_t *g2_bytes, size_t n2) {
    uint64_t i;
    blst_p2_affine g2_affine;
    g1_t *g1_projective = NULL;
    C_KZG_RET ret;

    out->fs = NULL;
    out->g1_values = NULL;
    out->g2_values = NULL;

    ret = new_g1_array(&out->g1_values, n1);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_g2_array(&out->g2_values, n2);
    if (ret != C_KZG_OK) goto out_error;
    ret = new_g1_array(&g1_projective, n1);
    if (ret != C_KZG_OK) goto out_error;

    for (i = 0; i < n1; i++) {
        ret = validate_kzg_g1(&g1_projective[i], (Bytes48 *)&g1_bytes[48 * i]);
        if (ret != C_KZG_OK) goto out_error;
    }

    for (i = 0; i < n2; i++) {
        blst_p2_uncompress(&g2_affine, &g2_bytes[96 * i]);
        blst_p2_from_affine(&out->g2_values[i], &g2_affine);
    }

    unsigned int max_scale = 0;
    while (((uint64_t)1 << max_scale) < n1) max_scale++;

    ret = c_kzg_malloc((void**)&out->fs, sizeof(FFTSettings));
    if (ret != C_KZG_OK) goto out_error;
    ret = new_fft_settings((FFTSettings*)out->fs, max_scale);
    if (ret != C_KZG_OK) goto out_error;
    ret = fft_g1(out->g1_values, g1_projective, true, n1, out->fs);
    if (ret != C_KZG_OK) goto out_error;
    ret = bit_reversal_permutation(out->g1_values, sizeof(g1_t), n1);
    if (ret != C_KZG_OK) goto out_error;

    goto out_success;

out_error:
    free((void *)out->fs);
    free(out->g1_values);
    free(out->g2_values);
out_success:
    free(g1_projective);
    return ret;
}

/*
 * Load trusted setup from a file.
 *
 * @remark The file format is n1 n2 g1_1 g1_2 ... g1_n1 g2_1 ... g2_n2
 * @remark where the first two numbers are in decimal and the remainder
 * @remark are hexstrings and any whitespace can be used as separators.
 * @remark See also #load_trusted_setup.
 *
 * @param[out] out Pointer to the loaded trusted setup data
 * @param[in]  in  File handle for input - will not be closed
 */
C_KZG_RET load_trusted_setup_file(KZGSettings *out, FILE *in) {
    uint64_t i;
    int num_matches;

    num_matches = fscanf(in, "%" SCNu64, &i);
    CHECK(num_matches == 1);
    CHECK(i == FIELD_ELEMENTS_PER_BLOB);
    num_matches = fscanf(in, "%" SCNu64, &i);
    CHECK(num_matches == 1);
    CHECK(i == 65);

    uint8_t g1_bytes[FIELD_ELEMENTS_PER_BLOB * 48];
    uint8_t g2_bytes[65 * 96];

    for (i = 0; i < FIELD_ELEMENTS_PER_BLOB * 48; i++) {
        num_matches = fscanf(in, "%2hhx", &g1_bytes[i]);
        CHECK(num_matches == 1);
    }

    for (i = 0; i < 65 * 96; i++) {
        num_matches = fscanf(in, "%2hhx", &g2_bytes[i]);
        CHECK(num_matches == 1);
    }

    return load_trusted_setup(out, g1_bytes, FIELD_ELEMENTS_PER_BLOB, g2_bytes, 65);
}

/*
 * Free a trusted setup (KZGSettings).
 */
void free_trusted_setup(KZGSettings *s) {
    free_fft_settings((FFTSettings*)s->fs);
    free_kzg_settings(s);
}
