#ifndef TRUSTED_SETUP_H
#define TRUSTED_SETUP_H

// Allow a library built from this code to be used from C++
#ifdef __cplusplus
extern "C" {
#endif


C_KZG_RET load_trusted_setup(KZGSettings *out,
                             const uint8_t *g1_bytes, /* n1 * 48 bytes */
                             size_t n1,
                             const uint8_t *g2_bytes, /* n2 * 96 bytes */
                             size_t n2);

C_KZG_RET load_trusted_setup_file(KZGSettings *out,
                                  FILE *in);

void free_trusted_setup(
    KZGSettings *s);

#ifdef __cplusplus
}
#endif

#endif // C_KZG_4844_H
