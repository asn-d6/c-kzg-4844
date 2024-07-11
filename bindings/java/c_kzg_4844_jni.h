/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class ethereum_ckzg4844_CKZG4844JNI */

#ifndef _Included_ethereum_ckzg4844_CKZG4844JNI
#define _Included_ethereum_ckzg4844_CKZG4844JNI
#ifdef __cplusplus
extern "C" {
#endif
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_G1
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_G1 48L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_G2
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_G2 96L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_FIELD_ELEMENT
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_FIELD_ELEMENT 32L
#undef ethereum_ckzg4844_CKZG4844JNI_BITS_PER_FIELD_ELEMENT
#define ethereum_ckzg4844_CKZG4844JNI_BITS_PER_FIELD_ELEMENT 255L
#undef ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_BLOB
#define ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_BLOB 4096L
#undef ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_EXT_BLOB
#define ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_EXT_BLOB 8192L
#undef ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_CELL
#define ethereum_ckzg4844_CKZG4844JNI_FIELD_ELEMENTS_PER_CELL 64L
#undef ethereum_ckzg4844_CKZG4844JNI_CELLS_PER_EXT_BLOB
#define ethereum_ckzg4844_CKZG4844JNI_CELLS_PER_EXT_BLOB 128L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_COMMITMENT
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_COMMITMENT 48L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_PROOF
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_PROOF 48L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_BLOB
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_BLOB 131072L
#undef ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_CELL
#define ethereum_ckzg4844_CKZG4844JNI_BYTES_PER_CELL 2048L
/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    loadTrustedSetup
 * Signature: (Ljava/lang/String;J)V
 */
JNIEXPORT void JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_loadTrustedSetup__Ljava_lang_String_2J
  (JNIEnv *, jclass, jstring, jlong);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    loadTrustedSetup
 * Signature: ([B[B[BJ)V
 */
JNIEXPORT void JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_loadTrustedSetup___3B_3B_3BJ
  (JNIEnv *, jclass, jbyteArray, jbyteArray, jbyteArray, jlong);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    freeTrustedSetup
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_freeTrustedSetup
  (JNIEnv *, jclass);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    blobToKzgCommitment
 * Signature: ([B)[B
 */
JNIEXPORT jbyteArray JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_blobToKzgCommitment
  (JNIEnv *, jclass, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    computeKzgProof
 * Signature: ([B[B)Lethereum/ckzg4844/ProofAndY;
 */
JNIEXPORT jobject JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_computeKzgProof
  (JNIEnv *, jclass, jbyteArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    computeBlobKzgProof
 * Signature: ([B[B)[B
 */
JNIEXPORT jbyteArray JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_computeBlobKzgProof
  (JNIEnv *, jclass, jbyteArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    verifyKzgProof
 * Signature: ([B[B[B[B)Z
 */
JNIEXPORT jboolean JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_verifyKzgProof
  (JNIEnv *, jclass, jbyteArray, jbyteArray, jbyteArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    verifyBlobKzgProof
 * Signature: ([B[B[B)Z
 */
JNIEXPORT jboolean JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_verifyBlobKzgProof
  (JNIEnv *, jclass, jbyteArray, jbyteArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    verifyBlobKzgProofBatch
 * Signature: ([B[B[BJ)Z
 */
JNIEXPORT jboolean JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_verifyBlobKzgProofBatch
  (JNIEnv *, jclass, jbyteArray, jbyteArray, jbyteArray, jlong);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    computeCellsAndKzgProofs
 * Signature: ([B)Lethereum/ckzg4844/CellsAndProofs;
 */
JNIEXPORT jobject JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_computeCellsAndKzgProofs
  (JNIEnv *, jclass, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    recoverCellsAndKzgProofs
 * Signature: ([J[B)Lethereum/ckzg4844/CellsAndProofs;
 */
JNIEXPORT jobject JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_recoverCellsAndKzgProofs
  (JNIEnv *, jclass, jlongArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    verifyCellKzgProof
 * Signature: ([BJ[B[B)Z
 */
JNIEXPORT jboolean JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_verifyCellKzgProof
  (JNIEnv *, jclass, jbyteArray, jlong, jbyteArray, jbyteArray);

/*
 * Class:     ethereum_ckzg4844_CKZG4844JNI
 * Method:    verifyCellKzgProofBatch
 * Signature: ([B[J[B[B)Z
 */
JNIEXPORT jboolean JNICALL Java_ethereum_ckzg4844_CKZG4844JNI_verifyCellKzgProofBatch
  (JNIEnv *, jclass, jbyteArray, jlongArray, jbyteArray, jbyteArray);

#ifdef __cplusplus
}
#endif
#endif
