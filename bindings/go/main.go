package ckzg4844

// #cgo CFLAGS: -I${SRCDIR}/../../src
// #cgo CFLAGS: -I${SRCDIR}/blst_headers
// #include "c_kzg_4844.c"
import "C"

import (
	"bytes"
	"encoding/hex"
	"errors"
	"fmt"
	"unsafe"

	// So its functions are available during compilation.
	_ "github.com/supranational/blst/bindings/go"
)

const (
	BytesPerBlob       = C.BYTES_PER_BLOB
	BytesPerCommitment = C.BYTES_PER_COMMITMENT
	BytesPerProof      = C.BYTES_PER_PROOF
	BytesPerCell       = C.BYTES_PER_CELL

	// Only used in testing, shouldn't be exposed.
	bytesPerFieldElement = C.BYTES_PER_FIELD_ELEMENT

	// Used to define return types, but clients should use the NUMBER_OF_COLUMNS
	// constant from the consensus specs. These will be the same.
	cellsPerExtBlob = C.CELLS_PER_EXT_BLOB
)

type (
	Bytes32       [32]byte
	Bytes48       [48]byte
	KZGCommitment Bytes48
	KZGProof      Bytes48
	Blob          [BytesPerBlob]byte
	Cell          [BytesPerCell]byte
)

var (
	loaded     = false
	settings   = C.KZGSettings{}
	ErrBadArgs = errors.New("bad arguments")
	ErrError   = errors.New("unexpected error")
	ErrMalloc  = errors.New("malloc failed")
)

///////////////////////////////////////////////////////////////////////////////
// Helper Functions
///////////////////////////////////////////////////////////////////////////////

// makeErrorFromRet translates an (integral) return value, as reported
// by the C library, into a proper Go error. This function should only be
// called when there is an error, not with C_KZG_OK.
func makeErrorFromRet(ret C.C_KZG_RET) error {
	switch ret {
	case C.C_KZG_BADARGS:
		return ErrBadArgs
	case C.C_KZG_ERROR:
		return ErrError
	case C.C_KZG_MALLOC:
		return ErrMalloc
	}
	return fmt.Errorf("unexpected error from c-library: %v", ret)
}

///////////////////////////////////////////////////////////////////////////////
// Unmarshal Functions
///////////////////////////////////////////////////////////////////////////////

func (b *Bytes32) UnmarshalText(input []byte) error {
	if bytes.HasPrefix(input, []byte("0x")) {
		input = input[2:]
	}
	if len(input) != 2*len(b) {
		return ErrBadArgs
	}
	l, err := hex.Decode(b[:], input)
	if err != nil {
		return err
	}
	if l != len(b) {
		return ErrBadArgs
	}
	return nil
}

func (b *Bytes48) UnmarshalText(input []byte) error {
	if bytes.HasPrefix(input, []byte("0x")) {
		input = input[2:]
	}
	if len(input) != 2*len(b) {
		return ErrBadArgs
	}
	l, err := hex.Decode(b[:], input)
	if err != nil {
		return err
	}
	if l != len(b) {
		return ErrBadArgs
	}
	return nil
}

func (b *Blob) UnmarshalText(input []byte) error {
	if bytes.HasPrefix(input, []byte("0x")) {
		input = input[2:]
	}
	if len(input) != 2*len(b) {
		return ErrBadArgs
	}
	l, err := hex.Decode(b[:], input)
	if err != nil {
		return err
	}
	if l != len(b) {
		return ErrBadArgs
	}
	return nil
}

func (c *Cell) UnmarshalText(input []byte) error {
	if bytes.HasPrefix(input, []byte("0x")) {
		input = input[2:]
	}
	if len(input) != 2*len(c) {
		return ErrBadArgs
	}
	l, err := hex.Decode(c[:], input)
	if err != nil {
		return err
	}
	if l != len(c) {
		return ErrBadArgs
	}
	return nil
}

///////////////////////////////////////////////////////////////////////////////
// Interface Functions
///////////////////////////////////////////////////////////////////////////////

/*
LoadTrustedSetup is the binding for:

	C_KZG_RET load_trusted_setup(
	    KZGSettings *out,
	    const uint8_t *g1_monomial_bytes,
	    const uint8_t *g1_lagrange_bytes,
	    size_t num_g1_points,
	    const uint8_t *g2_monomial_bytes,
	    size_t num_g2_points,
	    size_t precompute);
*/
func LoadTrustedSetup(g1MonomialBytes, g1LagrangeBytes, g2MonomialBytes []byte, precompute uint) error {
	if loaded {
		panic("trusted setup is already loaded")
	}
	if len(g1MonomialBytes)%C.BYTES_PER_G1 != 0 {
		panic(fmt.Sprintf("len(g1MonomialBytes) (%v) is not a multiple of %v", len(g1MonomialBytes), C.BYTES_PER_G1))
	}
	if len(g1LagrangeBytes)%C.BYTES_PER_G1 != 0 {
		panic(fmt.Sprintf("len(g1LagrangeBytes) (%v) is not a multiple of %v", len(g1LagrangeBytes), C.BYTES_PER_G1))
	}
	if len(g2MonomialBytes)%C.BYTES_PER_G2 != 0 {
		panic(fmt.Sprintf("len(g2MonomialBytes) (%v) is not a multiple of %v", len(g2MonomialBytes), C.BYTES_PER_G2))
	}
	numG1Monomial := len(g1MonomialBytes) / C.BYTES_PER_G1
	numG1Lagrange := len(g1LagrangeBytes) / C.BYTES_PER_G1
	numG2Monomial := len(g2MonomialBytes) / C.BYTES_PER_G2
	if numG1Monomial != numG1Lagrange {
		panic(fmt.Sprintf("numG1Monomial (%v) != numG1Lagrange (%v)", numG1Monomial, numG1Lagrange))
	}
	ret := C.load_trusted_setup(
		&settings,
		*(**C.uint8_t)(unsafe.Pointer(&g1MonomialBytes)),
		*(**C.uint8_t)(unsafe.Pointer(&g1LagrangeBytes)),
		(C.size_t)(numG1Monomial),
		*(**C.uint8_t)(unsafe.Pointer(&g2MonomialBytes)),
		(C.size_t)(numG2Monomial),
		(C.size_t)(precompute))
	if ret == C.C_KZG_OK {
		loaded = true
		return nil
	}
	return makeErrorFromRet(ret)
}

/*
LoadTrustedSetupFile is the binding for:

	C_KZG_RET load_trusted_setup_file(
	    KZGSettings *out,
	    FILE *in,
	    size_t precompute);
*/
func LoadTrustedSetupFile(trustedSetupFile string, precompute uint) error {
	if loaded {
		panic("trusted setup is already loaded")
	}
	cTrustedSetupFile := C.CString(trustedSetupFile)
	defer C.free(unsafe.Pointer(cTrustedSetupFile))
	cMode := C.CString("r")
	defer C.free(unsafe.Pointer(cMode))
	fp := C.fopen(cTrustedSetupFile, cMode)
	if fp == nil {
		panic("error reading trusted setup")
	}
	ret := C.load_trusted_setup_file(&settings, fp, (C.size_t)(precompute))
	C.fclose(fp)
	if ret == C.C_KZG_OK {
		loaded = true
		return nil
	}
	return makeErrorFromRet(ret)
}

/*
FreeTrustedSetup is the binding for:

	void free_trusted_setup(
	    KZGSettings *s);
*/
func FreeTrustedSetup() {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	C.free_trusted_setup(&settings)
	loaded = false
}

/*
BlobToKZGCommitment is the binding for:

	C_KZG_RET blob_to_kzg_commitment(
	    KZGCommitment *out,
	    const Blob *blob,
	    const KZGSettings *s);
*/
func BlobToKZGCommitment(blob *Blob) (KZGCommitment, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if blob == nil {
		return KZGCommitment{}, ErrBadArgs
	}

	var commitment KZGCommitment
	ret := C.blob_to_kzg_commitment(
		(*C.KZGCommitment)(unsafe.Pointer(&commitment)),
		(*C.Blob)(unsafe.Pointer(blob)),
		&settings)

	if ret != C.C_KZG_OK {
		return KZGCommitment{}, makeErrorFromRet(ret)
	}
	return commitment, nil
}

/*
ComputeKZGProof is the binding for:

	C_KZG_RET compute_kzg_proof(
	    KZGProof *proof_out,
	    Bytes32 *y_out,
	    const Blob *blob,
	    const Bytes32 *z_bytes,
	    const KZGSettings *s);
*/
func ComputeKZGProof(blob *Blob, zBytes Bytes32) (KZGProof, Bytes32, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if blob == nil {
		return KZGProof{}, Bytes32{}, ErrBadArgs
	}

	var proof, y = KZGProof{}, Bytes32{}
	ret := C.compute_kzg_proof(
		(*C.KZGProof)(unsafe.Pointer(&proof)),
		(*C.Bytes32)(unsafe.Pointer(&y)),
		(*C.Blob)(unsafe.Pointer(blob)),
		(*C.Bytes32)(unsafe.Pointer(&zBytes)),
		&settings)

	if ret != C.C_KZG_OK {
		return KZGProof{}, Bytes32{}, makeErrorFromRet(ret)
	}
	return proof, y, nil
}

/*
ComputeBlobKZGProof is the binding for:

	C_KZG_RET compute_blob_kzg_proof(
	    KZGProof *out,
	    const Blob *blob,
	    const Bytes48 *commitment_bytes,
	    const KZGSettings *s);
*/
func ComputeBlobKZGProof(blob *Blob, commitmentBytes Bytes48) (KZGProof, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if blob == nil {
		return KZGProof{}, ErrBadArgs
	}
	var proof KZGProof
	ret := C.compute_blob_kzg_proof(
		(*C.KZGProof)(unsafe.Pointer(&proof)),
		(*C.Blob)(unsafe.Pointer(blob)),
		(*C.Bytes48)(unsafe.Pointer(&commitmentBytes)),
		&settings)

	if ret != C.C_KZG_OK {
		return KZGProof{}, makeErrorFromRet(ret)
	}
	return proof, nil
}

/*
VerifyKZGProof is the binding for:

	C_KZG_RET verify_kzg_proof(
	    bool *out,
	    const Bytes48 *commitment_bytes,
	    const Bytes32 *z_bytes,
	    const Bytes32 *y_bytes,
	    const Bytes48 *proof_bytes,
	    const KZGSettings *s);
*/
func VerifyKZGProof(commitmentBytes Bytes48, zBytes, yBytes Bytes32, proofBytes Bytes48) (bool, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	var result C.bool
	ret := C.verify_kzg_proof(
		&result,
		(*C.Bytes48)(unsafe.Pointer(&commitmentBytes)),
		(*C.Bytes32)(unsafe.Pointer(&zBytes)),
		(*C.Bytes32)(unsafe.Pointer(&yBytes)),
		(*C.Bytes48)(unsafe.Pointer(&proofBytes)),
		&settings)

	if ret != C.C_KZG_OK {
		return false, makeErrorFromRet(ret)
	}
	return bool(result), nil
}

/*
VerifyBlobKZGProof is the binding for:

	C_KZG_RET verify_blob_kzg_proof(
	    bool *out,
	    const Blob *blob,
	    const Bytes48 *commitment_bytes,
	    const Bytes48 *proof_bytes,
	    const KZGSettings *s);
*/
func VerifyBlobKZGProof(blob *Blob, commitmentBytes, proofBytes Bytes48) (bool, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if blob == nil {
		return false, ErrBadArgs
	}

	var result C.bool
	ret := C.verify_blob_kzg_proof(
		&result,
		(*C.Blob)(unsafe.Pointer(blob)),
		(*C.Bytes48)(unsafe.Pointer(&commitmentBytes)),
		(*C.Bytes48)(unsafe.Pointer(&proofBytes)),
		&settings)

	if ret != C.C_KZG_OK {
		return false, makeErrorFromRet(ret)
	}
	return bool(result), nil
}

/*
VerifyBlobKZGProofBatch is the binding for:

	C_KZG_RET verify_blob_kzg_proof_batch(
	    bool *out,
	    const Blob *blobs,
	    const Bytes48 *commitments_bytes,
	    const Bytes48 *proofs_bytes,
	    const KZGSettings *s);
*/
func VerifyBlobKZGProofBatch(blobs []Blob, commitmentsBytes, proofsBytes []Bytes48) (bool, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if len(blobs) != len(commitmentsBytes) || len(blobs) != len(proofsBytes) {
		return false, ErrBadArgs
	}

	var result C.bool
	ret := C.verify_blob_kzg_proof_batch(
		&result,
		*(**C.Blob)(unsafe.Pointer(&blobs)),
		*(**C.Bytes48)(unsafe.Pointer(&commitmentsBytes)),
		*(**C.Bytes48)(unsafe.Pointer(&proofsBytes)),
		(C.size_t)(len(blobs)),
		&settings)

	if ret != C.C_KZG_OK {
		return false, makeErrorFromRet(ret)
	}
	return bool(result), nil
}

/*
ComputeCellsAndKZGProofs is the binding for:

	C_KZG_RET compute_cells_and_kzg_proofs(
	    Cell *cells,
	    KZGProof *proofs,
	    const Blob *blob,
	    const KZGSettings *s);
*/
func ComputeCellsAndKZGProofs(blob *Blob) ([cellsPerExtBlob]Cell, [cellsPerExtBlob]KZGProof, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}

	cells := [cellsPerExtBlob]Cell{}
	proofs := [cellsPerExtBlob]KZGProof{}
	ret := C.compute_cells_and_kzg_proofs(
		(*C.Cell)(unsafe.Pointer(&cells)),
		(*C.KZGProof)(unsafe.Pointer(&proofs)),
		(*C.Blob)(unsafe.Pointer(blob)),
		&settings)

	if ret != C.C_KZG_OK {
		return [cellsPerExtBlob]Cell{}, [cellsPerExtBlob]KZGProof{}, makeErrorFromRet(ret)
	}
	return cells, proofs, nil
}

/*
RecoverCellsAndKZGProofs is the binding for:

	C_KZG_RET recover_cells_and_kzg_proofs(
	    Cell *recovered_cells,
	    KZGProof *recovered_proofs,
	    const uint64_t *cell_indices,
	    const Cell *cells,
	    size_t num_cells,
	    const KZGSettings *s);
*/
func RecoverCellsAndKZGProofs(cellIndices []uint64, cells []Cell) ([cellsPerExtBlob]Cell, [cellsPerExtBlob]KZGProof, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	if len(cellIndices) != len(cells) {
		return [cellsPerExtBlob]Cell{}, [cellsPerExtBlob]KZGProof{}, ErrBadArgs
	}

	recoveredCells := [cellsPerExtBlob]Cell{}
	recoveredProofs := [cellsPerExtBlob]KZGProof{}
	ret := C.recover_cells_and_kzg_proofs(
		(*C.Cell)(unsafe.Pointer(&recoveredCells)),
		(*C.KZGProof)(unsafe.Pointer(&recoveredProofs)),
		*(**C.uint64_t)(unsafe.Pointer(&cellIndices)),
		*(**C.Cell)(unsafe.Pointer(&cells)),
		(C.size_t)(len(cells)),
		&settings)

	if ret != C.C_KZG_OK {
		return [cellsPerExtBlob]Cell{}, [cellsPerExtBlob]KZGProof{}, makeErrorFromRet(ret)
	}
	return recoveredCells, recoveredProofs, nil
}

/*
VerifyCellKZGProof is the binding for:

	C_KZG_RET verify_cell_kzg_proof(
	    bool *ok,
	    const Bytes48 *commitment_bytes,
	    uint64_t cell_index,
	    const Cell *cell,
	    const KZGProof *proof,
	    const KZGSettings *s);
*/
func VerifyCellKZGProof(commitmentBytes Bytes48, cellIndex uint64, cell Cell, proofBytes Bytes48) (bool, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}

	var result C.bool
	ret := C.verify_cell_kzg_proof(
		&result,
		(*C.Bytes48)(unsafe.Pointer(&commitmentBytes)),
		(C.uint64_t)(cellIndex),
		(*C.Cell)(unsafe.Pointer(&cell)),
		(*C.Bytes48)(unsafe.Pointer(&proofBytes)),
		&settings)

	if ret != C.C_KZG_OK {
		return false, makeErrorFromRet(ret)
	}
	return bool(result), nil
}

/*
VerifyCellKZGProofBatch is the binding for:

	C_KZG_RET verify_cell_kzg_proof_batch(
	    bool *ok,
	    const Bytes48 *commitments_bytes,
	    const uint64_t *column_indices,
	    const Cell *cells,
	    const Bytes48 *proofs_bytes,
	    size_t num_cells,
	    const KZGSettings *s);
*/
func VerifyCellKZGProofBatch(commitmentsBytes []Bytes48, columnIndices []uint64, cells []Cell, proofsBytes []Bytes48) (bool, error) {
	if !loaded {
		panic("trusted setup isn't loaded")
	}
	cellCount := len(cells)
	if len(commitmentsBytes) != cellCount || len(columnIndices) != cellCount || len(proofsBytes) != cellCount {
		return false, ErrBadArgs
	}

	var result C.bool
	ret := C.verify_cell_kzg_proof_batch(
		&result,
		*(**C.Bytes48)(unsafe.Pointer(&commitmentsBytes)),
		*(**C.uint64_t)(unsafe.Pointer(&columnIndices)),
		*(**C.Cell)(unsafe.Pointer(&cells)),
		*(**C.Bytes48)(unsafe.Pointer(&proofsBytes)),
		(C.size_t)(len(cells)),
		&settings)

	if ret != C.C_KZG_OK {
		return false, makeErrorFromRet(ret)
	}
	return bool(result), nil
}
