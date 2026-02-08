# CHOOT v0 File Format Specification

## Overview
CHOOT stores a continuous 2D radiance field as a list of anisotropic Gaussian splat atoms. There is no pixel grid in the file. Decoders reconstruct images by evaluating the field at the output resolution you want.

## Endianness
All integer and floating-point values are little-endian.

## Byte Layout Summary (v0)
The file is a fixed-size header followed by `atom_count` packed atom records.

```
Offset  Size  Description
0       8     Magic ASCII: 43 48 4F 4F 54 00 00 ("CHOOT\0\0")
8       2     version (u16)
10      2     flags (u16)
12      4     atom_count (u32)
16      4     header_size (u32)
20      4     reserved (u32)
24      N     atom records, each 20 bytes
```

## Numeric Types
- `u16`, `u32`: unsigned integers
- `f16`: IEEE 754 binary16 (half float)
- `f32`: IEEE 754 binary32 (single float) for runtime math only (not stored in v0)

## Color Model
- Color is stored in linear-light YCoCg.
- The transform is round-trip safe:
  $$Y = 0.25R + 0.5G + 0.25B$$
  $$Co = 0.5R - 0.5B$$
  $$Cg = -0.25R + 0.5G - 0.25B$$
- Inverse:
  $$R = Y + Co - Cg$$
  $$G = Y + Cg$$
  $$B = Y - Co - Cg$$
- $R$, $G$, $B$ are linear (not gamma-encoded).

## Header (fixed size, 24 bytes)
| Offset | Size | Type | Name | Description |
|---:|---:|---|---|---|
| 0 | 8 | bytes | magic | ASCII `"CHOOT\0\0"` |
| 8 | 2 | u16 | version | Must be `0` for v0 |
| 10 | 2 | u16 | flags | Reserved, must be `0` |
| 12 | 4 | u32 | atom_count | Number of atoms |
| 16 | 4 | u32 | header_size | Size of header in bytes (must be 24) |
| 20 | 4 | u32 | reserved | Must be `0` |

Total header size: 24 bytes.

## Atom Record (v0)
Each atom is exactly 20 bytes, packed without padding or alignment gaps.
Alignment for reading: byte-aligned; do not assume any padding. Read field-by-field.

| Offset | Size | Type | Name | Description |
|---:|---:|---|---|---|
| 0 | 2 | f16 | x | Center x in normalized space [0,1] |
| 2 | 2 | f16 | y | Center y in normalized space [0,1] |
| 4 | 2 | f16 | sxx | Covariance element $\sigma_{xx}$ |
| 6 | 2 | f16 | sxy | Covariance element $\sigma_{xy}$ |
| 8 | 2 | f16 | syy | Covariance element $\sigma_{yy}$ |
| 10 | 2 | f16 | alpha | Opacity $\alpha$ in [0,1] |
| 12 | 2 | f16 | Y | Luma |
| 14 | 2 | f16 | Co | Chroma Co |
| 16 | 2 | f16 | Cg | Chroma Cg |
| 18 | 2 | u16 | flags | Optional flags (v0: must be 0) |

### Covariance matrix
The covariance matrix is:
$$
\Sigma = \begin{bmatrix}
\sigma_{xx} & \sigma_{xy}\\
\sigma_{xy} & \sigma_{yy}
\end{bmatrix}
$$
It must be symmetric and positive semi-definite. Decoders should clamp or skip atoms with invalid covariance.

## Field Evaluation (normative)
Given a sample point $p = (x, y)$ and atom center $c$, define $d = p - c$. The Gaussian weight is:
$$
G = \exp\left(-\tfrac{1}{2} d^T \Sigma^{-1} d\right)
$$
Each atom contributes a premultiplied radiance weight:
$$
w = \alpha \cdot G
$$
Accumulation is order-independent:
$$
C_{sum} = \sum w \cdot C_{atom}, \quad A_{sum} = \sum w
$$
Final linear RGB is computed by converting YCoCg to RGB and normalizing:
$$
C_{final} = \frac{C_{sum}}{\max(A_{sum}, \epsilon)}
$$
with $\epsilon = 1e{-8}$.

## Determinism
All math is performed in 32-bit floats during decoding. The file itself uses f16 for compactness. Decoders must not reorder atoms or use stochastic sampling.

## Versioning Strategy
- `version` is a 16-bit integer.
- v0 decoders must reject any `version` != 0.
- Future versions may increase `header_size` and append additional header fields after offset 24.
- v0 decoders must ignore any extra bytes between `header_size` and the first atom record.

## Validation Rules (v0)
Decoders must validate the following and fail fast on violation:
1. `magic` matches `"CHOOT\0\0"`.
2. `version` == 0.
3. `flags` == 0.
4. `header_size` == 24.
5. File size >= `header_size + atom_count * 20`.
6. Atom `flags` == 0.
7. Covariance matrix is symmetric by construction and must be positive semi-definite: $\det(\Sigma) = sxx\cdot syy - sxy^2 \ge 0$. If invalid, the atom must be skipped.
8. `x`, `y`, `alpha` are clamped to [0,1] on decode; NaNs or infinities must cause the atom to be skipped.

## Example Hex Dump (minimal valid file)
This example has `atom_count = 0` and no atom records.

```
00000000  43 48 4F 4F 54 00 00 00  00 00 00 00 00 00 00 00
00000010  18 00 00 00 00 00 00 00
```

Breakdown:
- Magic: 43 48 4F 4F 54 00 00 00
- version: 00 00
- flags: 00 00
- atom_count: 00 00 00 00
- header_size: 18 00 00 00
- reserved: 00 00 00 00

## Future-Proofing
- `flags` fields are reserved.
- `header_size` allows optional extension headers in future versions.

