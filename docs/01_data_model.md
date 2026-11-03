# Data model and trajectory conventions

## Required LAMMPS dump columns

Each custom dump frame must contain:

```text
ITEM: TIMESTEP
<integer step>
ITEM: NUMBER OF ATOMS
<integer atom count>
ITEM: BOX BOUNDS ...
...
ITEM: ATOMS mol type x y z
...
```

Atom type 1 is oxygen and type 2 is hydrogen.  Every molecule identifier must
map to exactly one oxygen and two hydrogens.  Frames violating this topology
are rejected, not repaired.

## Molecular reconstruction

For each water molecule,

\[
\mathbf R_{\mathrm{COM}} =
\frac{m_O\mathbf R_O+m_H(\mathbf R_{H_1}+\mathbf R_{H_2})}
{m_O+2m_H},
\]

with `m_O = 15.999` and `m_H = 1.008`.  NPBC radial coordinates and shell
assignments use this molecular COM.  O--O RDF pair distances use oxygen
coordinates.

## R30 trajectory assembly

R30 contains two dump segments.  A 550000-step offset is applied to the second
segment before sorting and duplicate removal.  The global step is the only
time coordinate used after assembly.

## Sampling windows

| Analysis | R15 | R20 | R30 | PBC40 |
|---|---:|---:|---:|---:|
| Whole radius | 100--1100 ps, inclusive | 0--1000 ps, inclusive | 100--1100 ps, inclusive | not used |
| RDF extension | 41 frames | 41 frames | 41 frames | 33 frames, 200--1000 ps |
| Radial layers | `[100,1100)` ps | `[0,1000)` ps | `[100,1100)` ps | `[200,1000)` ps |
| Radial-layer frames | 200 | 200 | 200 | 160 |

The half-open radial-layer windows create exactly equal 100 ps blocks.  The
endpoint-inclusive whole-radius windows intentionally retain 201 frames.

## Input provenance

Every production kernel writes a SHA-256 inventory of the trajectories and
other inputs it used.  Treat those inventories as part of the scientific
result.  A file with the same basename but a different checksum is a different
input.
