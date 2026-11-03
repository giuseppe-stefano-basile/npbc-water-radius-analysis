#!/usr/bin/env python3
"""Write a tiny valid LAMMPS water dump for parser demonstrations only.

The geometry is deliberately artificial and must not be used as scientific
reference data.
"""

from __future__ import annotations

import argparse
from pathlib import Path

WATERS = (
    ((0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (-0.24, 0.93, 0.0)),
    ((3.0, 0.0, 0.0), (3.96, 0.0, 0.0), (2.76, 0.93, 0.0)),
    ((0.0, 3.0, 0.0), (0.96, 3.0, 0.0), (-0.24, 3.93, 0.0)),
    ((0.0, 0.0, 3.0), (0.96, 0.0, 3.0), (-0.24, 0.93, 3.0)),
)


def write_frame(handle, step: int, shift: float) -> None:
    handle.write("ITEM: TIMESTEP\n")
    handle.write(f"{step}\n")
    handle.write("ITEM: NUMBER OF ATOMS\n12\n")
    handle.write("ITEM: BOX BOUNDS pp pp pp\n-10 10\n-10 10\n-10 10\n")
    handle.write("ITEM: ATOMS id mol type x y z\n")
    atom_id = 1
    for molecule, sites in enumerate(WATERS, start=1):
        for site_index, xyz in enumerate(sites):
            atom_type = 1 if site_index == 0 else 2
            x, y, z = xyz
            handle.write(
                f"{atom_id} {molecule} {atom_type} "
                f"{x + shift:.6f} {y:.6f} {z:.6f}\n"
            )
            atom_id += 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, nargs="?", default=Path("toy_water.dump"))
    args = parser.parse_args()
    with args.output.open("w") as handle:
        write_frame(handle, 0, 0.0)
        write_frame(handle, 5000, 0.05)
    print(f"Wrote parser demonstration to {args.output}")


if __name__ == "__main__":
    main()
