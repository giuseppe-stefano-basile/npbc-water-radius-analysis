#!/usr/bin/env python3

src = "traj_pbc40_prod_1ns.dump"
dst = "traj_pbc40_prod_stationary_200_1000ps.dump"

STEP_MIN = 200000
STEP_MAX = 1000000

kept = 0
seen = 0

with open(src) as fin, open(dst, "w") as fout:
    while True:
        line = fin.readline()
        if not line:
            break

        if not line.startswith("ITEM: TIMESTEP"):
            raise RuntimeError(
                f"Unexpected dump structure: {line!r}"
            )

        step_line = fin.readline()
        step = int(step_line.strip())
        seen += 1

        nat_header = fin.readline()
        nat_line = fin.readline()
        nat = int(nat_line.strip())

        box_header = fin.readline()
        box_lines = [fin.readline() for _ in range(3)]

        atoms_header = fin.readline()
        atom_lines = [fin.readline() for _ in range(nat)]

        if STEP_MIN <= step <= STEP_MAX:
            fout.write(line)
            fout.write(step_line)
            fout.write(nat_header)
            fout.write(nat_line)
            fout.write(box_header)

            for x in box_lines:
                fout.write(x)

            fout.write(atoms_header)

            for x in atom_lines:
                fout.write(x)

            kept += 1

print("frames seen =", seen)
print("frames kept =", kept)
print("output      =", dst)

if kept != 801:
    raise SystemExit(
        f"ERROR: expected 801 stationary frames, got {kept}"
    )
