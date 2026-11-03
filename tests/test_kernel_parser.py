from pathlib import Path

import numpy as np

from npbc_water_analysis.kernels import load_kernel


def test_dump_parser_and_strict_water_reconstruction(tmp_path: Path):
    dump = tmp_path / "water.dump"
    dump.write_text(
        "ITEM: TIMESTEP\n0\n"
        "ITEM: NUMBER OF ATOMS\n6\n"
        "ITEM: BOX BOUNDS pp pp pp\n-5 5\n-5 5\n-5 5\n"
        "ITEM: ATOMS id mol type x y z\n"
        "1 1 1 0 0 0\n2 1 2 0.96 0 0\n3 1 2 -0.24 0.93 0\n"
        "4 2 1 3 0 0\n5 2 2 3.96 0 0\n6 2 2 2.76 0.93 0\n"
    )
    kernel = load_kernel("whole")
    frames = list(kernel.dump_frames(dump))
    assert len(frames) == 1
    step, mol, typ, xyz = frames[0]
    assert step == 0
    oxygen, h1, h2 = kernel.water_arrays(mol, typ, xyz)
    assert oxygen.shape == h1.shape == h2.shape == (2, 3)
    np.testing.assert_allclose(oxygen[:, 0], [0.0, 3.0])
