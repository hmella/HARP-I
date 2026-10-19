"""Command line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def _cmd_process(args):
    import harpi

    ds = harpi.load_mat_2d(args.input)
    print(f"loaded {args.input}: {ds.image.shape}, {ds.n_frames} frames, "
          f"{ds.pixels_per_wavelength:.2f} px/wavelength")

    flt = harpi.HarpFilter2D(ds.image, ds.ke, kind=args.filter)
    Ih = flt(ds.image)

    a = args.a if args.a is not None else ds.suggested_a
    print(f"motion: kernel={args.kernel} a={a:.3f} eta={args.eta}")

    def progress(i, n, stage):
        if args.verbose:
            print(f"  {stage} {i + 1}/{n}", end="\r", file=sys.stderr)

    res = harpi.estimate_motion_2d(
        Ih, ds.mask, kernel=args.kernel, a=a, eta=args.eta,
        smooth_reference=args.smooth_reference, temporal_fit=args.temporal_fit,
        temporal_order=args.temporal_order, downsample=args.downsample,
        progress=progress, phase_jump_rule=args.phase_jump_rule,
    )
    if not res.unwrap_complete.all():
        bad = np.flatnonzero(~res.unwrap_complete)
        print(f"WARNING: unwrapping incomplete on frame(s) {bad.tolist()}", file=sys.stderr)

    st = harpi.pixel_strain_2d(res.best, ds.mask[..., 0])
    gcs = st.global_curve("CC")
    grs = st.global_curve("RR")
    print("\nglobal circumferential strain:", np.array2string(gcs, precision=4))
    print("global radial strain:        ", np.array2string(grs, precision=4))

    if args.out:
        harpi.save_mat(args.out, motion=res, strain=st)
        print(f"wrote {args.out}")
    return 0


def _cmd_info(args):
    import harpi

    ds = harpi.load_mat_2d(args.input)
    print(f"image           {ds.image.shape}")
    print(f"mask            {ds.mask.shape} ({ds.mask[..., 0].sum()} px in frame 0)")
    print(f"pixel size      {ds.pixel_size}")
    print(f"tag spacing     {ds.tag_spacing}")
    print(f"px/wavelength   {ds.pixels_per_wavelength:.3f}")
    print(f"ke (rad/px)     {ds.ke}")
    print(f"suggested a     {ds.suggested_a:.3f}")
    return 0


def _cmd_gui(args):
    try:
        from .gui.app import main as gui_main
    except ImportError as exc:  # pragma: no cover
        print(f"the GUI needs PySide6: pip install 'harpi[gui]'  ({exc})", file=sys.stderr)
        return 1
    return gui_main(args.project)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="harpi", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("process", help="run the 2-D pipeline on a .mat dataset")
    pr.add_argument("input", type=Path)
    pr.add_argument("--out", type=Path, help="write results as a MATLAB .mat")
    pr.add_argument("--kernel", default="multiquadric3")
    pr.add_argument("--a", type=float, default=None,
                    help="RBF shape parameter (default: 150 / px-per-wavelength)")
    pr.add_argument("--eta", type=float, default=15.0)
    pr.add_argument("--filter", default="transmission",
                    choices=["transmission", "butterworth"])
    pr.add_argument("--downsample", type=int, default=1)
    pr.add_argument("--smooth-reference", action="store_true")
    pr.add_argument("--temporal-fit", action="store_true")
    pr.add_argument("--temporal-order", type=int, default=10)
    pr.add_argument("--phase-jump-rule", default="published",
                    choices=["published", "nearest"],
                    help="rule for the 2*pi offset of each frame")
    pr.add_argument("-v", "--verbose", action="store_true")
    pr.set_defaults(func=_cmd_process)

    pi = sub.add_parser("info", help="summarise a dataset")
    pi.add_argument("input", type=Path)
    pi.set_defaults(func=_cmd_info)

    pg = sub.add_parser("gui", help="launch the graphical interface")
    pg.add_argument("project", nargs="?", default=None)
    pg.set_defaults(func=_cmd_gui)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
