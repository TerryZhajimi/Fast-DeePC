# Fast-DeePC

Fast-DeePC moves the data-length-dependent work in Data-Enabled Predictive
Control (DeePC) offline. A Woodbury/K-form reduction eliminates the behavioral
coordinate `g`, leaving a fixed-dimensional box-constrained quadratic program
in the future input and output trajectory. The online problem is solved by a
dense C predictor-corrector interior-point method.

## Benchmarks

Mean online solve times from the final paper experiments are shown below.
Times are in seconds and depend on the machine and solver versions.

| Benchmark | `T` | Original baseline | Full soft-`L2` / OSQP | Reduced / OSQP | Fast-DeePC / C |
| --- | ---: | ---: | ---: | ---: | ---: |
| Two-tank | 1500 | 1.3506 | 0.2303 | 0.2372 | **0.0121** |
| Rocket | 225 | 0.0227 | 0.0058 | **0.00058** | 0.00068 |
| Mixed traffic | 1200 | 22.0861 | 4.7027 | 0.2203 | **0.0330** |

The original controllers are control-performance references and do not solve
the same optimization problem as the soft quadratic methods. Full soft OSQP,
reduced OSQP, and Fast-DeePC/C solve the identical soft quadratic problem and
separate the benefits of eliminating `g` from those of the specialized solver.

The repository contains three executable notebooks:

- `examples/two_tank.ipynb` compares the original hard-`L1` controller with
  Fast-DeePC and reports the full and reduced soft-QP baselines.
- `examples/rocket_video.ipynb` reproduces the fixed DeePC-Hunt landing case
  from initial condition `(0.65, 0.92, 0.05)`, compares the identical full and
  reduced soft problems, and records both landing videos.
- `examples/mixed_traffic.ipynb` runs the eight-vehicle DeeP-LCC benchmark.
  Its AEB safety filter is applied only to the original hard-`L1` baseline;
  Fast-DeePC uses its optimized input without AEB modification.

## Install and run

Install Python 3.10 or newer, Git, and a Windows GCC toolchain such as MSYS2
UCRT64, then run:

```powershell
python -m pip install -r requirements.txt
pwsh -ExecutionPolicy Bypass -File .\scripts\build_windows.ps1
jupyter notebook
```

The two-tank hard-`L1` comparison additionally requires PyDeePC. Set
`PYDEEPC_REPO` when that repository is not in the notebook's default search
locations. The rocket notebook requires DeePC-Hunt and its rocket environment;
set `DEEPC_HUNT_DIR` when it is stored elsewhere.

The solver source is in `src/`. Generated DLLs are written to the ignored
`build/` directory.

## Media

- `media/mixed_traffic.gif`
- `media/mixed_traffic_cav_spacing.png`
- `media/mixed_traffic_timing_vs_T_measured.png`
- `media/rocket_original_deepc_hunt.mp4`
- `media/rocket_fast_deepc.mp4`

MIT License.
