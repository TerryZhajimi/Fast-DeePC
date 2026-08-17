# Fast-DeePC

Dense C BoxQP implementation of Fast Data-Enabled Predictive Control. The
offline phase eliminates the Hankel coefficient `g`; online, the solver only
optimizes the bounded input/output trajectory.

## Three benchmark cases

| Benchmark | `T` | `L1` DeePC (s) | `L2` OSQP (s) | Fast-DeePC (s) |
| --- | ---: | ---: | ---: | ---: |
| Two-tank | 1500 | 1.7225 | 0.0430 | 0.005237 |
| Mixed traffic | 1200 | 28.5917 | 2.3059 | 0.04774 |
| Rocket | 225 | 0.02914 | 0.00656 | `8e-4` |

Each notebook runs one fixed case only:

- `examples/two_tank.ipynb`: prints the three timings and draws one figure of
  the two tank-level outputs for original `L1` DeePC and Fast-DeePC.
- `examples/mixed_traffic.ipynb`: runs the one emergency-braking case and
  writes `media/mixed_traffic.gif`.
- `examples/rocket_video.ipynb`: uses the fixed initial condition
  `(0.6, 0.9, 0.1)` and records the original and Fast-DeePC landing videos.

## Run

```powershell
python -m pip install -r requirements.txt
pwsh -ExecutionPolicy Bypass -File .\scripts\build_windows.ps1
jupyter notebook
```

The two-tank hard-`L1` comparison additionally needs PyDeePC. Set
`PYDEEPC_REPO` to its local directory. The rocket notebook additionally needs
DeePC-Hunt and its rocket environment; set `DEEPC_HUNT_DIR` to that repository.

The C source is in `src/`. The portable wrapper lets the examples compile
without a separate BLAS/LAPACK installation.

## Media

- `media/mixed_traffic.gif`
- `media/rocket_original_deepc_hunt.mp4`
- `media/rocket_fast_deepc.mp4`

MIT License.
