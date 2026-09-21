# Fast-DeePC

Dense C BoxQP implementation of Fast Data-Enabled Predictive Control. The
offline phase eliminates the Hankel coefficient `g`; online, the solver uses
only fixed-size maps and optimizes the bounded input/output trajectory.

## Benchmarks

| Benchmark | `T` | `L1` DeePC (s) | `L2` OSQP (s) | Fast-DeePC (s) |
| --- | ---: | ---: | ---: | ---: |
| Two-tank | 1500 | 1.7225 | 0.0430 | 0.005237 |
| Mixed traffic | 1200 | 28.5917 | 2.3059 | 0.04774 |
| Rocket | 225 | 0.02914 | 0.00656 | `8e-4` |

These are reference measurements; fresh timings depend on the machine and
solver versions. Each notebook runs one fixed case:

- `examples/two_tank.ipynb` prints the three timings and plots both tank levels.
- `examples/mixed_traffic.ipynb` retains the shared AEB safety filter, uses a
  genuinely `T`-independent Fast-DeePC online step, and writes the traffic GIF.
- `examples/rocket_video.ipynb` uses `T=225`, seed `42`, and initial condition
  `(0.65, 0.92, 0.05)`, then records both landing videos.

## Install and run

Install Python 3.10 or 3.11, Git, and MSYS2 UCRT64 GCC, then run:

```powershell
git clone https://github.com/TerryZhajimi/Fast-DeePC.git
cd Fast-DeePC
python -m pip install -r requirements.txt
pwsh -ExecutionPolicy Bypass -File .\scripts\build_windows.ps1
jupyter notebook
```

The two-tank and mixed-traffic notebooks need no external controller repository.
The mixed-traffic notebook builds its benchmark-specific DLL automatically.

The rocket notebook additionally needs DeePC-Hunt and the rocket environment:

```powershell
git clone https://github.com/michael-cummins/DeePC-Hunt.git external/DeePC-Hunt-main
python -m pip install -e .\external\DeePC-Hunt-main
python -m pip install "gymnasium==0.27.1" "pygame==2.1.3.dev8" "moviepy<2"
python -m pip install --only-binary=:all: "Box2D==2.3.10"
python -m pip install --no-deps "git+https://git.ethz.ch/bsaverio/coco-project.git"
```

Set `DEEPC_HUNT_DIR` only when DeePC-Hunt is stored somewhere else. The solver
source is in `src/`; generated DLLs remain in the ignored `build/` directory.

## Media

- `media/mixed_traffic.gif`
- `media/rocket_original_deepc_hunt.mp4`
- `media/rocket_fast_deepc.mp4`

MIT License.
