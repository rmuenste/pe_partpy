Project notes for agents working in this repo

Scope: applies to the entire repository rooted here.

Overview
- Purpose: Partition quadrilateral/hexahedral meshes for parallel computing using METIS and custom geometric methods.
- Entrypoint: `PyPartitioner.py` (CLI). Core logic in `partitioner/part_main.py` and `partitioner/part.py`.
- Mesh I/O and structures: `mesh/mesh_io.py`, `mesh/mesh.py`.
- Utilities and examples: `hex_mesh_2.py` (VTK/ParaView demo), `tri2vtk_converter.py`.

Unit Cube Workflow (VTK -> TRI -> PAR case)
- Goal: unit cube [0,1]^3 with nx=ny=nz=25, then generate .par files and a final case folder.
- Steps (uses `pe_partpy/unit_cube_vtk.py`):
  1) Create legacy ASCII VTK (run from repo root):
     - `./.venv/bin/python pe_partpy/unit_cube_vtk.py --out pe_partpy/mesh/unit_cube_25.vtk --nx 25 --ny 25 --nz 25`
  2) Convert VTK -> TRI (run from `pe_partpy/`):
     - `../.venv/bin/python tri2vtk_converter.py mesh/unit_cube_25.vtk mesh/unit_cube_25.tri`
  3) Generate .par files and case folder (run from repo root):
     - `./.venv/bin/python pe_partpy/gen_par_from_tri.py pe_partpy/mesh/unit_cube_25.tri --outdir pe_partpy/mesh/unit_cube_25_case`
- Output: `pe_partpy/mesh/unit_cube_25_case` contains `file.prj`, `unit_cube_25.tri`, and `xmin/xmax/ymin/ymax/zmin/zmax.par`.

How to Run (CLI)
- Command: `python PyPartitioner.py <n_part> <part_method> <n_sub_part> <mesh_name> <project_file>`
- Examples (do NOT pass `-f` flag; see gotchas):
  - METIS recursive: `python PyPartitioner.py 15 1 1 NEWCASE ./path/to/case.prj`
  - Axis-based (x-y-z): `python PyPartitioner.py 27 -4 x3-y3-z3 NEWCASE ./path/to/case3d.prj`
  - Plane-based ring: `python PyPartitioner.py 30 -7 x15-y2-z1 NEWCASE ./path/to/case.prj`
- Output: `_mesh/<mesh_name>/sub####/GRID.tri` plus per-sub `.par` files and a copied `.prj`. Format of `####` depends on formatting (see Formatting IDs).

Partition Methods (high level)
- `1`, `2`, `3`: METIS (Recursive, VKway, Kway) on subgrids.
- Negative methods: geometry-driven splits
  - `-4`: axis-based partitioning (binary splits per axis of `x`, `y`, `z`, `4` maps to xyz). Use `n_sub_part` as `xA-yB-zC`.
  - `-5`: single-plane based partitions (needs `single_plane.txt`).
  - `-6`: dual-plane based partitions (needs `x_planes_part34.txt` and `y_planes_part34.txt`).
  - `-7`: plane ring partitioning (needs `single_plane.txt`, interpreted as ordered ring).

Inputs and Case Layout
- `.prj`: Lists one `.tri` grid file and N `.par` boundary files; paths are relative to the project file’s folder.
- `.tri` expectations: sections contain `DCORVG` (coords), `KVERT` (connectivity), `KNPR` (boundary markers) in that order.
- Plane files (when required):
  - `single_plane.txt`: lines of `px py pz nx ny nz` in case folder.
  - `x_planes_part34.txt`, `y_planes_part34.txt`: same line structure; both must exist in the case folder.

Output Structure
- Work root: `_mesh/<mesh_name>/`.
- Subgrids: created under `sub####` where `####` is zero-padded (see Formatting IDs).
- For some flows (e.g., generating subgrids first), top-level `GRID.tri` and `GRID.prj` are copied alongside the `.par` files.
- Boundary nodes at partition interfaces are promoted in `KNPR`.

Formatting IDs
- `partitioner.part.get_formatted_value()` supports `v1` (3 digits) and `v2` (4 digits). Default is effectively `v1` unless overridden.
- The CLI parser in `PyPartitioner.py` defines `-f/--format` (default `v2`) and passes `args` to `main_process` so subdir IDs render with `v2` by default.

Important Gotchas
- Do not pass `-f` on the command line: `partitioner.part_main.checkParameters(sys.argv)` expects exactly 6 argv entries and will mis-parse when `-f` is present. Use the default formatting or let `PyPartitioner.py` set it.
- Programmatic API hazards:
  - `partitioner.part_main.main_process` expects an `all_args` object with `.format` when called from `PyPartitioner.py`. The helper `partitioner.part_main.partition()` omits this argument and is likely outdated; prefer the CLI or call `main_process` directly with a lightweight object providing `.format`.
- METIS library loading:
  - `partitioner.part._try_in_place_first` searches `./libmetis.so`, then `../lib64/libmetis.so`, then system paths.
  - This repo includes `libmetis.so` at project root and in `partitioner/`. Running from repo root works out-of-the-box; ensure CWD contains a usable lib or that system METIS is available.
- Axis-based `-4` expects `n_sub_part` as `xA-yB-zC` (e.g., `x3-y3-z3`). Other negatives (`-5`, `-6`, `-7`) also require `xA-yB-zC`. `checkParameters` validates that `product(A,B,C) == n_part`.
- Atomic splitting: If grid elements equal `n_part`, METIS stage uses atomic element-per-partition splitting.

Coding Conventions
- The repo has undergone modernization (see `modernization_plan.md`). Follow:
  - Type hints where practical.
  - f-strings for formatting.
  - `snake_case` naming.
  - Keep changes minimal and aligned with existing style.

Files of Interest
- `PyPartitioner.py`: CLI parsing, calls `partitioner.mkdir`, `checkParameters(sys.argv)`, and `main_process(...)`.
- `partitioner/part_main.py`: Parameter checks, plane loading, directory orchestration, sequencing through subgrid generation and METIS partitioning.
- `partitioner/part.py`: METIS bindings and algorithms: neighborhood construction, `GetParts`, axis and plane-based partitioning, subgrid output.
- `mesh/mesh_io.py` and `mesh/mesh.py`: Grid and boundary parsing/serialization; mesh data structures and quality utilities.
- `hex_mesh_2.py`: Standalone VTK/ParaView example to generate and render a regular hexahedral grid (requires `paraview.simple`, `vtk`, `numpy`). Not used by the CLI.

Common Tasks
- Add a new partitioning strategy: implement in `partitioner/part.py` (keep tuple-of-tuples structures), wire into `part_main.main_process` switch on `orig_method`.
- Change subdir numbering: extend `get_formatted_value` mapping; ensure the CLI passes desired `-f/--format` (but remember the gotcha: do not pass `-f` to the raw argv path).
- Debug METIS calls: verify adjacency arrays `MetisA/MetisB` and one-based indexing, and ensure `libmetis.so` resolves from the current CWD.

Validation Tips
- Start with small cases and `n_part=1` (bypasses METIS) to validate I/O.
- Use axis-based methods (`-4`) to exercise subgrid writing without METIS.
- For plane methods, ensure plane files live next to the `.prj` and normals point in the intended half-space.
