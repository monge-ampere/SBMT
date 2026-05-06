# Structured Bitmap-to-Mesh Triangulation (SBMT)

This repository contains a research prototype of **Structured Bitmap-to-Mesh Triangulation (SBMT)**, a structured triangulation framework for **geometry-aware discretization of bitmap- or segmentation-derived image domains**.

SBMT starts from a **regular triangular scaffold** and embeds polygonal boundaries through finite **lookup-table-driven local retriangulation templates**. The method is designed as a lightweight discretization front-end for image-derived domains, emphasizing exact boundary embedding, regular interior preservation, low-complexity local updates, parallel-friendly execution, and compatibility with downstream PDE / FEM / geometric analysis tasks.

> SBMT is not intended to replace general-purpose mesh generators in all scenarios.  
> It is designed for raster / segmentation-derived domains where boundary-exact embedding, structured interior preservation, deterministic local updates, and numerical-analysis-friendly discretization are especially important.

---

## Overview

Many geometry-processing and numerical-analysis pipelines start from image-derived domains, such as binary masks, bitmap regions, or segmentation outputs. Directly operating on raster grids is often inconvenient for PDE-based or FEM-like computations, while general-purpose meshing tools may destroy the regular interior structure or introduce global remeshing behavior that is unnecessary for this type of input.

SBMT addresses this setting by:

1. laying down a regular triangular background grid,
2. extracting / prescribing polygonal boundaries from raster domains,
3. classifying local segment–triangle intersection patterns,
4. applying finite lookup-table-based retriangulation templates,
5. preserving most of the regular interior structure while modifying only a thin boundary band.

The result is a boundary-conforming triangular mesh suitable for downstream geometry-aware numerical processing.

---

## Key Features

- **Bitmap / segmentation-native input**  
  Designed for image-derived domains rather than general CAD-style geometry.

- **Exact local boundary embedding**  
  Prescribed polygonal boundaries are embedded into the triangular grid through local retriangulation.

- **Regular triangular scaffold**  
  The interior is initialized as a regular triangular grid and is preserved as much as possible.

- **Finite lookup-table-driven retriangulation**  
  Local intersection configurations are classified into a finite set of symbolic cases, each associated with a predefined retriangulation template.

- **Deterministic local updates**  
  The method avoids global remeshing cascades and performs boundary-local modifications.

- **Low complexity and parallel-friendly structure**  
  Local templates and conflict-free updates make the framework naturally suitable for large raster domains and future parallel implementations.

- **PDE-ready discretization**  
  The generated meshes are designed to support interpolation, heat diffusion, elliptic PDE solvers, Hodge-type computations, and other geometry-aware numerical tasks.

---

## Publication

The corresponding paper has been published in **Graphical Models**:

> **Structured Bitmap-to-Mesh Triangulation for Geometry-Aware Discretization of Image-Derived Domains**  
> *Graphical Models*, Vol. 145, 2026, Article 101326  
> DOI: https://doi.org/10.1016/j.gmod.2026.101326

Preprint:

https://arxiv.org/abs/2602.19474

A Chinese technical introduction is available here:

https://zhuanlan.zhihu.com/p/2023809046206006841

---

## Repository Structure

```text
SBMT/
├── SBMT/                 # Main source files
│   ├── ImgProc class     # Image processing and boundary extraction
│   ├── KdTree/           # KD-tree code for nearest-neighbor queries
│   ├── Algo.h            # Core triangulation and lookup-table logic
│   ├── FileProc.h        # Raster and polygon I/O utilities
│   └── data/             # Example inputs
├── SBMT.sln              # Visual Studio solution file
└── README.md             # This file
```

---

## Build and Platform

This repository is currently a **Windows / Visual Studio research prototype**.

- **Tested with:** Microsoft Visual Studio 2015
- **Platform:** Windows
- **Current status:** Single-threaded prototype
- **Not yet supported:** Linux / macOS / GCC / Clang / CMake

> **Note:** Earlier internal project files may have been created under older Visual Studio settings.  
> If you encounter project-version issues, please retarget the solution in Visual Studio.

---

## Quick Start

At present, the codebase is intended primarily as a research prototype rather than a polished software package.

Typical usage:

1. Open `SBMT.sln` in Microsoft Visual Studio 2015.
2. Build the project in the appropriate configuration.
3. Use the sample raster / polygon inputs in the `data/` directory.
4. Run the executable to generate triangulation outputs or visualizations.

A more streamlined command-line interface, documentation, and cross-platform build system are planned for future versions.

---

## Current Implementation Status

This implementation focuses on demonstrating the core algorithmic principles of SBMT. It is not yet an optimized production implementation.

Current limitations include:

- No half-edge or fully modular mesh data structure;
- Triangles are mainly stored and accessed through relatively simple array-based structures;
- Some boundary-related and retriangulation logic remains prototype-oriented;
- Some local cases, especially certain `b`-type handling routines, are implemented in a slow, non-optimized manner;
- No multithreading or GPU acceleration;
- No adaptive refinement or local resolution control;
- No CMake or cross-platform build system.

Despite these limitations, the repository captures the main structural ideas of SBMT:

- Regular triangular scaffold construction;
- Boundary interaction classification;
- Lookup-table-based local retriangulation;
- Geometric preprocessing;
- Boundary-conforming mesh generation.

---

## Method Scope

SBMT is most suitable for problems where:

- The input domain is derived from raster images, bitmap masks, or segmentation outputs;
- Boundary conformity is required;
- A mostly regular interior triangulation is desirable;
- The mesh will be used for downstream geometry processing, interpolation, PDE, FEM, or related numerical tasks;
- Deterministic and local update behavior is preferred.

SBMT is not designed as a universal replacement for mature mesh generators such as Triangle or Gmsh. Those tools remain more appropriate for many general-purpose meshing problems, especially when global mesh-quality optimization, adaptive refinement, or broad CAD-style input support is required.

---

## Third-Party Code Notice

This repository includes KD-tree code derived from an implementation attributed to:

> Matthew B. Kennel, Institute for Nonlinear Science, UCSD, 2004.

It is used for nearest-neighbor queries in boundary-related processing.

Please review the original attribution and licensing status carefully before redistributing or reusing this component beyond research / prototype purposes.

If you know the canonical source or licensing details of this KD-tree implementation, contributions or corrections are welcome.

---

## Citation

If you use SBMT or refer to this repository in academic work, please cite the corresponding paper.

Please verify the author list against the official publication page before use.

```bibtex
@article{feng2026sbmt,
  title   = {Structured Bitmap-to-Mesh Triangulation for Geometry-Aware Discretization of Image-Derived Domains},
  author  = {Feng, Wei and Zheng, Haiyong},
  journal = {Graphical Models},
  volume  = {145},
  pages   = {101326},
  year    = {2026},
  doi     = {10.1016/j.gmod.2026.101326}
}
```
---

## Planned Future Work

- Modularization of retriangulation rules;
- Improved mesh data structures;
- Half-edge or adjacency-aware mesh representation;
- Performance profiling and benchmark scripts;
- OpenMP or other CPU-side parallel execution;
- Possible GPU acceleration;
- Adaptive refinement and local resolution control;
- CMake-based build system;
- Linux / macOS compatibility;
- Tighter integration with PDE / FEM solvers.

---

## License

A formal open-source license will be added after the third-party code status is fully clarified.

Until then, please treat this repository as a research prototype made publicly available for academic reference.

For redistribution, commercial use, or integration into other projects, please contact the author and review all third-party code notices carefully.

---

## Contact

For questions, comments, or collaboration related to SBMT, geometry processing, bitmap-to-mesh conversion, or PDE-ready discretization, please contact the author through GitHub.
