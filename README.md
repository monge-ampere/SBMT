# 🧩 Structured Bitmap-to-Mesh Triangulation (SBMT)

This repository contains a prototype implementation of **Structured Bitmap-to-Mesh Triangulation (SBMT)**, a rule-driven triangulation framework designed for **raster domains**.  
SBMT deterministically embeds polygonal boundaries into a regular triangular grid using a finite set of **lookup-table-based retriangulation templates**, supporting structure-aware numerical analysis on bitmap data.

---

## 🔧 Status

This is a **single-threaded prototype** focused on demonstrating the core algorithmic principles. Current limitations include:

- No half-edge data structure: triangles are stored via flat arrays and accessed through double loops.
- Brute-force handling of `$b$`-type retriangulation logic: currently implemented in a slow, non-optimized manner.
- No parallel execution: neither multithreading nor GPU acceleration is currently supported.
- No adaptive refinement or local resolution control.
- **Platform limitation:** This version is tested and built exclusively with **Microsoft Visual Studio 2013 (Windows)**.  
  It currently **does not support cross-platform compilation** on Linux or macOS.

Performance optimization, algorithm modularization, and platform portability are planned for future versions.

---

## 📘 Reference

If you use or reference this code in academic work, please cite the corresponding paper:

> **Structured Bitmap-to-Mesh Triangulation: A Boundary-Embedding Framework for Numerical Geometry on Raster Domains**  
> _Submitted to SIIMS, 2025._

---

## 📂 Structure

```
SBMT/
├── SBMT/                 # Main source files
│   ├── ImgProc class     # Image processing & boundary extraction
│   ├── KdTree/           # 3rd-party KD-tree for nearest neighbor queries
│   ├── Algo.h            # Core triangulation and lookup logic
│   ├── FileProc.h        # Raster and polygon I/O utilities
│   └── data/             # Example inputs (raster domain and polygon boundary)
├── SBMT.sln              # Visual Studio 2013 solution file
└── README.md             # This file
```

---

## 🔍 Third-Party Code Notice

This project includes an adapted version of the `kdtree` implementation by:

> **Matthew B. Kennel**, Institute for Nonlinear Science, UCSD (2004)

- Purpose: Used for nearest-neighbor queries during boundary-segment classification.
- Original source: [URL or citation if known]
- License: *License not specified by the original author*
- Status: Integrated in near-original form with minimal changes.

All rights belong to the original author. If you use this code, please also acknowledge the original source.

---

## 🖥️ Platform Compatibility

- ✔️ Tested with: **Windows 7+ / Visual Studio 2013**
- ❌ Not tested on: Linux / macOS / gcc / clang
- No makefile or CMake script provided yet.
- Some C++11 features may require adjustment for compatibility with modern compilers.

Future versions will aim to provide:
- A modular build system (e.g., CMake)
- POSIX compatibility for Linux/macOS
- Support for modern compilers (MSVC ≥ 2019, GCC, Clang)

---

## 🚧 Planned Features

- ✅ Modularization of retriangulation rules
- ✅ Refactor `$b$`-type logic for efficiency
- 🚀 Parallel execution (OpenMP or GPU)
- 📈 Performance profiling and benchmarks
- 🎯 Adaptive refinement for geometric detail
- 🛠️ CMake-based cross-platform build
- 🧪 Integration with numerical solvers (PDE, FEM)

---

Feel free to fork, modify, or contact the author if you're interested in collaboration or contributing to this project.
