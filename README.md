# Structured Bitmap-to-Mesh Triangulation (SBMT)

This repository contains a prototype implementation of **Structured Bitmap-to-Mesh Triangulation (SBMT)**, a rule-driven triangulation framework designed for raster domains. SBMT deterministically embeds polygonal boundaries into a regular triangular grid using a finite set of lookup templates.

---

## 🔧 Status

This is a **single-threaded prototype** implementation focused on algorithmic correctness. Current limitations include:

- No half-edge structure (triangles are stored and traversed via double loops)
- Brute-force handling of `$b$`-type retriangulation logic (sub-optimal and slow)
- No parallelism or GPU acceleration
- No support for adaptive refinement

Performance and data structure optimizations are planned for future versions.

---

## 📘 Reference

If you use or reference this code in academic work, please cite the corresponding paper:

> **Structured Bitmap-to-Mesh Triangulation: A Boundary-Embedding Framework for Numerical Geometry on Raster Domains**  
> _Submitted to SIIMS, 2025._

---

## 📂 Structure

```bash
SBMT/SBMT/            # Main source code
  KdTree              # KdTree code
  ImgProc             # Image analysis code
  Algo.h              # general algorithms code
  FileProc.h          # reading and writing files code
  data/               # Example input raster and boundary

## 🔍 Third-Party Code Notice

This project includes an adapted version of the `kdtree` implementation by **Matthew B. Kennel**, from the **Institute for Nonlinear Science, UCSD (2004)**.

- Original code: [kdtree by Matthew B. Kennel (2004)]
- License: [Specify license, if known; otherwise write "license not specified by author"]
- Status: Used in original or slightly modified form for nearest-neighbor queries during segment-triangle classification.

All rights belong to the original author. If you use this code, please also acknowledge the original source.

