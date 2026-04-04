# Package index

## Visualization

Core functions for producing peak visualizations on gene and genomic
region structures.

- [`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
  : Plot RNA-Binding Protein Peaks on a Single Gene
- [`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
  : Plot RNA-Binding Protein Peaks on a Genomic Region

## Splicing and Motif Analysis

Functions for quantifying RNA-binding protein occupancy and sequence
motif frequency around exon-intron boundaries.

- [`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)
  : Create Splicing Map
- [`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
  : Analyzes the frequency of a target sequence motif across splicing
  junction regions. Compares motif frequency between Retained, Excluded,
  and Control splicing events to identify position-specific enrichment
  patterns.

## Utilities

Helper functions for loading annotations and preparing input data.

- [`LoadGTF()`](https://krushna-b.github.io/RNAPeaks/reference/LoadGTF.md)
  :

  Loads GTF gene annotation data from AnnotationHub for Human or Mouse.
  This function can be called once to load the annotation, which can
  then be passed to other functions like
  [`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
  or
  [`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
  to avoid repeated downloads.

- [`checkBed()`](https://krushna-b.github.io/RNAPeaks/reference/checkBed.md)
  : Validate and Normalize a BED-Format Data Frame

## Sample Data

Curated datasets bundled with the package for testing and demonstration.

- [`sample_bed`](https://krushna-b.github.io/RNAPeaks/reference/sample_bed.md)
  : Sample K562 RBP Binding Peaks
- [`sample_se.mats`](https://krushna-b.github.io/RNAPeaks/reference/sample_se.mats.md)
  : Sample SE.MATS Skipped-Exon Splicing Events
- [`gtf_human`](https://krushna-b.github.io/RNAPeaks/reference/gtf_human.md)
  : Human GTF Gene Annotation

## Package

Package-level documentation.

- [`RNAPeaks`](https://krushna-b.github.io/RNAPeaks/reference/RNAPeaks-package.md)
  [`RNAPeaks-package`](https://krushna-b.github.io/RNAPeaks/reference/RNAPeaks-package.md)
  : RNAPeaks: Visualize RNA-Binding Protein Peaks on Gene Structures
