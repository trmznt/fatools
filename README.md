
FATools - Fragment Analysis Tools
=================================


FATools is a Python library and a collection of command line tools to process and
analysis DNA fragments (STRs/microsatellites) on electropherograms from capillary
sequencing machine assay, an activity also known as Fragment Analysis (FA).

In short, this tool provides more or less similar function as ABI GeneMapper software
in command line interface (sans GUI, but an assay viewer can be built using this
library).

Features:
- baseline normalization (using white tophat morphology transformation)
- peak scanning (CWT-based algorithm)
- ladder-size assignment (mostly using greedy approach with DP algorithm)
- peak sizing (cubic spline, least square, and local southern method)
- peak binning
- peak annotation and filtering (for stutter, dye overlapping, and peak artifact)

