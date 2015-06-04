FATools - Fragment Analysis Tools

This is a collection of tools for performing fragment analysis on electropherograms from
capillary sequencing machine.

Features:
- peak finding & sizing (using cubicspline, leastsquare)
- peak binning
- peak filtering (stutter, dye overlapping, marker artefact)


Procedure:

import FSA into database & renormalize

scan each channel

preannotate all channels

align ladder

call peak for all non-ladder channels

bin channels

post-annotate channels




