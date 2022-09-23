# chromo-looper
Chromatin loops can connect genetic elements together, this package attempt to do the same at the data-level and builds on tools such as pybedtools. The main goal of this package is to annotate loops, not 1D data but loops. It can become easy to tangle yourself in cobbweb of  in the intersection madness, I just want to lay out that the best way to intersect loop data is to: 1) using the BEDPE fields as an ID and 2) for each intersection between a single loop file and multiple 1D files, perform each loop-1d intersection one at a time, then merge all intersection data using the loop ID.

Note: Package setup according to: https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html
