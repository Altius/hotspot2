hotspot2 is a program developed at the Altius Institute for Biomedical Sciences in Seattle, U.S.A.
for identifying genomic regions with statistically significant "hotspots," or enrichments,
of cleavage activity in DNase-seq experiments.  It is designed to run in a UNIX environment,
and to work with alignment files in BAM format.

hotspot2 requires each of the following programs to be installed and accessible from the user's path
before it can be used:
* g++ (or an equivalent C++ compiler; replace "g++" with its name in Makefile in this case)
* [bedops](https://github.com/bedops/bedops)
* [samtools](https://github.com/samtools)
* [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils)
* [modwt](https://github.com/StamLab/modwt)

hotspot2 tallies cleavages within a small region ("window") around each site.  It slides the window
across the genome, and statistically evaluates cleavage tallies within their local context, i.e.,
within a larger "background window" that itself slides across the genome.  To enable hotspot2
to distinguish between "no data" (e.g., an unmappable region or the end of a chromosome) and
"no cleavages" (a valid region in which no cleavages were observed in the DNase-seq experiment),
a file of chromosome sizes is required.  Ideally, a file containing all the mappable regions
in the genome (e.g., regions uniquely mappable by 36mers) should also be supplied by the user.
(If a "blacklist" of, e.g., mappable but problematic microsatellite regions is available
for the genome, it should be subtracted from the file of mappable regions before being supplied
to hotspot2.)

Note:  In the file of chromosome sizes, the "start" of each chromosome (column 2) must be 0.
This file must also be an uncompressed BED file; .starch format is currently disallowed.

Before hotspots can be identified, the set of viable positions that can serve as centers of
sliding windows must be determined.  The script `extractCenterSites.sh` in the scripts subdirectory
must be run to determine these positions.  This script requires a file of chromosome sizes
and the name of the file to which the positions should be written.  A file of mappable regions
can optionally be provided as well, and this is highly recommended.  The default radius around
each position is 100 bp, yielding a sliding window of width 201 bp; a different value can be
supplied if desired.  To see the usage information for this script, type

    scripts/extractCenterSites.sh -h

Note:  `extractCenterSites.sh` only needs to be run once per genome.  (If analyses with different
window sizes are desired, `extractCenterSites.sh` needs to be run once per window size per genome.)

Hotspots are called for an input alignment file in BAM format via a set of scripts that are
each executed by the `hotspot2.sh` script.  This script also executes two programs that need to
be compiled or "made" on the computer where hotspot2.sh will run.  To make these programs, and
two others that are optionally run by the `density-peaks.bash` script, simply type the command `make`
from within the hotspot2 directory. (This will place the programs in the subdirectory named "bin.")

Note:  After the programs are made, their location (subdirectory "bin") must be added to the user's PATH.

Once the hotspot2 programs have been compiled and the center sites file has been created
by running `extractCenterSites.sh`, `hotspot2.sh` will be ready to run.  To see the usage information
for this script, including descriptions of its various parameters and their default settings, type

    scripts/hotspot2.sh -h

To run hotspot2.sh with default values for its various parameters, type

    scripts/hotspot2.sh yourData.bam yourOutputDirectory

(`hotspot2.sh` will create the output directory `yourOutputDirectory` if it does not already exist.)
In this example, after `hotspot2.sh` completes its tasks, the following files will be located
in `yourOutputDirectory`:

* yourData.allcalls.starch
* yourData.cleavage.total
* yourData.cutcounts.starch
* yourData.density.bw
* yourData.density.starch
* yourData.fragments.sorted.starch
* yourData.hotspots.fdr0.05.starch
* yourData.peaks.narrowpeaks.starch
* yourData.peaks.starch
* yourData.SPOT.txt

In one typical use case, only two of these files might be of interest, `yourData.hotspots.fdr0.05.starch`
and `yourData.SPOT.txt`.  The former contains the hotspots called at the specified (in this case, default)
FDR threshold.  The latter contains the SPOT score, or Signal Portion Of Tags; it is a metric that
gives an indication of the quality of the sample and/or of the experiment.  The SPOT score is simply
the number of (mappable) cleavages observed in hotspots divided by the total number of (mappable) cleavages;
it therefore can range from 0 to 1.

After hotspots have been called at the specified threshold, the user might be interested in examining
hotspot calls at a different FDR threshold.  This can be done without re-running the entire hotspot2
pipeline, using the `yourData.allcalls.starch` file and the `hsmerge.sh` script.  In this example, to call
hotspots at FDR threshold 0.01, type

    scripts/hsmerge.sh -f 0.01 yourOutputDirectory/yourData.allcalls.starch yourOutputDirectory/yourData.hotspots.fdr0.01.starch

Note:  The `SITECALL_THRESHOLD` (`-F`) that was supplied to hotspot2.sh will be the upper limit at which
hotspots can be re-called using `hsmerge.sh`.  If it was set to, e.g., 0.05 for the sake of speed and
the size of the `yourData.allcalls.starch` file, and hotspots called at threshold 0.10 are then desired,
`hotspot2.sh` would need to be re-run with `-f 0.10` and `-F 0.10` (or a higher threshold for the latter).

In another typical use case, alignment files from many DNase-seq experiments are available
(they might correspond to a comprehensive set of cell and tissue types for an organism,
or to samples of a given cell type taken from patients that constitute "cases" and "controls"
in a study of a disease or trait), and it is desired to create a "master list" or "Index"
of the universe of consensus DHSs across the samples. Finer granularity is needed to produce
such an Index than hotspots typically provide, and for this purpose (and other select purposes),
narrower regions of highly focused sensitivity to DNase I, called "peaks" or "DHSs,"
can be called by hotspot2 via its script `density-peaks.bash` and then used as input
to the Index-building script that is part of our [Index](https://github.com/Altius/Index) project.
For flexibility and for historical reasons, `density-peaks.bash` can produce variable-width
or fixed-width (150-bp) peaks. To produce an Index for a set of samples, variable-width peaks
are required.

To produce variable-width peaks during a full run of `hotspot2.sh`, supply the option `-p varWidth_nn_ID`
as an argument to `hotspot2.sh`, with "nn" replaced by the minimum width (in bp) desired for a peak
(20 is suggested and commonly used for analyses at Altius) and "ID" replaced by a unique identifier
for the sample in question (e.g., a unique portion of the name of the sample's alignment file).

Variable-width peaks can also be produced via the script `density-peaks.bash` without running
hotspot2 in its entirety. To do so, the following input files are needed:
* the `yourData.cutcounts.starch` file produced by the initial run of `hotspot2.sh` (if the user
accidentally or purposely deletes this file, s/he can recreate it via the script `cutcounts.bash`)
* the `yourData.cleavage.total` file produced by `hotspot2.sh` (this file can likewise be recreated
via `cutcounts.bash`)
* the file of hotspots called at the user's desired threshold (hotspots can be called at any
threshold of interest via `hsmerge.sh`, as described above)
* an appropriate file of chromosome sizes, in BED (not .starch) format, with column 2 set to 0
in every row.
See the usage statement at the beginning of the code for `density-peaks.bash` for the full list
of arguments and the order in which they must be provided. The name of a temporary directory
must be given as the first argument, and `varWidth_nn_ID` (as described above, without the
"-p" preceding it) must be given as the second argument.

### Output file formats

All output files in BED or .starch format are 0-based. In the following, the first 3 columns
in such files are named "seqname", "beg", and "end", and the 4th column is a placeholder
("i", lowercase letter i) and/or identifying string ("ID").

* `allcalls.starch`: seqname, beg, end, "i", FDR
* `cleavage.total`: this one-line text file contains the total number of mapped cleavages
* `cutcounts.starch`: seqname, beg, end, "i", number of mapped cleavages at that position
* `density.bw`: this is a [bigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) version of density.starch
* `density.starch`: seqname, beg, end, ID, unnormalized density (the number of mapped cleavages
in a 150-bp window centered on the given 20-bp interval)
* `hotspots.starch`: seqname, beg, end, ID, -10\*log10(FDR) rounded to the nearest integer
and capped at 1000, ".", "-1", "-1", -log10(FDR) capped at 100
* `peaks.starch` (_variable-width option_): seqname, beg, end, ID, maximum normalized density within
this element, FWHM summit coordinate, wavelet summit coordinate (maximum normalized density =
(1000000/cleavage.total) \* (unnormalized density); see `cleavage.total` and `density.starch` above)
* `peaks.starch` (_all other options_): seqname, beg, end, "i", maximum unnormalized density within this element
* `narrowpeaks.starch`: seqname, beg, end, ".", "0", ".", column 5 from the corresponding peaks.starch file,
"-1", "-1", "75"
* `SPOT.txt`: this one-line text file contains the SPOT score (see above)


hotspot2 was developed by Eric Rynes, Jeff Vierstra, Jemma Nelson, Richard Sandstrom, Shane Neph,
and Audra Johnson.

Questions and feature requests are welcome, and may be e-mailed to erynes@altius.org.
