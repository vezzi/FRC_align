INSTALLATION
==============


```
git clone https://github.com/vezzi/FRC_align/
cd FRC_align
```

You may need to install `cmake`, e.g. `sudo apt install cmake`

From the FRCurve directory run:

```
mkdir build
cd build
cmake ..
make
```

You will find the binaries in the main directory under bin. In case of problems the majority of the times there is a problem
with the local installation of boost.


DESCRIPTION
==============
 FRCbam is a tool able to evaluate and analyze de novo assembly/assemblers. The tool has been already successfully 
applied in several de novo .
 
 FRCbam: tool to compute Feature Response Curves in order to validate and rank assemblies and assemblers


FRCbam
--------------

### USAGE: basic, no CE-stats tuning**
 
* `--pe-sam` A_tool1_PE_lib.bam`: sorted bam file obtained aligning PE library against assembly obtained with tool A;
* `--pe-min-insert MIN_PE_INS`: estimated min insert length
* `--pe-max-insert MAX_PE_INS` : estimated max insert length
* `--mp-sam A_tool1_MP_lib.bam`: sorted bam file obtained aligning MP library against assembly obtained with tool A;
* `--mp-min-insert MIN_MP_INS` : estimated min insert length
* `--mp-max-insert MAX_MP_INS` : estimated max insert length
* `--genome-size ESTIMATED_GENOME_SIZE`: estimated genome size;
* `--output OUTPUT_HEADER`: output header;
	
#### IMPORTANT:

If `--genome-size` is not specified the assembly length is used to compute FRCurve. In order to be able to compare FRCurves
obtained with different tools (and hence producing slightly different assembly sizes) the same `ESTIMATED_GENOME_SIZE`
must be specified.
		
#### OUTPUT:

* `OUTPUT_HEADER_Features.txt`: human readable description of features: contig start end feature_type
* `OUTPUT_HEADER_FRC.txt`: FRCurve computed with all the features (to be plotted)
* `OUTPUT_HEADER_FEATURE.txt`: FRCurve for the corresponding feature
* `OUTPUT_HEADER_featureType.txt`: for each featureType the specific FRCurve
* `Features.gff`: features description in GFF format (for visualization)
* `OUTPUT_HEADER_CEstats_PE.txt`: CEvalues distribution (for CE_stats tuning)
* `OUTPUT_HEADER_CEstats_MP.txt`: CEvalues distribution (for CE_stats tuning)
		
### USAGE: advanced, CE-stats tuning

CE-stats are able to identify the presence of insertion and deletion events. Different insert sizes give the possibility to
identify different events. In order to avoid too many False Positives (or too many False Negatives) a tuning phase is 
highly recommended.

Once step 3 of USAGE is done, the user can already plot the FRCurves (for all the features or for only some of them).
CE_stats based features have been computed with default (i.e., not optimal) parameters. Each run of FRCbam produces two
files: `OUTPUT_HEADER_CEstats_PE.txt` and `OUTPUT_HEADER_CEstats_MP.txt`. These files contain the distribution of the `CE_values`
on each assembly. These values must be plotted as suggested in page 3 of Supplementary Material 
(see http://www.nada.kth.se/~vezzi/publications/supplementary.pdf) to estimate the optimal `CE_min` and `CE_max` values for 
 the PE and MP library respectively.
 Once the optimal parameters are estimated FRCurves must be recomputed for all assemblies (only `COMPR` and `STRECH` features
 will change) specifying the following extra parameters:
 
* `--CEstats-PE-min CE_PE_MIN`: all position with CE values computed with PE library lower than this are considered compressions
* `--CEstats-PE-max CE_PE_MAX`: all position with CE values computed with PE library higher than this are considered expansions 
* `--CEstats-MP-min CE_MP_MIN`: all position with CE values computed with MP library lower than this are considered compressions
* `--CEstats-MP-max CE_MP_MAX`: all position with CE values computed with MP library higher than this are considered expansions
 





LICENCE
==============
All the tools distributed with this package are distributed under GNU General Public License version 3.0 (GPLv3). 



