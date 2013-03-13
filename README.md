mineral-deposit-segmentation-pipeline
=====================================

A pipeline for segmenting mineral deposits on clinical magnetic resonance imaging (MRI) volumes of the brain [1,2,3].

The pipeline was originally developed at the Brain Imaging Research Centre [4] as part of a PhD project supported by SINAPSE/Spirit [5] on Red Hat Enterprise Linux 5 [6] and Matlab 2011b [7].

Requirements:
- Recent Linux distribution that with Bash [8] and Perl [9] and
  utilities such as make, rsync and common command line utilities. 
- FMRIB FSL library [10]
- N4 intensity inhomogeneity correction tool [11] (e.g. from Slicer3 [12])
- Tools for NIfTI and ANALYZE image in MATLAB [13]
- LIBRA: a Matlab library for robust analysis [14]

References:
[1] http://www.biomedical-image-analysis.co.uk/images/stories/glatz-posters1-32.pdf
[2] http://www.neuroinformatics2012.org/abstracts/intuitive-and-efficient-deployment-of-neuroimaging-pipelines-in-clinical-research-with-bricpipe
[3] Digital object identifier (DOI): 10.1594/ecr2013/C-2293
[4] http://www.bric.ed.ac.uk
[5] http://www.sinapse.ac.uk
[6] http://www.redhat.com/products/enterprise-linux
[7] http://www.mathworks.com
[8] http://www.gnu.org/software/bash
[9] http://www.perl.org
[10] http://fsl.fmrib.ox.ac.uk/fsl
[11] http://www.insight-journal.org/browse/publication/640
[12] http://www.slicer.org/
[13] http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/
[14] http://wis.kuleuven.be/stat/robust/LIBRA/LIBRA-home

