mineral-deposit-segmentation-pipeline
=====================================

A pipeline for segmenting mineral deposits on clinical magnetic resonance imaging (MRI) volumes of the brain [1,2,3,4].

The pipeline was originally developed at the Brain Imaging Research Centre [5] as part of a PhD project supported by SINAPSE/Spirit [6] on Red Hat Enterprise Linux 5 [7] and Matlab 2011b [8].

Requirements:
- Recent Linux distribution (e.g. Neurodebian [9]) with Bash [10] and Perl [11] and common utilities, such as make and rsync. 
- FMRIB FSL library [12]
- N4 intensity inhomogeneity correction tool [13] e.g. from Slicer3 [14]
- Tools for writing and reading NIfTI and ANALYZE image in MATLAB [15]
- LIBRA: a Matlab library for robust analysis [16]

References:
[1] http://www.biomedical-image-analysis.co.uk/images/stories/glatz-posters1-32.pdf
[2] http://www.neuroinformatics2012.org/abstracts/intuitive-and-efficient-deployment-of-neuroimaging-pipelines-in-clinical-research-with-bricpipe
[3] Digital object identifier (DOI): 10.1594/ecr2013/C-2293
[4] https://www.ucl.ac.uk/cabi/PDF/PG-BCISMRM_2013_Abstract_booklet_final.pdf
[5] http://www.bric.ed.ac.uk
[6] http://www.sinapse.ac.uk
[7] http://www.redhat.com/products/enterprise-linux
[8] http://www.mathworks.com
[9] http://neuro.debian.net/
[10] http://www.gnu.org/software/bash
[11] http://www.perl.org
[12] http://fsl.fmrib.ox.ac.uk/fsl
[13] http://www.insight-journal.org/browse/publication/640
[14] http://www.slicer.org/
[15] http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/
[16] http://wis.kuleuven.be/stat/robust/LIBRA/LIBRA-home

