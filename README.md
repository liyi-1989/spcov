# sicb2013

## Background

This repository includes the source code for the analyses I presented 
at the [EDEN](http://edenrcn.com) sponsored symposium ["Understanding First 
Order Phenotypes: Transcriptomics for Emerging Model Systems"](http://www.sicb.org/meetings/2013/symposia/phenotypes.php) 
at the 2013 SICB conference, described in 
[abstract S4-2.1](http://www.sicb.org/meetings/2013/SICB%202013%20abstracts.pdf).

This study is a collaborative project by Brown University faculty members:

- [Casey Dunn](http://www.brown.edu/Faculty/Dunn_Lab/)

- [Xi Luo](http://www.stat.brown.edu/FacultyDisplay.aspx?id=1317070681)

- [Zhijin Wu](http://www.stat.brown.edu/FacultyDisplay.aspx?id=1128605312)

All collaborators contributed to the code presented here.

Additional details will be provided in a forthcoming publication. A preprint of 
this publication is available at:

Dunn, CW, X Luo, and Z Wu (2013) Phylogenetic analysis of gene expression. 
[arXiv:1302.2978](http://arxiv.org/abs/1302.2978)

Many thanks to [Joe Felsenstein](http://www.gs.washington.edu/faculty/felsenstein.htm) 
for helpful discussions about multivariate comparative analyses.

## Use

To run the simulation and analyses, execute the following shell command:

    R CMD BATCH simulation.r
    
This will regenerate the pdf file.

## Citations

This code implements the two regularization methods described in:

Bickel, P. J. & Levina, E. Covariance regularization by thresholding. 
Ann. Statist. 36, 2577–2604 (2008). [doi:10.1214/08-AOS600](http://dx.doi.org/10.1214/08-AOS600)

Luo, X. High Dimensional Low Rank and Sparse Covariance Matrix Estimation via 
Convex Minimization. arXiv.org (2011). [arxiv](http://arxiv.org/abs/1111.1133)


## License

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU General Public License as published by
the Free Software Foundation](http://www.gnu.org/licenses/), either 
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
