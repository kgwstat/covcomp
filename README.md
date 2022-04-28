# covcomp
This is an R package for implementing the covariance completion algorithm in [1].
The syntax is as follows:
'''
Compln(Cov, Omega, VectorOfTuning) % for manual choice of truncation parameters used in the algorithm.
ComplnFVE(Cov, Omega, FVE) % for automatic choice of truncation parameters in accordance with the fraction of variance explained (FVE) criterion.
'''

## References
<a id="1">[1]</a> 
Waghmare, K.G. and Panaretos, V.M., 2021. The Completion of Covariance Kernels. arXiv preprint arXiv:2107.07350. (https://arxiv.org/abs/2107.07350)
