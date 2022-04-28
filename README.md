# covcomp
This is an R package for implementing the covariance completion algorithm in [1].
The syntax is as follows:
1. For manual choice of truncation parameters used in the algorithm:
```
Compln(Cov, Omega, VectorOfTuning) 
```
2. For automatic choice of truncation parameters in accordance with the fraction of variance explained (FVE) criterion:
```
ComplnFVE(Cov, Omega, FVE) 
```
Here `Omega` is a matrix of 1s and 0s which describes the domain Ω of the partial covariance K_{Ω}.

<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">


## References
<a id="1">[1]</a> 
Waghmare, K.G. and Panaretos, V.M., 2021. The Completion of Covariance Kernels. arXiv preprint arXiv:2107.07350. (https://arxiv.org/abs/2107.07350)
