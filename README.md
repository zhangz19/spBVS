# spBVS
Bayesian variable selection for spatial areal data with Conditional Autoregressive (CAR) model. 
Please see [Introduction](http://stt.msu.edu/~zhangz19/spBVS.html).

* The program is developed for modeling spatial continuous/count data with Conditional Autoregressive (CAR) model on errors/random effects. The program allow flexible choices of L_0 and L_1/L_2 penalties. For Poisson model, it implements the truncated MALA for efficiently block-sampling the high-dimensional latent parameters and takes advantage of the sparse precision matrix under CAR model. It can reduce to handeling general Bayesian variable selectoin for non-spatial data.

* Please contact the authors if there are any questions or implementation issues: Zhen Zhang, zhangz19@galton.uchicago.edu (Date coded: 2013-11-19)

* Example publications:
  1. Woznicki, S. A., Nejadhashemi, A. P., Ross, D. M., Zhang, Z., Wang, L. and Esfahanian, A. (2015). Ecohydrological model parameter selection for stream health evaluation. _Science of the Total Environment_. 511, 341-353. [URL](http://www.sciencedirect.com/science/article/pii/S0048969714017689)
  
  2. Herman, M., Nejadhashemi, A. P., Daneshvar, F., Ross, D. M., Woznicki, S. A., Zhang, Z. and Esfahanian, A. (2015). Optimization of Conservation Practice Implementation Strategies in the Context of Stream Health. _Ecological Engineering_. 84, 1-12.
  
  3. Woznicki, S. A., Nejadhashemi, A. P., Abouali, M., Herman, H. R., Esfahanian, E., Hamaamin, A. Y. and Zhang, Z. (2016). Ecohydrological Modeling for Large-scale Environmental Impact Assessment. _Science of the Total Environment_. 543, 274-286.
  
  4. Herman, M., Nejadhashemi, A. P., Daneshvar, F., Abouali, M., Ross, D. M., Woznicki, S. A. and Zhang, Z. (2016). Optimization of bioenergy crop selection and placement based on a stream health indicator using an evolutionary algorithm. _Journal of Environmental Management_. 181, 413-424.
  
  5. Li, Y., Bence, J. R., Zhang, Z. and Ebener, M. P. (2017). Why do lake whitefish move long distances in Lake Huron? Bayesian variable selections of factors explaining fish movement distance. _Fisheries Research_. 195, 169-179. [URL](http://www.sciencedirect.com/science/article/pii/S0165783617301960)
