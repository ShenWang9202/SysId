# Ordinary least-squares (OLS) for estimating partially observed LTI systems

This repository contains MATLAB scripts that implement the numerical examples in the following paper:

1) Zheng, Y. and Li, N. (2020). [Non-asymptotic Identification of Linear Dynamical Systems Using Multiple Trajectories](https://arxiv.org/abs/2009.00739), preprint.

## Instructions

To run the scripts in this repository, you only need a working MATLAB installation. We implement four types of OLS methods and the celebrated Ho-Kalman algorithm Ref. \[1\]
* OLS using mutliple independent trajectories, where all data points are utilized (our method)
* OLS using multiple independent trajectories, where only the last data point of each trajectory is used; Ref \[2\]
* OLS using a single trajectory; Ref \[3\]
* OLS + prefilter using a single trajectory; Ref \[4\]

**References**

1) Ho, Β. L., and Rudolf E. Kálmán. "Effective construction of linear state-variable models from input/output functions." at-*Automatisierungstechnik* 14.1-12 (1966): 545-548
2) Sun, Y., Oymak, S., & Fazel, M. (2020, July). [Finite sample system identification: Optimal rates and the role of regularization](http://proceedings.mlr.press/v120/sun20a/sun20a.pdf). In Learning for Dynamics and Control (pp. 16-25). PMLR.
3) Oymak, S., & Ozay, N. (2019, July). [Non-asymptotic identification of lti systems from a single trajectory](https://arxiv.org/pdf/1806.05722.pdf). In 2019 American Control Conference (ACC) (pp. 5655-5661). IEEE.
4) Simchowitz, M., Boczar, R., & Recht, B. (2019). [Learning linear dynamical systems with semi-parametric least squares](https://arxiv.org/pdf/1902.00768.pdf). arXiv preprint arXiv:1902.00768.


## Troubleshooting
If you have any trouble running the scripts in this repository, please email [Yang Zheng](mailto:zhengy@g.harvard.edu?Subject=SOS-csp).
