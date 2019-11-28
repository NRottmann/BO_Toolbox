This folder contains the implementation of the Gaussian Process Regression, used as a model in Baysiean Optimization.

## How to use
The GP Regression is provided as a function. To use the GP Regression implementation, you have to provide the training points ```x```, test points ```s``` (where you want to sample from the GP Regression) and target values ```y```, for the training data.
The funtion returns the the predictive mean and variance at the sample points ```s```. In the following code example ```D``` is the number of dimensions and ```n``` the number of training data.
```matlab
x = training points; % Size: D x n
s = sample points;   % Size: D x l
y = target values;   % Size: n x 1 
[mu,sigma] = GP(x,s,y);
```

Next to the required input, one can provide additional parameters, like the kernel function, or pretrained parameters of the kernel function.
```matalb
[mu,sigma] = GP(x,s,y, 'CovFunc', 'se_kernel_var', 'CovParam', [0.1; 1]);
```
For additional parameters, please check out the documentation in [GP.m](https://github.com/NRottmann/BO_Toolbox/blob/master/GP/GP.m).
