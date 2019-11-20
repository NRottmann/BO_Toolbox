# BO_Toolbox

This repository provides the implementation of Bayesian Optimization in Matlab.

## Requirements
Only the Matlab enviroment is required. All required code is provided by this repository

## How to use
In [Functions](https://github.com/NRottmann/BO_Toolbox/tree/master/Functions) different benchmark functions are implemented, while the folder [BO](https://github.com/NRottmann/BO_Toolbox/tree/master/BO) contains the implementation of different algorithms in the context of Bayesian Optimization.

In order to use this toolbox the script [initToolbox.m](https://github.com/NRottmann/BO_Toolbox/blob/master/initToolbox.m) has to be run initially, to add the requiered subfolder to the matlab path.

To perform classic Bayesian Optimization on a 30-dimensioanl Rosenbrock function over 50 iterations, one could use the following example.

```matlab
intiToolbox
f = Rosenbrock(30);
results=BO(@f.call,f.vars, 'minimize', f.minimize, 'maxIter', 50)
```

Another example is the optimization vor the Dropout for Bayesian Optimization, at the Ackley function.
```matlab
f = Ackley(40);
results=DropoutBO(@f.call,f.vars, 'minimize', f.minimize, 'maxIter', 50)
```
All provided BO algorithms, retrun a matlab struct with data about the optimization process:
* results.valueHistory: the values received from function evaluation
* results.maxValueHistory: for each iteration the best function value so far
* results.paramHistory: for each iteration, the parameters used for evaluation
* results.nextFeature: for each iteration the created feature (not provided by BO since no feature is generated)
* results.bestValue - best seen function value
* results.bestParams - parameters for the best function value
