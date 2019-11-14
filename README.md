# BO_Toolbox

This repository provides the implementation of Bayesian Optimization in Matlab.

## Requirements
Only the Matlab enviroment is required. All required code is provided by this repository

## How to use
In [Functions](https://github.com/NRottmann/BO_Toolbox/tree/master/Functions) different benchmark functions are implemented, while the folder [BO](https://github.com/NRottmann/BO_Toolbox/tree/master/BO) contains the implementation of different algorithms in the context of Bayesian Optimization.

To perform classic Bayesian Optimization on a 30-dimensioanl Rosenbrock function over 50 iterations, one could use the following example.

```matlab
f = Rosenbrock(30);
results=BO(@f.call,f.vars, 'minimize', f.minimize, 'maxIter', 50)
```

Another example is the optimization vor the Dropout for Bayesian Optimization, at the Ackley function.
```matlab
f = Ackley(40);
results=DropoutBO(@f.call,f.vars, 'minimize', f.minimize, 'maxIter', 50)
```
