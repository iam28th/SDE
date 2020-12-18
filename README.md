# Numerical Monte-Carlo solution to the system of stochastic differential equations of population dynamics

## Goals and objectives
The final aim of this project was to apply It̂o's stochastic differential equations (SDEs) approach to the field of population dynamics and to compare a stochastic population dynamics model with a deterministic one. In order to use this approach, several tasks should be covered:
1. Stochastic process modelling
2. It̂o's stochastic integral modelling and approximate calculation of its expectations
3. Obtaining numerical solutions to stochastic differential equation (SDE) using Monte-Carlo simulation

## Methods
### It̂o's lemma
Let F be a function F(t, X(t)) and X be a stochastic process such that

<a href="https://www.codecogs.com/eqnedit.php?latex=dX&space;=&space;\mu&space;dt&space;&plus;&space;\sigma&space;dW" target="_blank"><img src="https://latex.codecogs.com/gif.latex?dX&space;=&space;\mu&space;dt&space;&plus;&space;\sigma&space;dW" title="dX = \mu_t dt + \sigma dW" /></a>

then (assuming the necessary assumptions) function F satisfies the following SDE:

<a href="https://www.codecogs.com/eqnedit.php?latex=dF&space;=&space;(\frac{\partial&space;F}{\partial&space;t}&space;&plus;&space;\mu&space;\frac{\partial&space;F}{\partial&space;x}&space;&plus;&space;\frac{1}{2}\sigma^2&space;\frac{\partial^2&space;F}{\partial&space;x^2})dt&space;&plus;&space;\sigma&space;\frac{\partial&space;F}{\partial&space;x}&space;dW" target="_blank"><img src="https://latex.codecogs.com/gif.latex?dF&space;=&space;(\frac{\partial&space;F}{\partial&space;t}&space;&plus;&space;\mu&space;\frac{\partial&space;F}{\partial&space;x}&space;&plus;&space;\frac{1}{2}\sigma^2&space;\frac{\partial^2&space;F}{\partial&space;x^2})dt&space;&plus;&space;\sigma&space;\frac{\partial&space;F}{\partial&space;x}&space;dW" title="dF = (\frac{\partial F}{\partial t} + \mu \frac{\partial F}{\partial x} + \frac{1}{2}\sigma^2 \frac{\partial^2 F}{\partial x^2})dt + \sigma \frac{\partial F}{\partial x} dW" /></a>

It̂o's lemma allows us to find closed-form solution to some SDEs and exactly calculate some stochastic integrals.

### Monte-Carlo simulation 
Monte-Carlo simulation is a computational algorithm that relies on random sampling to obtain numerical result. Pseude-random numbers obtained from Numpy generator were used to simulate randomness

