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

## System requirements
Most of the project (except population dynamics modelling) was implemented in Python using laptop with the following resources:
Intel(R) Core(TM) i5-8300H CPU @ 2.30Ghz, 8GB RAM, Windows 10

Python and libraries:
- Python 3.8.5, 
- Numpy 1.19.1
- Scipy 1.5.0
- Matplotlib 3.3.1
- Seaborn 0.11.0

Population dynamics modelling was implemented in Google colab (with free access). Language and libraries versions are the following:
- Python 3.6.9
- Numpy 1.19.4
- Scipy 1.4.1
- Matplotlib 3.2.2
- Seaborn 0.11.0
- Numba 0.48.0

## API

## Results

## References 

1. Allen, Edward. (2007). Modeling with Itô Stochastic Differential Equations. 22. 10.1007/978-1-4020-5953-7. 
2. Choongbum Lee, Peter Kempthorne, Vasily Strela, Jake Xia. 18.S096 Topics in Mathematics with Applications in Finance. Fall 2013. Massachusetts Institute of Technology: MIT OpenCourseWare, https://ocw.mit.edu/. License: Creative Commons BY-NC-SA.
