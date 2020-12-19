from scipy.stats import norm, poisson
import numpy as np


def weiner_process(m, N, t0, T):
	"""
	Function for modelling a Weiner process

		Parameters:

			m (int): number of sample paths
			N (int): number of approximation points
			t0 (float): simulation start time
			T (float): simulation end time

		Returns:

			trajectories (numpy.ndarray): a matrix with m lines and N columns where each line corresponds to a different
			sample path of a Weiner process from t0 to T, approximated at N points
			time (numpy.ndarray): an array equal to numpy.linspace(t0, T, N)
	"""

	dt = (T - t0) / (N - 1)
	trajectories = np.zeros((m, N))

	distr = norm(0, 1)

	for j in range(1, N):
		trajectories[:, j] = trajectories[:, j-1] + distr.rvs(size=m) * np.sqrt(dt)

	time = np.linspace(t0, T, N)
	return time, trajectories


def brownian_bridge(existent_points, new_time_points):
	"""
	Function for adding new approximation points to existing trajectory via the brownian bridge

		Parameters:

			existent_points (2-dim matrix (list or numpy.ndarray)): a matrix with the first line corresponding to
			time points and the second line corresponding to the process values at these time points
			new_time_points (iterable): list of new approximation points

		Returns:

			Two-dimensional numpy.ndarray with the first line corresponding to both existing and added time points
			and the second line corresponding to the process values at these time points
	"""
	if isinstance(existent_points, list):
		existent_points = np.array(existent_points)

	ep = existent_points[:, existent_points[0].argsort()]  # just in case
	new_time_points.sort()

	tn = ep[0][-1]
	for t in new_time_points:
		if t > tn:
			# append this point using the Weiner process definition

			Wtn = ep[1][-1]

			xi = norm(0, t - tn).rvs()
			Wt = Wtn + xi

			ep = np.c_[ep, [t, Wt]]
			tn = t
		else:
			# construct Brownian Bridge between a = W(t_{i - 1}) and b = W(t_i)

			i = np.searchsorted(ep[0], t)

			t1 = ep[0][i - 1]
			t2 = ep[0][i]

			if t == t1 or t == t2:
				continue

			a = ep[1][i - 1]
			b = ep[1][i]

			mu = a + (b - a) * (t - t1) / (t2 - t1)
			sigma = (t - t1) * (t2 - t) / (t2 - t1)

			Wt = norm(mu, sigma).rvs()

			ep = np.c_[ep[:, :i], [t, Wt], ep[:, i:]]

	return ep


def poisson_process(m, N, t0, T, lam):
	"""
	Function for modelling a Poisson process

		Parameters:

			m (int): number of sample paths
			N (int): number of approximation points
			t0 (float): simulation start time
			T (float): simulation end time
			lam (float): lambda parameter for Poisson process. Should be > 0

		Returns:

			time (numpy.ndarray): an array equal to numpy.linspace(t0, T, N)
			trajectories (numpy.ndarray): a matrix with m lines and N columns where each line corresponds to a different
			sample path of a Poisson process from t0 to T, approximated at N points

	"""
	dt = (T - t0) / N
	trajectories = np.zeros((m, N))

	p = poisson(dt * lam)

	for j in range(1, N):
		trajectories[:, j] = trajectories[:, j-1] + p.rvs(size=m)

	time = np.linspace(t0, T, N)
	return time, trajectories


def _ito_integration(f, a, b, k=2, step=0.1, m=1000):
	"""
	Function for modelling Ito's integrals and approximate calculation of its expectations

		Parameters:

			f (callable): integrand
			a (float): lower integration limit
			b (float): upper integration limit
			k (int): order of expectation
			step (float): integration step
			m (int): number of sample paths in Monte-Carlo simulation

		Returns:

			approximate value of k-th order expectation and matrix of sample paths
	"""
	n = int((b - a) / step)  # number of intervals
	t = np.linspace(a, b, n+1)
	
	_, trajectories = weiner_process(m, n+1, a, b)
	F = np.zeros((m, n))  # matrix of function values

	for i in range(m):
		F[i] = f(t[:-1], trajectories[i][:-1])

	I_f = np.zeros((m, n))
	for i in range(n):
		I_f[:, i] = F[:, i] * (trajectories[:, i + 1] - trajectories[:, i])

	return (I_f.sum(axis=1) ** k).sum() / m, trajectories


def ito_int_expect(f, a, b, k=2, step=0.1, m=1000):
	"""
	Returns approximate value of k-th order expectation of Ito's integral f(t, W(t))dW from a to b 
	calculated using Monte-Carlo simulation.

		Parameters:

			f (callable): integrand
			a (float): lower integration limit
			b (float): upper integration limit
			k (int): order of expectation
			step (float): integration step
			m (int): number of sample paths in Monte-Carlo simulation

		Returns:

			expect(float): approximate value of k-th order expectation
	"""
	expect, _ = _ito_integration(f, a, b, k=k, step=step, m=m)
	return expect


def ito_int_paths(f, a, b, step=0.1, m=1000):
	"""
	Returns sample paths of Ito's integral f(t, W(t))dW

		Parameters:

			f (callable): integrand
			a (float): lower integration limit
			b (float): upper integration limit
			step (float): integration step
			m (int): number of sample paths to generate

		Returns:

			expect(float): array of time points and a matrix in which each line corres
	"""
	_, paths = _ito_integration(f, a, b, k=1, step=step, m=m)
	time = np.arange(a, b+step, step)
	return time, paths


def sde(dy, y0, t, step, nsamp=1000, method='euler'):
	"""
	Solves Ito's stochastic differential equations of the form "dy = f*dt + g*dW"
	using Monte-Carlo simulation and Milstein or Euler's method

		Parameters:

			dy (callable): function for calculation of dy; should take 3 arguments corresponding to y, dt and dW
			and return tuple with values of f and g for Euler's method or values of f, g, and partial derivative
			dg/dy for Milstein's method
			y0 (float): initial condition
			t (tuple <float, float>): starting and final time of the solution
			step (float): integration step
			nsamp (int), optional: number of sample paths used in Monte-Carlo simulation; default 1000
			method (str), optional: numerical method; either 'milstein' or 'euler' (default);
			Milstein's method is more precise, but requires additional derivation
			
		Returns:

			time: 1d numpy.ndarray with time points of approximation
			paths: matrix in which each row corresponds to a different solution sample path
	"""

	n = norm(0, step)
	N = int((t[1] - t[0]) / step)  # number of intervals

	paths = np.zeros((nsamp, N + 1))
	if method == 'euler':
		for j in range(nsamp):
			y = np.zeros(N + 1)
			y[0] = y0
			W = n.rvs(size=N)
			for i in range(N):
				f, g = dy(y[i], step, W[i])
				y[i + 1] = y[i] + f * step + g * W[i]
			paths[j] = y

	elif method == 'milstein':
		for j in range(nsamp):
			y = np.zeros(N + 1)
			y[0] = y0
			W = n.rvs(size=N)
			for i in range(N):
				f, g, dgdx = dy(y[i], step, W[i])
				y[i + 1] = y[i] + f * step + g * W[i] + 1 / 2 * g * dgdx * (W[i] - step) ** 2
			paths[j] = y
	time = np.arange(t[0], t[1] + step, step)
	return time, paths
