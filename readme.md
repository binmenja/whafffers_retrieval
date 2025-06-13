# Derivation of the Optimal Estimation Method (OEM) for Atmospheric Retrieval

**Author**: Benjamin Riot-Bretecher (based on code by Lei Liu)  
**Date**: June 13, 2025

## Abstract

This document provides a detailed derivation of the Optimal Estimation Method (OEM), a widely used statistical approach for retrieving atmospheric state parameters from remote sensing measurements. The method, based on Bayes' theorem, combines prior knowledge of the atmospheric state with observed measurements to obtain an optimal estimate of the true state. The derivation covers the cost function, the iterative solution (specifically the Levenberg-Marquardt approach), and the characterization of the retrieval through the averaging kernel and posterior covariance.

## Introduction

Atmospheric remote sensing aims to infer the state of the atmosphere (e.g., temperature, water vapor, trace gas profiles) from measurements of radiation. This is an ill-posed inverse problem, meaning multiple atmospheric states could produce similar measurements, or the solution might be highly sensitive to noise. The Optimal Estimation Method (OEM), pioneered by C.D. Rodgers [1], provides a robust statistical framework to address this challenge by optimally combining measurement information with prior knowledge about the atmospheric state.

The core idea is to find the atmospheric state vector $\mathbf{x}$ that minimizes a cost function, quantifying the discrepancy between simulated and actual measurements, and between the retrieved state and an *a priori* estimate.

## Forward Model and State Vectors

Let $\mathbf{y}$ be the vector of $m$ atmospheric measurements (e.g., radiance at different wavenumbers or channels).

Let $\mathbf{x}$ be the vector of $n$ atmospheric state parameters to retrieve (e.g., temperature and water vapor concentrations at different altitude levels).

The relationship between the atmospheric state and measurements is described by a non-linear forward model $F(\mathbf{x})$:

$$
\mathbf{y} = F(\mathbf{x}) + \varepsilon
$$

where $\varepsilon$ represents measurement noise and forward model errors, assumed to follow a Gaussian distribution with zero mean and a known measurement error covariance matrix $\mathbf{S}_e$:

$$
\varepsilon \sim \mathcal{N}(\mathbf{0}, \mathbf{S}_e)
$$

## Bayesian Formulation and Cost Function

OEM is rooted in Bayes' theorem, seeking the probability density function (PDF) of the atmospheric state $\mathbf{x}$ given measurements $\mathbf{y}$, denoted $P(\mathbf{x}|\mathbf{y})$:

$$
P(\mathbf{x}|\mathbf{y}) = \frac{P(\mathbf{y}|\mathbf{x}) P(\mathbf{x})}{P(\mathbf{y})}
$$

where:
- $P(\mathbf{y}|\mathbf{x})$ is the likelihood of observing measurements $\mathbf{y}$ given state $\mathbf{x}$.
- $P(\mathbf{x})$ is the *a priori* PDF of the state $\mathbf{x}$.
- $P(\mathbf{y})$ is the normalization constant.

Assuming the *a priori* state $\mathbf{x}_a$ is known with covariance matrix $\mathbf{S}_a$:

$$
\mathbf{x} \sim \mathcal{N}(\mathbf{x}_a, \mathbf{S}_a)
$$

The *a priori* term is:

$$
P(\mathbf{x}) \propto \exp\left(-\frac{1}{2}(\mathbf{x} - \mathbf{x}_a)^T \mathbf{S}_a^{-1} (\mathbf{x} - \mathbf{x}_a)\right)
$$

The likelihood term is:

$$
P(\mathbf{y}|\mathbf{x}) \propto \exp\left(-\frac{1}{2}(\mathbf{y} - F(\mathbf{x}))^T \mathbf{S}_e^{-1} (\mathbf{y} - F(\mathbf{x}))\right)
$$

Maximizing the posterior PDF $P(\mathbf{x}|\mathbf{y})$ is equivalent to minimizing the cost function $J(\mathbf{x})$:

$$
J(\mathbf{x}) = \frac{1}{2} (\mathbf{y} - F(\mathbf{x}))^T \mathbf{S}_e^{-1} (\mathbf{y} - F(\mathbf{x})) + \frac{1}{2} (\mathbf{x} - \mathbf{x}_a)^T \mathbf{S}_a^{-1} (\mathbf{x} - \mathbf{x}_a)
$$

## Iterative Solution: Levenberg-Marquardt Algorithm

Since $F(\mathbf{x})$ is typically non-linear, an iterative approach like the Levenberg-Marquardt algorithm is used. At iteration $i$, the forward model is linearized around the current state $\mathbf{x}_i$:

$$
F(\mathbf{x}) \approx F(\mathbf{x}_i) + \mathbf{K}_i (\mathbf{x} - \mathbf{x}_i)
$$

where $\mathbf{K}_i$ is the Jacobian matrix:

$$
\mathbf{K}_i = \frac{\partial F(\mathbf{x})}{\partial \mathbf{x}} \Big|_{\mathbf{x}=\mathbf{x}_i}
$$

The iterative update equation is:

$$
\mathbf{x}_{i+1} = \mathbf{x}_i + (\mathbf{K}_i^T \mathbf{S}_e^{-1} \mathbf{K}_i + (1 + \lambda_i)\mathbf{S}_a^{-1})^{-1} \left[ \mathbf{K}_i^T \mathbf{S}_e^{-1} (\mathbf{y} - F(\mathbf{x}_i)) - \mathbf{S}_a^{-1} (\mathbf{x}_i - \mathbf{x}_a) \right]
$$

Alternatively:

$$
\mathbf{x}_{i+1} = \mathbf{x}_i + (\mathbf{K}_i^T \mathbf{S}_e^{-1} \mathbf{K}_i + \mathbf{S}_a^{-1} + \lambda_i \mathbf{S}_a^{-1})^{-1} \left[ \mathbf{K}_i^T \mathbf{S}_e^{-1} (\mathbf{y} - F(\mathbf{x}_i)) - \mathbf{S}_a^{-1} (\mathbf{x}_i - \mathbf{x}_a) \right]
$$

This matches the MATLAB code:

```matlab
x(:, i+1) = gather(x(:,i) + inv((1+lambda)*inv(Sa) + K'*inv(Se)*K)*(K'*inv(Se)*(measurement - F)-inv(Sa)*(x(:,i)-xa)));
```

The damping parameter $\lambda_i$ controls the step:
- If $$J(\mathbf{x}_{i+1})$$ < $$J(\mathbf{x}_i)$$, accept the step and decrease $$\lambda_i$$ (e.g., $$\lambda_{i+1} = \lambda_i / 10$$).
- If $$J(\mathbf{x}_{i+1}) \ge J(\mathbf{x}_i)$$, reject the step, revert to $$x_i$$, and increase $$\lambda_i$$ (e.g., $$\lambda_{i+1} = \lambda_i \times 10$$).

Iteration continues until convergence, e.g., when the change in the state vector satisfies a threshold.

## Retrieval Characterization

### Posterior Covariance Matrix

The posterior covariance matrix $\hat{\mathbf{S}}$ describes the uncertainty of the retrieved state $$\hat{\mathbf{x}}$$:

$$
\hat{\mathbf{S}} = (\mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1})^{-1}
$$

For Levenberg-Marquardt, accounting for $\lambda$:

$$
\mathbf{S}_{pos} = (\lambda \mathbf{S}_a^{-1} + \mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1})^{-1} (\lambda^2 \mathbf{S}_a^{-1} + \mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1}) (\lambda \mathbf{S}_a^{-1} + \mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1})^{-1}
$$

This matches the MATLAB code:

```matlab
inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K)*((lambda_output(i)+1).^2*inv(Sa)+K'*inv(Se)*K)*inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K);
```

### Averaging Kernel Matrix

The averaging kernel matrix $$\mathbf{A}$$ relates the retrieval $$\hat{\mathbf{x}}$$ to the true state $$\mathbf{x}_{t}$$:

$$
\hat{\mathbf{x}} - \mathbf{x}_a = \mathbf{A} (\mathbf{x}_{t} - \mathbf{x}_a) + \mathbf{G} \epsilon
$$

where $\mathbf{G}$ is the gain matrix. For a linear problem:

$$
\mathbf{A} = (\mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1})^{-1} \mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K}
$$

Or:

$$
\mathbf{A} = \mathbf{G K}, \quad \mathbf{G} = (\mathbf{K}^T \mathbf{S}_e^{-1} \mathbf{K} + \mathbf{S}_a^{-1})^{-1} \mathbf{K}^T \mathbf{S}_e^{-1}
$$

This matches the MATLAB code:

```matlab
A = gather(inv(K'*inv(Se)*K+(lambda_output(i)+1).*inv(Sa)))*gather(K'*inv(Se)*K);
```

### Degrees of Freedom for Signal (DFS)

The DFS, the trace of $\mathbf{A}$, quantifies the number of independent pieces of information:

$$
\text{DFS} = \text{trace}(\mathbf{A})
$$

A high DFS indicates measurement-driven retrieval, while a low DFS suggests *a priori* dominance.

## Conclusion

The OEM provides a statistically rigorous framework for solving ill-posed inverse problems in atmospheric remote sensing. By minimizing a cost function balancing measurement fidelity and prior knowledge, OEM yields an optimal state estimate. Diagnostics like the posterior covariance and averaging kernel matrices characterize the uncertainty and information content, essential for understanding the reliability of retrieved atmospheric products.

## References

[1] Rodgers, C. D. (2000). *Inverse Methods for Atmospheric Sounding: Theory and Practice*. World Scientific.
