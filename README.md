# Measuring the M-$\sigma$ relation using measurements of galactic bulge velocity dispersion and supermassive black holes masses

---

By Kyle Soni

08 December 2023

## Summary

Here is the result:

## Introduction

Previous studies of elliptical galaxies and bulges in spiral galaxies have shown evidence that they contain supermassive blackholes at their center. Moreover, an important relationship which correlates the mass of the supermassive black holes to the velocity dispersion of stars in the surrounding galactic bulge has been discovered, often known as the M-$\sigma$ relation [1]. This relation is significant because it provides insights into the co-evolution of galaxies and their central black holes. It also helps with galaxy classification and cosmology, fitting our understanding of galaxies into a larger cosmological framework.

The relation has been estimated many times, structured as a linear relationship between $\text{log}(M_{BH})$ and $\text{log}(\sigma)$ of the form [1]
$$\text{log}(M_{BH}/M_\odot) = b + m \text{log}(\sigma / 200 \text{ km s$^{-1}$})$$

In this report, I will perform Bayesian analysis on 51 mass and $\sigma$ measurements to compare my calculation for this relation with previous work.

## Data

The data was compiled by a previous paper in which they examine various properties of galactic bulges [1]. Specifically, they surveyed the literature and collected mass measurements ($M_{BH}$) for dynamically detected central black holes, along with the velocity dispersion of the bulge in the host galaxies. They used the effective velocity dispersion ($\sigma_e$) if available, and the central stellar velocity dispersion ($\sigma_c$) found in HyperLEDA (an astronomical database) if not. The data consists of 51 galaxies and corresponding black holes, with uncertainties for both $\sigma$ and $M_{BH}$. A visualization of the data and its uncertainties is shown in Figure 1.

One important note is that the original data had some asymmetric uncertainties for the mass measurmenets. For the purposes of this project, I averaged the uncertainty for each mass and then converted to log space.

[1] Kayhan G ̈ultekin, Douglas O. Richstone, Karl Gebhardt, Tod R. Lauer, Scott Tremaine, M. C. Aller, Ralf Bender, Alan Dressler, S. M. Faber, Alexei V. Filippenko, Richard Green, Luis C. Ho, John Kormendy, John Magorrian, Jason Pinkney, and Christos Siopis. The m– and m–l relations in galactic bulges, and determinations of their intrinsic scatter. The Astrophysical Journal, 698(1):198, may 2009.

## Analysis/Results

The first step to the analysis was formulating the likelihood. I made a couple of key assumptions: the error distributions for both the mass and $\sigma$ are approximately Gaussian and the measurements of $M_{BH}$ and $\sigma$ are indepedent. Using these assumptions, I can project the data to be orthogonal to a best fit line and account equally for both uncertainties. I'll start by expressing the slope $m$ in terms of a normal vector

$\vec{n} =
\begin{pmatrix}
-\sin(\alpha)\\
\cos(\alpha)
\end{pmatrix}$

where $m = \arctan(\alpha)$. Since $M_{BH}$ and $\sigma$ are independent, our covariance matrix is

$\Sigma_i = 
\begin{pmatrix}
\sigma_{x_i}^2 & 0\\
0 & \sigma_{y_i}^2
\end{pmatrix}$

Continuing the derivation, the distance to the line will be 
$$\Delta_i = \vec{n}^T \begin{pmatrix} x_i \\ y_i \end{pmatrix} - b \cos(\alpha) = \begin{pmatrix} -\sin(\alpha) & \cos(\alpha) \end{pmatrix} \begin{pmatrix} x_i \\ y_i \end{pmatrix} - b \cos(\alpha) = y_i \cos(\alpha) - x_i \sin(\alpha) - b \cos(\alpha)$$

To calculate the standard deviation estimate, we follow

$$S_i^2 = \vec{n}^T \Sigma_i \vec{n} = \begin{pmatrix} -\sin(\alpha) & \cos(\alpha) \end{pmatrix} 
\begin{pmatrix}
\sigma_{x_i}^2 & 0\\
0 & \sigma_{y_i}^2
\end{pmatrix} 
\begin{pmatrix} -\sin(\alpha) \\ \cos(\alpha) \end{pmatrix} = \sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2$$

This means our pdf is 
$$p(y_i \mid m, b, x_i, \sigma_{y_i}, \sigma_{x_i}) = \frac{1}{\sqrt{2\pi(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}}\exp{- \frac{(y_i \cos(\alpha) - x_i \sin(\alpha) - b \cos(\alpha))^2}{2(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}}$$

However, this model has no way of accounting for outliers, so I'll introduce a mixture model to compensate for this. This gives a pdf
$$p(y_i\mid m, b, x_i, \sigma_{y_i}, \sigma_{x_i}, \omega, \mu_b, \sigma_b) = \frac{(1-\omega) }{\sqrt{2\pi(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}}\exp{- \frac{(y_i \cos(\alpha) - x_i \sin(\alpha) - b \cos(\alpha))^2}{2(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}} + \omega \mathcal{N}(\mu_b, \sigma_b)$$

where the model parameters are defined as follows: $\omega$ is the proportion of outliers in the data, $\mu_b$ is the mean of the background (causing the outliers), and $\sigma_b$ is the standard deviation of the background. With this pdf, we can now write the full likelihood.

$$\mathcal{L} = \prod_i p(y_i \mid m, b, x_i, \sigma_{y_i}, \sigma_{x_i}, \omega, \mu_b, \sigma_b) = \frac{(1-\omega) }{\sqrt{2\pi(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}}\exp{- \frac{(y_i \cos(\alpha) - x_i \sin(\alpha) - b \cos(\alpha))^2}{2(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}} + \omega \mathcal{N}(\mu_b, \sigma_b)$$

$$\log \mathcal{L} = \sum_i \log( \frac{(1-\omega) }{\sqrt{2\pi(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}}\exp{- \frac{(y_i \cos(\alpha) - x_i \sin(\alpha) - b \cos(\alpha))^2}{2(\sin^2(\alpha) \sigma_{y_i}^2 + \cos^2(\alpha) \sigma_{x_i}^2)}} + \omega \mathcal{N}(\mu_b, \sigma_b) )$$

This likelihood is implemented in code in the cell below.