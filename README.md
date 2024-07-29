# Measuring the M-U+03C3 relation using measurements of galactic bulge velocity dispersion and supermassive black holes masses

---

By Kyle Soni

08 December 2023

## Summary

There is a well-known relation between the mass of the supermassive black holes to the velocity dispersion of stars in the surrounding galactic bulge, often known as the M-$\sigma$ relation. However, the data for these objects has uncertainties for both the mass and the sigma (on both axes). I sought to fit a model to this relationship with a Bayesian analysis, accounting for both uncertainties and outliers with a mixture model. The results are shown in this figure:

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig5.png?raw=true)

Data with the best fit line, following the equation $\text{log}(M_{BH}/M_\odot) = (8.20 \pm 0.16) + (4.09 \pm 0.01) * \text{log}(\sigma / 200 km*s^{-1})$.

## Introduction

Previous studies of elliptical galaxies and bulges in spiral galaxies have shown evidence that they contain supermassive blackholes at their center. Moreover, an important relationship which correlates the mass of the supermassive black holes to the velocity dispersion of stars in the surrounding galactic bulge has been discovered, often known as the M-$\sigma$ relation [1]. This relation is significant because it provides insights into the co-evolution of galaxies and their central black holes. It also helps with galaxy classification and cosmology, fitting our understanding of galaxies into a larger cosmological framework.

The relation has been estimated many times, structured as a linear relationship between $\text{log}(M_{BH})$ and $\text{log}(\sigma)$ of the form [1]
$$\text{log}(M_{BH}/M_\odot) = b + m \text{log}(\sigma / 200 \text{ km s$^{-1}$})$$

In this report, I will perform Bayesian analysis on 51 mass and $\sigma$ measurements to compare my calculation for this relation with previous work.

## Data

The data was compiled by a previous paper in which they examine various properties of galactic bulges [1]. Specifically, they surveyed the literature and collected mass measurements ($M_{BH}$) for dynamically detected central black holes, along with the velocity dispersion of the bulge in the host galaxies. They used the effective velocity dispersion ($\sigma_e$) if available, and the central stellar velocity dispersion ($\sigma_c$) found in HyperLEDA (an astronomical database) if not. The data consists of 51 galaxies and corresponding black holes, with uncertainties for both $\sigma$ and $M_{BH}$. A visualization of the data and its uncertainties is shown in Figure 1.

One important note is that the original data had some asymmetric uncertainties for the mass measurmenets. For the purposes of this project, I averaged the uncertainty for each mass and then converted to log space.

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig1.png?raw=true)

**Figure 1.** Plot of the $M-\sigma$ relation in log space. The data is composed of 51 measurements with uncertainties for both variables.

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

For the prior, I chose to have a flexible prior to prevent heavy influence of the result. Specifically, I chose a uniform distribution for each parameter with these ranges: $m \in [0, 300]$, $b \in [-100, 100]$, $\omega \in [0, 1]$, $\mu_b \in [-100, 100]$, and $\sigma_b \in [0, 500]$. I also assumed that the probabilitiy distribution for each parameter was independent.

Lastly, I calculated the log of the posterior, which is simply the sum of the log likelihood and log prior.

With the log posterior, I could now perform Markov Chain Monte Carlo via the Metropolis-Hastings algorithm. A few diagnostic plots are shown in Figure 2. As shown, we have high density in the high probability region, trace plots with no significant structure for $m$ and $b$, and an acceptance ratio of a little under 0.2. The acceptance ratio may be lower than ideal, but these suggest that the MH sampler had reasonable success.

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig2.png?raw=true)

**Figure 2.** Various diagnostic plots for the MH sampler. Top left is $m$ vs. $b$ with lighter colors having higher probability. Bottom left is the dataset with the best fit-lines overplotted. The 4 plots on the right are trace plots for the quantities $m$, $b$, log of the posterior, and acceptance ratio.

Looking at the outlier parameters, their trace plots are shown in Figure 3. We also see that $\omega \approx 0.02$, which is a low fraction of the data. This may explain why $\mu_b$ and $\sigma_b$ were quite unstable, continously fluctuating between the ranges of the prior. However, after testing with many different configurations, the variance in the background mean and standard deviation did not seem to have a measurable effect on the best fit for $m$ and $b$, so that does strengthen my confidence in the results.

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig3.png?raw=true)

**Figure 3.** Trace plots for the outlier parameters $\omega$, $\mu_b$, and $\sigma_b$.

With the Metropolis-Hastings sampling completed, the only step left is to marginalize over the different parameters to get marginal distributions for $m$ and $b$. The marginal distributions are shown in Figure 4, 95% credible regions are shown in the cell below, and the best-fit line overplotted on the data is shown in Figure 5. The credible regions were calculated by making the reasonable assumption that the distributions were Gaussian and using the standard deviation.

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig4.png?raw=true)

**Figure 4.** Marginal distributions for the parameters $m$ (left) and $b$ (right).

The 95% credible region for m is: 4.09 +/- 0.01.
The 95% credible region for b is: 8.20 +/- 0.16.

![alt text](https://github.com/kylesoni/Measuring-the-M-sigma-relation/blob/main/figures/fig5.png?raw=true)

**Figure 5.** Data with the best fit line, following the equation $\text{log}(M_{BH}/M_\odot) = (8.20 \pm 0.16) + (4.09 \pm 0.01) * \text{log}(\sigma / 200 km*s^{-1})$.

## Conclusions

After attempting to fit a line to 51 observations and measure the $M-\sigma$ relation, I found
$$\text{log}(M_{BH}/M_\odot) = (8.20 \pm 0.16) + (4.09 \pm 0.01) * \text{log}(\sigma / 200 \text{ km s$^{-1}$})$$
Previous studies found $m$ to be around 3.5-5 and $b$ to be around 8, so these values are not unreasonable given previous knowledge. However, the uncertainty for $m$ is suspicously low, so it may be possible that my outlier correction was too harsh. Future work may seek to use a more robust likelihood for the outliers or investigate the most inconsistent galaxies directly. Regardless, this measurement is still useful for solidying our understanding of the $M-\sigma$ relation, and it can be used for comparison in future studies.