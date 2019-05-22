# psiDMLab

## What is psiDMLab

psiDMLab is meant as a numerical playground for the Schroedinger-Poisson equations, 
a model that utilises quantum mechanical techniques to approximate the numerically 
intractable Vlasov-Poisson system for the evolution of dark matter. 
Currently, our focus lies on modularity not utmost performance. The idea is to draft 
a numerically stable and accurate scheme that can later on be optimized to run 
on large scale systems and hopefully is capable of simulating 3+1 dimensional 
cosmological setups comparable to what is known for N-Body simulations.

## Theory behind psiDM
psiDM aims to model the dynamics of (cold) dark by employing a time evolution
described by the Schroedinger-Poisson (SP) PDE system:
\begin{equation}
    i \partial_t \psi(x,t) = \left[-\frac{\hbar^2}{2m} \nabla^2 +
    V(x,t)\right]\psi(x,t)
\end{equation}
From a cosmological point of view, two distinct applications of wavelike dark
matter can be identified.
1. Ultralight Axion Dark Matter
Axions are considered a viable dark matter particle candidate. In this
framework, axions (spin 0 bosons) form a cosmic Bose-Einstein condensate that
obeys SP, i.e. evolves quantum mechanically under its self gravity. For
sufficiently small axion masses, structure formation on the kpc-scale is controlled
by quantum mechanical effects such as Heissenberg's uncertainty princlple acting
as a "regularizer" of the emerging density structures. This solves long standing
difficulites of the standard $\Lambda$CDM model, such as the core cusp problem
or the predicition of diverging densities in cosmic halos. On scales vastly
larger than the axion de-Broglie wavelength $\psi$DM recovers the standard CDM
description.

2. Approximation to the Vlasov-Poisson system (VP)
Standard CDM obeys a collisionless Boltzmann equation (also known as Vlasov
equation in plasma physics) coupled to Poisson's
equation. Assuming a collisionless setting, this can be shown to be true even
for the one-particle distribution function $f(x,p,t)$. Unfortunately, solving
VP is hard both numerically and analytically. Therefore, people (almost always)
rely on N-Body simulations yielding a noisy sample from the true phase space
distributionn. Partial analytical results, however, show a good correspondence
between SP and VP which is not fully understood. Numerical results support those
analytical findings. For instance, the associated potential of the wavefunction
was shown to converge to the classical potential $\propto{\left(\frac{\hbar}{m}\right)^2}$.
Moreover, choosing Husimi's distribution to compute the associated phase
space distribution of $\psi$ to obtain a "smoothed eveolution" produces
distribution moments (such as matter densities or veloctiy dispersions) that
resemble closely what we expect in the classical case.
### 

## What dependenicies does the simulation have
TODO

## How to compile
TODO

## How to run
TODO

## Data Post Processing
TODO

## psiDM in the literature
TODO
