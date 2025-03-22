---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Atmospheric pressure and altitude

## Introduction

In this chapter, we will try to find in an ab initio manner how atmospheric pressure and altitude are related to each other. 

## The barometric formula

Let’s start from considering an ideal gas. This approximation works well with the Earth’s atmosphere, because it has a sufficiently low density and high temperature.  From the ideal gas law, the pressure is given by:

$$
\begin{equation}
    p=nk_BT
\end{equation}
$$ (idealgas)

where $n$ the numerical density of the gas, and $T$ its temperature. To model Earth’s atmosphere, we imagine an infinitely high column of gas that is subjected to an uniform gravitational field. This approximation stems from the assumptions that (i) Earth is spherical and that (ii) the vast majority of the gas mass stays below an altitude that is much smaller than Earth’s radius. Both approximations are very good. We will remove approximation (ii) in a later section. 

The pressure of the gas at certain altitude $h$ must be equivalent to the pressure exerted by the column of gas above it, due to the gravitational field, i.e.:

$$
\begin{equation}
    p(h)=\dfrac{M(h)g}{A}
\end{equation}
$$ (pressure_gravfield)

with $M$ the mass of the gas column, and $A$ the area. We can calculate the mass per unit of area, $M/A$, as the integral of the mass density $\rho(z)$ of the gas along the altitude:

$$
\begin{equation}
p(h) = g\int_h^{+\infty}\rho(z) dz  = g\int_h^{+\infty}m_0(z)n(z) dz 
\end{equation}
$$

where we expressed the mass density $\rho$ as the product between the average mass of the gas particle, $m_0$, and the numerical density, $n$, both depending on the altitude. We can rearrange and convert the equation in its differential form, knowing that $n(+\infty)=0$:

$$
\begin{equation}
\dfrac{d}{dh}p(h)=-g\,m_0(h)n(h) = -g\,m_0(h)\dfrac{p(h)}{k_BT(h)}
\end{equation}
$$(diff_pressure)

Where we used equation {eq}'idealgas' to express the numerical density in terms of density and temperature. It is important to note that equation (6) now only depends on variables at the altitude $h$ and its proximity (the derivative), and not on the whole gas column that sits above. This mathematical step permits us to reach for a solution, but at the cost of not being able to find the absolute value of the pressure, expect for getting it empirically. Rearranging equation {eq}`diff_pressure` we have:

$$
\begin{equation}
\dfrac{\dfrac{d}{dh}p(h)}{p(h)}=-\dfrac{g}{k_B}\dfrac{m_0(h)}{T(h)}
\end{equation}
$$

By integrating in $dz$ from $0$ to $h$ we eventually obtain the general solution of how the pressure varies with the altitude:

$$
\begin{equation}
p(h)=p(0)\exp\bigg[-\dfrac{g}{k_B}\int_0^h\dfrac{m_0(z)}{T(z)}dz\bigg]
\end{equation}
$$

Where $p(0)$, the pressure at $h=0$, is not known, and must be obtained from the observed data. If we choose $h=0$ to be the sea level, then $p(0)$ is the barometric pressure, which in standard conditions is 1013.25 hPa. We will understand how temperature and average mass vary with altitude in the following section. For now, let’s assume that are both constants and see how the pressure varies in such simple case. The resolution of the integral becomes trivial, as we have:

$$
\begin{equation}
p(h)=p(0)\exp\bigg[-\dfrac{m_0gh}{k_BT}\bigg]
\end{equation}
$$(barometric_formula)

We recognize $m_0gh$ as the potential energy of a single gas particle, and $k_BT$ its thermal energy. Wait, what? Well, I rediscovered hot water of course, as equation (9) is called barometric formula. But it is nonetheless interesting to notice that the exponential term is the Boltzmann factor ($e^{-E/k_BT}$), which, in a canonical ensemble (NVT, our case), represents the probability of system to be in a state with energy $E$. In our case, $E$ is the potential energy of a mass in an uniform gravitational field, and the Boltzmann factor represents the probability of a particle to be at that altitude. Macroscopically, this becomes the actual pressure of the gas.