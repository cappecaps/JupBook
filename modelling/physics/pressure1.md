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
$$

Where we used equation $\ref{idealgas}$ to express the numerical density in terms of density and temperature. It is important to note that equation (6) now only depends on variables at the altitude $h$ and its proximity (the derivative), and not on the whole gas column that sits above. This mathematical step permits us to reach for a solution, but at the cost of not being able to find the absolute value of the pressure, expect for getting it empirically. Rearranging equation (6) we have: