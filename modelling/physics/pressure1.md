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

We start from considering an ideal gas. This approximation works well with the Earth’s atmosphere, because it has a sufficiently low density and high temperature. From the ideal gas law, the pressure is given by:

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

Where we used equation {eq}`idealgas` to express the numerical density in terms of density and temperature. It is important to note that equation (6) now only depends on variables at the altitude $h$ and its proximity (the derivative), and not on the whole gas column that sits above. This mathematical step permits us to reach for a solution, but at the cost of not being able to find the absolute value of the pressure, expect for getting it empirically. Rearranging equation {eq}`diff_pressure` we have:

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

## International Standard Atmosphere
The International Standard Atmosphere (ISA) is model that describes how the atmospheric parameters change with altitude. It assumes a constant gravitational field, dry air, and it divides the atmosphere is various layers, with different characteristics.

### Geopotential altitude

Before moving on to our modelling, we need to first understand how temperature and atmospheric composition vary with altitude. Unfortunately, things are not as easy as one would tell, and digging deeper made me realize that even the concept of altitude is not trivial. The vertical distance from the Earth's mean sea level (MSL) is called the **geometric altitude**. In aviation and meteorology, the **geopotential altitude** is used instead, and it is defined as:

Geopotential altitude
: vertical coordinate referenced to Earth's MSL that represents the work performed when lifting one unit of mass over one unit of length through a hypothetical space in which the acceleration of gravity is assumed constant.

As you may know, Earth's gravitational field changes not only with altitude, but also with the latitude and, to a minor extent, longitude. The geopotential altitude arises when assuming a constant gravitational field with $g_0=9.80665\ \mathrm{m/s^2}$, the standard gravity at MSL, and it is related to potential energy, $E=mg_0h$. Specifically, a geopotential difference of $1\ \mathrm{m}$ corresponds to a potential energy difference of $9.80665\ \mathrm{J}$. For example, on the North Pole, where $g_{np}>g_0$, a geometric height of $1\ \mathrm{m}$ corresponds to a geopotential altitude that is larger than one meter:
$$
mg_0h_{geop} = mg_{np}h_{geom} \implies h_{geop} = \dfrac{g_{np}}{g_{0}}h_{geom} > h_{geom}
$$ 
More rigorously, denoting Earth's radius with $R$, the geopotential height $h$ is related to the geometric height $z$ according to the formula:
$$
h = \dfrac{R}{R+z}\,z
$$
The geopotential altitude is the one that we used in the section [](#the-barometric-formula), where we assumed $g$ constant. We will get rid of this approximation in a later section. 

### Temperature
The variation of temperature with altitude is called **lapse rate**, $\Gamma$:
$$
\Gamma = -\dfrac{dT}{dh}
$$
Lapse rates at the different atmospheric layers are usually given in °C/km. The dry air approximation gives the dry adiabatic lapse rate (DALR. $\Gamma_d$):
$$
\Gamma_d = \dfrac{g}{c_p} = 9.8\ ^\circ\mathrm{C/km}
$$
However, this is not a formula I obtained nor I understand, so that it is outside the scope of this blog. The key point here is that the ISA defined a set of lapse rates for each atmospheric layer. Starting from a temperature of 15°C at the MSL:
1. **Troposphere** (0-11 km): 6.5 °C/km
2. **Tropopause** (11-20 km): 0.0 °C/km
3. **Stratosphere** (20-32 km): 1.0 °C/km
4. **Stratosphere** (32-47 km): -2.8 °C/km
5. **Stratopause** (47-51 km): 0.0 °C/km
6. **Mesosphere** (51-71 km): 2.8 °C/km
7. **Mesosphere** (71-84.852 km): 2.0 °C/km

There are other layers above, but can be ignored due to the atmosphere being extremely rarefied. 

### Chemical composition