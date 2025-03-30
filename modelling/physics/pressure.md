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
    p=nk_BT
$$ (idealgas)

where $n$ the numerical density of the gas, and $T$ its temperature. To model Earth’s atmosphere, we imagine an infinitely high column of gas subjected to Earth's gravitational field. For now, we assume that Earth is spherically symmetric. The pressure of the gas at certain altitude $h$ must be equivalent to the pressure exerted by the column of gas above it, due to the gravitational field, i.e.:

$$
    nk_BT(h) = p(h)=\dfrac{F(h)}{A}
$$ (pressure_gravfield)

with $F$ the weigth of the gas column above $h$, which is given by the integral along the vertical direction $z$ of the mass density $\rho(z)$ of the gas, times the gravity $g(z)$:

$$
p(h) = \int_h^{+\infty}g(z)\rho(z) dz  = \int_h^{+\infty}g(z)m_0(z)n(z) dz 
$$

where we expressed the mass density $\rho$ as the product between the average mass of the gas particle, $m_0$, and the numerical density, $n$. We can rearrange and convert the equation in its differential form, knowing that $n(+\infty)=0$:

$$
\dfrac{d}{dh}p(h)=-g(h)\,m_0(h)n(h) = -g(h)\,m_0(h)\dfrac{p(h)}{k_BT(h)}
$$(diff_pressure)

Where we used equation {eq}`idealgas` to express the numerical density in terms of density and temperature. It is important to note that equation {eq}`diff_pressure` only depends on variables at the altitude $h$ and its proximity (the derivative), and not on the whole gas column that sits above. This mathematical step permits us to reach for a solution, but at the cost of not being able to find the absolute value of the pressure.  Rearranging equation {eq}`diff_pressure` we have:

$$
\dfrac{\dfrac{d}{dh}p(h)}{p(h)}=-\dfrac{1}{k_B}\dfrac{g(h)m_0(h)}{T(h)}
$$

By integrating in $dz$ from $0$ to $h$ we eventually obtain the general solution of how the pressure varies with the altitude:

$$
p(h)=p(0)\exp\bigg[-\dfrac{1}{k_B}\int_0^h\dfrac{g(z)m_0(z)}{T(z)}dz\bigg]
$$(general_formula)

Where $p(0)$, the pressure at $h=0$, is not known, and must be obtained from the observed data. If we choose $h=0$ to be the sea level, then $p(0)$ is the barometric pressure, which in standard conditions is $1013.25\ \mathrm{hPa}$. We will understand how temperature and average mass vary with altitude in the following section. For now, let’s find the simplest solution by assuming that are all variables inside the integral, i.e. composition, temperature, and gravity, are constants. Such approximation is valid close to the Earth's surface level. We then obtain:

$$
p(h)=p(0)\exp\bigg[-\dfrac{m_0gh}{k_BT}\bigg]
$$(barometric_formula)

We recognize $m_0gh$ as the potential energy of a single gas particle, and $k_BT$ its thermal energy. Wait, what? Equation {eq}`barometric_formula` is called barometric formula. The exponential term is the Boltzmann factor ($e^{-E/k_BT}$), which, in a canonical ensemble (NVT, our case), represents the probability of system to be in a state with energy $E$. In our case, $E$ is the potential energy of a mass in an uniform gravitational field, and the Boltzmann factor represents the probability of a particle to be at that altitude. Macroscopically, this becomes the actual pressure of the gas.


## Adding lapse rate

:::{note}
:class: dropdown

## International Standard Atmosphere

The International Standard Atmosphere (ISA) is a model that describes how the atmospheric parameters change with altitude. It assumes a constant gravitational field, dry air, and it divides the atmosphere is various layers, with different characteristics.

### Geopotential altitude

The vertical distance from the Earth's mean sea level (MSL) is called the **geometric altitude**. In aviation and meteorology, the **geopotential altitude** is used instead, and it is defined as:

Geopotential altitude
: The vertical coordinate referenced to Earth's MSL that represents the work performed when lifting one unit of mass over one unit of length through a hypothetical space in which the acceleration of gravity is assumed constant.

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
The atmospheric temperature depends on many factors, such as irradiation from Earth's surface, convection, chemical reactions, and interaction with high-energy photons from the Sun. The variation of temperature with altitude is called **lapse rate**, $\Gamma$:
$$
\Gamma = -\dfrac{dT}{dh}
$$
The dry air approximation gives the dry adiabatic lapse rate (DALR. $\Gamma_d$):
$$
\Gamma_d = \dfrac{g}{c_p} = 9.8\ ^\circ\mathrm{C/km}
$$
which is valid only at the vicinity of Earth's surface. The ISA provides a set of empirical lapse rates for each atmospheric layer. Starting from a temperature of 15°C (288.15 K) at MSL:
1. **Troposphere** (0-11 km): 6.5 °C/km
2. **Tropopause** (11-20 km): 0.0 °C/km
3. **Stratosphere** (20-32 km): 1.0 °C/km
4. **Stratosphere** (32-47 km): -2.8 °C/km
5. **Stratopause** (47-51 km): 0.0 °C/km
6. **Mesosphere** (51-71 km): 2.8 °C/km
7. **Mesosphere** (71-84.852 km): 2.0 °C/km

There are other layers above, but can be ignored since the atmosphere is extremely rarefied. The ranges are given in geopotential altitude.

### Chemical composition
According to the NRLMSIS empirical model, Earth's atmospheric composition remains rather constant up to $h=80\ \mathrm{km}$. 

:::

Let's start from the barometric formula {eq}`barometric_formula` and remove approximations one by one to finally arrive at the general formula {eq}`general_formula`. First, we introduce the empirical lapse rates that we learned above. We thus leave only the temperature term in the integral:

$$
p(h)=p(0)\exp\bigg[-\dfrac{gm_0}{k_B}\int_0^h\dfrac{1}{T(z)}dz\bigg]
$$(with_lapse)

Now the temperature can be written as the general expression:

$$
T(h)= T(h_{i}) - \Gamma_i (h-h_{i})
$$

with $i$ the atmospheric layer in which $h$ lies, and $h_{i}$ the base altitude of the layer $i$, and that:

$$
T(h_{i}) = T(0) - \sum_{j=0}^{i-1} \Gamma_j (h_{j+1}-h_{j})
$$


This can be easily calcuated with a simple code. Let's import the packages.

```{literalinclude} pressure.py
:start-line: 0
:end-line: 1
:lineno-match:
```

```{literalinclude} pressure.py
:start-line: 3
:end-line: 19
:lineno-match:
```


