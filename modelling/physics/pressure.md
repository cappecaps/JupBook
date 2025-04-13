---
title: The atmospheric pressure
subtitle: and how it varies with altitude
short_title: Pressure and altitude
subject: Physics - ab-initio modelling
tags: 
    - modelling
    - atmosphere
    - earth
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


:::{warning} Beware scientists!
:label: abinitio-warning
Here reported is an **_ab initio_ modelling** project. You will find many mistakes, errors, and inaccuracies because I like to do things by myself, often without looking for solutions or explanations. My projects are the opposite of standing on the shoulders of giants. They are like trying to build my own giant, and he's barely alive.
:::

%```{embed} #abinitio-warning
%```

:::{seealso} Preface
:icon: false
I did this project in Desmos back in 2021. For some reason, I've always been obsessed with air pressure, perhaps because I discovered that it can be measured outstandingly precisely with a smartphone. It can even detect altitude changes on the order of meters. So I asked myself: can I derive how the atmospheric pressure varies with altitude, and can I do it as accurately as possible? This way, I could estimate a change in altitude by measuring the change in pressure with my phone. And it seems to work rather nicely. Follow me through it and you'll see.
:::

# Introduction

In this chapter, we will try to find, in an _ab initio_ manner, how atmospheric pressure and altitude are related to each other. 


# Theoretical model


(heading-barometric-formula)=
## The barometric formula

We start from considering an ideal gas. This approximation works well with the Earth’s atmosphere, because it has a sufficiently low density and high temperature. From the ideal gas law, the pressure is given by:

$$
    p=nk_BT
$$ (idealgas)

where $n$ the numerical density of the gas, and $T$ its temperature. To model Earth’s atmosphere, we imagine an infinitely high column of gas subjected to Earth's gravitational field. For now, we assume that Earth is spherically symmetric. The pressure of the gas at a certain altitude $h$ must be equivalent to the pressure exerted by the column of gas above it, due to the gravitational field, i.e.:

$$
    nk_BT(h) = p(h)=\dfrac{F(h)}{A}
$$ (pressure_gravfield)

with $F$ the weight of the gas column above $h$, which is given by the integral along the vertical direction $z$ of the mass density $\rho(z)$ of the gas, times the gravity $g(z)$:

$$
p(h) = \int_h^{+\infty}g(z)\rho(z) dz  = \int_h^{+\infty}g(z)m_0(z)n(z) dz 
$$

where we expressed the mass density $\rho$ as the product of the average mass of the gas particle, $m_0$, and the numerical density, $n$. We can rearrange and convert the equation in its differential form, knowing that $n(+\infty)=0$:

$$
\dfrac{d}{dh}p(h)=-g(h)\,m_0(h)n(h) = -g(h)\,m_0(h)\dfrac{p(h)}{k_BT(h)}
$$(diff_pressure)

Where we used equation {eq}`idealgas` to express the numerical density in terms of density and temperature. Note that equation {eq}`diff_pressure` only depends on the altitude $h$ and its proximity (the derivative), and not on the whole gas column that sits above. This mathematical step permits us to reach for a solution, but at the cost of not being able to find the absolute value of the pressure. Rearranging equation {eq}`diff_pressure`, we have:

$$
\dfrac{\dfrac{d}{dh}p(h)}{p(h)}=-\dfrac{1}{k_B}\dfrac{g(h)m_0(h)}{T(h)}
$$

By integrating in $dz$ from $0$ to $h$ we eventually obtain the general solution of how the pressure varies with the altitude:

$$
p(h)=p(0)\exp\bigg[-\dfrac{1}{k_B}\int_0^h\dfrac{g(z)m_0(z)}{T(z)}dz\bigg]
$$(general_formula)

Where $p(0)$, the pressure at $h=0$, is not known and must be obtained from the observed data. If we choose $h=0$ to be the sea level, then $p(0)$ is the barometric pressure, which in standard conditions is $1013.25\ \mathrm{hPa}$. We will understand how temperature and average mass vary with altitude in the following section. For now, let’s find the simplest solution by assuming that are all variables inside the integral, i.e. composition, temperature, and gravity, are constants. Such approximation is valid close to the Earth's surface level. We then obtain:

$$
p(h)=p(0)\exp\bigg[-\dfrac{m_0gh}{k_BT}\bigg]
$$(barometric_formula)

We recognize $m_0gh$ as the potential energy of a single gas particle, and $k_BT$ as its thermal energy. Wait, what? Equation {eq}`barometric_formula` is called the barometric formula. The exponential term is the Boltzmann factor ($e^{-E/k_BT}$), which, in a canonical ensemble (NVT, our case), represents the probability of the system to be in a state with energy $E$. In our case, $E$ is the potential energy of a mass in a uniform gravitational field, and the Boltzmann factor represents the probability of a particle to be at that altitude. Macroscopically, this becomes the actual pressure of the gas.


## Empirical interlude: Earth's atmosphere 

:::{note}

## Reference atmospheric models

The International Standard Atmosphere (ISA) is a model that describes how the atmospheric parameters change with altitude. It assumes a constant gravitational field, dry air, and it divides the atmosphere in various layers, with different characteristics.

### Geopotential altitude

The vertical distance from the Earth's mean sea level (MSL) is called the **geometric altitude**. In aviation and meteorology, the **geopotential altitude** is used instead, and it is defined as:

Geopotential altitude
: The vertical coordinate referenced to Earth's MSL that represents the work performed when lifting one unit of mass over one unit of length through a hypothetical space in which the acceleration of gravity is assumed to be constant.

As you may know, Earth's gravitational field changes not only with altitude but also with the latitude and, to a minor extent, longitude. The geopotential altitude arises when assuming a constant gravitational field with $g_0=9.80665\ \mathrm{m/s^2}$, the standard gravity at MSL, and it is related to potential energy, $E=mg_0h$. Specifically, a geopotential difference of $1\ \mathrm{m}$ corresponds to a potential energy difference of $9.80665\ \mathrm{J}$. For example, on the North Pole, where $g_{np}>g_0$, a geometric height of $1\ \mathrm{m}$ corresponds to a geopotential altitude that is larger than one meter:
$$
mg_0h_{geop} = mg_{np}h_{geom} \implies h_{geop} = \dfrac{g_{np}}{g_{0}}h_{geom} > h_{geom}
$$ 
More rigorously, denoting Earth's radius with $R$, the geopotential height $h$ is related to the geometric height $z$ according to the formula:
$$
h = \dfrac{R}{R+z}\,z
$$
The geopotential altitude is the one that we used in the section [](#heading-barometric-formula), where we assumed $g$ constant. We will get rid of this approximation in a later section. 

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


(subheading-chemical-composition)=
### Chemical composition
The composition of the dry atmosphere is kindly provided by the [NOAA](https://www.noaa.gov/jetstream/atmosphere). However, the molar fractions they provide sum to more than one, due to experimental error (as also [wikipedia](https://en.wikipedia.org/wiki/Atmosphere_of_Earth#Composition) reports). While looking for the most precise values, I noticed an incongruence in the reported amount of $\mathrm{CO_2}$. I quickly realized a shocking fact: the atmospheric concentration of $\mathrm{CO_2}$ is rising so quickly that most values are now outdated. For example, [engineering toolbox](https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html) uses a value of $f_{CO_2}=0.033\%$, which is the fraction from circa 50 years ago (in the '70s). The value is now (2025) $0.042\%$, giving an outstanding 27% increase. This makes me wonder whether NOAA takes into account the change in fractional concentrations due to $\mathrm{CO_2}$ emissions and $\mathrm{O_2}$ depletion. A strong hint is that by substituting the present $\mathrm{CO_2}$ concentration with $0.033\%$ in NOAA's value, the sum magically becomes $1$. For this reason, I took NOAA's value and assumed that the $0.011\%$ increase is due to combustion, and it substitutes $\mathrm{O_2}$ molecules
$$
f_{CO_2}(2004) = 0.033\% &\longrightarrow f_{CO_2}(2025) = 0.042\% \\
f_{O_2}(2004) = 20.946\% &\longrightarrow f_{O_2}(2025) = 20.937\%
$$
This, however, results from the crude approximation that the additional carbon dioxide has been produced by a stoichiometric reaction between atmospheric oxygen and carbon, while all other species stay constant. In reality, $\mathrm{CO}_2 production is less "efficient", because of the formation of water, among other compounds. 

```{table} Chemical composition of Earth's dry atmosphere, modified data from [NOAA](https://www.noaa.gov/jetstream/atmosphere) to sum to one. In the ideal gas approximation, mole fractions and volume fractions are equivalent. 
:label: composition
:align: center
| Element | mole fraction |
| --- | --- |
| N{sub}`2` | 78.084% |
| O{sub}`2` | 20.937% |
| Ar | 0.934% |
| CO{sub}`2` | 0.042% |
| Ne | 18.182 ppm |
| He | 5.24 ppm |
| CH{sub}`4` | 1.92 ppm |
| Kr | 1.14 ppm |
| H{sub}`2` | 0.55 ppm |
| N{sub}`2`O | 0.33 ppm |
| CO | 0.10 ppm |
| Xe | 0.09 ppm |
| O{sub}`3`| 0.07 ppm |
| NO{sub}`2` | 0.02 ppm |
| I{sub}`2` | 0.01 ppm |
| other | traces | 
```
Which gives an average molar mass $m = 28.9656\ \mathrm{g/mol}$. The value is higher compared to $28.9647\ \mathrm{g/mol}$ encountered online, obtained from a lower level of atmospheric $\mathrm{CO_2}$ of 332 ppm. However, my calculation might be too rough, or blatantly wrong; I'm not too sure. See the [Scripps FAQ page](https://scrippso2.ucsd.edu/faq.html) for additional information.

Of course, the atmosphere is never dry. Local concentrations of water vapor range from 0% to 4%. The molar fractions, including humidity, are simply given by

$$
f_{A}^{(hum)}=\big(1-f_{H_2O}\big)\cdot f_{A}^{(dry)}
$$(f_dry_to_humid)

with $A$ any species.
```{tip} Quick proof
:class: dropdown
In the dry case it holds
$$
\sum_A f_{A}^{(dry)} = 1
$$
while for humid air
$$
f_{H_2O} + \sum_A f_{A}^{(hum)} = 1
$$
When water vapor is added to the dry air, all $f_{A}^{(dry)}$ must decrease by the same multiplying factor $k$:
$$
f_{H_2O} + k\sum_A f_{A}^{(dry)} &= 1 \\
f_{H_2O} + k &= 1  \\
k &= 1-f_{H_2O} 
$$
which gives equation {eq}`f_dry_to_humid`.
```

According to the [NRLMSIS empirical model](https://swx-trec.com/msis/?lz=N4Igtg9gJgpgNiAXCYAdEUCGAXG2CWYM6iAjAOykAsArAEykBsADK6wL4gA0ImcB2AK6wkKdHBz4hsEqWZdxEAHYBzKcOJI6zTjwDOggE4AzTAGMYotL37qZSOTu4gAbpkP5MAIziXk6ADkAeXRnNw9vXz1RAG10AEFDdAUQAAlk9FTNFICMkAC6POC8kO50IMKykABZTD09PIAVGDAABxhDHCNNAF0QdiA), Earth's atmospheric composition remains rather constant up to $h\approx 80\ \mathrm{km}$. This is due to diffusion and turbulent convection, resulting in a well-mixing of the various gases. This is true up to a certain altitude

:::


## Temperature change

Let's start from the barometric formula {eq}`barometric_formula` and remove approximations one by one to finally arrive at the general formula {eq}`general_formula`. First, we introduce the empirical lapse rates that we learned above. We thus leave only the temperature term in the integral:

$$
p(h)=p(0)\exp\bigg[-\dfrac{gm}{R}\int_0^h\dfrac{1}{T(z)}dz\bigg]
$$(with_lapse)
Where we used the average molar mass $m$ and the gas constant $R$ instead of $m_0$ and $k_B$.
Now the temperature can be written as the general expression:

$$
T(h)= T(h_{i}) - \Gamma_i (h-h_{i})
$$(T_fromlapse)

with $i$ the atmospheric layer in which $h$ lies, and $h_{i}$ the base altitude of the layer $i$. This however can be easily calculated with a simple code. Let's import the packages.

```{code-cell} ipython
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
```

and define a function that returns the temperature at a certain altitude, using equation {eq}`T_fromlapse` for each layer.

```{code-cell} ipython
def temperature(altitude):
    # Define base altitudes and temperatures for each layer
    T0= 288.15  # MSL standard temperature in Kelvin
    base_altitudes = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852]
    lapse_rates = [0.0065, 0.0, -0.001, -0.0028, 0.0, 0.0028, 0.002, 0.0]

    temperature = T0  

    # Find the layer corresponding to the altitude
    for i in range(len(base_altitudes) - 1):
        if altitude <= base_altitudes[i+1]:
            return temperature - lapse_rates[i] * (altitude - base_altitudes[i])
        else:
            temperature = temperature - lapse_rates[i] * (base_altitudes[i+1]-base_altitudes[i])

    # If altitude is above the highest defined layer
    return temperature
```


```{code-cell} ipython
:tags: ["hide-input"]
# Generate altitudes from 0 to 85 km and calculate temperature
altitudes = np.linspace(0, 85, 500)
temperatures = [temperature(1000*alt) for alt in altitudes]

plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(7, 3))
plt.plot(altitudes, temperatures, color="darkred",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Temperature (K)")
plt.grid(True)
plt.show()
```

Now we can define a function to calculate the pressure, which contains the evaluation of the integral of $1/T(h)$ in $dh$.


```{code-cell} ipython
def pressure(altitude):
    R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
    g = 9.80665  # Standard gravity in m/s^2
    m = 28.9656e-3  # Molar mass of air in kg/mol
    P0 = 101325  # MSL standard atmospheric pressure in Pa

    integral, err = quad(lambda h: 1 / temperature(h), 0, altitude, limit=100, points=[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    pressure = P0 * np.exp(-m*g/R * integral)

    return pressure
```

This results in a very good estimation of the pressure with altitude, especially at lower altitudes. 

```{code-cell} ipython
:tags: ["hide-input"]
altitudes = np.linspace(0, 85, 500)
pressures = [pressure(1000*alt) for alt in altitudes]
plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, pressures, color="black",lw=2)
plt.xlabel("Altitude (km)")
plt.ylabel("Pressure (Pa)")
plt.grid(True)
plt.show()
```

We can also substitute the lapse rate with the actual data from [atmospheric sounding](wiki:Atmospheric_sounding), if we want to calculate the change in pressure more accurately. For now, we are happy with the standard lapse rate. 

## Humidity



## Gravitational field


## Earth as a spheroid


## Centrifugal force



## Altitude from pressure

%Take equation T_from lapse and solve it with T = T(0)-Lz or T(h)-Lh-Lz and rearrange



# Empirical model

## The NRLMSIS model

The value that we calculated in [](subheading-chemical-composition) refers to the global average of the atmospheric composition. We are now interested in how such composition changes with altitude. We know that turbulence and diffusion make the atmospheric composition rather constant up to 80 km ([](#composition-altitude)). The vertical profile of carbon dioxide, one of the heaviest molecules in the air, starts decreasing from an altitude of 60 km. 

:::{figure} https://upload.wikimedia.org/wikipedia/commons/b/bd/Chemical_composition_of_atmosphere_accordig_to_altitude.png
:label: composition-altitude
:align: center
:w: 400px

Atmospheric composition vs altitude

:::

:::{figure} ../../images/CO2_altitude.jpg
:label: CO2-altitude
:align: center
:w: 400px

Altitude profile of $\mathrm{CO}_2$ molar fraction, in ppm. From [Brown et al. (2024)](https://doi.org/10.1029/2024JA032659).

:::


The [NRLMSIS empirical model](https://swx-trec.com/msis/?lz=N4Igtg9gJgpgNiAXCYAdEUCGAXG2CWYM6iAjAOykAsArAEykBsADK6wL4gA0ImcB2AK6wkKdHBz4hsEqWZdxEAHYBzKcOJI6zTjwDOggE4AzTAGMYotL37qZSOTu4gAbpkP5MAIziXk6ADkAeXRnNw9vXz1RAG10AEFDdAUQAAlk9FTNFICMkAC6POC8kO50IMKykABZTD09PIAVGDAABxhDHCNNAF0QdiA) provides an enormous quantity of data , that ...

%The model that we have built so far is actually rather accurate, especially if we want to use it at the surface of the Earth, where the atmospheric composition, the gravitational field, and 


## Using sounding data


## The tool