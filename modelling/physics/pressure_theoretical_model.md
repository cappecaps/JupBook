---
title: The atmospheric pressure
subtitle: and how it varies with altitude
short_title: Pressure and altitude
subject: Physics - ab-initio modelling
tags: 
    - modelling
    - atmosphere
    - earth
numbering: 
    headings: true
    heading_4: false
    figure: true
    table: true
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


%TODO:
% - Rewrite "pressure at sea level" and redefine \chi to also include gravity change and such
%   - I need to find f_w such that avg f_w = 0.4%. But considering that f_w(h)=0 for h>20km
% - Take care that I moved some sections 
% - Rewrite almost everything to have it more fluid
%


:::{warning} Beware scientists!
:label: abinitio-warning
Here reported is an **_ab initio_ modelling** project. You will find many mistakes, errors, and inaccuracies because I like to do things by myself, often without looking for solutions or explanations. My projects are the opposite of standing on the shoulders of giants. They are like trying to build my own giant, and he's barely alive.
:::

%```{embed} #abinitio-warning
%```

:::{seealso} Preface
:icon: false
I did this project in Desmos back in 2021. For some reason, I've always been obsessed with air pressure, perhaps because I discovered that it can be measured outstandingly precisely with a smartphone. It can even detect altitude changes on the order of meters. So I asked myself: can I derive how the atmospheric pressure varies with altitude, and can I do it as accurately as possible? This way, I could estimate a change in altitude by measuring the change in pressure with my phone. And it works rather nicely!
:::

## Introduction

In this chapter, we will try to find, in an _ab initio_ manner, how atmospheric pressure and altitude are related to each other. 


## Theoretical model


(heading-barometric-formula)=
### The barometric formula

We start from considering an ideal gas. This approximation works well with the Earth’s atmosphere, because it has a sufficiently low density and high temperature. From the ideal gas law, the pressure is given by:

$$
    p=nk_BT
$$ (idealgas)

where $n$ the numerical density of the gas, and $T$ its temperature. To model Earth’s atmosphere, we imagine an infinitely high column of gas subjected to Earth's gravitational field. For now, we assume that Earth is spherically symmetric. We also assume that the mass of the atmosphere is negligible, so that there is no self-interaction, which would make the problem not analytically solvable. The pressure of the gas at a certain altitude $h$ must be equivalent to the pressure exerted by the column of gas above it, due to the gravitational field, i.e.:

$$
    nk_BT(h) = p(h)=\dfrac{F(h)}{A}
$$ (pressure_gravfield)

with $F$ the weight of the gas column above $h$. The pressure of the column can also be calculated as the integral along the vertical direction $z$ of the mass density $\rho(z)$ of the gas, times the gravitational acceleration $g(z)$:

$$
p(h) = \int_h^{+\infty}g(z)\rho(z) dz  = \int_h^{+\infty}g(z)m_0(z)n(z) dz 
$$(pressure_h)

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

Interestingly, to compute the pressure at a certain altitude $h$ we just need to know how the integrand varies below that point, and not above. This is because the formula assumes we know $p(0)$, the pressure at $h=0$, so that the information of the air column above is implicitly contained there. If we choose $h=0$ to be the sea level, then $p(0)$ is the barometric pressure, which in standard conditions is $1013.25\ \mathrm{hPa}$. In reality, the sea-level pressure varies continuously and it must be measured. We will understand how temperature and average mass vary with altitude in the following section. For now, let’s find the simplest solution by assuming that are all variables inside the integral, i.e. composition, temperature, and gravity, are constants. Such approximation is valid close to the Earth's surface level. We then obtain:

$$
p_{bar}(h)=p(0)\exp\bigg[-\dfrac{m_0g_0h}{k_BT}\bigg]
$$(barometric_formula)

We recognize $m_0g_0h$ as the potential energy of a single gas particle, and $k_BT$ as its thermal energy. Wait, what? Equation {eq}`barometric_formula` is called the **barometric formula**. The exponential term is the Boltzmann factor ($e^{-E/k_BT}$), which, in a canonical ensemble (NVT, our case), represents the probability of the system to be in a state with energy $E$. In our case, $E$ is the potential energy of a mass in a uniform gravitational field, and the Boltzmann factor represents the probability of a particle to be at that altitude. Macroscopically, this becomes the actual pressure of the gas.

From the formula we can infer how the pressure profile changes with different parameters. A smaller particle mass makes the gas less attracted to the surface, and thus more spread toward higher altitudes. The same effect is achieved with higher temperatures, because of the higher kinetic energy of the molecules. Vice-versa, larger masses and lower temperatures make the gas more "compressed" at the surface.


### Empirical interlude: Earth's atmosphere 

:::{note}

### Reference atmospheric models

The International Standard Atmosphere (ISA) is a model that describes how the atmospheric parameters change with altitude. It assumes a constant gravitational field, dry air, and it divides the atmosphere in various layers, with different characteristics.

#### Geopotential altitude

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

(ISA_temp)=
#### Temperature
The atmospheric temperature depends on many factors, such as irradiation from Earth's surface, convection, chemical reactions, and interaction with high-energy photons from the Sun. The variation of temperature with altitude is called **lapse rate**, $\Gamma$:
$$
\Gamma = -\dfrac{dT}{dh}
$$
The dry air approximation gives the dry adiabatic lapse rate (DALR. $\Gamma_d$):
$$
\Gamma_d = \dfrac{g_0}{c_p} = 9.8\ ^\circ\mathrm{C/km}
$$
which is valid only at the vicinity of Earth's surface. The ISA provides a set of empirical lapse rates for each atmospheric layer. 

```{table} Atmospheric layer data provided by ISA.
:label: tab:lapse-rates
:align: center
| Layer | h range (km) | $\Gamma$ (K/km) | Base T (K) |
| --- | --- | --- | --- |
| **Troposphere** | 0 - 11 | +6.5 | 288.15 |
| **Tropopause** | 11 - 20 | 0.0 | 216.65 |
| **Stratosphere** | 20 - 32 | -1.0 | 216.65 |
| **Stratosphere** | 32 - 47 | -2.8 | 228.65 |
| **Stratopause** | 47 - 51 | 0.0 | 270.65 |
| **Mesosphere** | 51 - 71 | +2.8 | 270.65 |
| **Mesosphere** | 71 - 84.852 | +2.0 | 214.65 |
| **Mesopause** | 84.852 -  | 0.0 | 186.946 |
```

There are other layers above, but can be ignored for now since the atmosphere is extremely rarefied there. The ranges are given in geopotential altitude.


(subheading-chemical-composition)=
#### Chemical composition
The composition of dry atmosphere is kindly provided by [NOAA](https://www.noaa.gov/jetstream/atmosphere). Unfortunately, the molar fractions they provide sum to a number greater than one, due to round-off and experimental errors (see [wikipedia](https://en.wikipedia.org/wiki/Atmosphere_of_Earth#Composition) reports). While I was looking for more accurate values, I noticed an incongruence in the reported amount of $\mathrm{CO_2}$. I quickly realized a shocking fact: the atmospheric concentration of $\mathrm{CO_2}$ is rising so quickly that most values are now outdated. For example, [engineering toolbox](https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html) uses a value of $f_{CO_2}=0.033\%$, which is the fraction from circa 50 years ago (in the '70s). The value is now (2025) $0.042\%$, giving an outstanding 27% increase. This makes me wonder whether NOAA takes into account the change in fractional concentrations due to $\mathrm{CO_2}$ emissions and $\mathrm{O_2}$ depletion. A strong hint is that by substituting the present $\mathrm{CO_2}$ concentration with $0.033\%$ in NOAA's value, the sum magically becomes $1$. For this reason, I took NOAA's value and assumed that the $0.011\%$ increase is due to combustion, and it substitutes $\mathrm{O_2}$ molecules

$$
\begin{aligned}
f_{CO_2}(2004) = 0.033\% &\longrightarrow f_{CO_2}(2025) = 0.042\% \\
f_{O_2}(2004) = 20.946\% &\longrightarrow f_{O_2}(2025) = 20.937\%
\end{aligned}
$$

This, however, results from the crude approximation that the additional carbon dioxide has been produced by a stoichiometric reaction between atmospheric oxygen and carbon, while all other species stay constant. In reality, $\mathrm{CO}_2$ production is less "efficient", because of the formation of water, among other compounds. 

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
Which gives an average molar mass for dry air $m_d = 28.9656\ \mathrm{g/mol}$. The value is higher compared to $28.9647\ \mathrm{g/mol}$ encountered online, obtained from a lower level of atmospheric $\mathrm{CO_2}$ of 332 ppm. However, my calculation might be too rough, or blatantly wrong; I'm not too sure. See the [Scripps FAQ page](https://scrippso2.ucsd.edu/faq.html) for additional information.

Of course, the atmosphere is never dry. Water vapor makes circa 0.25\% of the atmospheric mass (i.e. 0.40\% by volume) and its local concentration ranges within 0-4%. The molar fractions, including humidity, are simply given by

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
\begin{aligned}
f_{H_2O} + k\sum_A f_{A}^{(dry)} &= 1 \\
f_{H_2O} + k &= 1  \\
k &= 1-f_{H_2O} 
\end{aligned}
$$
which gives equation {eq}`f_dry_to_humid`.
```

According to the [NRLMSIS empirical model](https://swx-trec.com/msis/?lz=N4Igtg9gJgpgNiAXCYAdEUCGAXG2CWYM6iAjAOykAsArAEykBsADK6wL4gA0ImcB2AK6wkKdHBz4hsEqWZdxEAHYBzKcOJI6zTjwDOggE4AzTAGMYotL37qZSOTu4gAbpkP5MAIziXk6ADkAeXRnNw9vXz1RAG10AEFDdAUQAAlk9FTNFICMkAC6POC8kO50IMKykABZTD09PIAVGDAABxhDHCNNAF0QdiA), Earth's atmospheric composition remains rather constant up to $h\approx 80\ \mathrm{km}$. This is due to diffusion and turbulent convection, resulting in a well-mixing of the various gases. This is true up to a certain altitude

:::


### Lapse rates

Let's start from the barometric formula {eq}`barometric_formula` and remove approximations one by one to finally arrive at the general formula {eq}`general_formula`. First, we introduce the empirical lapse rates that we learned above. We thus leave only the temperature term in the integral:

$$
p_{dry}(h)=p(0)\exp\bigg[-\dfrac{g_0m_d}{R}\int_0^h\dfrac{1}{T(z)}dz\bigg]
$$(with_lapse)
Where we used the average molar mass of dry air $m_d$ and the gas constant $R$ instead of $m_0$ and $k_B$.
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

define the global variables that we're going to use across this article

```{code-cell} ipython
    R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
    g0 = 9.80665  # Standard gravity in m/s^2
    m_dry = 28.9656e-3  # Molar mass of air in kg/mol
    m_water = 18.01528e-3  # Molar mass of water in kg/mol
    T0 = 288.15  # MSL standard temperature in Kelvin
    p0 = 101325  # MSL standard atmospheric pressure in Pa
```

and define a function that returns the temperature at a certain altitude, using equation {eq}`T_fromlapse` for each layer.

```{code-cell} ipython
def ISA_temperature(h):
    # Define base altitudes and temperatures for each layer
    base_altitudes = [0, 11, 20, 32, 47, 51, 71, 84.852]
    lapse_rates = [6.5, 0.0, -1.0, -2.8, 0.0, 2.8, 2.0, 0.0]

    temperature = T0

    # Find the layer corresponding to the altitude
    for i in range(len(base_altitudes) - 1):
        if h <= base_altitudes[i+1]:
            return temperature - lapse_rates[i] * (h - base_altitudes[i])
        else:
            temperature = temperature - lapse_rates[i] * (base_altitudes[i+1]-base_altitudes[i])

    # If altitude is above the highest defined layer
    return temperature
```


```{code-cell} ipython
:tags: ["hide-input"]
# Generate altitudes from 0 to 85 km and calculate temperature
altitudes = np.linspace(0, 85, 100)
temperatures = [ISA_temperature(alt) for alt in altitudes]

plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(7, 3))
plt.plot(altitudes, temperatures, color="darkred",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Temperature (K)")
plt.grid(True)
plt.show()
```

We can define a function that calculates the pressure in dry air, and compare it with the barometric approximation.


```{code-cell} ipython
def pressure_barometric(h):
    pressure = p0 * np.exp(-m_dry * g0 * h / (R * T0))

    return pressure

def pressure_dry(h):
    integral, err = quad(lambda h: 1 / ISA_temperature(h), 0, h, limit=100, points=[0, 11, 20, 32, 47, 51, 71, 84.852])
    pressure = p0 * np.exp(-m_dry * g0 / R * 1000 * integral )

    return pressure
```


```{code-cell} ipython
:tags: ["hide-input"]
altitudes = np.linspace(0, 85, 500)     # altitude array in km
pressure_dry_arr = np.array([pressure_dry(alt)/100 for alt in altitudes])   # divided by 100 to return hPa
pressure_barom_arr = np.array([pressure_barometric(alt)/100 for alt in altitudes])
plt.rcParams.update({'font.size': 10})
fig, ax1 = plt.subplots(figsize=(7, 3))
ax2 = ax1.twinx()
ax1.plot(altitudes, pressure_dry_arr, color='darkred',lw=2,label='with Γ ($p_{dry}$)')
ax1.plot(altitudes, pressure_barom_arr, color='c',lw=2, label='barometric approx. ($p_{bar}$)')
ax2.plot(altitudes, pressure_barom_arr - pressure_dry_arr,linestyle='dashed', color='r', lw=2, label=r'$p_{bar}-p_{dry}$')
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel("Pressure (hPa)")
ax2.set_ylabel("Pressure difference (hPa)")
ax1.set_xlim(0,60)
ax1.legend()
ax2.legend(loc=7)
plt.show()
```

Our current model already results in a very good estimation of the pressure, especially close to the surface. As a recap, our model has so far been built within the following approximations:
- Ideal gas law
- Spherical Earth, no spin
- Constant gravitational field
- Constant atmospheric composition, dry air
- Empirical lapse rates

Let's keep improving our model by removing most of them one by one. The most impactful one is arguably dry air: on Earth, air is never completely dry, but some water vapor is mixed with the other gases. In practice, addition of water vapor affects the composition, namely the molar mass $m$ in our model.



:::{note}
#### Sea-level pressure reduction

Atmospheric pressure is always reported at the MSL. When a weather station at a certain altitude measures the local pressure, that value is then reduced to the sea level. In other words, the station must estimate the pressure that would be measured if someone digged down to the sea level. This is called **sea-level pressure reduction**, and since no air exists below ground, it is purely hypothetical, so that many assumptions need to be made. For example, how temperature and humidity would vary going down cannot be properly defined.

The sea-level pressure reduction is carried out by means of the <wiki:hypsometric_equation>:

$$
p_{MSL} = p_{obs}\cdot \exp\bigg[\frac{m_mg_0h}{R\overline{T}}\bigg]
$$

where $h$ is the altitude in which the pressure $p_{obs}$ is measured, $p_{MSL}$ is the pressure at the MSL, and $\overline{T}$ is the mean temperature of the (moist) air with molar mass $m_m$ (I omitted to explain the <wiki:virtual_temperature>, but the formula is equivalent). The hypsometric equation is basically the barometric formula (equation {eq}`barometric_formula`), but with an average temperature across the vertical distance. The average temperature can simply be computed using the standard lapse rate $\Gamma$, and a 12-hour average surface temperature, an attempt to exclude the effect of the irradiation of surface of the Earth:
$$
    \overline{T} \approx \dfrac{T_{obs}(t) + T_{obs}(t-12h)}{2} - \Gamma \cdot \dfrac{h}{2}
$$
Additional refinements can be employed, but they are often empirical and specific to the station site. Mind that this pressure reduction is fictitious, and it is better to be consistent (across the world) rather than accurate. And I don't like this.

:::


### Humidity
The amount of water vapor in the air is usually measured in relative humidity (RH or $\phi$), which is the fraction of the water vapor in the air relative to the "maximum" potential at that temperature. 


#### From relative humidity

:::{note}

#### Relative humidity and water vapor pressure

**Relative humidity** $\phi$ is defined as the ratio between the measured partial pressure of water $p_w$ and its equilibrium (saturation) vapor pressure $p_{vap,w}$:

$$
\varphi = \dfrac{p_w}{p_{vap,w}}
$$(relative_humidity)

Note that these terms assume a different meaning in physics and in meteorology. In physics, the **vapor pressure** of a substance is the partial pressure of the gas phase in *equilibrium* with the liquid phase. In meteorology, however, **vapor pressure** refers to the measured partial pressure of water vapor, even when not in equilibrium, whereas the **saturation vapor pressure** is the actual vapor pressure of water, i.e., in equilibrium conditions. This nomenclature comes from the erroneous idea of air dissolving water vapor, eventually reaching a saturation limit. In reality, the vapor pressure of a substance only depends on the liquid-phase temperature, in first approximation. Here, we will try to use the correct physical definitions, albeit minding such incongruences.

The dependence of the vapor pressure of a substance on temperature can be estimated from the Clausius-Clapeyron equation knowing its boiling point $T_b$ at standard pressure $p^\circ$  ($p^\circ$ = 1 atm = 101325 Pa, and $T_b$ = 99.97°C for water)

$$
p_{vap}(T) = p^\circ \cdot \exp\bigg[{-\dfrac{\Delta_{vap} H}{R}\bigg(\dfrac{1}{T}-\dfrac{1}{T_b}\bigg)}\bigg]
$$

With $\Delta_{vap} H$ the enthalpy of vaporization. However, this formula lacks of the desired accuracy, because of the numerous approximations that led to it: ideal gas, constant $\Delta_{vap} H$ with temperature, no volume change, etc. This is particularly true for water, which deviates from the ideality due to the strong intermolecular interactions it can establish. For this reason, it is common to use empirical formula. Here, we are going to use the <wiki:Tetens_equation>

$$
p_{vap,w}(T) = a \exp\bigg[\dfrac{bT}{T+c}\bigg]
$$(tetens)

with $a=610.78,b=17.27,c=237.3$ for $p_{vap,w}:$ °C $\rightarrow$ Pa

```{code-cell} ipython
def vapor_pressure(T):
    # Tetens equation
    a = 610.78
    b = 17.27
    c = 237.3  
    p_vap = a * np.exp((b * (T - 273.15)) / (T - 273.15 + c))   
    return p_vap
```

```{code-cell} ipython
:tags: ["hide-input"]
T_array = np.linspace(0, 100, 100)  # Temperature range from 0 to 100 °C
vapor_pressure_array = np.array([vapor_pressure(T+273.15)/100 for T in T_array])
p0_array = [p0/100 for T in T_array]

plt.figure(figsize=(7, 3))
plt.plot(T_array, vapor_pressure_array, color="darkblue",lw=2,label=r'$p_{vap,w}$')
plt.plot(T_array, p0_array,linestyle='dashed',color="k",lw=2,label=r'$p\!^\circ$')
plt.xlim(0,100)
plt.xlabel("Temperature (°C)")
plt.ylabel("Pressure (hPa)")
plt.legend()
plt.show()
```

:::

In order to account for humidity in our model, we need to calculate $f_{H_2O}$, the molar fraction of water in the air, from the relative humidity. Given equation {eq}`relative_humidity`:

$$
f_{H_2O}(p,T) = \dfrac{p_w}{p} = \dfrac{\varphi \cdot p_{vap,w}(T)}{p} 
$$(water_molar_frac)

with $p$ the atmospheric pressure. The average mass of a mole of humid air $m_{m}$ (subscript "m" from moist) now includes the molar mass of water $m_w$:

$$
\begin{align}
m_m(p,T) &= m_d\big( 1 - f_{H_2O}(p,T) \big) + m_w f_{H_2O}(p,T) \\[5pt]
         &=  m_d  - (m_d - m_w) \cdot f_{H_2O}(p,T)
\end{align}
$$(moist_molar_mass)

Water has a smaller mass compared to the other major species in the air, therefore humidity reduces the average molar mass of a parcel of air, making the atmospheric pressure smaller.

Notice from {eq}`water_molar_frac` that the molar fraction of water, needed to compute the average molar mass of moist air $m_m$, depends on the atmospheric pressure itself. This causes a problem, since the pressure is our sought variable. We can get around this through a trick, that is using a "first-order" dry-air pressure profile $p_{dry}(h)$ from the (equation {eq}`with_lapse`) instead of the real $p_{moist}(h)$. 

$$
\begin{cases}
f_{H_2O}(h) \approx  \varphi(h) \cdot \dfrac{p_{vap,w}(T(h))}{p_{dry}(h)} \\[15pt]
p_{vap,w}(T(h)) = 610.78\cdot e^{17.27(T(h)+273.15)/(T(h)+35.85)} \\[10pt]
p_{dry}(h) = p(0)e^{-g_0m_d/R\int_0^{h}dz/T(z)}
\end{cases}
$$(water_molar_frac_dry)

where we allowed the relative humidity to vary with altitude. Such approximation must not worry us, as the fraction of water vapor in the air never exceeds 5\% and its effect on the atmospheric pressure is very small. I was actually surprised when I discovered this, and I thought that clouds were the real water carriers. Turns out not to be true. We will return on clouds soon. We can do a quick calculation to give ourself an estimate of the water content in the atmosphere. A MSL atmospheric pressure of 1008 hPa is often considered the threshold for a low-pressure area. Assuming that the change in pressure is only attributed to a change in water content, we can use equation {eq}`water_molar_frac` to roughly evaluate an average value of the water fraction, $\overline{f_{H_2O}}$:
- 1010 hPa: $\overline{f_{H_2O}} \approx 0.85\%$
- 1008 hPa: $\overline{f_{H_2O}} \approx 1.37\%$
- 1000 hPa: $\overline{f_{H_2O}} \approx 3.46\%$

The lowest pressure ever recorded is 870 hPa, which gives an average water fraction of 38%, but our assumption is too crude for that case. Atmospheric pressure is also influenced by temperature and air circulation, which play a greater role in stormy weathers.

Let's built a function to compute $f_{H_2O}$ depending on relative humidity and either temperature or altitude (from the barometric formula). 

```{code-cell} ipython
def water_molar_fraction(RH,T=None,h=None,p=p0):
    if h is not None:
        if h > 20:
            return 0.0
        T = ISA_temperature(h)
        p_dry = pressure_dry(h)
    
    else:
        p_dry = p

    p_vap = vapor_pressure(T)

    # partial pressure of water vapor
    p_water = RH * p_vap

    # molar fraction of water vapor
    f_water = p_water / p_dry

    return f_water
```

Assuming a constant relative humidity with altitude, we can observe how the molar fraction of water in the air varies with temperature:

```{code-cell} ipython
:tags: ["hide-input"]

RHs = [0.25, 0.50, 0.75, 1.0]
f_water_RHs = [[water_molar_fraction(RH=RH,T=T+273.15)*100 for T in T_array] for RH in RHs]
plt.figure(figsize=(7, 3))
cmap = plt.get_cmap('coolwarm')
colors = cmap(np.linspace(0, 1, len(RHs)))
for idx,RH in enumerate(RHs):
    plt.plot(T_array, f_water_RHs[:][idx],lw=2,label=f'RH = {RH}',color=colors[len(RHs)-idx-1])
plt.xlim(0,40)
plt.ylim(0,8)
plt.xlabel("Temperature (°C)")
plt.ylabel("Water vapor molar fraction (%)")
plt.legend()
plt.show()
```

and it makes sense that $f_{H_2O}$ never exceeds 5\% on Earth. Taking the temperatures from [ISA](#ISA_temp) (15°C at the surface), we can see that the fraction of water decreases exponentially with altitude:


```{code-cell} ipython
:tags: ["hide-input"]

f_water_RHs = [[water_molar_fraction(RH=RH,h=alt)*100 for alt in altitudes] for RH in RHs]
plt.figure(figsize=(7, 3))
for idx,RH in enumerate(RHs):
    plt.plot(altitudes, f_water_RHs[:][idx],lw=2,label=f'RH = {RH}',color=colors[len(RHs)-idx-1])
plt.xlim(0,15)
plt.xlabel("Altitude (km)")
plt.ylabel("Water vapor molar fraction (%)")
plt.legend()
plt.show()
```

We can finally reach an expression for the atmospheric pressure with moist air

$$
p_{moist}(h)= p(0)\exp\bigg[-\dfrac{g_0}{R}\int_0^h\dfrac{m_m(z)}{T(z)}dz\bigg]
$$

Notice that $m_m(h)$ can be split into a $m_d$ term, and a term that depends on $z$ (equation {eq}`moist_molar_mass`). We can therefore split the integral into two terms, and we find back the expression for $p_{dry}(h)$, equation {eq}`with_lapse`:

$$
\begin{align}
p_{moist}(h) &= p_{dry}(h)\cdot\exp\bigg[\dfrac{g_0}{R}(m_d - m_w)\int_0^h \dfrac{f_{H_2O}(z)}{T(z)}dz\bigg] \\[10pt]
             &\approx p_{dry}(h) \cdot \exp\bigg[\dfrac{g_0}{R}(m_d - m_w)\int_0^h \varphi(z) \dfrac{p_{vap,w}(z)}{p_{dry}(z)T(z)}dz\bigg]
\end{align}
$$(p_moist)

Since $(m_d - m_w) > 0$, the exponential term in equation {eq}`p_moist` is greater than one. Humidity thus seems to effectively increase the pressure, which is the opposite of what we would expect! This is indeed not true. The effect of a smaller air molar mass is twofold. First, it reduces the slope of the pressure vertical profile, because less mass "pushes down" the air column. Secondly, a lighter air column produces a smaller pressure at the surface, $p(0)$. Here, we fixed $p(0)\equiv p^\circ$, so that we are violating the law of conservation of mass. We will return on this topic in the following section.


#### From dew point 

The data avaiable from [atmospheric soundings](wiki:atmospheric_sounding) does not usually provide the relative humidity, but the **dew point** ([](#fig:sounding)). 

:::{figure}https://www.greenskychaser.com/blog/wp-content/uploads/2011/05/OUN.gif
:width: 350px
:label: fig:sounding
Atmospheric sounding chart showing the temperature (red line) and the dew point (green line), among other parameters, measured at different altitudes.
:::

:::{note} Dew point
The <wiki:dew_point> is the temperature the air needs to be cooled to at constant pressure in order to reach a relative humidity of 100\%. 

The vapor pressure of water is the partial pressure of water $p_w$ in the air, and is given by

$$
p_w(T) = \varphi\cdot p_{vap,w}(T) 
$$

Mathematically speaking, the dew point is the temperature $T_{dew}$ at which $p_{vap,w}(T_{dew}) \equiv p_w(T)$, i.e. the equilibrium (saturated) vapor pressure is equivalent to the observed partial pressure of water:

$$
\begin{align}
p_{vap,w}(T_{dew}) &\equiv \varphi\cdot p_{vap,w}(T)\\[10pt]
A \exp\bigg[\dfrac{bT_{dew}}{T_{dew}+c}\bigg] &\equiv \varphi \cdot A \exp\bigg[\dfrac{bT}{T+c}\bigg] \\[10pt]
\dfrac{bT_{dew}}{T_{dew}+c} &\equiv \ln(\varphi) + \dfrac{bT}{T+c}
\end{align}
$$(dew_point)
where we used Tetens equation {eq}`tetens`. Solving for $T_{dew}$ we obtain:

$$
\begin{cases}
T_{dew} =\dfrac{c\gamma(T,\varphi)}{b-\gamma(T,\varphi)} \\[15pt]
\gamma(T,\varphi) = \dfrac{bT}{T+c}+\ln\left(\varphi\right)
\end{cases}
$$
:::


The relative humidity at a dry-bulb temperature $T$ (as measured by a conventional thermometer) can be obtained from the dew point $T_{dew}$ by (cfr. {eq}`dew_point`):

$$
\varphi(T,T_{dew}) = \dfrac{p_{vap,w}(T_{dew})}{p_{vap,w}(T)}
$$

Thus, the fraction of water vapor in the air (equation {eq}`water_molar_frac_dry`) can be rewritten as

$$
f_{H_2O}(h) \approx \dfrac{p_{vap,w}(T_{dew}(h))}{p_{dry}(h)}
$$




#### Clouds

Clouds are aerosols of liquid droplets or crystals, which are mainly water. They form when the relative humidity exceeds 100\%, or, equivalently, when the (dry-bulb) temperature reaches the dew point. The amount of water in clouds is measured by the <wiki:liquid_water_content> (LWC), which depends on the type of the cloud. Contrary to what one might expect, only a tiny fraction of the cloud volume is occupied by liquid water. LWC ranges within 0.03-3.0 g/m{sup}`3`, i.e. grams of liquid water per cubic meter of air. Considering that 1 m{sup}`3` of air at the limit of the troposhere, which is roughly the upper limit for clouds, weigths circa 364 grams:

$$
M(11 \, \mathrm{km}) = \dfrac{m_d\,p_{dry}(11\, \mathrm{km})\cdot 1\, \mathrm{m}^3}{RT(11\, \mathrm{km})} \approx 364 \, \mathrm{g}
$$

we deduce that the weight of clouds can be safely ignored.

#### Average water fraction 

The global average water fraction if 0.40\%. However, it is essentially zero above the troposphere. 


### Earth as a spinning spheroid

Let's now further improve our model by considering Earth as a spheroid (or ellipsoid) that spins and generates a gravitational field. Our description will then depend not only on altitude $h$, but also on the geographic latitude $\varphi$ (distinct from $\phi$ for relative humidity).

#### Gravity with altitude

Our first step is to include the variation of gravity with altitude. According to Newton's law of gravitation:

$$
g_0(h) = G\dfrac{M}{(R+h)^2}  
$$

Where the subscript 0 refers to the spherical symmetry of Earth. The formula can be expanded in powers of $h$ around the point $h=0$:

$$
g_0(h) = g_0\sum_{n=0}^{\infty} (-1)^{n}\left(n+1\right) \dfrac{h^n}{R^n}= \dfrac{g_0}{\left(1+\dfrac{h}{R}\right)^2}
$$(g_h)

where we used the power series. We call the term $\left(1+h/R\right)^{-2}$ the altitude factor, which makes the gravitational acceleration decrease with $h$ squared, as expected.


#### Reference ellipsoid

We now more accurately account for Earth's shape. We will use the World Geodetic System 1984 (WGS 84), which is used by the GPS system and suggested by the <wiki:International_Civil_Aviation_Organization>. 

:::{note}

#### World Geodetic System 1984

World Geodetic System 1984 describes Earth as a reference ellipsoid with the following parameters

```{table} Earth's parameters as defined by the WGS84 model.
:label: WGS84_data
:align: center
| Parameter | Symbol | Value |
| --- | --- | --- |
| equatorial semi-axis | $a$ or $R_e$ | $6378137.0 \ \mathrm{m}$ |
| polar semi-axis | $b$ or $R_p$ | $\approx 6356752.314140 \ \mathrm{m}$ |
| gravitational constant | $GM$ | $3.986004418\cdot 10^{-14}\ \mathrm{m^3/s^2}$ |
| angular velocity | $\omega$ | $72.92115\cdot10^{-6}\ \mathrm{rad/s}$ |
```

From these parameters are derived:
- eccentricity $e = \sqrt{1-b^2/a^2} \approx 0.0818$
- equatorial gravity $g_e = 9.7803253359\ \mathrm{m/s}^2$ 
- polar gravity $g_p = 9.8321849378\ \mathrm{m/s}^2$

The Ellipsoidal Gravity Formula (see <wiki:Theoretical_gravity#Somigliana_equation>) gives the gravitational acceleration depending on the latitude $\varphi$:

$$
g_{WGS84}(\phi) = g_e  \dfrac{1+k\sin^2(\phi)}{\sqrt{1-e^2\sin^2(\phi)}} 
$$(WGS84)

with $k$ a constant

$$
k = \dfrac{R_pg_p-R_eg_e}{R_eg_e} 
$$

:::

As much as I would love to, modeling Earth's gravitational field _ab initio_ is not worth it. Not because the mathematical modeling is too complicated (albeit cumbersome), but because I suppose that WGS84 relies on empirical data, and it's of course better than any model I could ever devise. 

We now want to combine equation {eq}`WGS84` with {eq}`g_h` to obtain a general expression for the gravitational acceleration for any $\varphi$ and $h$. But first, we need to find how Earth's radius varies with the latitude. The radius of the WGS84 ellipsoid with altitude is:

$$
R(\phi) = \dfrac{R_eR_p}{\sqrt{\left(R_p\cos(\phi)\right)^2+\left(R_e\sin(\phi)\right)^2}}
$$

so that we can combine it with the altitude factor to achieve a general expression for the gravitational acceleration depending on both altitude and latitude

$$
g(h,\phi) = g_e  \dfrac{1+k\sin^2(\phi)}{\sqrt{1-e^2\sin^2(\phi)}}\cdot  \dfrac{1}{\left(1+\dfrac{h}{R(\phi)}\right)^2}
$$(WGS84_h)

Note that simply applying the altitude factor is an approximate solution since, in an ellipsoid, the center of gravity is not intersected by the normal of the surface, except at the poles and at the equator ([](#fig:vertical-deflection)). Moreover, the altitude factor does not take into account the increase in centrifugal force with altitude. We can consider such correction negligible.


:::{figure} https://upload.wikimedia.org/wikipedia/commons/8/8f/Geodetic_coordinates.svg
:label: fig:vertical-deflection
Representation of an ellipsoid and the <wiki:vertical_deflection>.
:::


A python function that returns $g(h,\phi)$ may be:

```{code-cell} ipython
def g_WGS84_altitude(h,latitude):
    # WGS84 model for gravity times the altitude factor. Latitude in degrees, altitude in km
    a = 6378137.0
    b = 6356752.314140 
    g_e = 9.7803253359
    g_p = 9.8321849378
    e2 = 1 - (b**2 / a**2)
    k = (b * g_p - a * g_e) / (a * g_e)
    h = h * 1000.0 # Convert altitude to meters

    latitude = latitude * np.pi / 180.0  # Convert latitude to radians

    g_lat = g0 * (1 + k * np.sin(latitude)**2) / (1 - e2 * np.sin(latitude)**2)**0.5
    R_lat = a*b / (a * np.sin(latitude)**2 + b * np.cos(latitude)**2)
    g = g_lat / ( 1 + (h / R_lat) )**2
    return g
```

```{code-cell} ipython
:tags: ["hide-input"]
altitudes = np.linspace(0, 800, 100)
fig, ax = plt.subplots(figsize=(7, 3))
gravity_eq = [g_WGS84_altitude(alt,0) for alt in altitudes]
gravity_pol = [g_WGS84_altitude(alt,90) for alt in altitudes]
ax.plot(altitudes, gravity_eq,lw=2,color='darkred',label='Equator')
ax.plot(altitudes, gravity_pol,lw=2,color='lightblue',label='Pole')
ax.set_xlabel("Altitude (km)")
ax.set_ylabel("Gravity (m/s²)")
ax.set_xlim(0)
ax.legend()
plt.show()
```


### The thermosphere

The thermosphere is the outer layer of the atmosphere, above 80 km of altitude. The name stems from the high temperatures that are reached due to the ionizing radiation from the sun. 

#### Temperature of the thermosphere

The ISA model does not include the thermosphere, but reaches a maximum altitude of 86 km, where the temperature is 186.946 K. The data for the thermosphere is provided by the [NRLMSIS empirical model](https://swx-trec.com/msis/). The temperature profile of the two models combined is shown in [](#fig:thermo-T-altitude), where I extended the constant temperature of the stratopause (186.946 K) up to 107.41 km, and used an exponential regression upward (see red line). The empirical fitting gives the function:

$$
T(h) =  -9799 e^{-0.0238x} + 947.23, \quad 107.41\,\mathrm{km} \le h \le 1000 \,\mathrm{km}
$$

:::{figure} ../../images/thermosphere_T_h.png
:label: fig:thermo-T-altitude
:align: center
:w: 400px

Temperature profile with altitude (black dots). The red line is the union between the ISA model (below 86 km), and my fitting. Data from the [NRLMSIS empirical model](https://swx-trec.com/msis/) at 2024-05-01 00:00 UTC over 0°N 50°E. Graph from my [Desmos](https://www.desmos.com/calculator/pnt2qmypuf).

:::

We can include the thermosphere in a new function that returns the temperature up to 1000 km above sea level, and extend the functions to calculate $\chi$ and $p_{dry}$. Let's also add the possibility to change the surface temperature and the lapse rate, as long as it is positive (temperature decreases with altitude).

```{code-cell} ipython
def ISA_temperature_1000km(h,T_surf=None,L0=None):
    # Define base altitudes and temperatures for each layer
    base_altitudes = [0, 11, 20, 32, 47, 51, 71, 84.852, 107.41] 
    lapse_rates = [6.5, 0.0, -1.0, -2.8, 0.0, 2.8, 2.0, 0.0]
    base_temperatures = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946]

    if h > 1000 or h < 0:
        raise TypeError("Altitude range is 0-1000 km")
    elif h >= base_altitudes[-1]:   # Thermosphere
        return -9799 * np.exp(-0.0238 * h) + 947.23
    elif h < base_altitudes[2]:     # If in troposphere, check for custom parameters
        if T_surf is not None:
            temperature = T_surf
        else:    
            temperature = base_temperatures[0]

        if L0 is not None and L0 <= 0:
            raise TypeError("Troposphere lapse rate must be positive and non-zero")
        elif L0 is not None:
            lapse_rates[0] = L0

        if T_surf is not None or L0 is not None:
            if ( temperature - base_temperatures[1] ) / lapse_rates[0] > base_altitudes[2]:
                raise TypeError("Surface temperature is too high or lapse rate is too small")
            
            base_altitudes[1] = ( temperature - base_temperatures[1] ) / lapse_rates[0]


    # Find the layer corresponding to the altitude
    for i in range(len(base_altitudes) - 1):
        if h <= base_altitudes[i+1]:
            return temperature - lapse_rates[i] * (h - base_altitudes[i])
        else:
            temperature = base_temperatures[i+1]

```

#### Composition in the thermosphere

The value that we calculated in [](#subheading-chemical-composition) refers to the global average of the atmospheric composition. We are now interested in how such composition changes with altitude. We know that turbulence and diffusion make the atmospheric composition constant up to 85 km ([](#fig:composition-altitude)). The vertical profile of carbon dioxide, one of the heaviest molecules in the air, starts decreasing from an altitude of 60 km ([](#fig:CO2-altitude)). 

:::{figure} https://upload.wikimedia.org/wikipedia/commons/b/bd/Chemical_composition_of_atmosphere_accordig_to_altitude.png
:label: fig:composition-altitude
:align: center
:w: 400px

Atmospheric composition versus altitude

:::

:::{figure} ../../images/CO2_altitude.jpg
:label: fig:CO2-altitude
:align: center
:w: 400px

Altitude profile of $\mathrm{CO}_2$ molar fraction, in ppm. From [Brown et al. (2024)](https://doi.org/10.1029/2024JA032659). Graph from my [Desmos](https://www.desmos.com/calculator/ohzv68xy5w).

:::

Since the atmospheric pressure at 85 km is $p_{dry}(85\,\mathrm{km})\approx 0.5\,\mathrm{Pa}$. We are therefore talking about minuscole changes, so that we are allowed to be rough. [](#fig:avgm-altitude) shows the vertical profile of the average molar mass of air, computed from the [NRLMSIS empirical model](https://swx-trec.com/msis/) data. A simple exponential regression from an altitude of 85 km can be made (blue dashed line). We can use the following expression for the average molar mass of dry air, in g/mol:

$$
m_d(h)=\begin{cases}
28.9656,\quad h\le 85\,\mathrm{km}\\
28.9656\cdot e^{-0.002(h-85\,\mathrm{km})}, \quad h\gt 85\,\mathrm{km}
\end{cases}
$$


:::{figure} ../../images/avgm_altitude.png
:label: fig:avgm-altitude
:align: center
:w: 400px

Average molar mass of air (black dots) versus altitude. The blue dashed line is the (exponential) regression. Data from the [NRLMSIS empirical model](https://swx-trec.com/msis/) at 2024-05-01 00:00 UTC over 0°N 50°E. 

:::

```{code-cell} ipython
def air_avg_molar_mass(h,RH):
    #altitude in km
    if h < 20:
        f_water = water_molar_fraction(RH,h=h)
        return (1 - f_water) * m_dry + f_water * m_water
    elif h < 85:
        return m_dry
    else:
        return m_dry * np.exp(-0.002*(h-85))
```


### Pressure at sea level

Until now, we blindly used the value of 1013.25 hPa as surface pressure. Such value is the average pressure at MSL which then includes the average molar fraction of water, $\overline{f_w}\approx 0.40\%$. We thus need to consider that dry air would produce a larger pressure at sea level, and vice-versa, moist air with a total $f_w > 0.40\%$ produces a smaller surface pressure. 

Pressure at MSL can be evaluated as total mass of the air column (cfr. equation {eq}`pressure_h`):

$$
p(0) = \dfrac{F}{A} = \int_{0}^{\infty} n(z)m_{m}(z)g(z) dz = \int_{0}^{\infty} \dfrac{p(z)m_{m}(z)g(z)}{RT(z)}dz
$$

Which can be computed only if we already know the vertical profile of the pressure. We can of course use the dry-air approximation and estimate $p_{moist}(0)$ by means of a multiplying factor $\chi_m$, given by the ratio between the pressures exerted by the moist air column and the standard air (DEFINE STANDARD) column, i.e. $p_{moist}(0) = \chi_m \, p_{std}(0)$:

$$
\chi_m = \dfrac{p_{moist}(0)}{p_{std}(0)} \approx \dfrac{\int_{0}^{\infty} \frac{p_{dry}(z)m_{m}(z)g(z)}{T(z)}dz}{\int_{0}^{\infty} \frac{p_{dry}(z)m_{sdt}(z)g(z)}{T(z)}dz} 
$$

where we approximated the moist air pressure with the dry air pressure. By explicitating $m_m(z)$ according to equation {eq}`water_molar_frac` we obtain

$$
\chi_m \approx 1 - \dfrac{m_d-m_w}{m_d} \dfrac{\int_{0}^{\infty} f_{H_2O}(z)\frac{p_{dry}(z)}{T(z)}dz}{\int_{0}^{\infty} \frac{p_{dry}(z)}{T(z)}dz} 
$$


With the ISA temperature profile and $\varphi = 1$ up to 20 km I compute a value of 0.9983, namely an average water fraction of $\overline{f_{w}}\approx 0.17\%$, and $p_{moist}(0) \approx 1011.5 \ \mathrm{hPa}$. 

Let's then build a function to calculate the factor $\chi$ and then include it in a function that evaluates the pressure with moist air. 


```{code-cell} ipython
def calc_chi(RH,maxh=84.852):
    integral_frac, _ = quad(lambda z: water_molar_fraction(RH,h=z) * pressure_dry(z) / ISA_temperature(z), 0, maxh, limit=100, points=[0, 11, 20, 32, 47, 51, 71, 84.852])
    integral_tot, _ = quad(lambda z: pressure_dry(z) / ISA_temperature(z), 0, maxh, limit=100, points=[0, 11, 20, 32, 47, 51, 71, 84.852])
    
    chi = 1 - ((m_dry - m_water) / m_dry) * (integral_frac / integral_tot)

    return chi

def pressure_moist(h,RH,chi):
    integral, _ = quad(lambda z: water_molar_fraction(RH,h=z) / ISA_temperature(z), 0, h, limit=100, points=[0, 11, 20, 32, 47, 51, 71, 84.852])
    p_moist = chi * pressure_dry(h) * np.exp(g0 / R * (m_dry - m_water) * 1000 * integral)

    return p_moist
```

```{code-cell} ipython
:tags: ["hide-input"]
chi = calc_chi(RH)
pressure_moist_arr = np.array([pressure_moist(alt,1.0,chi)/100 for alt in altitudes])   # divided by 100 to return hPa
diff_moist_dry = pressure_moist_arr - pressure_dry_arr
fig, ax1 = plt.subplots(figsize=(7, 3))
ax2 = ax1.twinx()
ax2.axhline(0,color='r',alpha=0.5,linestyle='dashed')
ax1.plot(altitudes, pressure_moist_arr,lw=2,label='$p_{moist}$ (RH=1)',color='b')
ax1.plot(altitudes, pressure_dry_arr,lw=2,label='$p_{dry}$',color='darkred',linestyle='dashed')
ax2.plot(altitudes, diff_moist_dry,lw=2,label=r'$p_{moist}-p_{dry}$',color='r',linestyle='dashed')
plt.xlim(0,35)
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel("Pressure (hPa)")
ax2.set_ylabel("Pressure difference (hPa)")
ax1.set_ylim(0)
ax1.legend(loc=(0.65,0.25))
ax2.legend(loc=(0.65,0.50))
plt.show()
```

We can indeed see that the effect of water vapor is small, and the two pressure profiles very close to each other. To appreciate the deviation, also plotted is the relative difference, in red. The difference is in the order of 1 hPa (~ 0.1\%). We are therefore very happy with the dry-air approximation that we made for the molar fraction of water, equation {eq}`water_molar_frac_dry`. The inclusion of the second-order term will provide negligible improvement over substantially higher computational costs.


### The final model

We thus have built the best model we can concieve for atmospheric model. If we allow temperature and humidity as well to change with latitude, we can write

$$
p(h,\phi)=p(0)\exp\bigg[-\dfrac{1}{R}\int_0^h\dfrac{g(z,\phi)m_m(z,\phi)}{T(z,\phi)}dz\bigg]
$$

Our python function will then include the recently defined `g_WGS84_altitude`, `air_avg_molar_mass`, and `ISA_temperature_1000km` functions to evaluate the integral. We don't need to include our recent enhancements of the model in the calculation of $p_{dry}$, since that is needed only to obtain $f_{H_2O}$.

```{code-cell} ipython
def pressure_dry_1000km(h,T_surf=288.15,L0=6.5):  
    H1 = ( T_surf - 216.65 ) / L0
    if H1 > 11:
        raise TypeError("Surface temperature is too high or lapse rate is too small")
    
    integral, _ = quad(lambda h: air_avg_molar_mass(z,RH=0.0) / ISA_temperature_1000km(h,T_surf,L0), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    pressure = p0 * np.exp(-g0 / R * 1000 * integral)

    return pressure

def water_molar_fraction_1000km(RH,h,T_surf=288.15,L0=6.5):
    if h > 20 or RH == 0:
        return 0.0
    
    H1 = ( T_surf - 216.65 ) / L0
    if H1 > 11:
        raise TypeError("Surface temperature is too high or lapse rate is too small")

    T = ISA_temperature_1000km(h,T_surf,L0)
    p_dry = pressure_dry_1000km(h,T_surf,L0)
    p_vap = vapor_pressure(T)
    f_water = RH * p_vap / p_dry

    return f_water

def calc_chi_1000km(RH,maxh=1000,T_surf=288.15,L0=6.5):
    H1 = ( T_surf - 216.65 ) / L0
    if H1 > 11:
        raise TypeError("Surface temperature is too high or lapse rate is too small")
    
    integral_frac, _ = quad(lambda z: water_molar_fraction_1000km(RH,z,T_surf,L0) * pressure_dry_1000km(z,T_surf,L0) / ISA_temperature_1000km(z,T_surf,L0), 0, 20, limit=100, points=[0, H1, 20])  #up to 20 km since above that the water fraction is 0 
    integral_norm, _ = quad(lambda z: pressure_dry_1000km(z,T_surf,L0) / ISA_temperature_1000km(z,T_surf,L0), 0, maxh, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    
    chi = 1 - ((m_dry - m_water) / m_dry) * (integral_frac / integral_norm)

    return chi

def pressure_moist_WGS84_1000km(h,lat,RH=0.0,chi=1.0,T_surf=288.15,L0=6.5):
    H1 = ( T_surf - 216.65 ) / L0
    integral, _ = quad(lambda z: g_WGS84_altitude(z,lat) * air_avg_molar_mass(z,RH) / ISA_temperature_1000km(z,T_surf,L0), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])

    p_moist = chi * p0 * np.exp(- 1000 * integral / R)

    return p_moist
```

#### Total atmospheric mass

With our model, we can calculate the total mass of our atmosphere, assuming a MSL standard pressure of 1013.25 hPa. What we need to do is just integrate the density through the entire atmospheric spherical shell. The density $\rho_d$ of an ideal gas is given by

$$
\rho = \dfrac{mp}{RT}
$$

with $m$ the molar mass. We can then calculate the total atmospheric mass $M_{atm}$ as

$$
\begin{align}
 M_{atm} &= \int_{0}^{2\pi} \int_{0}^{\pi} \int_{R(\phi)}^{\infty} \rho(r,\phi)\, d\theta\, \sin(\phi) d\phi\, r^2 dr \\[15pt]
 &= \int_{0}^{2\pi} \int_{0}^{\pi} \int_{R(\phi)}^{\infty} \dfrac{m_m(r-R(\phi),\phi)\,p(r-R(\phi),\phi)}{RT(r-R(\phi))} d\theta\, \sin(\phi) d\phi\, r^2 dr 
\end{align}
$$

The radial integration should start from the Earth's surface ($R(\phi)$), and since our functions are defined with respect to the altitude $h$, we need to shift them by $R(\phi)$. We can then substitute the integrand $r$ with $h=r-R(\phi)$, so that the integral in $dh=dr$ will go from 0 to infinity. Moreover, we can integrate in $d\theta$, as no functions depend on the longitude $\theta$:

$$
\begin{align}
 M_{atm} &= \dfrac{2\pi}{R} \int_{0}^{\pi} \int_{0}^{\infty} \dfrac{m_m(h,\phi)\,p(h,\phi)}{T(h,\phi)}  \sin(\phi) d\phi\, \left(h+R(\phi)\right)^2 dh \\[15pt]
&= \dfrac{2\pi}{R} p(0) \int_{0}^{\pi} \int_{0}^{\infty} \dfrac{m_m(h,\phi)}{T(h,\phi)}\exp\bigg[-\dfrac{1}{R}\int_0^{h} \dfrac{g(z,\phi)m_m(z,\phi)}{T(z,\phi)}dz\bigg]  \sin(\phi) d\phi\, \left(h+R(\phi)\right)^2 dh
\end{align}
$$

where we explicitated $p(h,\phi)$. To calculate $M_{atm}$ thus requires to evaluate three integrals. It may make your laptop hot.




### Altitude from pressure

%Take equation T_from lapse and solve it with T = T(0)-Lz or T(h)-Lh-Lz and rearrange



## Empirical model

### The NRLMSIS model



The [NRLMSIS empirical model](https://swx-trec.com/msis/?lz=N4Igtg9gJgpgNiAXCYAdEUCGAXG2CWYM6iAjAOykAsArAEykBsADK6wL4gA0ImcB2AK6wkKdHBz4hsEqWZdxEAHYBzKcOJI6zTjwDOggE4AzTAGMYotL37qZSOTu4gAbpkP5MAIziXk6ADkAeXRnNw9vXz1RAG10AEFDdAUQAAlk9FTNFICMkAC6POC8kO50IMKykABZTD09PIAVGDAABxhDHCNNAF0QdiA) provides an enormous quantity of data , that ...

%The model that we have built so far is actually rather accurate, especially if we want to use it at the surface of the Earth, where the atmospheric composition, the gravitational field, and 

%We can also substitute the lapse rate with the actual data from [atmospheric sounding](wiki:Atmospheric_sounding), if we want to calculate the change in pressure more accurately. For now, we are happy with the standard lapse rate. 


### Using sounding data


### The tool