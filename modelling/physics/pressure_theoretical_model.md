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
I did this project in Desmos back in 2021. For some reason, I've always been obsessed with air pressure, perhaps because I discovered that it can be measured outstandingly precisely with a smartphone. It can even detect altitude changes on the order of meters. So I asked myself: can I derive how the atmospheric pressure varies with altitude, and can I do it as accurately as possible? This way, I could estimate a change in altitude by measuring the change in pressure with my phone. And it works rather nicely!
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

where $n$ the numerical density of the gas, and $T$ its temperature. To model Earth’s atmosphere, we imagine an infinitely high column of gas subjected to Earth's gravitational field. For now, we assume that Earth is spherically symmetric. We also assume that the mass of the atmosphere is negligible, so that there is no self-interaction, which would make the problem not analytically solvable. The pressure of the gas at a certain altitude $h$ must be equivalent to the pressure exerted by the column of gas above it, due to the gravitational field, i.e.:

$$
    nk_BT(h) = p(h)=\dfrac{F(h)}{A}
$$ (pressure_gravfield)

with $F$ the weight of the gas column above $h$. The pressure of the column can also be calculated as the integral along the vertical direction $z$ of the mass density $\rho(z)$ of the gas, times the gravitational acceleration $g(z)$:

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

Interestingly, to compute the pressure at a certain altitude $h$ we just need to know how the integrand varies below that point, and not above. This is because the formula assumes we know $p(0)$, the pressure at $h=0$, so that the information of the air column above is implicitly contained there. If we choose $h=0$ to be the sea level, then $p(0)$ is the barometric pressure, which in standard conditions is $1013.25\ \mathrm{hPa}$. In reality, the sea-level pressure varies continuously and it must be measured. We will understand how temperature and average mass vary with altitude in the following section. For now, let’s find the simplest solution by assuming that are all variables inside the integral, i.e. composition, temperature, and gravity, are constants. Such approximation is valid close to the Earth's surface level. We then obtain:

$$
p_{bar}(h)=p(0)\exp\bigg[-\dfrac{m_0gh}{k_BT}\bigg]
$$(barometric_formula)

We recognize $m_0gh$ as the potential energy of a single gas particle, and $k_BT$ as its thermal energy. Wait, what? Equation {eq}`barometric_formula` is called the **barometric formula**. The exponential term is the Boltzmann factor ($e^{-E/k_BT}$), which, in a canonical ensemble (NVT, our case), represents the probability of the system to be in a state with energy $E$. In our case, $E$ is the potential energy of a mass in a uniform gravitational field, and the Boltzmann factor represents the probability of a particle to be at that altitude. Macroscopically, this becomes the actual pressure of the gas.


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

(ISA_temp)=
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
3. **Stratosphere** (20-32 km): -1.0 °C/km
4. **Stratosphere** (32-47 km): -2.8 °C/km
5. **Stratopause** (47-51 km): 0.0 °C/km
6. **Mesosphere** (51-71 km): 2.8 °C/km
7. **Mesosphere** (71-84.852 km): 2.0 °C/km

There are other layers above, but can be ignored since the atmosphere is extremely rarefied. The ranges are given in geopotential altitude.


(subheading-chemical-composition)=
### Chemical composition
The composition of the dry atmosphere is kindly provided by the [NOAA](https://www.noaa.gov/jetstream/atmosphere). However, the molar fractions they provide sum to more than one, due to experimental error (as also [wikipedia](https://en.wikipedia.org/wiki/Atmosphere_of_Earth#Composition) reports). While looking for the most precise values, I noticed an incongruence in the reported amount of $\mathrm{CO_2}$. I quickly realized a shocking fact: the atmospheric concentration of $\mathrm{CO_2}$ is rising so quickly that most values are now outdated. For example, [engineering toolbox](https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html) uses a value of $f_{CO_2}=0.033\%$, which is the fraction from circa 50 years ago (in the '70s). The value is now (2025) $0.042\%$, giving an outstanding 27% increase. This makes me wonder whether NOAA takes into account the change in fractional concentrations due to $\mathrm{CO_2}$ emissions and $\mathrm{O_2}$ depletion. A strong hint is that by substituting the present $\mathrm{CO_2}$ concentration with $0.033\%$ in NOAA's value, the sum magically becomes $1$. For this reason, I took NOAA's value and assumed that the $0.011\%$ increase is due to combustion, and it substitutes $\mathrm{O_2}$ molecules

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


## Lapse rates

Let's start from the barometric formula {eq}`barometric_formula` and remove approximations one by one to finally arrive at the general formula {eq}`general_formula`. First, we introduce the empirical lapse rates that we learned above. We thus leave only the temperature term in the integral:

$$
p_{dry}(h)=p(0)\exp\bigg[-\dfrac{gm_d}{R}\int_0^h\dfrac{1}{T(z)}dz\bigg]
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
    T0 = 288.15  # MSL standard temperature in Kelvin
    p0 = 101325  # MSL standard atmospheric pressure in Pa
```

and define a function that returns the temperature at a certain altitude, using equation {eq}`T_fromlapse` for each layer.

```{code-cell} ipython
def temperature(altitude):
    # Define base altitudes and temperatures for each layer
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
def pressure_barometric(altitude):
    pressure = p0 * np.exp(-m_dry * g0 * altitude / (R * T0))

    return pressure

def pressure_dry(altitude):
    integral, err = quad(lambda h: 1 / temperature(h), 0, altitude, limit=100, points=[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    pressure = p0 * np.exp(-m_dry * g0 / R * integral)

    return pressure
```


```{code-cell} ipython
:tags: ["hide-input"]
altitudes = np.linspace(0, 85, 500)     # altitude array in km
pressure_dry_arr = [pressure_dry(1000*alt)/1000 for alt in altitudes]   # divided by 1000 to return hPa
pressure_barom_arr = [pressure_barometric(1000*alt)/1000 for alt in altitudes]
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, pressure_dry_arr, color='k',lw=2,label='with lapse rates')
plt.plot(altitudes, pressure_barom_arr, linestyle='dashed', color='c',lw=2, label='barometric approx.')
plt.xlabel("Altitude (km)")
plt.ylabel("Pressure (hPa)")
plt.legend()
plt.show()
```

This results in a very good estimation of the pressure, especially at lower altitudes. As a recap, our model has so far been built within the following approximations:
- Ideal gas law
- Spherical Earth, no spin
- Constant gravitational field
- Constant atmospheric composition, dry air
- Empirical lapse rates

Let's keep improving our model by removing most of them one by one. The most impactful one is arguably dry air: on Earth, air is never completely dry, but some water vapor is mixed with the other gases. In practice, addition of water vapor affects the composition, namely the molar mass $m$ in our model.
%Furthermore, water droplets and ice crystals can also be encountered, which are clouds. 


## Humidity
The amount of water vapor in the air is usually measured in relative humidity (RH or $\phi$), which is the fraction of the water vapor in the air relative to the "maximum" potential at that temperature. 


### From relative humidity

:::{note}

### Relative humidity and water vapor pressure

**Relative humidity** $\phi$ is defined as the ratio between the measured partial pressure of water $p_w$ and its equilibrium (saturation) vapor pressure $p_{vap,w}$:

$$
\varphi = \dfrac{p_w}{p_{vap,w}}
$$(relative_humidity)

Note that these terms assume a different meaning in physics and in meteorology. In physics, the **vapor pressure** of a substance is the partial pressure of the gas phase in *equilibrium* with the liquid phase. In meteorology, however, **vapor pressure** refers to the measured partial pressure of water vapor, even when not in equilibrium, whereas the **saturation vapor pressure** is the actual vapor pressure of water, i.e., in equilibrium conditions. This nomenclature comes from the erroneous idea of air dissolving water vapor, eventually reaching a saturation limit. In reality, the vapor pressure of a substance only depends on the liquid-phase temperature, in first approximation. Here, we will try to use the correct physical definitions, albeit minding such incongruences.

The dependence of the vapor pressure of a substance on temperature can be estimated from the Clausius-Clapeyron equation knowing its boiling point $T_b$ at standard pressure $p^\circ$  ($p^\circ$ = 1 atm = 101325 Pa, and $T_b$ = 99.97°C for water)

$$
p_{vap}(T) = p^\circ \cdot \exp\bigg[{-\dfrac{\Delta_{vap} H}{R}\bigg(\dfrac{1}{T}-\dfrac{1}{T_b}\bigg)}\bigg]
$$

With $\Delta_{vap} H$ the enthalpy of vaporization. However, this formula lacks of the desired accuracy, because of the numerous approximations that led to it: ideal gas, constant $\Delta_{vap} H$ with temperature, no volume change, etc. This is particularly true for water, which deviates from the ideality due to the strong intermolecular interactions it can establish. For this reason, it is common to use empirical formula. Here, we are going to use the <wiki:Tetens_equation> (where T is in °C, and it returns Pa)

$$
p_{vap,w}(T) = 610.78\cdot \exp\bigg[\dfrac{17.27T}{T+237.3}\bigg]
$$

```{code-cell} ipython
def vapor_pressure(T):
    # Tetens equation
    A = 610.78
    B = 17.27
    C = 237.3  
    p_vap = A * np.exp((B * (T - 273.15)) / (T - 273.15 + C))   
    return p_vap
```

```{code-cell} ipython
:tags: ["hide-input"]
T_array = np.linspace(0, 100, 100)  # Temperature range from 0 to 100 °C
vapor_pressure_array = [vapor_pressure(T+273.15)/100 for T in T_array]
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
         &=  m_d  + (m_w - m_d) \cdot f_{H_2O}(p,T)
\end{align}
$$(moist_molar_mass)

Water has a smaller mass compared to the other major species in the air, therefore humidity reduces the average molar mass of a parcel of air, making the atmospheric pressure smaller.

Notice from {eq}`water_molar_frac` that the molar fraction of water, needed to compute the average molar mass of moist air $m_m$, depends on the atmospheric pressure itself. This causes a problem, since the pressure is our sought variable. We can get around this through a trick, that is using a "zeroth-order" dry-air pressure profile $p_{bar}(h)$ from the **barometric formula** (equation {eq}`barometric_formula`) instead of the real $p_{moist}(h)$. 

$$
\begin{cases}
f_{H_2O}(h) \approx  \varphi(h) \cdot \dfrac{p_{vap,w}(T(h))}{p_{bar}(h)} \\[15pt]
p_{vap,w}(T(h)) = 610.78\cdot e^{17.27(T(h)+273.15)/(T(h)+35.85)} \\[10pt]
p_{bar}(h) = p(0)e^{-gm_dh/(RT_0)}
\end{cases}
$$(water_molar_frac_baro)

where we allowed the relative humidity to vary with altitude. This approximation must not worry us, as the fraction of water vapor in the air never exceeds 5\% and its effect on the atmospheric pressure is minimal. If such information doesn't sound right, note that this doesn't mean that water does not alter the pressure at all. It is the presence of condensed-phase water (clouds) that greatly affects the atmospheric pressure, and it does so by locally occupying a much greater fraction of volume of atmosphere. We will later return on clouds. 

Let's built a function to compute $f_{H_2O}$ depending on relative humidity and either temperature or altitude (from the barometric formula). 

```{code-cell} ipython
def water_molar_fraction(RH,T=None,h=None):
    if h is not None:
        if h > 20000:
            return 0.0
        T = temperature(h)
        p_barom = pressure_barometric(h)
    
    else:
        p_barom = p0

    p_vap = vapor_pressure(T)

    # partial pressure of water vapor
    p_water = RH * p_vap

    # molar fraction of water vapor
    f_water = p_water / p_barom

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

f_water_RHs = [[water_molar_fraction(RH=RH,h=1000*alt)*100 for alt in altitudes] for RH in RHs]
plt.figure(figsize=(7, 3))
for idx,RH in enumerate(RHs):
    plt.plot(altitudes, f_water_RHs[:][idx],lw=2,label=f'RH = {RH}',color=colors[len(RHs)-idx-1])
plt.xlim(0,15)
plt.xlabel("Altitude (km)")
plt.ylabel("Water vapor molar fraction (%)")
plt.legend()
plt.show()
```

We finally reach an expression for the atmospheric pressure with moist air

$$
p_{moist}(h)\approx p(0)\exp\bigg[-\dfrac{g}{R}\int_0^h\dfrac{m_m(z)}{T(z)}dz\bigg]
$$(p_moist)

Notice that $m_m(h)$ can be split into a $m_d$ term, and a term that depends on $z$ (equation {eq}`moist_molar_mass`). We can therefore split the integral into two terms, and we find back the expression for $p_{dry}(h)$, equation {eq}`with_lapse`:

$$
p_{moist}(h)\approx p_{dry}(h)\cdot\exp\bigg[-\dfrac{g}{R}(m_w - m_d)\int_0^h \dfrac{f_{H_2O}(z)}{T(z)}dz\bigg]
$$(p_moist)


### From dew point 

aa

:::{note}
### Sea-level pressure reduction

Atmospheric pressure is always reported at the MSL. When a weather station at a certain altitude measures the local pressure, that value is then reduced to the sea level. In other words, the station must estimate the pressure that would be measured if someone digged down to the sea level. This is called **sea-level pressure reduction**, and since no air exists below ground, it is purely hypothetical, so that many assumptions needs to be made. For example, how temperature and humidity would vary going down cannot be properly defined.

The sea-level pressure reduction is carried out by means of the <wiki:hypsometric_equation>:

$$
p_{MSL} = p_{obs}\cdot \exp\bigg[\frac{m_mgh}{R\overline{T}}\bigg]
$$

where $h$ is the altitude in which the pressure $p_{obs}$ is measured, $p_{MSL}$ is the pressure at the MSL, and $\overline{T}$ is the mean temperature of the (moist) air with molar mass $m_m$ (I omitted to explain the <wiki:virtual_temperature>, but the formula is equivalent). The hypsometric equation is basically the barometric formula (equation {eq}`barometric_formula`), but with an average temperature across the vertical distance. The average temperature can simply be computed using the standard lapse rate $\Gamma$, and a 12-hour average surface temperature, an attempt to exclude the effect of the irradiation of surface of the Earth:
$$
    \overline{T} \approx \dfrac{T_{obs}(t) + T_{obs}(t-12h)}{2} - \Gamma \cdot \dfrac{h}{2}
$$
Additional refinements can be employed, but they are often empirical and specific to the station site. Mind that this pressure reduction is fictitious, and it is better to be consistent (across the world) rather than accurate. And I don't like this.

:::


## Earth as a spinning spheroid

Let's now further improve our model by considering Earth as a spheroid (or ellipsoid) that spins and generates a gravitational field. Our description will then depend not only on altitude $h$, but also on the geographic latitude $\varphi$.

### Gravity with altitude

Our first step is to include the variation of gravity with altitude. According to Newton's law of gravitation:

$$
g_0(h) = G\dfrac{M}{(R+h)^2}  
$$

Where the subscript 0 refers to the spherical symmetry of Earth. The formula can be expanded at $h=0$:

$$
g_0(h) = g_0\sum_{n=0}^{\infty} (-1)^{n}\left(n+1\right) \dfrac{h^n}{R^n}= \dfrac{g_0}{\left(1+\dfrac{h}{R}\right)^2}
$$(g_h)

where we used the power series. We can call the term $\left(1+h/R\right)^{-2}$ the altitude factor.


### Reference ellipsoid

We now account more accurately for Earth's shape. We will use the World Geodetic System 1984 (WGS 84), which is used by the GPS system and suggested by the <wiki:International_Civil_Aviation_Organization>. 

:::{note}

### World Geodetic System 1984

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
- polar gravity $g_e = 9.8321849378\ \mathrm{m/s}^2$

The Ellipsoidal Gravity Formula (see <wiki:Theoretical_gravity#Somigliana_equation>) gives the gravitational acceleration depending on the latitude $\varphi$:

$$
g(\phi) = g_e  \dfrac{1+k\sin^2(\phi)}{\sqrt{1-e^2\sin^2(\phi)}} 
$$(WGS84)

with $k$ a constant

$$
k = \dfrac{R_pg_p-R_eg_e}{R_eg_e} 
$$

:::

As much as I would love to, modeling Earth's gravitational field _ab initio_ is not worth it. Not because the mathematical modelling is too complicate (albeit cumbersome), but because I suppose that WGS84 relies on empirical data, and it's of course better than any model I could ever devise. 

We now want to combine equation {eq}`WGS84` with {eq}`g_h` to obtain a general expression for the gravitational acceleration for any $\varphi$ and $h$. But first, we need to find how Earth's radius varies with the latitude. The radius the WGS84 ellipsoid is:

$$
R(\phi) = \dfrac{R_eR_p}{\sqrt{\left(R_p\cos(\phi)\right)^2+\left(R_e\sin(\phi)\right)^2}}
$$

so that we obtain

$$
g(h,\phi) = g_e  \dfrac{1+k\sin^2(\phi)}{\sqrt{1-e^2\sin^2(\phi)}}\cdot  \dfrac{1}{\left(1+\dfrac{h}{R(\phi)}\right)^2}
$$(WGS84_h)

Note that the altitude factor is approximate, since, in an ellipsoid, the center of gravity is not intersected by the normal of the surface, except at the poles and at the equator. Moreover, the factor does not take into account the increase of the centrifugal force with altitude. We can consider such correction negligible.


:::{figure} https://upload.wikimedia.org/wikipedia/commons/8/8f/Geodetic_coordinates.svg
Exaggerated representation of an ellipsoid and the <wiki:vertical_deflection>.
:::
 



## Altitude from pressure

%Take equation T_from lapse and solve it with T = T(0)-Lz or T(h)-Lh-Lz and rearrange



# Empirical model

## The NRLMSIS model

The value that we calculated in [](#subheading-chemical-composition) refers to the global average of the atmospheric composition. We are now interested in how such composition changes with altitude. We know that turbulence and diffusion make the atmospheric composition rather constant up to 80 km ([](#composition-altitude)). The vertical profile of carbon dioxide, one of the heaviest molecules in the air, starts decreasing from an altitude of 60 km. 

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

%We can also substitute the lapse rate with the actual data from [atmospheric sounding](wiki:Atmospheric_sounding), if we want to calculate the change in pressure more accurately. For now, we are happy with the standard lapse rate. 


## Using sounding data


## The tool