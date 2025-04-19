import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
g = 9.80665  # Standard gravity in m/s^2
m_dry = 28.9647e-3  # Molar mass of dry air in kg/mol
m_water = 18.01528e-3  # Molar mass of water in kg/mol
p0 = 101325  # Sea level standard atmospheric pressure in Pa
T0 = 288.15  # MSL standard temperature in Kelvin

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


"""# Plot the graph
plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, temperatures, color="red",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Temperature (K)")
plt.grid(True)
plt.show()"""


def barometric_formula(altitude):
    p_bar = p0 * np.exp(-m_dry * g * altitude / (R * T0))

    return p_bar

def pressure_dry(altitude):
    integral, _ = quad(lambda h: 1 / temperature(h), 0, altitude, limit=100, points=[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    p_dry = p0 * np.exp(-m_dry * g / R * integral)

    return p_dry


def sat_vapor_pressure(T):
    # Constants for the Tetens equation
    A = 610.78  # Pa
    B = 17.27
    C = 237.3  # °C

    p_vap = A * np.exp((B * (T - 273.15)) / (T - 273.15 + C))     # It converts temperature to Celsius

    return p_vap


def water_molar_fraction(altitude,RH):
    # Calculate the vapor pressure

    if altitude > 20000:
        return 0.0

    T = temperature(altitude)
    p_vap = sat_vapor_pressure(T)

    # Calculate the total pressure
    p_total = barometric_formula(altitude)

    # Calculate the partial pressure of water vapor
    p_water = RH * p_vap

    # Calculate the molar fraction of water vapor
    f_water = p_water / p_total

    return f_water


def moist_molar_mass(altitude,RH):
    f_water = water_molar_fraction(altitude,RH)

    # Calculate the molar mass of moist air
    m_moist = (1 - f_water) * m_dry + f_water * m_water

    return m_moist


def pressure_moist(altitude,RH):
    integral, _ = quad(lambda h: moist_molar_mass(h,RH) / temperature(h), 0, altitude, limit=100, points=[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    p_moist = p0 * np.exp(-g / R * integral)

    return p_moist

def water_fraction_T(T,RH):
    # Calculate the vapor pressure
    p_vap = sat_vapor_pressure(T)

    # Calculate the total pressure
    p_total = p0

    # Calculate the partial pressure of water vapor
    p_water = RH * p_vap

    # Calculate the molar fraction of water vapor
    f_water = p_water / p_total

    return f_water

def moist_molar_mass_T(T,RH):
    f_water = water_fraction_T(T,RH)

    # Calculate the molar mass of moist air
    m_moist = (1 - f_water) * m_dry + f_water * m_water

    return m_moist


RH = 1.0  # Relative humidity

# Generate altitudes from 0 to 85,000 meters
altitudes = np.linspace(0, 80000, 50)

# Calculate temperatures for each altitude
temperatures = [temperature(alt) for alt in altitudes]

barometric_pressure = [barometric_formula(alt) for alt in altitudes]
pressures_dry = [pressure_dry(alt) for alt in altitudes]
pressures_moist = [pressure_moist(alt,RH) for alt in altitudes]
perc_diff = [(pressures_moist[i] - pressures_dry[i]) / pressures_dry[i] * 100 for i in range(len(pressures_dry))]

moist_molar_mass_alt = [moist_molar_mass(alt,RH) for alt in altitudes]
water_fractions_alt = [water_molar_fraction(alt,RH) for alt in altitudes]

Ts = np.linspace(273.15, 373.15, 100)  # Temperature range from 0 to 100 °C
water_fractions = [water_fraction_T(T, RH) for T in Ts]
molar_mass = [moist_molar_mass_T(T, RH) for T in Ts]
sat_vapor_pressures_T = [sat_vapor_pressure(T) for T in Ts]
sat_vapor_pressures_alt = [sat_vapor_pressure(temperature(alt)) for alt in altitudes]

# Plot the graph
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(8, 4))
plt.yscale("log")
#plt.plot(altitudes, water_fractions_alt, color="black", lw=2, label="water fraction")
#plt.plot(altitudes, pressures_dry, color="red", lw=2, label="dry pressure")
#plt.plot(altitudes, pressures_moist, color="green", lw=2, label="moist pressure")
#plt.plot(altitudes, barometric_pressure, color="blue", lw=2, label="barometric pressure")
plt.plot(altitudes, perc_diff, color="orange", lw=2, label="pressure difference")
#plt.plot(altitudes, sat_vapor_pressures_alt, color="red", lw=2, label="saturation vapor pressure")
plt.plot(altitudes, moist_molar_mass_alt, color="blue", lw=2, label="molar mass")
#plt.plot(Ts, water_fractions, color="black", lw=2, label="water fraction")
#plt.plot(Ts, sat_vapor_pressures, color="blue", lw=2, label="molar mass")
plt.xlabel("altitude (m)")
plt.ylabel("Molar mass (kg/mol)")
plt.legend()
plt.grid(True)
plt.show()