import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

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


# Generate altitudes from 0 to 85,000 meters
altitudes = np.linspace(0, 80000, 500)

# Calculate temperatures for each altitude
temperatures = [temperature(alt) for alt in altitudes]

"""# Plot the graph
plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, temperatures, color="red",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Temperature (K)")
plt.grid(True)
plt.show()"""


def pressure(altitude):
    R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
    g = 9.80665  # Standard gravity in m/s^2
    m = 28.9647e-3  # Molar mass of air in kg/mol
    P0 = 101325  # Sea level standard atmospheric pressure in Pa

    integral, _ = quad(lambda h: 1 / temperature(h), 0, altitude, limit=100, points=[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852])
    pressure = P0 * np.exp(-m*g/R * integral)

    return pressure

pressures = [pressure(alt) for alt in altitudes]

# Plot the graph
plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, pressures, color="black",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Pressure (Pa)")
plt.grid(True)
plt.show()