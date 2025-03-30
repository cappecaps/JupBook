import matplotlib.pyplot as plt
import numpy as np

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
altitudes = np.linspace(0, 85000, 500)

# Calculate temperatures for each altitude
temperatures = [temperature(alt) for alt in altitudes]

# Plot the graph
plt.rcParams.update({'font.size': 9})
plt.figure(figsize=(8, 4))
plt.plot(altitudes, temperatures, color="red",lw=2)
plt.xlabel("Altitude (m)")
plt.ylabel("Temperature (K)")
plt.grid(True)
plt.show()