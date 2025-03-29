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
            return temperature - lapse_rates[i] * altitude
        else:
            temperature = temperature - lapse_rates[i] * altitude

    # If altitude is above the highest defined layer
    return temperature

print(temperature(0))
print(temperature(11000))
    

