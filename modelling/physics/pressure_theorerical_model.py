import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
g0 = 9.80665  # Standard gravity in m/s^2
m_dry = 28.9656e-3  # Molar mass of air in kg/mol
m_water = 18.01528e-3  # Molar mass of water in kg/mol
T0 = 288.15  # MSL standard temperature in Kelvin
T1 = 298.15
H1 = ( T1 - 216.65 ) / 0.0065

p0 = 101325  # MSL standard atmospheric pressure in Pa

def temperature(altitude):
    # Define base altitudes and temperatures for each layer
    base_altitudes = [0, H1, 20000, 32000, 47000, 51000, 71000, 84852]
    lapse_rates = [0.0065, 0.0, -0.001, -0.0028, 0.0, 0.0028, 0.002, 0.0]

    temperature = T1

    # Find the layer corresponding to the altitude
    for i in range(len(base_altitudes) - 1):
        if altitude <= base_altitudes[i+1]:
            return temperature - lapse_rates[i] * (altitude - base_altitudes[i])
        else:
            temperature = temperature - lapse_rates[i] * (base_altitudes[i+1]-base_altitudes[i])

    # If altitude is above the highest defined layer
    return temperature


def pressure_barometric(h):
    pressure = p0 * np.exp(-m_dry * g0 * h / (R * T0))

    return pressure

def pressure_dry(h):
    integral, err = quad(lambda h: 1 / temperature(h), 0, h, limit=100, points=[0, H1, 20000, 32000, 47000, 51000, 71000, 84852])
    pressure = p0 * np.exp(-m_dry * g0 / R * integral)

    return pressure


def vapor_pressure(T):
    # Tetens equation
    a = 610.78
    b = 17.27
    c = 237.3  
    p_vap = a * np.exp((b * (T - 273.15)) / (T - 273.15 + c))   
    return p_vap


def water_molar_fraction(RH,T=None,h=None,p=p0):
    if h is not None:
        if h > 20000:
            return 0.0
        T = temperature(h)
        p_dry = pressure_dry(h)
    
    else:
        p_dry = p

    p_vap = vapor_pressure(T)

    # partial pressure of water vapor
    p_water = RH * p_vap

    # molar fraction of water vapor
    f_water = p_water / p_dry

    return f_water



def moist_molar_mass(altitude,RH):
    f_water = water_molar_fraction(altitude,RH)

    # Calculate the molar mass of moist air
    m_moist = (1 - f_water) * m_dry + f_water * m_water

    return m_moist


def calc_chi(RH):
    integral_frac, _ = quad(lambda z: water_molar_fraction(RH,h=z) * pressure_dry(z) / temperature(z), 0, 84852, limit=100, points=[0, H1, 20000, 32000, 47000, 51000, 71000, 84852])
    integral_norm, _ = quad(lambda z: pressure_dry(z) / temperature(z), 0, 84852, limit=100, points=[0, H1, 20000, 32000, 47000, 51000, 71000, 84852])
    

    chi = 1 - ((m_dry - m_water) / m_dry) * (integral_frac / integral_norm)

    return chi


def pressure_moist(h,RH,chi):
    integral, _ = quad(lambda z: water_molar_fraction(RH,h=z) / temperature(z), 0, h, limit=100, points=[0, H1, 20000, 32000, 47000, 51000, 71000, 84852])

    p_moist = chi * pressure_dry(h) * np.exp(g0 / R * (m_dry - m_water) * integral)

    return p_moist



def water_fraction_T(T,RH):
    # Calculate the vapor pressure
    p_vap = vapor_pressure(T)

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

# Generate altitudes from 0 to 85 km
altitudes = np.linspace(0, 85, 500)

# Calculate temperatures for each altitude
temperatures = [temperature(alt) for alt in altitudes]

pressure_dry_arr = np.array([pressure_dry(1000*alt)/100 for alt in altitudes])   # divided by 100 to return hPa

"""moist_molar_mass_alt = [moist_molar_mass(alt,RH) for alt in altitudes]
water_fractions_alt = [water_molar_fraction(alt,RH) for alt in altitudes]

Ts = np.linspace(273.15, 373.15, 100)  # Temperature range from 0 to 100 °C
water_fractions = [water_fraction_T(T, RH) for T in Ts]
molar_mass = [moist_molar_mass_T(T, RH) for T in Ts]
sat_vapor_pressures_T = [vapor_pressure(T) for T in Ts]
sat_vapor_pressures_alt = [vapor_pressure(temperature(alt)) for alt in altitudes]"""

chi = calc_chi(RH)
print(f'chi: {chi}')
pressure_moist_arr = np.array([pressure_moist(1000*alt,1.0,chi)/100 for alt in altitudes])   # divided by 100 to return hPa
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