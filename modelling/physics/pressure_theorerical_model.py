import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

R = 8.31446  # Specific gas constant for dry air in J/(mol·K)
g0 = 9.80665  # Standard gravity in m/s^2
m_dry = 28.9656e-3  # Molar mass of air in kg/mol
m_water = 18.01528e-3  # Molar mass of water in kg/mol
T0 = 288.15  # MSL standard temperature in Kelvin
T1 = 298.15

p0 = 101325  # MSL standard atmospheric pressure in Pa

def ISA_temperature(altitude):
    # Define base altitudes and temperatures for each layer
    base_altitudes = [0, H1, 20, 32, 47, 51, 71, 84.852]
    lapse_rates = [6.5, 0.0, -1.0, -2.8, 0.0, 2.8, 2.0, 0.0]

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
    integral, err = quad(lambda h: 1 / ISA_temperature_up_to_1000(h), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    pressure = p0 * np.exp(-m_dry * g0 / R * 1000 * integral)

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
        if h > 20:
            return 0.0
        T = ISA_temperature_up_to_1000(h)
        p_dry = pressure_dry(h)
    
    else:
        p_dry = p

    p_vap = vapor_pressure(T)

    # partial pressure of water vapor
    p_water = RH * p_vap

    # molar fraction of water vapor
    f_water = p_water / p_dry

    return f_water



def air_avg_molar_mass(h,RH):
    #altitude in km
    if h < 20:
        f_water = water_molar_fraction(RH,h=h)
        return (1 - f_water) * m_dry + f_water * m_water
    elif h < 85:
        return m_dry
    else:
        return m_dry * np.exp(-0.002*(h-85))



def calc_chi(RH,maxh=1000):
    if T_surf is None:
        T_surf = 288.15
    if L0 is None:
        L0 = 6.5

    H1 = ( T_surf - 216.65 ) / L0
    integral_frac, _ = quad(lambda z: water_molar_fraction(RH,h=z) * pressure_dry(z) / ISA_temperature_up_to_1000(z), 0, 20, limit=100, points=[0, H1, 20]) #up to 20 km since above that the water fraction is 0
    integral_norm, _ = quad(lambda z: pressure_dry(z) / ISA_temperature_up_to_1000(z), 0, maxh, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    
    chi = 1 - ((m_dry - m_water) / m_dry) * (integral_frac / integral_norm)

    return chi


def pressure_moist(h,RH,chi):
    integral, _ = quad(lambda z: water_molar_fraction(RH,h=z) / ISA_temperature(z), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852])
    p_moist = chi * pressure_dry(h) * np.exp(g0 / R * (m_dry - m_water) * 1000 * integral)

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

def pressure_dry_up_to_1000(h,T_surf=288.15,L0=6.5):
    H1 = ( T_surf - 216.65 ) / L0
    if H1 > 11:
        raise TypeError("Surface temperature is too high or lapse rate is too small")
        
    integral, _ = quad(lambda h: 1 / ISA_temperature_up_to_1000(h,T_surf,L0), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    pressure = p0 * np.exp(-m_dry * g0 / R * 1000 * integral)

    return pressure

def moist_molar_mass_T(T,RH):
    f_water = water_fraction_T(T,RH)

    # Calculate the molar mass of moist air
    m_moist = (1 - f_water) * m_dry + f_water * m_water

    return m_moist

def ISA_temperature_up_to_1000(h,T_surf=None,L0=None):
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


def pressure_moist_WGS84(h,lat,RH=0.0,chi=1.0,T_surf=None,L0=None):
    if T_surf is None:
        T_surf = 288.15
    if L0 is None:
        L0 = 6.5

    H1 = ( T_surf - 216.65 ) / L0
    integral, _ = quad(lambda z: g_WGS84_altitude(z,lat) * air_avg_molar_mass(z,RH) / ISA_temperature_up_to_1000(z,T_surf,L0), 0, h, limit=100, points=[0, H1, 20, 32, 47, 51, 71, 84.852, 107.41])
    p_moist = chi * p0 * np.exp(- 1000 * integral / R)

    return p_moist


RH = 1.0  # Relative humidity

# Generate altitudes from 0 to 800 km
altitudes = np.linspace(0, 150, 500)

# Calculate temperatures for each altitude
#temperatures = [ISA_temperature_up_to_1000(alt) for alt in altitudes]

pressure_dry_arr = np.array([pressure_dry(alt)/100 for alt in altitudes])   # divided by 100 to return hPa

"""moist_molar_mass_alt = [moist_molar_mass(alt,RH) for alt in altitudes]
water_fractions_alt = [water_molar_fraction(alt,RH) for alt in altitudes]

Ts = np.linspace(273.15, 373.15, 100)  # Temperature range from 0 to 100 °C
water_fractions = [water_fraction_T(T, RH) for T in Ts]
molar_mass = [moist_molar_mass_T(T, RH) for T in Ts]
sat_vapor_pressures_T = [vapor_pressure(T) for T in Ts]
sat_vapor_pressures_alt = [vapor_pressure(temperature(alt)) for alt in altitudes]"""

chi = calc_chi(RH,altitudes[-1])
print(f'chi: {chi}')
pressure_final_poles_arr = np.array([pressure_moist_WGS84(alt,0,RH=1.0,chi=chi)/100 for alt in altitudes])   # divided by 100 to return hPa
pressure_final_dry_poles_arr = np.array([pressure_moist_WGS84(alt,0,RH=0.0,chi=1.0)/100 for alt in altitudes])   # divided by 100 to return hPa
#pressure_final_equat_arr = np.array([pressure_moist_WGS84(alt,90,RH=1.0,chi=1.0)/100 for alt in altitudes])   # divided by 100 to return hPa
diff_final_poles_dry = pressure_final_poles_arr - pressure_dry_arr
diff_final_dry_poler_dry = pressure_final_dry_poles_arr - pressure_dry_arr
fig, ax1 = plt.subplots(figsize=(7, 3))
ax2 = ax1.twinx()
ax2.axhline(0,color='r',alpha=0.5,linestyle='dashed')
ax1.plot(altitudes, pressure_dry_arr,lw=2,label='$p_{dry}$',color='darkred')
ax1.plot(altitudes, pressure_final_poles_arr,lw=2,label='$p_{final}$ (RH=1, poles)',color='lightblue')
ax1.plot(altitudes, pressure_final_dry_poles_arr,lw=2,label='$p_{final}$ (RH=0, poles)',color='red')
ax2.plot(altitudes, diff_final_poles_dry,lw=2,label=r'$p_{final}-p_{dry}$',color='lightblue',linestyle='dashed')
ax2.plot(altitudes, diff_final_dry_poler_dry,lw=2,label=r'$p_{final}-p_{dry}$ (RH=0)',color='red',linestyle='dashed')
plt.xlim(0)
ax1.set_yscale('log')
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel("Pressure (hPa)")
ax2.set_ylabel("Pressure difference (hPa)")
ax1.legend(loc=(0.65,0.2))
ax2.legend(loc=(0.65,0.60))
plt.show()

"""#Plot single axis
fig, ax1 = plt.subplots(figsize=(7, 4))
altitudes = np.linspace(0,800,50)
#temperatures_smaller = [ISA_temperature_up_to_1000(alt) for alt in altitudes]
#temperatures_smaller_mod = [ISA_temperature_up_to_1000(alt,T_surf=298.15,L0=5.0) for alt in altitudes]
gravity_eq = [g_WGS84_altitude(alt,0) for alt in altitudes]
gravity_pol = [g_WGS84_altitude(alt,90) for alt in altitudes]
ax1.plot(altitudes, gravity_eq,lw=2,color='darkred',label='Equator')
ax1.plot(altitudes, gravity_pol,lw=2,color='lightblue',label='Pole')
ax1.set_xlabel("Altitude (km)")
ax1.set_ylabel("Gravity (m/s²)")
ax1.set_xlim(0)
plt.show()"""