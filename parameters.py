Parameters = {
    'V1': 2.72,  # [mM/s], Pettersson 1988
    'V6': 1.6,  # [mM/s], Pettersson 1988
    'V9': 0.32,  # [mM/s], Pettersson 1988
    'V13': 8.0,  # [mM/s], Pettersson 1988
    'V16': 2.8,  # [mM/s], Pettersson 1988
    'Vst': 0.32,  # [mM/s], Pettersson 1988
    'Vx': 2.0,  # [mM/s], Pettersson 1988
    'Km1': 0.02,  # [mM], Pettersson 1988
    'Km6': 0.03,  # [mM], Pettersson 1988
    'Km9': 0.013,  # [mM], Pettersson 1988
    'Km131': 0.05,  # [mM], Pettersson 1988
    'Km132': 0.05,  # [mM], Pettersson 1988
    'Km161': 0.014,  # [mM], Pettersson 1988
    'Km162': 0.3,  # [mM], Pettersson 1988
    'Kmst1': 0.08,  # [mM], Pettersson 1988
    'Kmst2': 0.08,  # [mM], Pettersson 1988
    'Kpga': 0.25,  # [mM], Pettersson 1988
    'Kgap': 0.075,  # [mM], Pettersson 1988
    'Kdhap': 0.077,  # [mM], Pettersson 1988
    'Kpi': 0.63,  # [mM], Pettersson 1988
    'Kpxt': 0.74,  # [mM], Pettersson 1988
    'Ki11': 0.04,  # [mM], Pettersson 1988
    'Ki12': 0.04,  # [mM], Pettersson 1988
    'Ki13': 0.075,  # [mM], Pettersson 1988
    'Ki14': 0.9,  # [mM], Pettersson 1988
    'Ki15': 0.07,  # [mM], Pettersson 1988
    'Ki61': 0.7,  # [mM], Pettersson 1988
    'Ki62': 12.0,  # [mM], Pettersson 1988
    'Ki9': 12.0,  # [mM], Pettersson 1988
    'Ki131': 2.0,  # [mM], Pettersson 1988
    'Ki132': 0.7,  # [mM], Pettersson 1988
    'Ki133': 4.0,  # [mM], Pettersson 1988
    'Ki134': 2.5,  # [mM], Pettersson 1988
    'Ki135': 0.4,  # [mM], Pettersson 1988
    'Kist': 10.0,  # [mM], Pettersson 1988
    'Kast1': 0.1,  # [mM], Pettersson 1988
    'Kast2': 0.02,  # [mM], Pettersson 1988
    'Kast3': 0.02,  # [mM], Pettersson 1988
    'kRE': 800000000.0,  # Rapid Equilibrium speed
    'q2': 0.00031,  # [], Pettersson 1988
    'q3': 16000000.0,  # [], Pettersson 1988
    'q4': 22.0,  # [], Pettersson 1988
    'q5': 7.1,  # [1/mM]], Pettersson 1988
    'q7': 0.084,  # [], Pettersson 1988
    'q8': 13.0,  # [1/mM]], Pettersson 1988
    'q10': 0.85,  # [], Pettersson 1988
    'q11': 0.4,  # [], Pettersson 1988
    'q12': 0.67,  # [], Pettersson 1988
    'q14': 2.3,  # [], Pettersson 1988
    'q15': 0.058,  # [], Pettersson 1988
    'CO2': 0.2,  # [mM], Pettersson 1988
    'Cp': 15.0,  # [mM], Pettersson 1988
    'Ca': 0.5,  # [mM], Pettersson 1988
    'CN': 0.5,  # [mM], Pettersson 1988
    'Pext': 0.5,  # [mM], Pettersson 1988
    'pHmedium': 7.6,  # [], Pettersson 1988
    'pHstroma': 7.9,  # [], Pettersson 1988
    'protonsStroma': 1.2589254117941661e-05,  # [mM], Pettersson 1988
    'NADPH': 0.21,  # [mM], Pettersson 1988
    'NADP': 0.29}  # [mM], Pettersson 1988

ParametersNADPH = Parameters.copy()
ParametersNADPH.pop("NADPH")
ParametersNADPH.pop("NADP")
ParametersNADPH["Vnadph"] = 2.816
ParametersNADPH["Kmnadph"] = 0.19
