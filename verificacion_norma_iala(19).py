def calcular_boya_iala(D1, H1, D2, H2, Dp, t, rho_agua, m_estructura, m_lastre, m_equipos, rho_pu, V_viento, V_corriente):
    """
    Backend de cálculo hidrostático, estabilidad y verificación IALA G1006 / G1099.
    """
    R1 = D1 / 2.0
    R2 = D2 / 2.0
    
    # Dimensiones internas del cilindro superior
    Hi = H1 - t
    Di = D1 - (2 * t)
    
    # Volúmenes del cilindro superior
    Vext_cil = (np.pi * (D1**2) / 4.0) * H1
    Vint_cil = (np.pi * (Di**2) / 4.0) * Hi
    V_cil_neto = Vext_cil - Vint_cil
    
    # Dimensiones y volúmenes del tronco de cono inferior
    Re = R1
    Ri = R2
    re = Di / 2.0
    Dii = D2 - (2 * t)
    ri = Dii / 2.0
    H2i = H2 - t
    
    Vext_cono = (1.0 / 3.0) * np.pi * H2 * (Re**2 + Ri**2 + Re * Ri)
    Vint_cono = (1.0 / 3.0) * np.pi * H2i * (re**2 + ri**2 + re * ri)
    V_cono_neto = Vext_cono - Vint_cono
    
    # Volumen del pasador externo e interno
    Ht = H1 + H2
    Vpe = (np.pi * (Dp**2) / 2.0) * t
    Vpi = (np.pi * (Dp**2) / 4.0) * (Ht - (2 * t))
    
    # Volumen total del casco HDPE
    V_hdpe = V_cil_neto + V_cono_neto - Vpe
    m_hdpe = V_hdpe * 950.0  
    
    # Volumen y masa de poliuretano (PU)
    V_pu = Vint_cil + Vint_cono - Vpi
    m_pu = V_pu * rho_pu
    
    # Volumen total del casco
    V_total_casco = V_hdpe + V_pu
    
    # Masa Total actual de la boya
    m_total = m_hdpe + m_pu + m_estructura + m_lastre + m_equipos
    
    # Volumen sumergido real
    v_sum = m_total / rho_agua
    
    # 1. Criterio de Reserva de Flotabilidad IALA G1099 (2.0 <= V_total / V_sum <= 2.5)
    rf = V_total_casco / v_sum
    cumple_rf = 2.0 <= rf <= 2.5
    
    # Cálculo exacto para el Punto Óptimo IALA (RF = 2.25)
    v_sum_optimo_225 = V_total_casco / 2.25
    m_total_optima_225 = v_sum_optimo_225 * rho_agua
    lastre_adicional_optimo = max(0.0, m_total_optima_225 - m_total)

    # 2. Determinación del Calado T mediante root_scalar
    resultado_raiz = root_scalar(
        ecuacion_equilibrio_arquimedes,
        args=(m_total, rho_agua, D1, D2, H2),
        bracket=[0.01, Ht]
    )
    
    if resultado_raiz.converged:
        T = resultado_raiz.root
    else:
        T = H2 * (v_sum / Vext_cono) if v_sum <= Vext_cono else H2 + ((v_sum - Vext_cono) / (np.pi * (R1**2)))

    if T <= H2:
        caso_calado = f"Caso A: El calado se encuentra dentro del tronco de cono (T = {T:.4f} m <= H2 = {H2} m)"
    else:
        caso_calado = f"Caso B: El calado supera el cono y entra en el cilindro (T = {T:.4f} m > H2 = {H2} m)"
        
    # 3. Estabilidad Metacéntrica (GM) y Criterio IALA
    kg = 0.5 * Ht * (1.0 - (m_lastre / (m_total + 1e-6))) 
    kb = (0.5 * T) if (T <= H2) else (H2 * 0.4 + (T - H2) * 0.5)
    
    R_wl = R1 if T > H2 else (R2 + (R1 - R2) * (T / H2))
    iwl = (np.pi * (R_wl**4)) / 4.0
    bm = iwl / (v_sum + 1e-6)
    km = kb + bm
    gm = km - kg
    
    cumple_gm = 0.6 <= gm <= 1.2 and gm > 0
    
    # 4. Análisis Riguroso de Cargas Ambientales y Ángulo de Escora Estática
    # Altura de la torre expuesta al viento estimada en 2.0 m sobre el casco
    h_torre = 2.0 
    area_viento = D1 * h_torre  
    area_corriente = Vext_cono + (D1 * T * 0.5)
    
    # Velocidades en m/s (V_viento en km/h convertida a m/s; V_corriente ya ingresa en m/s)
    v_viento_ms = V_viento * (1000.0 / 3600.0)
    v_corr_ms = V_corriente  
    
    F_viento = 0.5 * 1.225 * (v_viento_ms**2) * 1.2 * area_viento
    F_corriente = 0.5 * rho_agua * (v_corr_ms**2) * 1.0 * area_corriente
    F_total = F_viento + F_corriente
    
    # Brazo de momento del viento (desde la línea de flotación hasta la mitad de la torre expuesta)
    z_centro_viento = Ht + (h_torre / 2.0)
    z_centro_corriente = T / 2.0
    brazo_momento = (F_viento * z_centro_viento + F_corriente * z_centro_corriente) / (F_total + 1e-6)
    
    # Momento escorante total
    momento_escorante = F_total * brazo_momento
    
    # Desplazamiento en Newtons
    desplazamiento_N = m_total * 9.81
    
    # Cálculo estricto del ángulo de escora (tan(theta) = Momento / (Desplazamiento * GM))
    tan_theta = momento_escorante / (desplazamiento_N * gm + 1e-6)
    theta_escora = np.degrees(np.arctan(tan_theta))
    cumple_escora = theta_escora <= 15.0

    resultados = {
        "m_hdpe": m_hdpe,
        "m_pu": m_pu,
        "m_total": m_total,
        "v_sum": v_sum,
        "V_total_casco": V_total_casco,
        "rf": rf,
        "cumple_rf": cumple_rf,
        "T": T,
        "caso_calado": caso_calado,
        "kg": kg,
        "km": km,
        "gm": gm,
        "cumple_gm": cumple_gm,
        "theta_escora": theta_escora,
        "cumple_escora": cumple_escora,
        "F_total": F_total,
        "v_sum_optimo_225": v_sum_optimo_225,
        "m_total_optima_225": m_total_optima_225,
        "lastre_adicional_optimo": lastre_adicional_optimo
    }
    return resultados
