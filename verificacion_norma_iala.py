import streamlit as st
import numpy as np

# ==========================================
# 1. MÓDULO BACKEND: LÓGICA DE CÁLCULO E IALA
# ==========================================

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
    m_hdpe = V_hdpe * 950.0  # Densidad HDPE ~950 kg/m3
    
    # Volumen y masa de poliuretano (PU)
    V_pu = Vint_cil + Vint_cono - Vpi
    m_pu = V_pu * rho_pu
    
    # Volumen total del casco
    V_total_casco = V_hdpe + V_pu
    
    # Masa Total de la boya
    m_total = m_hdpe + m_pu + m_estructura + m_lastre + m_equipos
    
    # Volumen sumergido (Desplazamiento)
    v_sum = m_total / rho_agua
    
    # 1. Criterio de Reserva de Flotabilidad IALA G1099 (2.0 <= V_total / V_sum <= 2.5)
    rf = V_total_casco / v_sum
    cumple_rf = 2.0 <= rf <= 2.5
    
    # 2. Determinación del Calado T (Iterativo / Condicional)
    # Volumen del cono exterior
    if v_sum <= Vext_cono:
        # El calado está dentro del cono
        # Aproximación por método iterativo o fórmula cúbica simplificada
        T = (3.0 * v_sum / (np.pi * (Re**2 + Ri*Re + Ri**2))) ** (1/3.0) # Simplificado para cono invertido/frustum
        caso_calado = "Caso A: El calado se encuentra dentro del tronco de cono (T <= H2)"
    else:
        # El calado supera el cono y entra al cilindro
        v_rem = v_sum - Vext_cono
        Area_cil = np.pi * (R1**2)
        T = H2 + (v_rem / Area_cil)
        caso_calado = "Caso B: El calado supera el cono y entra en el cilindro (T > H2)"
        
    # 3. Estabilidad Metacéntrica (GM) y Criterio IALA (0.6 <= GM <= 1.2 m, GM > 0)
    # Cálculo simplificado de altura metacéntrica KM y KG
    # Asumiendo centro de gravedad KG estimado según lastre y equipos
    kg = 0.5 * Ht * (1.0 - (m_lastre / (m_total + 1e-6))) # Estimación adaptativa
    kb = 0.6 * T # Centro de carena aproximado
    
    # Inercia del plano de flotación Iwl = pi * R^4 / 4
    R_wl = R1 if T > H2 else (R2 + (R1 - R2) * (T / H2))
    iwl = (np.pi * (R_wl**4)) / 4.0
    bm = iwl / (v_sum + 1e-6)
    km = kb + bm
    gm = km - kg
    
    cumple_gm = 0.6 <= gm <= 1.2 and gm > 0
    
    # 4. Análisis de Cargas Ambientales y Escora (IALA G1099: theta <= 10° - 15°)
    # Área expuesta al viento y corriente
    area_viento = D1 * 2.0  # Altura de torre estimada 2m
    area_corriente = Vext_cono + (D1 * T * 0.5)
    
    # Fuerzas de arrastre (Drag)
    F_viento = 0.5 * 1.225 * (V_viento * (1000/3600))**2 * 1.2 * area_viento
    F_corriente = 0.5 * 1025 * (V_corriente * (1000/3600))**2 * 1.0 * area_corriente
    F_total = F_viento + F_corriente
    
    # Momento escorante y ángulo de escora (theta en grados)
    py_escora_rad = np.arctan(F_total * (Ht * 0.5) / (m_total * 9.81 * gm + 1e-6))
    theta_escora = np.degrees(py_escora_rad)
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
        "F_total": F_total
    }
    return resultados

# ==========================================
# 2. MÓDULO FRONTEND: INTERFAZ CON STREAMLIT
# ==========================================

st.set_page_config(page_title="Validador IALA G1006/G1099 - Boyas Oceanográficas", layout="wide")

st.title("🌊 Validador de Normativa IALA para Boyas Oceanográficas")
st.markdown("""
Esta aplicación evalúa los parámetros de diseño geométrico, masa, **reserva de flotabilidad**, **calado ($T$)**, **estabilidad metacéntrica ($GM$)** y **ángulo de escora** conforme a las recomendaciones de la **normativa IALA G1006 / G1099**.
""")

st.sidebar.header("⚙️ Parámetros de Diseño y Geometría")
D1 = st.sidebar.number_input("Diámetro Cilindro Superior D1 (m)", value=2.0, step=0.1)
H1 = st.sidebar.number_input("Altura Cilindro Superior H1 (m)", value=1.04, step=0.01)
D2 = st.sidebar.number_input("Diámetro Menor Cono D2 (m)", value=1.2, step=0.1)
H2 = st.sidebar.number_input("Altura Cono Inferior H2 (m)", value=0.6, step=0.01)
Dp = st.sidebar.number_input("Diámetro Tubo Pasante Dp (m)", value=0.1, step=0.01)
t = st.sidebar.number_input("Espesor de pared HDPE t (m)", value=0.012, step=0.001)

st.sidebar.header("⚖️ Masas y Condiciones Ambientales")
rho_agua = st.sidebar.number_input("Densidad Agua de Mar (kg/m³)", value=1025.0, step=5.0)
m_estructura = st.sidebar.number_input("Masa Estructura/Torre (kg)", value=100.0, step=10.0)
m_lastre = st.sidebar.number_input("Masa Lastre Contrapeso (kg)", value=400.0, step=50.0)
m_equipos = st.sidebar.number_input("Masa Sensores y Equipos (kg)", value=99.63, step=5.0)
rho_pu = st.sidebar.number_input("Densidad Poliuretano (kg/m³)", value=50.0, step=5.0)

V_viento = st.sidebar.number_input("Velocidad del Viento (km/h)", value=100.4, step=5.0)
V_corriente = st.sidebar.number_input("Velocidad de la Corriente (km/h)", value=7.17, step=0.5)

if st.button("🚀 Ejecutar Análisis y Validación IALA", type="primary"):
    res = calcular_boya_iala(D1, H1, D2, H2, Dp, t, rho_agua, m_estructura, m_lastre, m_equipos, rho_pu, V_viento, V_corriente)
    
    st.header("📊 Resultados del Análisis Hidrostático")
    
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Masa Total ($M_{total}$)", f"{res['m_total']:.2f} kg")
    col2.metric("Volumen Sumergido ($V_{sum}$)", f"{res['v_sum']:.4f} m³")
    col3.metric("Calado Estimado ($T$)", f"{res['T']:.4f} m")
    col4.metric("Factor Reserva ($RF$)", f"{res['rf']:.2f}")
    
    st.info(f"📍 **Condición Operativa de Calado:** {res['caso_calado']}")
    
    st.markdown("---")
    st.subheader("🚨 Alertas de Cumplimiento Normativo IALA G1099")
    
    # 1. Alerta Reserva de Flotabilidad
    if res['cumple_rf']:
        st.success(f"✅ **Reserva de Flotabilidad CUMPLE IALA:** El factor de reserva es de **{res['rf']:.2f}** (Rango aceptable IALA: 2.0 – 2.5).")
    else:
        if res['rf'] > 2.5:
            st.error(f"❌ **NO CUMPLE (Exceso de flotabilidad / Boyante en exceso):** El factor actual es **{res['rf']:.2f}**. La masa es insuficiente para el volumen del casco. Se recomienda **aumentar el lastre inferior** para descender el factor al rango óptimo (2.0 - 2.5).")
        else:
            st.error(f"❌ **NO CUMPLE (Insuficiencia de flotabilidad):** El factor de reserva **{res['rf']:.2f}** está por debajo del mínimo de 2.0 requerido por IALA.")

    # 2. Alerta Estabilidad Metacéntrica (GM)
    if res['cumple_gm']:
        st.success(f"✅ **Estabilidad Metacéntrica (GM) CUMPLE IALA:** $GM$ = **{res['gm']:.3f} m** (Rango recomendado: 0.6 m – 1.2 m). La boya es estable y auto-adrizante sin exceso de rigidez.")
    else:
        st.warning(f"⚠️ **ADVERTENCIA GM FUERA DE RANGO ÓPTIMO:** $GM$ = **{res['gm']:.3f} m**. Un valor muy bajo compromete el adrizamiento, y un valor muy alto (>1.2 m) genera movimientos violentos que fatigan la torre y dañan equipos.")

    # 3. Alerta Ángulo de Escora
    if res['cumple_escora']:
        st.success(f"✅ **Ángulo de Escora CUMPLE IALA:** $\\theta$ = **{res['theta_escora']:.2f}°** (Límite operativo seguro: $\\le 15°$).")
    else:
        st.error(f"❌ **NO CUMPLE (Escora excesiva):** El ángulo de inclinación calculado es de **{res['theta_escora']:.2f}°**, superando el límite máximo recomendado de $15^\circ$ por la normativa IALA ante tormenta de diseño. Esto afectará los paneles solares y sensores.")