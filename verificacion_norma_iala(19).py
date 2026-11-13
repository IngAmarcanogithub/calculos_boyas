import streamlit as st
import numpy as np
from scipy.optimize import root_scalar
import plotly.graph_objects as go

# ==========================================
# 1. MÓDULO BACKEND: LÓGICA DE CÁLCULO E IALA
# ==========================================

def ecuacion_equilibrio_arquimedes(T, m_total, rho_agua, D1, D2, H2):
    R1 = D1 / 2.0
    R2 = D2 / 2.0
    Vext_cono = (1.0 / 3.0) * np.pi * H2 * (R1**2 + R2**2 + R1 * R2)
    
    if T <= H2:
        R_t = R2 + (R1 - R2) * (T / H2)
        V_sum = (1.0 / 3.0) * np.pi * T * (R2**2 + R_t**2 + R2 * R_t)
    else:
        v_rem = (T - H2) * np.pi * (R1**2)
        V_sum = Vext_cono + v_rem
        
    empuje = V_sum * rho_agua
    return empuje - m_total

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
    h_torre = 2.0 
    area_viento = D1 * h_torre  
    area_corriente = Vext_cono + (D1 * T * 0.5)
    
    v_viento_ms = V_viento * (1000.0 / 3600.0)
    v_corr_ms = V_corriente  
    
    F_viento = 0.5 * 1.225 * (v_viento_ms**2) * 1.2 * area_viento
    F_corriente = 0.5 * rho_agua * (v_corr_ms**2) * 1.0 * area_corriente
    F_total = F_viento + F_corriente
    
    z_centro_viento = Ht + (h_torre / 2.0)
    z_centro_corriente = T / 2.0
    brazo_momento = (F_viento * z_centro_viento + F_corriente * z_centro_corriente) / (F_total + 1e-6)
    
    momento_escorante = F_total * brazo_momento
    desplazamiento_N = m_total * 9.81
    
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

# ==========================================
# 2. MÓDULO VISUALIZACIÓN 3D (PLOTLY)
# ==========================================

def generar_grafico_boya(D1, H1, D2, H2, Dp, T):
    """
    Genera la visualización 3D paramétrica de la boya con sus cotas y el plano de flotación dinámico.
    """
    R1 = D1 / 2.0
    R2 = D2 / 2.0
    Rp = Dp / 2.0
    Ht = H1 + H2

    theta = np.linspace(0, 2 * np.pi, 72)
    z_steps = np.linspace(0, 1, 20)
    THETA, Z_NORM = np.meshgrid(theta, z_steps)

    fig = go.Figure()

    # 1. Tronco de Cono Inferior (z de 0 a H2)
    Z_cono = Z_NORM * H2
    R_cono = R2 + (R1 - R2) * Z_NORM
    X_cono = R_cono * np.cos(THETA)
    Y_cono = R_cono * np.sin(THETA)

    fig.add_trace(go.Surface(
        x=X_cono, y=Y_cono, z=Z_cono,
        colorscale=[[0, '#FFEE00'], [1, '#FFD700']],
        showscale=False, opacity=0.82,
        name='Tronco de Cono (H2)', hoverinfo='skip'
    ))

    # 2. Cilindro Superior (z de H2 a Ht)
    Z_cil = H2 + Z_NORM * H1
    X_cil = R1 * np.cos(THETA)
    Y_cil = R1 * np.sin(THETA)

    fig.add_trace(go.Surface(
        x=X_cil, y=Y_cil, z=Z_cil,
        colorscale=[[0, '#FFF200'], [1, '#FFE600']],
        showscale=False, opacity=0.75,
        name='Cilindro Superior (H1)', hoverinfo='skip'
    ))

    # 3. Tubo Pasante Central (z de 0 a Ht)
    Z_tubo = Z_NORM * Ht
    X_tubo = Rp * np.cos(THETA)
    Y_tubo = Rp * np.sin(THETA)

    fig.add_trace(go.Surface(
        x=X_tubo, y=Y_tubo, z=Z_tubo,
        colorscale=[[0, '#D4AF37'], [1, '#C5A000']],
        showscale=False, opacity=0.95,
        name=f'Tubo Pasante (Dp={Dp}m)', hoverinfo='skip'
    ))

    # 4. Aristas Estructurales (Bordes negros técnicos)
    for r_ring, z_ring in [(R2, 0), (R1, H2), (R1, Ht), (Rp, Ht), (Rp, 0)]:
        fig.add_trace(go.Scatter3d(
            x=r_ring * np.cos(theta),
            y=r_ring * np.sin(theta),
            z=np.full_like(theta, z_ring),
            mode='lines',
            line=dict(color='#333333', width=3),
            showlegend=False, hoverinfo='skip'
        ))

    # 5. Línea de Flotación Real (Intersección en Cono o en Cilindro)
    T_visual = min(T, Ht) # Evita desbordamiento gráfico si sobrepasa Ht
    if T_visual <= H2:
        R_wl = R2 + (R1 - R2) * (T_visual / H2)
        zona_txt = "Zona Cónica"
    else:
        R_wl = R1
        zona_txt = "Zona Cilíndrica"

    fig.add_trace(go.Scatter3d(
        x=R_wl * np.cos(theta),
        y=R_wl * np.sin(theta),
        z=np.full_like(theta, T_visual),
        mode='lines',
        line=dict(color='#006400', width=7),
        name=f'Línea de Flotación (T = {T:.4f} m)'
    ))

    # 6. Plano de Agua Translúcido
    lim = R1 * 1.45
    x_w = np.linspace(-lim, lim, 10)
    y_w = np.linspace(-lim, lim, 10)
    XW, YW = np.meshgrid(x_w, y_w)
    ZW = np.full_like(XW, T_visual)

    fig.add_trace(go.Surface(
        x=XW, y=YW, z=ZW,
        colorscale=[[0, '#2ECC71'], [1, '#27AE60']],
        showscale=False, opacity=0.38,
        name='Superficie del Agua', hoverinfo='skip'
    ))

    # 7. Líneas de Cota y Etiquetas de Parámetros en el Espacio 3D
    # Cota D1 (Superior)
    fig.add_trace(go.Scatter3d(
        x=[-R1, R1], y=[0, 0], z=[Ht + 0.12, Ht + 0.12],
        mode='lines+text',
        line=dict(color='black', width=3),
        text=['', f'D1 = {D1:.2f} m'],
        textposition='top center',
        textfont=dict(size=12, color='black'),
        showlegend=False
    ))

    # Cota D2 (Inferior)
    fig.add_trace(go.Scatter3d(
        x=[-R2, R2], y=[0, 0], z=[-0.12, -0.12],
        mode='lines+text',
        line=dict(color='black', width=3),
        text=['', f'D2 = {D2:.2f} m'],
        textposition='bottom center',
        textfont=dict(size=12, color='black'),
        showlegend=False
    ))

    # Cota H1 y H2 (Lateral Izquierdo)
    x_cota = -R1 * 1.25
    fig.add_trace(go.Scatter3d(
        x=[x_cota, x_cota, x_cota], y=[0, 0, 0], z=[0, H2, Ht],
        mode='lines+markers+text',
        line=dict(color='black', width=3),
        marker=dict(size=4, color='black'),
        text=['', f'H2 = {H2:.2f} m', f'H1 = {H1:.2f} m'],
        textposition='middle left',
        textfont=dict(size=11, color='black'),
        showlegend=False
    ))

    # Flecha Vertical de Calado Resultante T (Lateral Derecho)
    x_calado = R1 * 1.25
    fig.add_trace(go.Scatter3d(
        x=[x_calado, x_calado], y=[0, 0], z=[0, T_visual],
        mode='lines+markers+text',
        line=dict(color='#006400', width=5),
        marker=dict(size=5, color='#006400'),
        text=['Base (z=0)', f'<b>Calado (T) = {T:.4f} m</b><br>({zona_txt})'],
        textposition='middle right',
        textfont=dict(size=12, color='#004d00'),
        showlegend=False
    ))

    # Configuración de Cámara, Proporción Real (aspectmode='data') y Cuadro de Datos
    fig.update_layout(
        title=dict(
            text="Visualización Geométrica y de Calado de la Boya (IALA G1006)",
            x=0.5, font=dict(size=18)
        ),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(title='Altura (m)', showbackground=False, gridcolor='lightgray'),
            aspectmode='data',
            camera=dict(eye=dict(x=0.15, y=-2.1, z=0.45))
        ),
        annotations=[
            dict(
                x=0.01, y=0.02, xref='paper', yref='paper',
                text=(
                    f"<b>Resumen Geométrico:</b><br>"
                    f"• D1: {D1:.2f} m | H1: {H1:.2f} m<br>"
                    f"• D2: {D2:.2f} m | H2: {H2:.2f} m<br>"
                    f"• Dp: {Dp:.2f} m<br>"
                    f"• <b>Calado T: {T:.4f} m</b> ({zona_txt})"
                ),
                showarrow=False, align='left',
                bgcolor='rgba(245, 247, 249, 0.92)',
                bordercolor='#CCCCCC', borderwidth=1, borderpad=8,
                font=dict(size=12, color='#111111')
            )
        ],
        legend=dict(x=0.78, y=0.05, bgcolor='rgba(255,255,255,0.85)', bordercolor='#CCCCCC', borderwidth=1),
        margin=dict(l=10, r=10, t=50, b=10),
        height=620
    )

    return fig

# ==========================================
# 3. MÓDULO FRONTEND: INTERFAZ CON STREAMLIT
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
V_corriente = st.sidebar.number_input("Velocidad de la Corriente (m/s)", value=1.99, step=0.1)

if st.button("🚀 Ejecutar Análisis y Validación IALA", type="primary"):
    res = calcular_boya_iala(D1, H1, D2, H2, Dp, t, rho_agua, m_estructura, m_lastre, m_equipos, rho_pu, V_viento, V_corriente)
    
    # 1. PRIMERO: Resultados del Análisis Hidrostático
    st.header("📊 Resultados del Análisis Hidrostático")
    
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Masa Total ($M_{total}$)", f"{res['m_total']:.2f} kg")
    col2.metric("Volumen Sumergido ($V_{sum}$)", f"{res['v_sum']:.4f} m³")
    col3.metric("Calado Estimado ($T$)", f"{res['T']:.4f} m")
    col4.metric("Factor Reserva ($RF$)", f"{res['rf']:.2f}")
    
    st.info(f"📍 **Condición Operativa de Calado:** {res['caso_calado']}")

    # 2. DESPUÉS: Visualización 3D Interactiva de la Boya
    fig_boya = generar_grafico_boya(D1, H1, D2, H2, Dp, res['T'])
    st.plotly_chart(fig_boya, use_container_width=True)
    
    # ALERTA DINÁMICA DEL PUNTO ÓPTIMO IALA (RF = 2.25)
    st.warning("### 🎯 Alerta de Referencia: Punto Óptimo IALA (Factor de Reserva $RF = 2.25$)")
    st.markdown(f"""
    Para cumplir exactamente con el **Punto Óptimo de Reserva de Flotabilidad ($RF = 2.25$)** de la norma IALA G1099:
    - **Volumen Sumergido Requerido ($V_{{sum}}$):** `{res['v_sum_optimo_225']:.4f} m³`
    - **Masa Total Requerida ($M_{{total}}$):** `{res['m_total_optima_225']:.1f} kg`
    - **Contrapeso / Lastre Adicional a Añadir:** **`+{res['lastre_adicional_optimo']:.1f} kg`** 
      *(calculado desde su masa total actual de {res['m_total']:.1f} kg)*
    """)

    st.markdown("---")
    st.subheader("🚨 Alertas de Cumplimiento Normativo IALA G1099")
    
    # 1. Alerta Reserva de Flotabilidad
    if res['cumple_rf']:
        st.success(f"✅ **Reserva de Flotabilidad CUMPLE IALA:** El factor de reserva es de **{res['rf']:.2f}** (Rango aceptable IALA: 2.0 – 2.5).")
    else:
        if res['rf'] > 2.5:
            st.error(f"❌ **NO CUMPLE (Exceso de flotabilidad / Boyante en exceso):** El factor actual es **{res['rf']:.2f}**. La masa es insuficiente para el volumen del casco. Se recomienda **aumentar el lastre en +{res['lastre_adicional_optimo']:.1f} kg** para alcanzar el óptimo IALA (RF = 2.25).")
        else:
            st.error(f"❌ **NO CUMPLE (Insuficiencia de flotabilidad):** El factor de reserva **{res['rf']:.2f}** está por debajo del mínimo de 2.0 requerido por IALA.")

    # 2. Alerta Estabilidad Metacéntrica (GM)
    if res['cumple_gm']:
        st.success(f"✅ **Estabilidad Metacéntrica (GM) CUMPLE IALA:** $GM$ = **{res['gm']:.3f} m** (Rango recomendado: 0.6 m – 1.2 m). La boya es estable y auto-adrizante sin exceso de rigidez.")
    else:
        st.warning(f"⚠️ **ADVERTENCIA GM FUERA DE RANGO ÓPTIMO:** $GM$ = **{res['gm']:.3f} m**. Un valor muy bajo compromete el adrizamiento, y un valor muy alto (>1.2 m) genera movimientos violentos que fatigas la torre y dañan equipos.")

    # 3. Alerta Ángulo de Escora
    if res['cumple_escora']:
        st.success(f"✅ **Ángulo de Escora CUMPLE IALA:** $\\theta$ = **{res['theta_escora']:.2f}°** (Límite operativo seguro: $\\le 15°$).")
    else:
        st.error(f"❌ **NO CUMPLE (Escora excesiva):** El ángulo de inclinación calculado es de **{res['theta_escora']:.2f}°**, superando el límite máximo recomendado de $15^\\circ$ por la normativa IALA ante tormenta de diseño.")

    st.markdown("---")
    
 
