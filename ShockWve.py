import sys
import math
import csv

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLabel, QPushButton, QLineEdit, QTextEdit, QTabWidget,
    QFileDialog, QMessageBox, QScrollArea, QStackedWidget
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QGuiApplication
from PyQt6.QtWidgets import QCheckBox
from PyQt6.QtWidgets import (
    QTableWidget,
    QTableWidgetItem,
    QComboBox
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

USE_SCALING = True  # Global toggle: scale divergent to match Pe
CYBER_FONT = "Consolas, Monospace"
NEON_CYAN = "#00f5ff"   # C+
NEON_PINK = "#ff40c0"   # C-


# ============================================================
#  Small helper: safe root finder (bisection)
# ============================================================

def bisection(func, a, b, tol=1e-6, max_iter=100):
    fa = func(a)
    fb = func(b)
    if math.isnan(fa) or math.isnan(fb):
        raise RuntimeError("Function returned NaN in bisection.")
    if fa * fb > 0:
        raise RuntimeError("Bisection interval does not bracket a root.")

    for _ in range(max_iter):
        c = 0.5 * (a + b)
        fc = func(c)
        if math.isnan(fc):
            raise RuntimeError("Function returned NaN in bisection.")
        if abs(fc) < tol or 0.5 * (b - a) < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


# ============================================================
#  Core nozzle physics (translation of your MATLAB logic)
# ============================================================

def moc_bell_nozzle_original(P0, g, TR_mm, Me):
    """
    Direct translation of your MATLAB moc_bell_nozzle_original(P0,g,TR)
    Returns xw (mm), yw (mm) for the nozzle wall from MOC.
    No plotting inside, just geometry.
    """

    DtoR = math.pi / 180.0
    RtoD = 180.0 / math.pi
    P_list = []

    A = math.sqrt((g + 1.0) / (g - 1.0))
    B = (g - 1.0) / (g + 1.0)

    def V_PM(M):
        return A * math.atan(math.sqrt(B * (M**2 - 1.0))) - math.atan(
            math.sqrt(M**2 - 1.0)
        )

    # Total PM angle for Me
    T_max_deg = 0.5 * V_PM(Me) * RtoD
    DT = (90.0 - T_max_deg) - math.floor(90.0 - T_max_deg)

    T = [DT * DtoR]
    n = int(round(T_max_deg * 2.0))

    M_list = [None] * (n + 1)
    RR = [None] * (n + 1)
    LR = [None] * (n + 1)
    SL = [None] * (n + 1)

    TR = TR_mm  # mm

    for m in range(2, n + 2):
        theta_m = (DT + (m - 1)) * DtoR
        T.append(theta_m)

        def func_for_M(M):
            if M <= 1.0:
                return 1e6
            try:
                return theta_m - V_PM(M)
            except ValueError:
                return 1e6

        Mguess_low = 1.0
        Mguess_high = max(1.2, 1.5 * Me)

        f1 = func_for_M(Mguess_low)
        f2 = func_for_M(Mguess_high)
        if f1 * f2 > 0:
            Mguess_high = 5.0 * Me

        M_root = bisection(func_for_M, Mguess_low, Mguess_high, tol=1e-6, max_iter=200)
        M_list[m - 1] = M_root

        P_val = TR * math.tan(theta_m)
        P_list.append(P_val)
        RR[m - 1] = -TR / P_val
        LR[m - 1] = math.tan(theta_m + math.asin(1.0 / M_root))
        SL[m - 1] = -RR[m - 1]

    C_plus_lines = []  # right-running
    C_minus_lines = []  # left-running

    P_arr = np.array(P_list, dtype=float)
    RR_arr = np.array(RR[1:], dtype=float)
    LR_arr = np.array(LR[1:], dtype=float)
    SL_arr = np.array(SL[1:], dtype=float)

    Fline = RR_arr[-1]
    numP = len(P_arr)
    x_center = []
    y_center = []
    for c in range(numP - 1):
        x_c = (TR + SL_arr[c] * P_arr[c]) / (SL_arr[c] - Fline)
        y_c = Fline * x_c + TR

        x_center.append(x_c)
        y_center.append(y_c)

        # build dense C+ characteristic
        x_line = np.linspace(0.0, x_c, 50)
        slope = (y_c - TR) / x_c
        y_line = TR + slope * x_line

        C_plus_lines.append((x_line, y_line))

    TM = T_max_deg * DtoR
    xw = [(TR + SL_arr[0] * P_arr[0]) / (SL_arr[0] - math.tan(TM))]
    yw = [math.tan(TM) * xw[0] + TR]





    DTW = math.tan(TM) / (numP - 1)
    s = [math.tan(TM)]
    b_list = [TR]

    for k in range(2, numP):
        s_k = math.tan(TM) - (k - 1) * DTW
        s.append(s_k)
        b_k = yw[k - 2] - s_k * xw[k - 2]
        b_list.append(b_k)

        xw_k = (b_k + SL_arr[k - 1] * P_arr[k - 1]) / (SL_arr[k - 1] - s_k)
        yw_k = s_k * xw_k + b_k
        xw.append(xw_k)
        yw.append(yw_k)

    # C− characteristics: centerline → wall (DENSE)
    for i in range(len(x_center)):
        x0 = x_center[i]
        y0 = y_center[i]
        x1 = xw[i]
        y1 = yw[i]

        # build dense line
        x_line = np.linspace(x0, x1, 50)
        slope = (y1 - y0) / (x1 - x0)
        y_line = y0 + slope * (x_line - x0)

        C_minus_lines.append((x_line, y_line))

    return (
        np.array(xw),
        np.array(yw),
        C_plus_lines,
        C_minus_lines
    )

def rao_bell_nozzle(
    P0, T0, gamma, Pa, Pe, R, F,
    use_cea=False, Cf=None, eps_user=None,
    CR=None, Lstar_m=None,
    include_chamber=False,          # 🔥 ADD
    bell_percent=1.0, N=1500
):
    # ===== HARD INPUT GUARDS =====
    x_ch = np.array([])
    r_ch = np.array([])

    Lc = None
    L_conv = None

    if use_cea and not include_chamber:
        # CEA performance only — no chamber geometry
        Lc = None
        L_conv = None

    x_blend = np.array([])
    r_blend = np.array([])


    if use_cea:
        if Cf is None or eps_user is None or CR is None or Lstar_m is None:
            raise ValueError(
                "CEA mode requires Cf, ε, CR, and L* to be set before running."
            )


    import math
    from scipy.interpolate import interp1d, PchipInterpolator
    from scipy.optimize import brentq

    # ===== Exit Mach =====
    f = lambda M: (1+(gamma-1)/2*M**2)**(-gamma/(gamma-1)) - Pe/P0
    Me = brentq(f, 1.01, 10)

    # ===== Expansion ratio =====
    eps = (1/Me)*((2/(gamma+1))*(1+(gamma-1)/2*Me**2))**((gamma+1)/(2*(gamma-1)))


    # ===== Throat / exit =====
    # ===== Throat / exit =====
    if use_cea:
        # ===== CEA sizing (USER PROVIDED) =====
        At = F / (Cf * P0)
        Rt = math.sqrt(At / math.pi)

        eps = eps_user
        Re = Rt * math.sqrt(eps)
        Ae = math.pi * Re ** 2

        # ===== CEA CHAMBER SIZING =====
        if Lstar_m is None or CR is None:
            raise ValueError(
                "CEA mode is enabled, but L* or CR was not provided."
            )

        Ac = CR * At
        Rc = math.sqrt(Ac / math.pi)  # 🔥 THIS WAS MISSING
        Vc = Lstar_m * At


        # Compute Me for reporting only
        f = lambda M: (1 + (gamma - 1) / 2 * M ** 2) ** (-gamma / (gamma - 1)) - Pe / P0
        Me = brentq(f, 1.01, 10)

    else:
        # ===== CLASSIC ISENTROPIC SIZING =====
        f = lambda M: (1 + (gamma - 1) / 2 * M ** 2) ** (-gamma / (gamma - 1)) - Pe / P0
        Me = brentq(f, 1.01, 10)

        eps = (1 / Me) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * Me ** 2)) ** ((gamma + 1) / (2 * (gamma - 1)))

        Cf = math.sqrt(
            (2 * gamma ** 2 / (gamma - 1)) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)) *
            (1 - (101325 / P0) ** ((gamma - 1) / gamma))
        ) + (101325 - Pa) / P0 * eps

        At = F / (Cf * P0)
        Rt = math.sqrt(At / math.pi)
        Re = Rt * math.sqrt(eps)
        Ae = math.pi * Re ** 2

    # ============================================================
    #  CHAMBER DEFINITION (NOW Rt IS GUARANTEED TO EXIST)
    # ============================================================

    if use_cea:
        # already defined by CEA
        pass
    else:
        # Rao-only geometric chamber (visual + smoothness)
        Rc = 2.5 * Rt  # CR ≈ 6.25
        Lc = 2.0 * Rt
        Ac = math.pi * Rc ** 2

    # ===== Rao tables =====
    ar = np.array([3.5,4,5,10,20,30,40,50,100])

    theta_n = {
        0.6:[25.5,26.5,28,32,35,36.3,37.5,38,40.25],
        0.7:[23,24,25,28.5,31.3,32.6,33.6,34.5,36.5],
        0.8:[21,21.6,23,26.3,27.9,30,31,32.2,33.75],
        0.9:[19,20,21,24.25,27,28.5,29.3,30,32.5],
        1.0:[18.5,19,20,22.5,25.5,27,28,29,32]
    }

    theta_e = {
        0.6:[21.8,21,19,16,14.5,14,13.5,13,12],
        0.7:[18,17,16,13,12,11.2,10.8,10.5,9.75],
        0.8:[14.8,14,13,10.5,9,8.25,8,7.75,7],
        0.9:[12,11,10,8,7,6.5,6.25,6,6],
        1.0:[10,9,8,6,5.25,4.9,4.75,4.5,4.25]
    }

    # ===== Length (CRITICAL) =====
    f1 = ((math.sqrt(eps)-1)*Rt)/math.tan(math.radians(15))
    LN = bell_percent * f1

    # ===== Log-log interpolation =====
    thn = np.deg2rad(np.interp(np.log10(eps), np.log10(ar), theta_n[bell_percent]))
    the = np.deg2rad(np.interp(np.log10(eps), np.log10(ar), theta_e[bell_percent]))

    # ===== Throat arc =====
    th2 = np.linspace(-math.pi/2, thn-math.pi/2, 50)
    x2 = 0.382*Rt*np.cos(th2)
    y2 = 0.382*Rt*np.sin(th2) + 0.382*Rt + Rt

    Nx, Ny = x2[-1], y2[-1]

    # ===== Bezier =====
    m1, m2 = math.tan(thn), math.tan(the)
    C1 = Ny - m1*Nx
    C2 = Re - m2*LN

    Qx = (C2-C1)/(m1-m2)
    Qy = (m1*C2 - m2*C1)/(m1-m2)

    t = np.linspace(0,1,300)
    xb = (1-t)**2*Nx + 2*(1-t)*t*Qx + t**2*LN
    yb = (1-t)**2*Ny + 2*(1-t)*t*Qy + t**2*Re

    # ===== Convergent =====
    # ============================================================
    #  TRUE LINEAR CONVERGENT (CHAMBER → THROAT)
    # ============================================================
    # ============================================================
    #  CONVERGENT SELECTION LOGIC
    # ============================================================

    if use_cea and include_chamber:
        # --------------------------------------------------------
        # S-CURVE CONVERGENT (CEA + chamber ONLY)
        # --------------------------------------------------------
        theta_throat = math.radians(25.0)
        L_conv = (Rc - Rt) / math.tan(theta_throat)

        x_conv = np.linspace(-L_conv, 0.0, 200)
        s = (x_conv + L_conv) / L_conv

        r_conv = (
                Rc
                + (Rt - Rc) * (
                        10 * s ** 3
                        - 15 * s ** 4
                        + 6 * s ** 5
                )
        )


        # ============================================================
        #  CORRECT L* CHAMBER LENGTH (CEA + chamber ONLY)
        # ============================================================

        V_required = Lstar_m * At
        V_conv = np.trapezoid(np.pi * r_conv**2, x_conv)
        V_cyl = V_required - V_conv

        if V_cyl <= 0:
            raise ValueError("L* too small for selected convergent geometry")

        Lc = V_cyl / Ac





        x_conv_lin = x_conv
        r_conv_lin = r_conv

    else:
        # --------------------------------------------------------
        # ORIGINAL RAO CONVERGENT (MATCHES SHOCKSVE.py)
        # --------------------------------------------------------
        if use_cea:
            Rc_geom = Rc
        else:
            Rc_geom = 1.5 * Rt

        tc = np.linspace(-3 * math.pi / 4, -math.pi / 2, 80)
        x_conv_lin = Rc_geom * np.cos(tc)
        r_conv_lin = Rc_geom * np.sin(tc) + Rc_geom + Rt

    # ===== FORCE connection to chamber radius (CEA + chamber ONLY) =====
    if use_cea and include_chamber and Lc is not None and L_conv is not None:


        Rc_use = Rc
        Lc_use = Lc

        # ============================================================
        #  SMOOTH BLEND: CHAMBER → LINEAR CONVERGENT (C¹)
        # ============================================================

        L_blend = 0.25 * (Rc_use - Rt)

        x_blend = np.linspace(-L_conv - L_blend, -L_conv, 80)
        t = (x_blend - x_blend[0]) / L_blend

        drdx_end = -(Rc_use - Rt) / L_conv

        r_blend = (
                Rc_use * (2 * t ** 3 - 3 * t ** 2 + 1)
                + r_conv_lin[0] * (-2 * t ** 3 + 3 * t ** 2)
                + drdx_end * L_blend * (t ** 3 - t ** 2)
        )

        # Chamber cylinder
        x_ch = np.linspace(x_blend[0] - Lc_use, x_blend[0], 200)
        r_ch = np.full_like(x_ch, Rc_use)

    else:
        x_ch = np.array([])
        r_ch = np.array([])

    # ===== Combine + resample =====
    x = np.concatenate([
        x_ch,
        x_blend,
        x_conv_lin,
        x2,
        xb
    ])

    r = np.concatenate([
        r_ch,
        r_blend,
        r_conv_lin,
        y2,
        yb
    ])

    # ===== CLEANUP: enforce unique throat radius (prevents 3D spike) =====
    # Find throat location (x ≈ 0)


    # --- MATLAB: unique(x,'stable') ---
    x_unique, idx = np.unique(x, return_index=True)
    r_unique = r[idx]

    # --- Ensure strictly increasing ---
    order = np.argsort(x_unique)
    x_unique = x_unique[order]
    r_unique = r_unique[order]

    # --- PCHIP resampling ---
    # ===== SPLINE ONLY THE NOZZLE (NOT THE CHAMBER) =====
    if include_chamber and len(x_ch) > 0:
        x_nozzle = np.concatenate([x_conv_lin, x2, xb])
        r_nozzle = np.concatenate([r_conv_lin, y2, yb])

        # unique + sort nozzle only
        xu, idx = np.unique(x_nozzle, return_index=True)
        ru = r_nozzle[idx]
        order = np.argsort(xu)
        xu = xu[order]
        ru = ru[order]

        # resample nozzle only
        xi_n = np.linspace(xu[0], xu[-1], N)
        ri_n = PchipInterpolator(xu, ru)(xi_n)

        # FINAL contour = chamber + nozzle
        xi = np.concatenate([x_ch, xi_n])
        ri = np.concatenate([r_ch, ri_n])

    else:
        # original behavior (no chamber)
        xi = np.linspace(x_unique[0], x_unique[-1], N)
        ri = PchipInterpolator(x_unique, r_unique)(xi)

    # ===== Performance (NO calibrated pressure) =====
    Ve = math.sqrt(
        (2 * gamma * R * T0) / (gamma - 1.0)
        * (1.0 - (Pe / P0) ** ((gamma - 1.0) / gamma))
    )

    phi = math.sqrt(gamma / R) * (
            ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    mdot = P0 / math.sqrt(T0) * phi * At

    F_actual = mdot * Ve + (Pe - Pa) * Ae


    result = {
       # full contour (resampled)
       "x": xi * 1000,
       "r": ri * 1000,

    # throat arc (for continuity)
    "x_throat": x2 * 1000,
    "r_throat": y2 * 1000,

    # convergent (original resolution)
    "x_conv": x_conv_lin * 1000,
    "r_conv": r_conv_lin * 1000,


    # divergent bell (Bezier only)
    "x_div": xb * 1000,
    "r_div": yb * 1000,

    # length
    "LN_m": LN,

    # performance
    "Ve": Ve,
    "mdot": mdot,
    "F_actual": F_actual,
    "At": At,
    "Ae": Ae,
    "Rt_mm": Rt * 1000,
    "Re_mm": Re * 1000,
    "eps": eps,
    "Me": Me,
}

# ===== ADD CEA CHAMBER DATA ONLY IF USED =====
    # ===== ADD CHAMBER DATA IF INCLUDED =====
    if include_chamber:
        result.update({
            "x_chamber": x_ch * 1000,
            "r_chamber": r_ch * 1000,
        })

    # ===== ADD CEA META DATA (OPTIONAL) =====
    if use_cea:
        result.update({
            "CR": CR,
            "Lstar_m": Lstar_m,
            "Rc_m": Rc,
            "Lc_m": Lc,
            "Ac": Ac,
        })

    return result


def shockwve_compute(b, P0, T0, gamma, Pa, Pe, R, F):
    """
    Core physics + geometry calculation, translated from your MATLAB script.
    Returns a dictionary with all main values and arrays.
    Units:
      - Lengths in mm for x, r
      - Areas in m^2
    """
    # === CFD-calibrated exit pressure ===
    Pe_user = Pe
    Pe_cal = Pe_user - 25543.0

    # Safety clamp (important)
    Pe_cal = max(Pe_cal, 1.0)

    # ========== Ideal nozzle performance ==========
    Ve = math.sqrt(
        (2 * gamma * R * T0) / (gamma - 1.0)
        * (1.0 - (Pa / P0) ** ((gamma - 1.0) / gamma))
    )

    Me_geom = math.sqrt(
        (2.0 / (gamma - 1.0)) * ((P0 / (Pe_cal)) ** ((gamma - 1.0) / gamma) - 1.0)
    )

    Ae_At = (1.0 / Me_geom) * (
        (2.0 / (gamma + 1.0) * (1.0 + (gamma - 1.0) / 2.0 * Me_geom**2))
        ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    Me2 = math.sqrt(
        (2.0 / (gamma - 1.0)) * ((P0 / (Pe)) ** ((gamma - 1.0) / gamma) - 1.0)
    )

    Te = T0 / (1.0 + (gamma - 1.0) / 2.0 *Me2**2)

    phi = math.sqrt(gamma / R) * (
        ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    mdot_ideal = F / Ve
    At = mdot_ideal / (P0 / math.sqrt(T0) * phi)
    Ae = Ae_At * At

    At_mm2 = At * 1.0e6
    Ae_mm2 = Ae * 1.0e6
    R_throat_mm = math.sqrt(At_mm2 / math.pi)
    R_exit_mm = math.sqrt(Ae_mm2 / math.pi)

    dN = 2.0 * R_throat_mm

    # ========== Chamber / convergent ==========
    CR = 2.5
    Ac = CR * At
    Rc_m = math.sqrt(Ac / math.pi)
    dc_mm = 2.0 * Rc_m * 1e3

    theta_c = math.radians(20.0)
    Rt_m = math.sqrt(At / math.pi)
    Rc2_m = math.sqrt(Ac / math.pi)
    lN_m = (Rc2_m - Rt_m) / math.tan(theta_c)
    lN_mm = lN_m * 1e3

    # Convergent contour
    x_conv = np.linspace(0.0, lN_mm, 1000)

    d_conv = dN + (dc_mm - dN) * (1.0 - (x_conv / lN_mm)) ** b
    r_conv = d_conv / 2.0
    x_conv = x_conv - lN_mm  # throat at x=0
    R_throat_mm = r_conv[-1]

    # ========== Throat → Point N circular blend ==========

    Apm = math.sqrt((gamma + 1.0) / (gamma - 1.0))
    Bpm = (gamma - 1.0) / (gamma + 1.0)

    def V_PM(Me_geom):
        return Apm * math.atan(math.sqrt(Bpm * (Me_geom**2 - 1.0))) - math.atan(
            math.sqrt(Me_geom**2 - 1.0)
        )

    Tmax = 0.5 * V_PM(Me_geom)  # radians

    Rt = R_throat_mm
    Nt = 0.382 * Rt * math.sin(Tmax)
    Na = Rt + 0.382 * Rt * (1.0 - math.cos(Tmax))
    x_N = Nt
    r_N = Na

    # Parabolic connector
    x_throat = 0.0
    r_throat = R_throat_mm
    slope_start = 0.0

    A_mat = np.array(
        [
            [0.0**2, 0.0, 1.0],
            [x_N**2, x_N, 1.0],
            [2.0 * 0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    rhs_par = np.array([r_throat, r_N, slope_start], dtype=float)
    a_par, b_par, c_par = np.linalg.solve(A_mat, rhs_par)

    x_conn = np.linspace(x_throat, x_N, 10)

    r_conn = a_par * x_conn**2 + b_par * x_conn + c_par

    # Circular Rao blend
    Rb = 1.5 * Rt
    xc = -Rb
    rc = Rt
    phi_T = 0.0
    phi_N = math.atan2(r_N - rc, x_N - xc)

    phi_arc_full = np.linspace(phi_T, phi_N, 200)
    x_arc_full = xc + Rb * np.cos(phi_arc_full)
    r_arc_full = rc + Rb * np.sin(phi_arc_full)

    keep = r_arc_full < Rt
    x_arc = x_arc_full[keep]
    r_arc = r_arc_full[keep]

    # ========== MOC divergent section ==========
    TR_mm = r_N
    xw_moc, yw_moc, C_plus, C_minus = moc_bell_nozzle_original(
        P0, gamma, TR_mm, Me_geom
    )

    x_div = xw_moc + x_N
    y_div = yw_moc.copy()
    x_div[0] = x_N
    y_div[0] = r_N

    # ========== Optional scaling of divergent section to hit target Pe ==========
    def exit_pressure_for_scale(s):
        # Scale y_div about r_N
        y_scaled = r_N + s * (y_div - r_N)
        # Build temporary full radius array (mm)
        r_full_tmp = np.concatenate([r_conv, r_conn, r_arc, y_scaled])
        Rt_tmp_mm = float(np.min(r_full_tmp))
        At_tmp = math.pi * (Rt_tmp_mm * 1e-3)**2  # m^2
        r_exit_tmp_m = float(y_scaled[-1]) * 1e-3
        Ae_tmp = math.pi * r_exit_tmp_m**2
        AR_e = Ae_tmp / At_tmp

        def fAM(M):
            if M <= 1.0:
                return 1e6
            return (1.0 / M) * (
                (2.0 / (gamma + 1.0) * (1.0 + (gamma - 1.0) / 2.0 * M**2))
                ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
            ) - AR_e

        try:
            M_e = bisection(fAM, 1.01, 50.0, tol=1e-6, max_iter=100)
        except RuntimeError:
            M_e = 5.0

        pqP0 = (1.0 + (gamma - 1.0) / 2.0 * M_e**2) ** (-gamma / (gamma - 1.0))
        return pqP0 * P0

    if USE_SCALING:
        s_low, s_high = 0.1, 2.5
        f_low = exit_pressure_for_scale(s_low) - Pe
        f_high = exit_pressure_for_scale(s_high) - Pe
        expand_count = 0
        while f_low * f_high > 0 and expand_count < 20:
            s_low *= 0.5
            s_high *= 1.5
            f_low = exit_pressure_for_scale(s_low) - Pe
            f_high = exit_pressure_for_scale(s_high) - Pe
            expand_count += 1
        if f_low * f_high > 0:
            s_opt = 1.0
        else:
            def f_root(s):
                return exit_pressure_for_scale(s) - Pe_cal
            s_opt = bisection(f_root, s_low, s_high, tol=1e-4, max_iter=60)
    else:
        s_opt = 1.0

    # Apply final scaling about point N
    y_div = r_N + s_opt * (y_div - r_N)

    # ========== Use MOC exit radius for isentropic sizing ==========
    r_exit_m = float(y_div[-1]) * 1e-3
    Me_exit = math.sqrt(
        (2.0 / (gamma - 1.0)) * ((P0 / Pe_cal) ** ((gamma - 1.0) / gamma) - 1.0)
    )

    Ae_moc = math.pi * r_exit_m**2

    # Recompute throat area from geometry
    r_full_tmp_final = np.concatenate([r_conv, r_conn, r_arc, y_div])
    Rt_geom_mm = float(np.min(r_full_tmp_final))
    At_moc = math.pi * (Rt_geom_mm * 1e-3)**2
    Rt_m_moc = Rt_geom_mm * 1e-3

    Ae_At_moc = (1.0 / Me_exit) * (
            (2.0 / (gamma + 1.0) * (1.0 + (gamma - 1.0) / 2.0 * Me_exit ** 2))
            ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    At_moc = Ae_moc / Ae_At_moc

    Ve_moc = math.sqrt(
        (2 * gamma * R * T0) / (gamma - 1.0)
        * (1.0 - (Pe / P0) ** ((gamma - 1.0) / gamma))
    )
    mdot_moc = (F - (Pe - Pa) * Ae_moc) / Ve_moc

    # Back-computed exit pressure for final contour
    def exit_pressure_final():
        AR_e = Ae_moc / At_moc

        def fAM2(M):
            if M <= 1.0:
                return 1e6
            return (1.0 / M) * (
                (2.0 / (gamma + 1.0) * (1.0 + (gamma - 1.0) / 2.0 * M**2))
                ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
            ) - AR_e

        try:
            M_e = bisection(fAM2, 1.01, 50.0, tol=1e-6, max_iter=100)
        except RuntimeError:
            M_e = 5.0

        pqP0 = (1.0 + (gamma - 1.0) / 2.0 * M_e**2) ** (-gamma / (gamma - 1.0))
        return pqP0 * P0

    p_exit = Pe

    # ============================================================
    #  Actual thrust for user-chosen over/under-expanded Pe
    # ============================================================

    # Fixed mass flow from choked flow at the throat
    mdot_fixed = P0 / math.sqrt(T0) * phi * At

    # Actual thrust from isentropic thrust equation
    F_actual = mdot_fixed * Ve + (Pe - Pa) * Ae

    # ========== Combine full 2D contour ==========
    x_full = np.concatenate([x_conv, x_conn, x_arc, x_div])
    r_full = np.concatenate([r_conv, r_conn, r_arc, y_div])

    # === GLOBAL C² SPLINE SMOOTHING (MATLAB-EQUIVALENT) ===
    x_raw = x_full.copy()
    r_raw = r_full.copy()

    xu, ia = np.unique(x_raw, return_index=True)
    ru = r_raw[ia]

    Nsmooth = 2000
    x_full = np.linspace(xu[0], xu[-1], Nsmooth)

    from scipy.interpolate import CubicSpline
    cs = CubicSpline(xu, ru)
    r_full = cs(x_full)

    r_full[r_full < 0.0] = 0.0

    # 3D surface revolve
    r_safe = np.maximum(r_full, 1e-6)
    theta = np.linspace(0.0, 2.0 * math.pi, 100)
    X, TH = np.meshgrid(x_full, theta)
    Rg = np.tile(r_safe, (theta.size, 1))
    Y = Rg * np.cos(TH)
    Z = Rg * np.sin(TH)

    result = {
        # scalar performance
        "F_actual": F_actual,
        "mdot_fixed": mdot_fixed,
        "Ve": Ve,
        "Me": Me_geom,
        "Ae_At": Ae_At,
        "Te": Te,
        "phi": phi,
        "mdot_ideal": mdot_ideal,
        "At": At,
        "Ae": Ae,
        "R_throat_mm": R_throat_mm,
        "R_exit_mm": R_exit_mm,
        "CR": CR,
        "Ac": Ac,
        "Rc_m": Rc_m,
        "dc_mm": dc_mm,
        "lN_mm": lN_mm,
        # MOC results
        "r_exit_m_moc": r_exit_m,
        "Ae_moc": Ae_moc,
        "Ae_At_moc": Ae_At_moc,
        "C_plus": C_plus,
        "C_minus": C_minus,
        "At_moc": At_moc,
        "Rt_m_moc": Rt_m_moc,
        "Ve_moc": Ve_moc,
        "mdot_moc": mdot_moc,
        "s_opt": s_opt,
        "p_exit": p_exit,
        # geometry arrays
        "x_conv": x_conv,
        "r_conv": r_conv,
        "x_conn": x_conn,
        "r_conn": r_conn,
        "x_arc": x_arc,
        "r_arc": r_arc,
        "x_div": x_div,
        "y_div": y_div,
        "x_full": x_full,
        "r_full": r_full,
        # 3D mesh
        "X3D": X,
        "Y3D": Y,
        "Z3D": Z,
    }

    return result


# ============================================================
#  Matplotlib canvas wrapper
# ============================================================

class MplCanvas(FigureCanvas):
    def __init__(self, is_3d=False):
        fig = Figure(facecolor="black")
        if is_3d:
            ax = fig.add_subplot(111, projection="3d")
        else:
            ax = fig.add_subplot(111)

        fig.patch.set_facecolor("black")
        ax.set_facecolor("black")
        ax.tick_params(colors="white")

        # 2D axes have spines, 3D doesn't; guard it
        for spine in getattr(ax, "spines", {}).values():
            spine.set_color("white")

        super().__init__(fig)
        self.fig = fig
        self.ax = ax
        self.setStyleSheet("background-color: #000000;")


# ============================================================
#  PyQt6 GUI Application with QStackedWidget
# ============================================================


class GlassPanel(QWidget):
    def __init__(self):
        super().__init__()
        self.setStyleSheet("""
            QWidget {
                background-color: rgba(10,10,10,230);
                border: 1px solid #00f5ff;
                border-radius: 14px;
            }
        """)




class ShockWveApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Colors
        self.bg_color = "#000000"
        self.cyan = "#00f5ff"
        self.magenta = "#ff00ff"
        self.text = "#ffffff"
        self.accent = "#ff8c00"

        self.setWindowTitle("SHOCKWVE - Nozzle Designer")
        # self.setGeometry(100, 100, 1400, 800)
        #self.setGeometry(300, 200, 1100, 550)
        self.resize(300, 500)
        self.center_on_screen()

        # === Central widget is a QStackedWidget ===
        self.stack = QStackedWidget()
        self.setCentralWidget(self.stack)

        # Build both pages
        self.build_main_page()
        self.build_splash_page()

        # Add pages to stack
        self.stack.addWidget(self.splash_page)
        self.stack.addWidget(self.main_page)

        # Show splash first
        self.stack.setCurrentWidget(self.splash_page)

        self.apply_dark_style()

        self.refresh_engine_list()
    # -------------------------------------------
    # Build MAIN UI page (no setCentralWidget here)
    # -------------------------------------------
    def build_main_page(self):
        self.main_page = QWidget()
        root = QHBoxLayout(self.main_page)
        root.setContentsMargins(14, 14, 14, 14)
        root.setSpacing(14)

        # ==========================
        # LEFT CONTROL HUD
        # ==========================
        left_panel = GlassPanel()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setSpacing(12)

        title = QLabel("SHOCKWVE     CONTROL CORE")
        title.setStyleSheet("""
            QLabel {
                color: #00f5ff;
                font-size: 18px;
                font-weight: bold;
                font-family: Consolas;
            }
        """)
        left_layout.addWidget(title)

        # ==========================
        # ENGINE PROFILE CONTROLS
        # ==========================

        self.engine_selector = QComboBox()
        self.engine_selector.setStyleSheet("""
            QComboBox {
                background-color: black;
                color: #00f5ff;
                border: 1px solid #00f5ff;
                padding: 4px;
                font-family: Consolas;
            }
        """)
        self.engine_selector.currentIndexChanged.connect(self.load_selected_engine)
        left_layout.addWidget(self.engine_selector)

        save_engine_btn = QPushButton("SAVE ENGINE")
        save_engine_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #00f5ff;
                border: 1px solid #00f5ff;
                padding: 6px;
                font-family: Consolas;
                border-radius: 6px;
            }
            QPushButton:hover {
                background-color: #ff40c0;
                color: black;
                border: 1px solid #ff40c0;
            }
        """)
        save_engine_btn.clicked.connect(self.save_current_engine)
        left_layout.addWidget(save_engine_btn)

        delete_engine_btn = QPushButton("DELETE ENGINE")
        delete_engine_btn.setStyleSheet("""
            QPushButton {
                background-color: #111111;
                color: #00f5ff;
                border: 1px solid #00f5ff;
                padding: 6px;
                font-family: Consolas;
                border-radius: 6px;
            }
            QPushButton:hover {
                background-color: #ff40c0;
                color: black;
                border: 1px solid #ff40c0;
            }
        """)
        delete_engine_btn.clicked.connect(self.delete_selected_engine)
        left_layout.addWidget(delete_engine_btn)



        form = QFormLayout()
        form.setSpacing(6)

        def neon_edit(default):
            w = QLineEdit(default)
            w.setStyleSheet("""
                QLineEdit {
                    background-color: #000000;
                    color: #ff40c0;
                    border: 1px solid #ff40c0;
                    padding: 4px;
                    font-family: Consolas;
                }

                /* 🔥 HOVER = CYAN */
                QLineEdit:hover {
                    border: 1px solid #00f5ff;
                }

                /* 🔥 CLICK / FOCUS = CYAN (keep this) */
                QLineEdit:focus {
                    border: 1px solid #00f5ff;
                    color: #00f5ff;
                }
            """)
            return w

        self.b_edit = neon_edit("1.5")
        self.P0_edit = neon_edit("3000000")
        self.T0_edit = neon_edit("3600")
        self.gamma_edit = neon_edit("1.4")
        self.Pa_edit = neon_edit("101325")
        self.Pe_edit = neon_edit("101325")
        self.R_edit = neon_edit("287")
        self.F_edit = neon_edit("1500")

        for lab, w in [
            ("Shape exponent b", self.b_edit),
            ("P0 [Pa]", self.P0_edit),
            ("T0 [K]", self.T0_edit),
            ("Gamma", self.gamma_edit),
            ("Pa [Pa]", self.Pa_edit),
            ("Pe [Pa]", self.Pe_edit),
            ("R [J/Kg*K]", self.R_edit),
            ("F [N]", self.F_edit),
        ]:
            lbl = QLabel(lab)
            lbl.setStyleSheet("color:white; font-family:Consolas;")
            form.addRow(lbl, w)

            # ==========================
            # CEA-only inputs (hidden by default)
            # ==========================
        # ==========================
        # CEA-only inputs (hidden by default)
        # ==========================

        self.Cf_edit = neon_edit("1.5")
        self.eps_edit = neon_edit("10")
        self.Lstar_edit = neon_edit("1.0")
        self.CR_edit = neon_edit("2.5")

        self.cea_rows = []

        for label_text, widget in [
            ("Cf (sea-level)", self.Cf_edit),
            ("Expansion Ratio ε", self.eps_edit),
            ("L* [m]", self.Lstar_edit),
            ("Contraction Ratio", self.CR_edit),
        ]:
            lbl = QLabel(label_text)
            lbl.setStyleSheet("color:#ffaa00; font-family:Consolas;")
            lbl.setVisible(False)
            widget.setVisible(False)

            form.addRow(lbl, widget)
            self.cea_rows.append((lbl, widget))

        left_layout.addLayout(form)

        # RUN BUTTON
        self.run_button = QPushButton("▶ RUN SHOCKWVE")
        self.run_button.setStyleSheet("""
            QPushButton {
                background-color: #ff40c0;
                color: black;
                font-size: 14px;
                font-weight: bold;
                padding: 10px;
                border-radius: 10px;
            }
            QPushButton:hover {
                background-color: #00f5ff;
            }
        """)
        self.run_button.clicked.connect(self.run_computation)
        left_layout.addWidget(self.run_button)

        self.show_char_cb = QCheckBox("Show Characteristics")
        self.show_char_cb.setChecked(False)
        self.show_char_cb.setStyleSheet("color:#00f5ff;")
        self.show_char_cb.stateChanged.connect(self.update_2d_plot)
        left_layout.addWidget(self.show_char_cb)

        self.show_legend_cb = QCheckBox("Show Legend")
        self.show_legend_cb.setChecked(True)  # legend ON by default
        self.show_legend_cb.setStyleSheet("color:#00f5ff;")
        self.show_legend_cb.stateChanged.connect(self.update_2d_plot)
        left_layout.addWidget(self.show_legend_cb)

        self.nozzle_mode_cb = QCheckBox("Use Rao Bell Nozzle")
        self.nozzle_mode_cb.setChecked(False)  # default = MOC
        self.nozzle_mode_cb.setStyleSheet("color:#00ff66;")
        self.nozzle_mode_cb.stateChanged.connect(self.on_nozzle_mode_changed)

        left_layout.addWidget(self.nozzle_mode_cb)

        self.cea_mode_cb = QCheckBox("Implement CEA")
        self.cea_mode_cb.setChecked(False)
        self.cea_mode_cb.setStyleSheet("color:#ffaa00;")
        self.cea_mode_cb.setEnabled(False)  # disabled unless Rao is ON
        self.cea_mode_cb.stateChanged.connect(self.update_cea_ui)

        left_layout.addWidget(self.cea_mode_cb)

        self.show_chamber_cb = QCheckBox("Include Chamber Length (Rao)")
        self.show_chamber_cb.setChecked(False)
        self.show_chamber_cb.setStyleSheet("color:#ffaa00;")
        self.show_chamber_cb.setEnabled(False)
        self.show_chamber_cb.stateChanged.connect(self.run_computation)

        left_layout.addWidget(self.show_chamber_cb)

        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setStyleSheet("""
            QTextEdit {
                background-color: #000;
                border: 1px solid #00f5ff;
                color: white;
                font-family: Consolas;
            }
        """)
        left_layout.addWidget(self.output_text, stretch=1)

        export_button = QPushButton("EXPORT 2D CONTOUR (.xlsx)")
        export_button.setStyleSheet("""
            QPushButton {
                background-color: #00f5ff;
                color: black;
                font-weight: bold;
                padding: 6px;
                border-radius: 8px;
            }
            QPushButton:hover {
                background-color: #ff40c0;
            }
        """)
        export_button.clicked.connect(self.export_contour_csv)
        left_layout.addWidget(export_button)

        root.addWidget(left_panel, 1)

    # ==========================
    # RIGHT VISUAL DECK
    # ==========================
        right_panel = GlassPanel()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setSpacing(8)

    # ==========================
    # TOP RIGHT DATA HUD
    # ==========================
        data_bar = QHBoxLayout()
        data_bar.setAlignment(Qt.AlignmentFlag.AlignRight)

        self.point_selector = QComboBox()
        self.point_selector.addItems(["MOC Points", "Rao Bell Points"])
        self.point_selector.setStyleSheet("""
        QComboBox {
            background-color: black;
            color: #00f5ff;
            border: 1px solid #00f5ff;
            padding: 4px;
            font-family: Consolas;
        }
        QComboBox QAbstractItemView {
            background-color: black;
            color: #00f5ff;
            selection-background-color: #ff40c0;
        }
    """)
        self.point_selector.currentIndexChanged.connect(self.update_point_table)

        data_bar.addWidget(QLabel("DATA:"))
        data_bar.addWidget(self.point_selector)
        right_layout.addLayout(data_bar)

    # ==========================
    # EXCEL-STYLE TABLE
    # ==========================
        self.point_table = QTableWidget()
        self.point_table.setColumnCount(2)
        self.point_table.setHorizontalHeaderLabels(["Axial (mm)", "Radius (mm)"])
        self.point_table.setMaximumHeight(220)
        self.point_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.point_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.point_table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self.point_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.point_table.customContextMenuRequested.connect(self.show_table_context_menu)

        right_layout.addWidget(self.point_table)

    # ==========================
    # TABS (THIS WAS MISSING)
    # ==========================
        self.tabs = QTabWidget()  # 🔥 REQUIRED LINE

        self.point_table.setStyleSheet("""
            QTableWidget {
                background-color: black;
                color: white;
                gridline-color: #00f5ff;
                font-family: Consolas;
                selection-background-color: black;   /* 🔥 black highlight */
                selection-color: #ff40c0;            /* 🔥 neon pink text */
            }

            QTableWidget::item:selected {
                background-color: black;             /* 🔥 black */
                color: #ff40c0;                      /* 🔥 neon pink */
            }

            QHeaderView::section {
                background-color: #111;
                color: #00f5ff;
                padding: 4px;
                border: 1px solid #00f5ff;
            }
        """)

        self.tabs.setStyleSheet("""
        QTabBar::tab {
            background:#111;
            color:white;
            padding:6px 16px;
            border-radius:4px;
        }
        QTabWidget::pane {
            border-top: 8px solid #111;
        }
        QTabBar::tab:selected {
            background:#00f5ff;
            color:black;
        }
    """)

        self.canvas_2d = MplCanvas(is_3d=False)
        self.canvas_3d = MplCanvas(is_3d=True)

        t1 = QWidget()
        QVBoxLayout(t1).addWidget(self.canvas_2d)

        t2 = QWidget()
        QVBoxLayout(t2).addWidget(self.canvas_3d)

        self.tabs.addTab(t1, "2D CONTOUR")
        self.tabs.addTab(t2, "3D WIREFRAME")

        right_layout.addWidget(self.tabs)
        root.addWidget(right_panel, 2)

    def update_cea_ui(self):
        show = self.nozzle_mode_cb.isChecked() and self.cea_mode_cb.isChecked()
        for lbl, widget in self.cea_rows:
            lbl.setVisible(show)
            widget.setVisible(show)

    def copy_selected_column(self, col):
        selected = self.point_table.selectionModel().selectedRows()
        if not selected:
            return

        lines = [
            self.point_table.item(idx.row(), col).text()
            for idx in selected
        ]

        QApplication.clipboard().setText("\n".join(lines))

    def copy_selected_rows(self):
        selected = self.point_table.selectionModel().selectedRows()
        if not selected:
            return

        lines = []
        for idx in selected:
            x = self.point_table.item(idx.row(), 0).text()
            r = self.point_table.item(idx.row(), 1).text()
            lines.append(f"{x}\t{r}")  # tab = Excel friendly

        QApplication.clipboard().setText("\n".join(lines))

    def show_table_context_menu(self, pos):
        menu = QMessageBox()  # dummy init to grab palette
        from PyQt6.QtWidgets import QMenu
        from PyQt6.QtGui import QAction

        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #050505;
                color: #00f5ff;
                border: 1px solid #00f5ff;
                font-family: Consolas;
            }
            QMenu::item:selected {
                background-color: #00f5ff;
                color: black;
            }
        """)

        copy_all = QAction("Copy Selected Rows", self)
        copy_x = QAction("Copy Axial (mm) Only", self)
        copy_r = QAction("Copy Radius (mm) Only", self)

        copy_all.triggered.connect(self.copy_selected_rows)
        copy_x.triggered.connect(lambda: self.copy_selected_column(0))
        copy_r.triggered.connect(lambda: self.copy_selected_column(1))

        menu.addAction(copy_all)
        menu.addSeparator()
        menu.addAction(copy_x)
        menu.addAction(copy_r)

        menu.exec(self.point_table.viewport().mapToGlobal(pos))

    def update_point_table(self):
        if not hasattr(self, "last_result"):
            return

        use_rao = self.point_selector.currentText().startswith("Rao")

        if use_rao:
            x = self.last_result.get("x", [])
            r = self.last_result.get("r", [])
        else:
            x = self.last_result.get("x_full", [])
            r = self.last_result.get("r_full", [])

        self.populate_table(x, r)

    def populate_table(self, x, r):
        self.point_table.setRowCount(len(x))

        for i, (xi, ri) in enumerate(zip(x, r)):
            self.point_table.setItem(i, 0, QTableWidgetItem(f"{xi:.8f}"))
            self.point_table.setItem(i, 1, QTableWidgetItem(f"{ri:.8f}"))

        self.point_table.resizeColumnsToContents()

    # -------------------------------------------
    # Splash screen page: full ASCII logo + slogan + BRIDGE
    # -------------------------------------------
    def build_splash_page(self):
        self.splash_page = QWidget()
        layout = QVBoxLayout(self.splash_page)
        layout.setContentsMargins(10, 10, 10, 10)

        # ===== Load ASCII from logo.txt =====
        try:
            with open("logo.txt", "r", encoding="utf-8") as f:
                raw_logo = f.read()
        except Exception:
            raw_logo = "SHOCKWVE"

        # ===== Apply neon cyan + neon pink coloring =====
        lines = raw_logo.split("\n")

        colored_lines = []

        for i, line in enumerate(lines):

            if 0 <= i <= 20:
                line_html = ""
                for ch in line:
                    if ch in ("Q", "S"):
                        line_html += f"<span style='color:#00eaff'>{ch}</span>"
                    elif ch in ("O", "L", "N", "J"):
                        line_html += f"<span style='color:#ff40c0'>{ch}</span>"
                    else:
                        line_html += ch
                colored_lines.append(line_html)

            elif 21 <= i <= 35:
                line_html = ""
                for ch in line:
                    if ch.strip() == "":
                        line_html += ch
                    elif ch in ("L", "N", "O"):
                        line_html += f"<span style='color:#ff40c0'>{ch}</span>"
                    else:
                        line_html += f"<span style='color:#00eaff'>{ch}</span>"
                colored_lines.append(line_html)

            elif 36 <= i <= 45:
                line_html = ""
                for ch in line:
                    if ch in ("O", "Q"):
                        line_html += f"<span style='color:#ff40c0'>{ch}</span>"
                    else:
                        line_html += f"<span style='color:#00eaff'>{ch}</span>"
                colored_lines.append(line_html)

            else:
                if "玉" in line or "兔" in line:
                    line_html = line.replace("玉", "<span style='color:#00eaff'>玉</span>") \
                        .replace("兔", "<span style='color:#00eaff'>兔</span>")
                    colored_lines.append(line_html)
                else:
                    line_html = ""
                    for ch in line:
                        if ch.strip():
                            line_html += f"<span style='color:#00eaff'>{ch}</span>"
                        else:
                            line_html += ch
                    colored_lines.append(line_html)

        colored_logo = "<pre style='font-family:Consolas; font-size:12px;'>" + \
                       "\n".join(colored_lines) + "</pre>"
        from PyQt6.QtGui import QGuiApplication



        # ===== Scroll area so logo can be huge =====
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setStyleSheet("background-color:black; border:none;")

        logo_container = QWidget()
        logo_layout = QVBoxLayout(logo_container)
        logo_layout.setContentsMargins(0, 0, 0, 0)

        self.logo_label = QLabel("")
        self.logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.logo_label.setStyleSheet("background-color:black; color:#00eaff; font-family:Consolas; font-size:13px;")

        logo_layout.addWidget(self.logo_label)

        self.splash_lines = colored_lines  # ← your processed ASCII lines
        self.line_index = 0
        self.char_index = 0
        self.typed_text = ""

        self.typing_timer = QTimer()
        self.typing_timer.timeout.connect(self.type_next_character)
        self.typing_timer.start(1)  # 12 ms per character – cyber speed

        scroll.setWidget(logo_container)

        # ===== Credits =====
        credits = QLabel(
            "<span style='color:#666666; font-size:14px;'>© 2025 Jason Da Silva — All Rights Reserved</span>"
        )
        credits.setAlignment(Qt.AlignmentFlag.AlignLeft)

        # ===== Slogan =====
        slogan = QLabel(
            "<span style='color:#00eaff; font-size:26px;'>THE FUTURE IS CHROME</span>"
        )
        slogan.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # ===== BRIDGE button =====
        bridge_btn = QPushButton("BRIDGE")
        bridge_btn.setStyleSheet(
            """
            QPushButton {
                background-color: #ff40c0;
                color: black;
                font-size: 22px;
                padding: 10px 20px;
                border-radius: 12px;
            }
            QPushButton:hover {
                background-color: #00eaff;  /* neon cyan hover */
                color: black;
            }
            """
        )

        bridge_btn.clicked.connect(self.load_main_ui)

        # ===== Assemble splash layout =====
        layout.addSpacing(10)
        layout.addWidget(credits)
        layout.addWidget(scroll, stretch=1)
        layout.addSpacing(15)
        layout.addWidget(slogan)
        layout.addSpacing(15)
        layout.addWidget(bridge_btn, alignment=Qt.AlignmentFlag.AlignCenter)

    def type_next_character(self):
        if self.line_index >= len(self.splash_lines):
            # When ASCII is done, start typing "the future is chrome"
            if not hasattr(self, "slogan_started"):
                self.slogan_started = True
                self.start_slogan_typing()
            return

        SPEED = 500  # number of characters per tick

        line = self.splash_lines[self.line_index]

        # Take characters chunk-by-chunk
        if self.char_index < len(line):
            end = min(self.char_index + SPEED, len(line))
            self.typed_text += line[self.char_index:end]
            self.char_index = end
        else:
            self.typed_text += "\n"
            self.line_index += 1
            self.char_index = 0

        # 🔥 Always update the label each tick
        self.logo_label.setText(
            f"<pre style='font-family:Consolas; font-size:5px; line-height:5px;'>{self.typed_text}</pre>"
        )

    def center_on_screen(self, y_offset=40):
        screen = QGuiApplication.primaryScreen().availableGeometry()
        frame = self.frameGeometry()
        frame.moveCenter(screen.center())

        # Move slightly UP to avoid taskbar overlap
        pos = frame.topLeft()
        pos.setY(pos.y() - y_offset)

        self.move(pos)

    def start_slogan_typing(self):
        self.slogan_text = "   "
        self.slogan_out = ""
        self.slogan_index = 0

        self.slogan_label = QLabel("")
        self.slogan_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.slogan_label.setStyleSheet("color:#ff40c0; font-size:26px; font-family:Consolas;")
        self.splash_page.layout().addWidget(self.slogan_label)

        self.slogan_timer = QTimer()
        self.slogan_timer.timeout.connect(self.type_slogan_character)
        self.slogan_timer.start(10)

    def type_slogan_character(self):
        if self.slogan_index < len(self.slogan_text):
            self.slogan_out += self.slogan_text[self.slogan_index]
            self.slogan_index += 1
            self.slogan_label.setText(self.slogan_out)
        else:
            self.slogan_timer.stop()

    def load_main_ui(self):
        self.stack.setCurrentWidget(self.main_page)
        self.resize(900, 450)
        self.center_on_screen(y_offset=20)  # higher than splash

        #self.stack.setCurrentWidget(self.main_page)
        #self.setGeometry(300, 200, 1000, 700)  # <-- smaller main UI window

    def apply_dark_style(self):
        self.setStyleSheet("""
            QMainWindow {
                background-color: black;
            }
            QLabel {
                color: white;
                font-family: Consolas;
            }
        """)

    # ========================================================
    #   Numerical UI logic
    # ========================================================

    def parse_inputs(self):
        try:
            b_val = float(self.b_edit.text())
            P0_val = float(self.P0_edit.text())
            T0_val = float(self.T0_edit.text())
            gamma_val = float(self.gamma_edit.text())
            Pa_val = float(self.Pa_edit.text())
            Pe_val = float(self.Pe_edit.text())
            R_val = float(self.R_edit.text())
            F_val = float(self.F_edit.text())
        except ValueError:
            QMessageBox.warning(self, "Input error", "Please enter valid numeric values.")
            return None
        return b_val, P0_val, T0_val, gamma_val, Pa_val, Pe_val, R_val, F_val

    def run_computation(self):
        parsed = self.parse_inputs()
        if parsed is None:
            return



        b_val, P0_val, T0_val, gamma_val, Pa_val, Pe_val, R_val, F_val = parsed

        try:
            # ===== RAO BELL PATH =====
            if self.nozzle_mode_cb.isChecked():

                Cf = None
                eps = None
                CR = None
                Lstar_m = None

                if self.cea_mode_cb.isChecked():
                    try:
                        Cf = float(self.Cf_edit.text())
                        eps = float(self.eps_edit.text())
                        CR = float(self.CR_edit.text())
                        Lstar_m = float(self.Lstar_edit.text())


                    except ValueError:
                        QMessageBox.warning(
                            self,
                            "CEA Input Error",
                            "CEA mode is enabled, but one or more CEA inputs are invalid."
                        )
                        return



                else:
                    CR = None
                    Lstar_mm = None

                res = rao_bell_nozzle(
                    P0_val, T0_val, gamma_val,
                    Pa_val, Pe_val, R_val, F_val,
                    use_cea=self.cea_mode_cb.isChecked(),
                    Cf=Cf,
                    eps_user=eps,
                    CR=CR,
                    Lstar_m=Lstar_m,
                    include_chamber=self.show_chamber_cb.isChecked()  # 🔥 ADD
                )

                self.last_result = res
                self.display_results_rao(res)
                self.plot_2d_rao(res)
                self.plot_3d_rao(res)
                self.populate_table(res["x"], res["r"])
                self.point_selector.setCurrentText("Rao Bell Points")

            # ===== MOC PATH =====
            else:
                res = shockwve_compute(
                    b_val, P0_val, T0_val,
                    gamma_val, Pa_val, Pe_val,
                    R_val, F_val
                )

                self.last_result = res
                self.display_results(res)
                self.plot_2d(res)
                self.plot_3d(res)
                self.populate_table(res["x_full"], res["r_full"])
                self.point_selector.setCurrentText("MOC Points")

        except Exception as e:
            QMessageBox.critical(self, "Computation error", str(e))

    def update_2d_plot(self):
        if not hasattr(self, "last_result"):
            return

        if self.nozzle_mode_cb.isChecked():
            self.plot_2d_rao(self.last_result)
        else:
            self.plot_2d(self.last_result)

    def on_nozzle_mode_changed(self):
        rao_on = self.nozzle_mode_cb.isChecked()

        # Enable/disable Rao-only options
        self.cea_mode_cb.setEnabled(rao_on)
        self.show_chamber_cb.setEnabled(rao_on)

        if not rao_on:
            self.cea_mode_cb.setChecked(False)
            self.show_chamber_cb.setChecked(False)

        # Clear previous results
        if hasattr(self, "last_result"):
            del self.last_result

        # Clear plots
        self.canvas_2d.ax.clear()
        self.canvas_2d.draw()

        self.canvas_3d.ax.clear()
        self.canvas_3d.draw()

    def display_results(self, res):
        self.output_text.clear()

        lines = []
        lines.append("<span style='color:#00eaff; font-weight:bold;'>MOC Nozzle Parameters</span>")
        lines.append(f"Exit Mach number (Me):        {res['Me']:.3f}")
        lines.append(f"Expansion Ratio (Ae/At):      {res['Ae_At']:.3f}")
        lines.append(f"Exhaust Velocity (Ve):        {res['Ve']:.2f} m/s")
        lines.append(f"Mass flow rate (mdot ideal):  {res['mdot_ideal']:.4f} kg/s")
        lines.append(f"Throat area (At):             {res['At']:.6e} m^2")
        lines.append(f"Exit area (Ae):               {res['Ae']:.6e} m^2")
        lines.append(f"Throat radius (R_throat):     {res['R_throat_mm']:.3f} mm")
        lines.append(f"Exit radius  (R_exit):        {res['R_exit_mm']:.3f} mm")
        lines.append(f"Exhaust Temperature (Te):        {res['Te']:.2f} K")
        lines.append("")
        lines.append("=== Chamber / Convergent ===")
        lines.append(f"Contraction ratio (Ac/At):    {res['CR']:.3f}")
        lines.append(f"Chamber area  Ac:             {res['Ac']:.6e} m^2")
        lines.append(f"Chamber radius Rc:            {res['Rc_m']:.6f} m")
        lines.append(f"Chamber diameter dc:          {res['dc_mm']:.3f} mm")
        lines.append(f"Convergent length lN:         {res['lN_mm']:.2f} mm")
        lines.append("")
        lines.append("=== MOC / Divergent Section ===")
        lines.append(f"MOC exit radius r_e:          {res['r_exit_m_moc']:.6f} m")
        lines.append(f"MOC exit area Ae_moc:         {res['Ae_moc']:.6e} m^2")
        lines.append(f"MOC Ae/At (geom):             {res['Ae_At_moc']:.3f}")
        lines.append(f"MOC throat area At_moc:       {res['At_moc']:.6e} m^2")
        lines.append(f"MOC throat radius Rt_moc:     {res['Rt_m_moc']:.6f} m")
        lines.append(f"MOC exhaust velocity Ve_moc:  {res['Ve_moc']:.2f} m/s")
        lines.append(f"MOC mass flow mdot_moc:       {res['mdot_moc']:.6f} kg/s")
        lines.append(f"Divergent scale factor s_opt: {res['s_opt']:.5f}")
        lines.append(f"Back-computed p_exit:         {res['p_exit']:.1f} Pa")
        lines.append("")
        lines.append("=== Thrust Consistency Warning ===")
        lines.append("Because the throat area is sized from the design thrust,")
        lines.append("the mass flow is fixed by choked flow and does NOT adjust")
        lines.append("when you choose an over/under-expanded exit pressure.")
        lines.append("")
        lines.append(f"Requested thrust (design):    {self.F_edit.text()} N")
        lines.append(f"Actual thrust delivered:       {res['F_actual']:.2f} N")
        lines.append(f"Fixed mass flow rate:          {res['mdot_fixed']:.5f} kg/s")

        # ===== Expansion Condition Analysis (COLOR CODED) =====
        Pe_val = float(self.Pe_edit.text())
        Pa_val = float(self.Pa_edit.text())

        if Pe_val > Pa_val:
            expand_state = "<span style='color:#00ff00; font-weight:bold;'>UNDER-EXPANDED</span>"
            expand_msg = "<span style='color:#00ff00;'>Pe &gt; Pa → Pressure thrust positive → Slightly higher thrust.</span>"
        elif Pe_val < Pa_val:
            expand_state = "<span style='color:#ff0066; font-weight:bold;'>OVER-EXPANDED</span>"
            expand_msg = "<span style='color:#ff0066;'>Pe &lt; Pa → Pressure thrust negative → Thrust reduction.</span>"
        else:
            expand_state = "<span style='color:#00eaff; font-weight:bold;'>PERFECT EXPANSION</span>"
            expand_msg = "<span style='color:#00eaff;'>Pe = Pa → Ideal thrust condition.</span>"

        lines.append("")
        lines.append("<b>=== Expansion Condition ===</b>")
        lines.append(f"Expansion regime: {expand_state}")
        lines.append(expand_msg)

        self.output_text.setHtml("<br>".join(lines))

    def display_results_rao(self, res):
        self.output_text.clear()

        lines = []
        lines.append("<span style='color:#00ff00; font-weight:bold;'>Rao Bell Nozzle Parameters</span>")
        lines.append(f"Exit Mach number (Me):        {res['Me']:.3f}")
        lines.append(f"Expansion Ratio (Ae/At):      {res['eps']:.3f}")
        lines.append(f"Exhaust Velocity (Ve):        {res['Ve']:.2f} m/s")
        lines.append(f"Mass Flow Rate (mdot):        {res['mdot']:.4f} kg/s")
        lines.append(f"Actual Thrust:               {res['F_actual']:.2f} N")
        lines.append("")
        lines.append("=== Geometry ===")
        lines.append(f"Throat Radius:               {res['Rt_mm']:.3f} mm")
        lines.append(f"Exit Radius:                 {res['Re_mm']:.3f} mm")
        lines.append(f"Nozzle Length:               {res['LN_m'] * 1000:.2f} mm")
        if "Rc_m" in res and res.get("Lc_m") is not None:
            lines.append("")
            lines.append("=== CEA Chamber ===")
            lines.append(f"Contraction Ratio (Ac/At):   {res['CR']:.3f}")
            lines.append(f"Characteristic Length L*:   {res['Lstar_m']:.1f} m")
            lines.append(f"Chamber Radius Rc:           {res['Rc_m'] * 1000:.2f} mm")
            lines.append(f"Chamber Length Lc:           {res['Lc_m'] * 1000:.2f} mm")

        # Expansion condition
        Pe = float(self.Pe_edit.text())
        Pa = float(self.Pa_edit.text())

        if Pe > Pa:
            state = "<span style='color:#00ff00; font-weight:bold;'>UNDER-EXPANDED</span>"
            msg = "<span style='color:#00ff00;'>Pe &gt; Pa → Pressure thrust positive → Slightly higher thrust.</span>"
        elif Pe < Pa:
            state = "<span style='color:#ff0066; font-weight:bold;'>OVER-EXPANDED</span>"
            msg = "<span style='color:#ff0066;'>Pe &lt; Pa → Pressure thrust negative → Thrust reduction.</span>"
        else:
            state = "<span style='color:#00eaff; font-weight:bold;'>PERFECT EXPANSION</span>"
            msg = "<span style='color:#00eaff;'>Pe = Pa → Ideal thrust condition.</span>"

        lines.append("")
        lines.append("<b>=== Expansion Condition ===</b>")
        lines.append(f"Expansion Regime: {state}")
        lines.append(msg)

        self.output_text.setHtml("<br>".join(lines))

    def plot_2d(self, res):
        ax = self.canvas_2d.ax
        ax.clear()
        ax.set_facecolor("black")

        x_conv = res["x_conv"]
        r_conv = res["r_conv"]
        x_conn = res["x_conn"]
        r_conn = res["r_conn"]
        x_arc = res["x_arc"]
        r_arc = res["r_arc"]
        x_div = res["x_div"]
        y_div = res["y_div"]

        from scipy.interpolate import interp1d

        # Build wall interpolator (positive radius only)
        wall_interp = interp1d(
            x_div,
            y_div,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan
        )
# helper function
        def clip_line_to_wall(x, y, wall_func):
            """
            Clips a line (x, y) so it stops at first intersection with wall.
            Assumes y >= 0.
            """
            r_wall = wall_func(x)
            valid = np.isfinite(r_wall)

            if not np.any(valid):
                return None, None

            x = x[valid]
            y = y[valid]
            r_wall = r_wall[valid]

            # Find where line crosses wall
            diff = y - r_wall
            idx = np.where(diff >= 0)[0]

            if len(idx) == 0:
                return None, None

            i = idx[0]
            return x[: i + 1], y[: i + 1]

        if self.show_char_cb.isChecked():

            x_shift = res["x_div"][0]

            # ---- C+ characteristics (NEON CYAN) ----
            for xp, yp in res["C_plus"]:
                xp = np.array(xp) + x_shift
                yp = np.array(yp)

                ax.plot(xp, yp, color=NEON_CYAN, lw=0.8)
                ax.plot(xp, -yp, color=NEON_CYAN, lw=0.8)

            # ---- C- characteristics (NEON PINK) ----
            for xm, ym in res["C_minus"]:
                xm = np.array(xm) + x_shift
                ym = np.array(ym)

                x_clip, y_clip = clip_line_to_wall(xm, ym, wall_interp)
                if x_clip is not None:
                    ax.plot(x_clip, y_clip, color="#00f5ff", lw=0.8)
                    ax.plot(x_clip, -y_clip, color="#ff40c0", lw=0.8)

        ax.plot(x_conv, r_conv, color="#00aaff", linewidth=2, label="Convergent")
        ax.plot(x_conv, -r_conv, color="#00aaff", linewidth=1)


        ax.plot(x_conn, r_conn, color="#ff00ff", linewidth=2, label="Connector")
        ax.plot(x_conn, -r_conn, color="#ff00ff", linewidth=1)

        #ax.plot(x_arc, r_arc, color="#ffff00", linewidth=2, label="Circular Blend")
        #ax.plot(x_arc, -r_arc, color="#ffff00", linewidth=1)

        ax.plot(x_div, y_div, color="#00ffff", linewidth=2, label="MOC Divergent")
        ax.plot(x_div, -y_div, color="#00ffff", linewidth=1)

        ax.plot([], [], color="#00f5ff", linewidth=1.2, label="C⁺ Characteristics")
        ax.plot([], [], color="#ff40c0", linewidth=1.2, label="C⁻ Characteristics")

        ax.set_xlabel("Axial Distance (mm)", color="white")
        ax.set_ylabel("Radius (mm)", color="white", labelpad=-1)
        ax.set_title("Method of Characteristics Nozzle (2D)", color="white")

        ax.relim()
        ax.autoscale_view()

        if self.show_legend_cb.isChecked():
            ax.legend(
                facecolor="black",
                edgecolor="white",
                labelcolor="white"
            )

        ax.grid(True, color="#444444")

        ax.tick_params(colors="white")
        for spine in getattr(ax, "spines", {}).values():
            spine.set_color("white")

        self.canvas_2d.draw()

    def plot_2d_rao(self, res):
        ax = self.canvas_2d.ax
        ax.clear()
        ax.set_facecolor("black")

        # 1) Plot FULL contour first (this guarantees the cylindrical section shows)
        x_full = np.asarray(res["x"])
        r_full = np.asarray(res["r"])

        ax.plot(x_full, r_full, color="lime", lw=2, label="Rao Bell (Full Contour)")
        ax.plot(x_full, -r_full, color="lime", lw=1)

        # 2) Optional: overlay segments (helps you visually debug the join)
        if "x_chamber" in res and "r_chamber" in res and len(res["x_chamber"]) > 0:
            ax.plot(res["x_chamber"], res["r_chamber"], color="#ffaa00", lw=2, label="Chamber")
            ax.plot(res["x_chamber"], -res["r_chamber"], color="#ffaa00", lw=1)

        ax.plot(res["x_conv"], res["r_conv"], color="yellow", lw=2, label="Convergent")
        ax.plot(res["x_conv"], -res["r_conv"], color="yellow", lw=1)

        ax.plot(res["x_throat"], res["r_throat"], color="#ff00ff", lw=2, label="Throat Arc")
        ax.plot(res["x_throat"], -res["r_throat"], color="#ff00ff", lw=1)

        # 3) Proper zoom MUST be set before draw
        xmin = float(np.min(x_full))
        xmax = float(np.max(x_full))
        rmax = float(np.max(np.abs(r_full)))

        x_pad = 0.1 * (xmax - xmin)
        r_pad = 0.2 * rmax

        ax.set_xlim(xmin - x_pad, xmax + x_pad)
        ax.set_ylim(-(rmax + r_pad), (rmax + r_pad))

        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Axial Distance (mm)", color="white")
        ax.set_ylabel("Radius (mm)", color="white", labelpad=-1)
        ax.set_title("Rao Bell Nozzle (2D)", color="white")
        ax.grid(True, color="#444")

        if self.show_legend_cb.isChecked():
            ax.legend(facecolor="black", edgecolor="white", labelcolor="white")

        ax.tick_params(colors="white")
        for spine in getattr(ax, "spines", {}).values():
            spine.set_color("white")

        self.canvas_2d.draw()

    def plot_3d(self, res):
        ax = self.canvas_3d.ax
        ax.clear()
        ax.set_facecolor("black")

        X = res["X3D"]
        Y = res["Y3D"]
        Z = res["Z3D"]
        ax.plot_wireframe(X, Y, Z, linewidth=0.5, color="cyan")
        # --- Neon Pink → Cyan Gradient Wireframe ---
        num_slices = X.shape[1]  # number of axial slices

        # Colors: start = neon pink (#ff00ff), end = cyan (#00ffff)
        pink = np.array([1.0, 0.0, 1.0])   # RGB
        cyan = np.array([0.0, 1.0, 1.0])   # RGB

        # Linear gradient along the nozzle axis
        t = np.linspace(0, 1, num_slices)
        gradient_colors = [pink * (1 - g) + cyan * g for g in t]

        # Draw each constant-x ring as its own colored wire
        stride = 26  # try 4, 5, or 6 for more speed
        for i in range(0, num_slices, stride):
            ax.plot(
                X[:, i], Y[:, i], Z[:, i],
                color=gradient_colors[i],
                linewidth=0.9
            )

            ax.plot(
                X[:, i], Y[:, i], Z[:, i],
                color=gradient_colors[i],
                linewidth=0.9
            )



        ax.set_xlabel("Axial [mm]", color="white")
        ax.set_ylabel("Y [mm]", color="white")
        ax.set_zlabel("Z [mm]", color="white")
        ax.set_title("Method of Characteristics Nozzle (3D)", color="white")

        ax.tick_params(colors="white")
        # ax.spines not present on 3D, so no spine coloring here

        self.canvas_3d.draw()

    def plot_3d_rao(self, res):
        ax = self.canvas_3d.ax
        ax.clear()
        ax.set_facecolor("black")

        x = np.asarray(res["x"])
        r = np.maximum(np.asarray(res["r"]), 1e-6)

        theta = np.linspace(0, 2 * np.pi, 120)
        X, T = np.meshgrid(x, theta)
        Rg = np.tile(r, (theta.size, 1))
        Y = Rg * np.cos(T)
        Z = Rg * np.sin(T)

        ax.plot_wireframe(X, Y, Z, color="lime", linewidth=0.6)

        # Force framing to include chamber + full bell
        xmin = float(np.min(x))
        xmax = float(np.max(x))
        rmax = float(np.max(r))

        x_pad = 0.1 * (xmax - xmin)
        r_pad = 0.2 * rmax

        ax.set_xlim(xmin - x_pad, xmax + x_pad)
        ax.set_ylim(-(rmax + r_pad), (rmax + r_pad))
        ax.set_zlim(-(rmax + r_pad), (rmax + r_pad))

        ax.set_title("Rao Bell Nozzle (3D)", color="white")
        ax.set_xlabel("Axial (mm)", color="white")
        ax.set_ylabel("Y (mm)", color="white")
        ax.set_zlabel("Z (mm)", color="white")

        self.canvas_3d.draw()

    def show_neon_warning(self, message: str):
        popup = QMessageBox(self)
        popup.setWindowTitle("SHOCKWVE – WARNING")
        popup.setWindowModality(Qt.WindowModality.ApplicationModal)
        popup.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint)

        popup.setText(
            f"<span style='color:#ff40c0; font-size:20px; font-family:Consolas;'>{message}</span>"
        )

        popup.setStyleSheet("""
            QMessageBox {
                background-color: #050505;
                border: 3px solid #ff40c0;
            }
            QMessageBox QLabel {
                color: #ff40c0;
                font-size: 18px;
                font-family: Consolas;
            }
            QPushButton {
                background-color: #ff40c0;
                color: black;
                font-size: 16px;
                font-family: Consolas;
                padding: 6px 20px;
                border-radius: 10px;
                border: 2px solid #00f5ff;
            }
            QPushButton:hover {
                background-color: #00f5ff;
                color: black;
                border: 2px solid #ff40c0;
            }
        """)

        popup.addButton("OK", QMessageBox.ButtonRole.AcceptRole)
        popup.exec()

    def show_neon_popup(self, message: str):
        popup = QMessageBox()
        popup.setWindowTitle("SHOCKWVE – EXPORT")
        popup.setWindowModality(Qt.WindowModality.ApplicationModal)
        popup.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint)

        popup.setText(
            f"<span style='color:#00f5ff; font-size:20px; font-family:Consolas;'>{message}</span>"
        )

        popup.setStyleSheet("""
            QMessageBox {
                background-color: #050505;
                border: 3px solid #ff00ff;
            }
            QMessageBox QLabel {
                color: #00f5ff;
                font-size: 18px;
                font-family: Consolas;
            }
            QPushButton {
                background-color: #ff00ff;
                color: black;
                font-size: 16px;
                font-family: Consolas;
                padding: 6px 18px;
                border-radius: 8px;
                border: 2px solid #00f5ff;
            }
            QPushButton:hover {
                background-color: #00f5ff;
                color: black;
                border: 2px solid #ff00ff;
            }
        """)

        popup.addButton("OK", QMessageBox.ButtonRole.AcceptRole)

        popup.exec()

    def export_contour_csv(self):

        # 🔥 HARD GUARD — must have run ShockWVE first
        if not hasattr(self, "last_result"):
            self.show_neon_warning(
                "NO DATA TO EXPORT<br><br>"
                "Run <b>SHOCKWVE</b> before exporting the nozzle contour."
            )
            return

        try:
            parsed = self.parse_inputs()
            if parsed is None:
                return
            b_val, P0_val, T0_val, gamma_val, Pa_val, Pe_val, R_val, F_val = parsed
            res = shockwve_compute(b_val, P0_val, T0_val, gamma_val, Pa_val, Pe_val, R_val, F_val)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to recompute contour:\n{e}")
            return

        if self.nozzle_mode_cb.isChecked():
            x_full = self.last_result["x"]
            r_full = self.last_result["r"]
        else:
            x_full = res["x_full"]
            r_full = res["r_full"]

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Contour XLSX",
            "Nozzle_2D_Contour.xlsx",
            "Excel Files (*.xlsx)"
        )
        if not path:
            return

        try:
            from openpyxl import Workbook
            wb = Workbook()
            ws = wb.active
            ws.title = "NozzleContour"
            ws.append(["Axial_mm", "Radius_mm"])

            for x, r in zip(x_full, r_full):
                ws.append([float(x), float(r)])

            wb.save(path)
            self.show_neon_popup("XLSX EXPORTED SUCCESSFULLY")

        except Exception as e:
            QMessageBox.critical(self, "Export error", f"REAL ERROR: {repr(e)}")
            print("EXPORT ERROR >>>", repr(e))





    def load_saved_engines(self):
        import json, os
        ENGINE_FILE = "saved_engines.json"

        if not os.path.exists(ENGINE_FILE):
            return {}

        with open(ENGINE_FILE, "r") as f:
            return json.load(f)

    def save_engines_to_file(self, engines):
        import json
        ENGINE_FILE = "saved_engines.json"

        with open(ENGINE_FILE, "w") as f:
            json.dump(engines, f, indent=4)

    def get_current_engine_data(self):
        return {
            "b": self.b_edit.text(),
            "P0": self.P0_edit.text(),
            "T0": self.T0_edit.text(),
            "gamma": self.gamma_edit.text(),
            "Pa": self.Pa_edit.text(),
            "Pe": self.Pe_edit.text(),
            "R": self.R_edit.text(),
            "F": self.F_edit.text(),
            "Cf": self.Cf_edit.text(),
            "eps": self.eps_edit.text(),
            "CR": self.CR_edit.text(),
            "Lstar": self.Lstar_edit.text(),
            "use_rao": self.nozzle_mode_cb.isChecked(),
            "use_cea": self.cea_mode_cb.isChecked(),
            "include_chamber": self.show_chamber_cb.isChecked()
        }

    def apply_engine_data(self, data):
        self.b_edit.setText(data.get("b", ""))
        self.P0_edit.setText(data.get("P0", ""))
        self.T0_edit.setText(data.get("T0", ""))
        self.gamma_edit.setText(data.get("gamma", ""))
        self.Pa_edit.setText(data.get("Pa", ""))
        self.Pe_edit.setText(data.get("Pe", ""))
        self.R_edit.setText(data.get("R", ""))
        self.F_edit.setText(data.get("F", ""))

        self.Cf_edit.setText(data.get("Cf", ""))
        self.eps_edit.setText(data.get("eps", ""))
        self.CR_edit.setText(data.get("CR", ""))
        self.Lstar_edit.setText(data.get("Lstar", ""))

        self.nozzle_mode_cb.setChecked(data.get("use_rao", False))
        self.cea_mode_cb.setChecked(data.get("use_cea", False))
        self.show_chamber_cb.setChecked(data.get("include_chamber", False))





    def refresh_engine_list(self):
        engines = self.load_saved_engines()
        self.engine_selector.clear()
        self.engine_selector.addItems(engines.keys())

    def save_current_engine(self):
        from PyQt6.QtWidgets import QInputDialog

        dialog = QInputDialog(self)
        dialog.setWindowTitle("SHOCKWVE – SAVE ENGINE")
        dialog.setLabelText("Engine Name:")
        dialog.setStyleSheet("""
            QInputDialog {
                background-color: #050505;
                border: 2px solid #ff40c0;
            }
            QLabel {
                color: #00f5ff;
                font-family: Consolas;
                font-size: 14px;
            }
            QLineEdit {
                background-color: black;
                color: #ff40c0;
                border: 1px solid #00f5ff;
                padding: 4px;
                font-family: Consolas;
            }
            QPushButton {
                background-color: #111111;
                color: #00f5ff;
                border: 1px solid #00f5ff;
                padding: 6px 14px;
                font-family: Consolas;
                border-radius: 6px;
            }
            QPushButton:hover {
                background-color: #ff40c0;
                color: black;
                border: 1px solid #ff40c0;
            }
        """)

        ok = dialog.exec()

        if ok:
            name = dialog.textValue()
            if not name:
                return

            engines = self.load_saved_engines()
            engines[name] = self.get_current_engine_data()
            self.save_engines_to_file(engines)
            self.refresh_engine_list()

    def load_selected_engine(self):
        name = self.engine_selector.currentText()
        engines = self.load_saved_engines()
        if name in engines:
            self.apply_engine_data(engines[name])

    def delete_selected_engine(self):
        name = self.engine_selector.currentText()
        engines = self.load_saved_engines()
        if name in engines:
            del engines[name]
            self.save_engines_to_file(engines)
            self.refresh_engine_list()


# ============================================================
#  Entry point
# ============================================================

def main():
    app = QApplication(sys.argv)
    from PyQt6.QtGui import QPalette, QColor

    #palette = QPalette()
    #palette.setColor(QPalette.ColorRole.Window, QColor(0, 0, 0))
    #palette.setColor(QPalette.ColorRole.WindowText, QColor(0, 245, 255))
    #palette.setColor(QPalette.ColorRole.Base, QColor(10, 10, 10))
    #palette.setColor(QPalette.ColorRole.Text, QColor(0, 245, 255))
    #palette.setColor(QPalette.ColorRole.Button, QColor(0, 0, 0))
    #palette.setColor(QPalette.ColorRole.ButtonText, QColor(0, 245, 255))
    #palette.setColor(QPalette.ColorRole.Highlight, QColor(0, 245, 255))

    #app.setPalette(palette)

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(0, 0, 0))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(0, 245, 255))
    palette.setColor(QPalette.ColorRole.Base, QColor(10, 10, 10))
    palette.setColor(QPalette.ColorRole.Text, QColor(0, 245, 255))
    palette.setColor(QPalette.ColorRole.Button, QColor(0, 245, 255))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(0, 245, 255))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(0, 245, 255))

    app.setPalette(palette)

    # 🔥 FORCE identical checkbox rendering
    app.setStyle("Fusion")
    win = ShockWveApp()
    win.show()
    QTimer.singleShot(0, win.center_on_screen)

    sys.exit(app.exec())




if __name__ == "__main__":
    main()
