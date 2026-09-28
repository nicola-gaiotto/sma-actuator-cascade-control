"""Independent Python check of the SMA cascade controller.

Re-implements the equilibrium computation, the linearised inner (thermal) and
outer (mechanical) plants and the loop-shaping design of
matlab/sma_cascade_control.m, then:
  * evaluates phase margins, step-response metrics and disturbance attenuation
    for the original exam tuning and for the revised tuning;
  * simulates the full NONLINEAR plant in closed loop.

Usage (from the repository root):
    pip install numpy scipy matplotlib
    python tools/verify_design.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sg
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.3})

# ---------------------------------------------------------------- parameters
T0 = 273.15
Re, Le = 4.8, 55e-3
c, l, d, Cth, Rth = 6.2, 15e-2, 0.2e-3, 25.1e-3, 16.7
Tp, h, Tamb = 15 + T0, 6.55, 24 + T0
r1, r2, J, m, beta = 0.45e-2, 4.5e-2, 2e-4, 15e-3, 1.5
Kmax, Tavg, Tdiff, z_eq = 3.92e3, 70.0, 20.0, 4e-2
m_eq, g = J / r2 + m * r2, 9.81
S = np.pi * d * l
Ramb = 1 / (h * S)
rho = r2 / r1
dl_max = 0.04 * l

# ---------------------------------------------------------------- equilibrium
K_star = m * g * r2 / (r1 * (dl_max - z_eq * r1 / r2))
T_eq_C = Tavg + Tdiff / c * np.log(K_star / (Kmax - K_star))
T_eq = T_eq_C + T0
I_eq = np.sqrt((h * S * (T_eq - Tamb) + (T_eq - Tp) / Rth) / Re)
V_eq = Re * I_eq


# ---------------------------------------------------------------- tf helpers
def tf(num, den):
    return (np.atleast_1d(np.asarray(num, float)), np.atleast_1d(np.asarray(den, float)))


def series(*sys):
    num, den = np.array([1.0]), np.array([1.0])
    for n, dd in sys:
        num, den = np.polymul(num, n), np.polymul(den, dd)
    return num, den


def feedback(L):
    return L[0], np.polyadd(L[1], L[0])


def freq(sys, w):
    return np.polyval(sys[0], 1j * w) / np.polyval(sys[1], 1j * w)


def ss_to_tf(A, B, C):
    num, den = sg.ss2tf(A, B, C, np.zeros((1, 1)))
    return tf(np.trim_zeros(np.round(num[0], 12), "f"), den)


# Inner plant: voltage -> wire temperature
At = np.array([[-Re / Le, 0], [2 * Re * I_eq / Cth, -1 / (Rth * Cth) - 1 / (Ramb * Cth)]])
Gt = ss_to_tf(At, np.array([[1 / Le], [0]]), np.array([[0, 1]]))
# Outer plant: temperature -> load position
n_ = np.exp(c / Tdiff * (T_eq_C - Tavg))
Am = np.array([[0, 1], [-Kmax * (1 - 1 / (1 + n_)) * r1 / (m_eq * rho), -beta * r2 / m_eq]])
Bm = np.array([[0], [Kmax * r1 / m_eq * (dl_max - z_eq / rho) * c / Tdiff * n_ / (1 + n_) ** 2]])
Gm = ss_to_tf(Am, Bm, np.array([[1, 0]]))


def phase_deg(sys, w):
    ang = np.degrees(np.angle(freq(sys, w)))
    return ang - 360 if ang > 0 else ang


def design(Mf1, Mf2, wc1=290.0, Ta=1.7):
    """Same algorithm as the MATLAB script (evaluated exactly at wc)."""
    Gts = series(Gt, tf([1], [1, 0]))
    mod = 20 * np.log10(abs(freq(Gts, wc1)))
    pd = -180 + Mf1 - phase_deg(Gts, wc1)
    r = np.tan(np.radians(pd / 2))
    tauz = r / wc1
    mu = 10 ** (-(2 * 20 * np.log10(np.sqrt(1 + r**2)) + mod) / 20)
    taup = 1 / (10 * wc1)
    R1 = tf(mu * np.polymul([tauz, 1], [tauz, 1]), np.polymul(np.polymul([1, 0], [taup, 1]), [1.2e-5, 1]))

    wc2 = 460 / (Ta * 65)
    Gms = series(Gm, tf([1], [1, 0]))
    mod = 20 * np.log10(abs(freq(Gms, wc2)))
    pd = -180 + Mf2 - phase_deg(Gms, wc2)
    r = np.tan(np.radians(pd))
    mu2 = 10 ** (-(20 * np.log10(np.sqrt(1 + r**2)) + mod) / 20)
    R2 = tf(mu2 * np.array([r / wc2, 1]), np.polymul(np.polymul([1, 0], [0.024, 1]), [1e-4, 1]))
    return R1, R2


def crossover(L):
    w = np.logspace(-2, 6, 400000)
    H = freq(L, w)
    mag = 20 * np.log10(abs(H))
    k = np.where(np.diff(np.sign(mag)) != 0)[0][0]
    ph = np.degrees(np.unwrap(np.angle(H)))
    return w[k], 180 + ph[k]


def step_metrics(F, t_end):
    t = np.linspace(0, t_end, 200001)
    _, y = sg.step(F, T=t)
    out = np.where(abs(y - y[-1]) > 0.01 * abs(y[-1]))[0]
    return (y.max() - y[-1]) / y[-1] * 100, t[out[-1]], t, y


def evaluate(R1, R2):
    L1 = series(R1, Gt)
    F1 = feedback(L1)
    L2 = series(R2, F1, Gm)
    wc1, pm1 = crossover(L1)
    wc2, pm2 = crossover(L2)
    os2, ts2, _, _ = step_metrics(feedback(L2), 5)
    os1, ts1, _, _ = step_metrics(F1, 0.2)
    return dict(wc1=wc1, pm1=pm1, wc2=wc2, pm2=pm2, os1=os1, ts1=ts1, os2=os2, ts2=ts2,
                L1_5=20 * np.log10(abs(freq(L1, 5))), L2_20=20 * np.log10(abs(freq(L2, 20))))


designs = {"original (exam)": design(65, 65), "revised": design(72, 75)}

print(f"Equilibrium: I* = {I_eq:.4f} A, T* = {T_eq_C:.2f} degC, z* = {z_eq*100:.1f} cm, V* = {V_eq:.4f} V")
print(f"Gt(s) = {Gt[0][-1]:.1f} / (s^2 + {Gt[1][1]:.2f} s + {Gt[1][2]:.2f}), poles {np.round(np.roots(Gt[1]), 3)}")
print(f"Gm(s) = {Gm[0][-1]:.4f} / (s^2 + {Gm[1][1]:.2f} s + {Gm[1][2]:.2f}), poles {np.round(np.roots(Gm[1]), 3)}")
for name, (R1, R2) in designs.items():
    e = evaluate(R1, R2)
    print(f"\n[{name}]")
    print(f"  inner: wc = {e['wc1']:.0f} rad/s, PM = {e['pm1']:.1f} deg, |L1(j5)| = {e['L1_5']:.1f} dB, "
          f"step overshoot = {e['os1']:.1f} %")
    print(f"  outer: wc = {e['wc2']:.2f} rad/s, PM = {e['pm2']:.1f} deg, |L2(j20)| = {e['L2_20']:.1f} dB, "
          f"overshoot = {e['os2']:.1f} %, Ts(1%) = {e['ts2']:.2f} s")

R1, R2 = designs["revised"]
R1o, R2o = designs["original (exam)"]

# ---------------------------------------------------------------- figures (linear)
L1, L2 = series(R1, Gt), series(R2, feedback(series(R1, Gt)), Gm)
fig, axes = plt.subplots(2, 2, figsize=(12, 6.5), sharex="col")
for col, (L, w, title, zones) in enumerate([
        (L1, np.logspace(-1, 4, 2000), "Inner loop $L_1 = R_1 G_T$", [((0.1, 5), (-60, 60))]),
        (L2, np.logspace(-1, 2, 2000), "Outer loop $L_2 = R_2 F_1 G_M$", [((0.1, 4), (-60, 0)), ((20, 100), (-20, 60))])]):
    H = freq(L, w)
    axes[0, col].semilogx(w, 20 * np.log10(abs(H)), lw=1.8)
    for (w0, w1), (a0, a1) in zones:
        axes[0, col].fill_between([w0, w1], a0, a1, color="red", alpha=0.15, lw=0)
    axes[0, col].axhline(0, color="k", lw=0.6)
    axes[0, col].set(title=title, ylabel="Magnitude [dB]")
    axes[1, col].semilogx(w, np.degrees(np.unwrap(np.angle(H))), lw=1.8)
    axes[1, col].axhline(-180 + 65, color="red", ls="--", lw=1, label="$-180° + 65°$")
    axes[1, col].set(xlabel="Frequency [rad/s]", ylabel="Phase [deg]")
    axes[1, col].legend(fontsize=8)
save = lambda f, n: (f.tight_layout(), f.savefig(OUT / n, bbox_inches="tight"), plt.close(f), print("saved", n))
save(fig, "loop_bode.png")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 3.8))
for (R1_, R2_), style, lab in [((R1o, R2o), "--", "original"), ((R1, R2), "-", "revised")]:
    F1 = feedback(series(R1_, Gt))
    _, _, t, y = step_metrics(F1, 0.08)
    a1.plot(t * 1000, y, style, lw=1.8, label=lab)
    _, _, t, y = step_metrics(feedback(series(R2_, F1, Gm)), 3)
    a2.plot(t, y, style, lw=1.8, label=lab)
a1.set(xlabel="Time [ms]", ylabel="Normalised temperature", title="Inner loop – closed-loop step response")
a2.axhspan(0.99, 1.01, color="green", alpha=0.12, label="±1 % band")
a2.axhline(1.05, color="red", ls=":", lw=1, label="5 % overshoot limit")
a2.axvline(1.8, color="gray", ls=":", lw=1, label="$T_{a,1\\%}$ = 1.8 s limit")
a2.set(xlabel="Time [s]", ylabel="Normalised position", title="Outer loop – closed-loop step response")
a1.legend(fontsize=8)
a2.legend(fontsize=8, loc="lower right")
save(fig, "step_responses.png")

# ---------------------------------------------------------------- nonlinear closed loop
A1c, B1c, C1c, D1c = sg.tf2ss(*R1)
A2c, B2c, C2c, D2c = sg.tf2ss(*R2)
n1, n2 = A1c.shape[0], A2c.shape[0]

w_step = lambda t: 5e-3 * (t >= 1.0)                 # +0.5 cm reference step at t = 1 s
noise = lambda t: 5e-4 * np.sin(20 * t)              # 0.05 cm position-sensor noise at 20 rad/s
dTp = lambda t: -5.0 * (t >= 4.0) + 5.0 * np.sin(5 * t) * (t >= 7.0)  # Peltier disturbances


def K_of_T(T):
    return Kmax * (1 - 1 / (1 + np.exp(c * (T - T0 - Tavg) / Tdiff)))


def rhs(t, x):
    I, T, z, zd = x[:4]
    xc2, xc1 = x[4:4 + n2], x[4 + n2:]
    e2 = (z_eq + w_step(t)) - (z + noise(t))
    dT_ref = (C2c @ xc2)[0] + D2c[0, 0] * e2
    e1 = (T_eq + dT_ref) - T
    V = V_eq + (C1c @ xc1)[0] + D1c[0, 0] * e1
    dI = (-Re * I + V) / Le
    dT = (Re * I**2 - h * S * (T - Tamb) - (T - (Tp + dTp(t))) / Rth) / Cth
    dzd = (r1 * K_of_T(T) * (dl_max - z / rho) - m * g * r2 - beta * r2 * zd) / m_eq
    return np.concatenate([[dI, dT, zd, dzd], A2c @ xc2 + B2c[:, 0] * e2, A1c @ xc1 + B1c[:, 0] * e1])


x0 = np.concatenate([[I_eq, T_eq, z_eq, 0.0], np.zeros(n2 + n1)])
sol = solve_ivp(rhs, (0, 10), x0, method="Radau", rtol=1e-7, atol=1e-10, max_step=2e-3, dense_output=True)
t = np.linspace(0, 10, 5001)
X = sol.sol(t)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True)
a1.plot(t, (z_eq + w_step(t)) * 100, "b", lw=1.8, label="Reference $z^\\star + w$")
a1.plot(t, X[2] * 100, "r", lw=1.4, label="Load position $z$ (nonlinear plant)")
a1.set(ylabel="Position [cm]", title="Nonlinear closed-loop simulation (revised tuning)")
a1.legend(fontsize=8, loc="lower right")
a2.plot(t, X[1] - T0, "tab:orange", lw=1.4, label="Wire temperature $T$")
a2.set(xlabel="Time [s]", ylabel="[°C]")
for tt, txt in [(1, "ref. step +0.5 cm"), (4, "$d_{Tp}$ step −5 °C"), (7, "$d_{Tp}$ = 5 sin(5t) °C")]:
    for a in (a1, a2):
        a.axvline(tt, color="gray", ls=":", lw=1)
    a2.text(tt + 0.05, a2.get_ylim()[1] - 0.3, txt, fontsize=7.5, va="top")
a2.legend(fontsize=8, loc="lower left")
save(fig, "nonlinear_simulation.png")
err = (z_eq + w_step(t)) - X[2]
print(f"\nNonlinear sim: final position {X[2, -1]*100:.3f} cm, "
      f"max |error| after settling (t>3 s) = {abs(err[t > 3]).max()*1e3:.3f} mm")
