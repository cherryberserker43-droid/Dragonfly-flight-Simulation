import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D

OUT = "/mnt/user-data/outputs"

def _sanity_check_rhat2_moment(n_strips=2000):
    r_hat = (np.arange(n_strips) + 0.5) / n_strips
    dr_hat = 1.0 / n_strips
    c_hat = np.ones_like(r_hat)
    rhat2_moment = np.sum((r_hat ** 2) * c_hat) * dr_hat
    return rhat2_moment

class _SmoothedTriangleGenerator:
    def __init__(self, period, sharpness=5.0, n_ref=4000):
        self.T = period
        self.omega = 2 * np.pi / period
        self.k = sharpness

        t_ref = np.linspace(0, period, n_ref, endpoint=False)
        g_ref = self._g(t_ref)
        phi_unscaled = np.concatenate([[0.0], np.cumsum((g_ref[:-1] + g_ref[1:]) / 2 * np.diff(t_ref))])
        self._unscaled_ptp = phi_unscaled.max() - phi_unscaled.min()
        self._t_ref = t_ref
        self._phi_unscaled_ref = phi_unscaled - phi_unscaled.mean()

    def _g(self, t):
        u = -np.cos(self.omega * t)
        return np.tanh(self.k * u)

    def _g_dot(self, t):
        u = -np.cos(self.omega * t)
        u_dot = self.omega * np.sin(self.omega * t)
        return self.k * u_dot * (1.0 - np.tanh(self.k * u) ** 2)

    def phi_dot(self, t, amplitude):
        vmax = amplitude / self._unscaled_ptp
        return vmax * self._g(t)

    def phi_ddot(self, t, amplitude):
        vmax = amplitude / self._unscaled_ptp
        return vmax * self._g_dot(t)

    def phi(self, t, amplitude):
        vmax = amplitude / self._unscaled_ptp
        t_mod = np.mod(t, self.T)
        return vmax * np.interp(t_mod, self._t_ref, self._phi_unscaled_ref)

def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)

class WingKinematics:
    def __init__(self, amplitude_deg, freq_hz,
                 aoa_mid_down_deg=45.0, aoa_mid_up_deg=45.0,
                 flip_start_frac=-0.05, flip_duration_frac=0.16,
                 deviation_amp_down_deg=0.0, deviation_amp_up_deg=0.0,
                 deviation_shape="oval",
                 phase_offset_frac=0.0):
        self.Phi = np.radians(amplitude_deg)
        self.freq = freq_hz
        self.T = 1.0 / freq_hz
        self.aoa_down = np.radians(aoa_mid_down_deg)
        self.aoa_up = np.radians(aoa_mid_up_deg)
        self.tau0 = flip_start_frac
        self.dtau = flip_duration_frac
        self.theta_down = np.radians(deviation_amp_down_deg)
        self.theta_up = np.radians(deviation_amp_up_deg)
        self.shape = deviation_shape
        self.phase_offset = phase_offset_frac * self.T
        self._tri = _SmoothedTriangleGenerator(self.T, sharpness=5.0)

    def phi(self, t):
        ts = t + self.phase_offset
        return self._tri.phi(ts, self.Phi)

    def phi_dot(self, t):
        ts = t + self.phase_offset
        return self._tri.phi_dot(ts, self.Phi)

    def phi_ddot(self, t):
        ts = t + self.phase_offset
        return self._tri.phi_ddot(ts, self.Phi)

    def reversal_phase(self, t):
        ts = (t + self.phase_offset + self.T / 4.0) % self.T
        return ts / self.T

    def _half_stroke_phase(self, t):
        cycle_frac = self.reversal_phase(t)
        is_down = cycle_frac < 0.5
        half_frac = np.where(is_down, cycle_frac / 0.5, (cycle_frac - 0.5) / 0.5)
        return half_frac, is_down

    def alpha(self, t):
        half_frac, is_down = self._half_stroke_phase(t)
        aoa_mid = np.where(is_down, self.aoa_down, self.aoa_up)
        aoa_next = np.where(is_down, self.aoa_up, self.aoa_down)

        flip_start = 1.0 + self.tau0
        flip_end = flip_start + self.dtau
        
        prev_flip_start = self.tau0
        prev_flip_end = self.tau0 + self.dtau

        out = np.array(aoa_mid, dtype=float, copy=True)

        in_prev_flip = (half_frac >= max(prev_flip_start, 0.0)) & (half_frac < max(prev_flip_end, 0.0)) & (prev_flip_end > 0.0)
        if np.any(in_prev_flip):
            local = (half_frac[in_prev_flip] - prev_flip_start) / self.dtau
            s = _smoothstep(local)
            prev_aoa_mid = np.where(is_down[in_prev_flip], self.aoa_up, self.aoa_down)
            out[in_prev_flip] = prev_aoa_mid + s * (aoa_mid[in_prev_flip] - prev_aoa_mid)

        in_next_flip = (half_frac >= min(flip_start, 1.0)) & (flip_start < 1.0)
        if np.any(in_next_flip):
            local = (half_frac[in_next_flip] - flip_start) / self.dtau
            s = _smoothstep(local)
            out[in_next_flip] = aoa_mid[in_next_flip] + s * (aoa_next[in_next_flip] - aoa_mid[in_next_flip])

        return out

    def theta(self, t):
        half_frac, is_down = self._half_stroke_phase(t)
        amp = np.where(is_down, self.theta_down, self.theta_up)
        if self.shape == "oval":
            return amp * np.sin(np.pi * half_frac)
        elif self.shape == "figure8":
            return amp * np.sin(2 * np.pi * half_frac)
        else:
            raise ValueError("deviation_shape must be 'oval' or 'figure8'")

    def state(self, t, dt=None):
        if dt is None:
            dt = self.T / 2000.0

        t = np.atleast_1d(np.asarray(t, dtype=float))

        def d1(f):
            return (f(t + dt) - f(t - dt)) / (2 * dt)

        def d2(f):
            return (f(t + dt) - 2 * f(t) + f(t - dt)) / (dt ** 2)

        phi = self.phi(t)
        alpha = self.alpha(t)
        theta = self.theta(t)

        return {
            "phi": phi, "phi_dot": self.phi_dot(t), "phi_ddot": self.phi_ddot(t),
            "alpha": alpha, "alpha_dot": d1(self.alpha), "alpha_ddot": d2(self.alpha),
            "theta": theta, "theta_dot": d1(self.theta),
        }

class WakeCaptureLiftModel:
    def __init__(self, AR_eff=4.0, switching_gain=6.0, t_cut=0.1,
                 cl_alpha_2d_per_deg=0.09):
        self.AReff = AR_eff
        self.n = switching_gain
        self.t_cut = t_cut
        self.cl_alpha_2d = cl_alpha_2d_per_deg * (180.0 / np.pi)

    def switching_function(self, half_frac):
        tc = self.t_cut
        n = self.n
        x = np.asarray(half_frac, dtype=float)

        s1 = 2.0 * np.exp(-n * (x / tc) ** n) - 1.0
        s2 = 2.0 * np.exp(-n * ((x - 0.5) / tc) ** n) - 1.0

        return np.where(x < 0.5, s1, s2)

    def coefficients(self, alpha, cycle_frac, CD0=0.4):
        S = self.switching_function(cycle_frac)
        cla2d = self.cl_alpha_2d
        denom = 1.0 - (cla2d / (np.pi * self.AReff)) * S
        CL_alpha = cla2d / denom

        CL = CL_alpha * np.sin(alpha) * np.cos(alpha)
        CD = CD0 + CL * np.tan(alpha)
        return CL, CD

class WingAeroModel:
    def __init__(self, span, mean_chord, n_strips=20,
                 rotation_axis_frac=0.25, rho=1.225,
                 AR_eff=None, switching_gain=6.0, t_cut=0.1, CD0=0.4):
        self.R = span
        self.c_bar = mean_chord
        self.rho = rho
        self.N = n_strips
        self.x0 = rotation_axis_frac
        self.CD0 = CD0

        self.r_hat = (np.arange(n_strips) + 0.5) / n_strips
        self.dr_hat = 1.0 / n_strips
        self.c_hat = np.ones_like(self.r_hat)
        self.S = self.R * self.c_bar

        AR_geom = self.R ** 2 / self.S
        self.AReff = AR_eff if AR_eff is not None else AR_geom

        self.wake = WakeCaptureLiftModel(AR_eff=self.AReff,
                                          switching_gain=switching_gain,
                                          t_cut=t_cut)

        self.Cr = np.pi * (0.75 - self.x0)
        self.I_c2 = np.sum(self.c_hat ** 2) * self.dr_hat

    def instantaneous_forces(self, state, cycle_frac):
        phi_dot = state["phi_dot"]
        phi_ddot = state["phi_ddot"]
        alpha = state["alpha"]
        alpha_dot = state["alpha_dot"]
        alpha_ddot = state["alpha_ddot"]

        rho = self.rho
        r = self.r_hat * self.R
        c = self.c_hat * self.c_bar

        phi_dot_b = np.asarray(phi_dot)[..., np.newaxis]
        alpha_b = np.asarray(alpha)[..., np.newaxis]
        alpha_dot_b = np.asarray(alpha_dot)[..., np.newaxis]
        cycle_frac_b = np.asarray(cycle_frac)[..., np.newaxis]

        CL, CD = self.wake.coefficients(alpha_b, cycle_frac_b, CD0=self.CD0)
        v_local = phi_dot_b * r
        q_local = 0.5 * rho * v_local ** 2

        dL_trans = CL * q_local * c * self.dr_hat * self.R
        dD_trans = CD * q_local * c * self.dr_hat * self.R

        L_trans = np.sum(dL_trans, axis=-1)
        D_trans = np.sum(dD_trans, axis=-1)

        phi_dot_arr = np.asarray(phi_dot)
        D_trans = D_trans * np.sign(phi_dot_arr + 1e-12) * -1.0

        dF_rot = rho * self.Cr * (c ** 2) * v_local * alpha_dot_b * (self.dr_hat * self.R)
        F_rot = np.sum(dF_rot, axis=-1)

        term1 = (np.pi / 16.0) * rho * (self.R ** 2) * (self.c_bar ** 2) * self.I_c2 \
                 * (phi_ddot * np.sin(alpha) + phi_dot * alpha_dot * np.cos(alpha))
        term2 = (np.pi / 64.0) * rho * self.R * (self.c_bar ** 3) * self.I_c2 * alpha_ddot
        F_added_mass = term1 - term2

        L_total = L_trans + F_rot + F_added_mass
        D_total = D_trans

        return {
            "lift": L_total,
            "drag": D_total,
            "lift_translation": L_trans,
            "drag_translation": D_trans,
            "lift_rotational": F_rot,
            "lift_added_mass": F_added_mass,
        }

class Wing:
    def __init__(self, name, root_pos, kinematics: WingKinematics, aero: WingAeroModel):
        self.name = name
        self.root_pos = np.array(root_pos, dtype=float)
        self.kin = kinematics
        self.aero = aero

    def forces(self, t):
        state = self.kin.state(t)
        cycle_frac = self.kin.reversal_phase(t)
        f = self.aero.instantaneous_forces(state, cycle_frac)
        return f

class Dragonfly:
    def __init__(self,
                 wing_span=0.045, wing_chord=0.008,
                 body_mass=0.0003,
                 flap_freq=30.0,
                 fore_amp_deg=60.0, hind_amp_deg=60.0,
                 fore_aoa_deg=45.0, hind_aoa_deg=45.0,
                 hind_phase_lag_frac=0.25,
                 rho=1.225, g=9.81,
                 rotation_axis_frac=0.25,
                 switching_gain=6.0):

        self.mass = body_mass
        self.g = g
        self.freq = flap_freq

        common_aero_kwargs = dict(
            span=wing_span, mean_chord=wing_chord, rho=rho,
            rotation_axis_frac=rotation_axis_frac, switching_gain=switching_gain,
        )

        def make_wing(name, root_pos, amp_deg, aoa_deg, phase_frac, mirror_theta=1.0):
            kin = WingKinematics(
                amplitude_deg=amp_deg, freq_hz=flap_freq,
                aoa_mid_down_deg=aoa_deg, aoa_mid_up_deg=aoa_deg,
                flip_start_frac=-0.05, flip_duration_frac=0.16,
                deviation_amp_down_deg=0.0, deviation_amp_up_deg=0.0,
                deviation_shape="oval",
                phase_offset_frac=phase_frac,
            )
            aero = WingAeroModel(**common_aero_kwargs)
            return Wing(name, root_pos, kin, aero)

        d_fore = 0.01
        d_hind = -0.012
        b_half = 0.006

        self.wings = {
            "fore_left": make_wing("fore_left", (d_fore, -b_half, 0.0),
                                    fore_amp_deg, fore_aoa_deg, phase_frac=0.0),
            "fore_right": make_wing("fore_right", (d_fore, +b_half, 0.0),
                                     fore_amp_deg, fore_aoa_deg, phase_frac=0.0),
            "hind_left": make_wing("hind_left", (d_hind, -b_half, 0.0),
                                    hind_amp_deg, hind_aoa_deg, phase_frac=hind_phase_lag_frac),
            "hind_right": make_wing("hind_right", (d_hind, +b_half, 0.0),
                                     hind_amp_deg, hind_aoa_deg, phase_frac=hind_phase_lag_frac),
        }

    def total_forces_and_moments(self, t):
        per_wing = {}
        Fz = 0.0
        roll_moment = 0.0
        pitch_moment = 0.0

        for name, wing in self.wings.items():
            f = wing.forces(np.array([t]))
            lift = float(np.squeeze(f["lift"]))
            per_wing[name] = f
            Fz += lift
            x, y, z = wing.root_pos
            roll_moment += lift * y
            pitch_moment += lift * x

        weight = self.mass * self.g
        return {
            "Fz_aero": Fz,
            "Fz_net": Fz - weight,
            "roll_moment": roll_moment,
            "pitch_moment": pitch_moment,
            "per_wing": per_wing,
        }

def simulate_hover(dragonfly: Dragonfly, t_end, dt=None, n_per_cycle=400):
    T = 1.0 / dragonfly.freq
    if dt is None:
        dt = T / n_per_cycle

    n_steps = int(t_end / dt)
    t_arr = np.zeros(n_steps)
    z_arr = np.zeros(n_steps)
    vz_arr = np.zeros(n_steps)
    roll_arr = np.zeros(n_steps)
    rollrate_arr = np.zeros(n_steps)
    Fz_history = np.zeros(n_steps)

    z, vz = 0.0, 0.0
    roll, rollrate = 0.0, 0.0
    I_roll = 1e-9

    for i in range(n_steps):
        t = i * dt
        out = dragonfly.total_forces_and_moments(t)
        az = out["Fz_net"] / dragonfly.mass
        vz += az * dt
        z += vz * dt

        roll_ang_acc = out["roll_moment"] / I_roll
        rollrate += roll_ang_acc * dt
        roll += rollrate * dt

        t_arr[i] = t
        z_arr[i] = z
        vz_arr[i] = vz
        roll_arr[i] = roll
        rollrate_arr[i] = rollrate
        Fz_history[i] = out["Fz_aero"]

    return {
        "t": t_arr, "z": z_arr, "vz": vz_arr,
        "roll": roll_arr, "rollrate": rollrate_arr,
        "Fz_aero": Fz_history,
    }

def demo_single_wing():
    rhat2 = _sanity_check_rhat2_moment()
    print(f"[Sanity check] blade-element r_hat_2^2 moment (rectangular wing): "
          f"{rhat2:.6f} (analytic value: {1.0/3.0:.6f})")

    freq = 30.0
    T = 1.0 / freq

    kin = WingKinematics(
        amplitude_deg=180.0, freq_hz=freq,
        aoa_mid_down_deg=45.0, aoa_mid_up_deg=45.0,
        flip_start_frac=-0.05, flip_duration_frac=0.16,
        deviation_amp_down_deg=0.0, deviation_amp_up_deg=0.0,
    )
    
    aero = WingAeroModel(span=0.045, mean_chord=0.008, rho=1.225,
                          rotation_axis_frac=0.25, switching_gain=6.0)

    t = np.linspace(0, 2 * T, 800)
    state = kin.state(t)
    cycle_frac = kin.reversal_phase(t)

    f = aero.instantaneous_forces(state, cycle_frac)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    axes[0].plot(t / T, np.degrees(state["phi"]), label="phi (stroke position)")
    axes[0].plot(t / T, np.degrees(state["alpha"]), label="alpha (angle of attack)")
    axes[0].set_ylabel("degrees")
    axes[0].legend()
    axes[0].set_title("Wing kinematics (2 cycles)")

    axes[1].plot(t / T, f["lift"] * 1e3, label="Total lift", color="tab:blue")
    axes[1].plot(t / T, f["lift_translation"] * 1e3, "--", label="Translation+WC only", color="tab:cyan")
    axes[1].plot(t / T, f["lift_rotational"] * 1e3, ":", label="Rotational", color="tab:orange")
    axes[1].plot(t / T, f["lift_added_mass"] * 1e3, ":", label="Added mass", color="tab:green")
    axes[1].axhline(0, color="k", linewidth=0.5)
    axes[1].set_ylabel("Force (mN)")
    axes[1].legend(fontsize=8)
    axes[1].set_title("Lift contributions by mechanism")

    axes[2].plot(t / T, f["drag"] * 1e3, label="Drag (translation+WC)", color="tab:red")
    axes[2].axhline(0, color="k", linewidth=0.5)
    axes[2].set_ylabel("Force (mN)")
    axes[2].set_xlabel("t / T (stroke cycles)")
    axes[2].legend(fontsize=8)
    axes[2].set_title("Drag")

    fig.tight_layout()
    fig.savefig(f"{OUT}/single_wing_force_traces.png", dpi=150)
    print("Saved single_wing_force_traces.png")
    
    plt.show()

def demo_hover():
    dfly = Dragonfly(
        wing_span=0.045, wing_chord=0.008,
        body_mass=0.0003, flap_freq=30.0,
        fore_amp_deg=90.0, hind_amp_deg=90.0,
        fore_aoa_deg=45.0, hind_aoa_deg=45.0,
        hind_phase_lag_frac=0.25,
    )

    weight_mN = dfly.mass * dfly.g * 1e3
    print(f"\n[Hover demo] Body weight: {weight_mN:.4f} mN")

    T = 1.0 / dfly.freq
    t_sample = np.linspace(0, T, 200)
    Fz_over_cycle = np.array([dfly.total_forces_and_moments(tt)["Fz_aero"] for tt in t_sample])
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_sample / T, Fz_over_cycle * 1e3, label="Net aero vertical force (4 wings)")
    ax.axhline(weight_mN, color="r", linestyle="--", label="Body weight")
    ax.set_xlabel("t / T")
    ax.set_ylabel("Force (mN)")
    ax.set_title("Four-wing net vertical aerodynamic force over one cycle")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUT}/hover_force_balance.png", dpi=150)
    print("Saved hover_force_balance.png")

    dfly_turn = Dragonfly(
        wing_span=0.045, wing_chord=0.008,
        body_mass=0.0003, flap_freq=30.0,
        fore_amp_deg=90.0, hind_amp_deg=90.0,
        fore_aoa_deg=45.0, hind_aoa_deg=45.0,
        hind_phase_lag_frac=0.25,
    )
    
    for side in ["fore_left", "hind_left"]:
        old = dfly_turn.wings[side]
        new_kin = WingKinematics(amplitude_deg=70.0, freq_hz=dfly_turn.freq,
                                  aoa_mid_down_deg=45.0, aoa_mid_up_deg=45.0,
                                  flip_start_frac=-0.05, flip_duration_frac=0.16,
                                  phase_offset_frac=old.kin.phase_offset / old.kin.T)
        dfly_turn.wings[side] = Wing(side, old.root_pos, new_kin, old.aero)
    
    plt.show()

def _rotate_about_axis(v, axis, angle):
    axis = axis / np.linalg.norm(axis)
    return (v * np.cos(angle)
            + np.cross(axis, v) * np.sin(angle)
            + axis * np.dot(axis, v) * (1.0 - np.cos(angle)))

def demo_animation(n_cycles=2, fps=25, frames_per_cycle=48, deviation_amp_deg=20.0,
                    deviation_shape="figure8"):
    freq = 30.0
    T = 1.0 / freq
    d_fore, d_hind, b_half = 0.01, -0.012, 0.006

    def make_wing(name, root_pos, amp_deg, aoa_deg, phase_frac):
        kin = WingKinematics(
            amplitude_deg=amp_deg, freq_hz=freq,
            aoa_mid_down_deg=aoa_deg, aoa_mid_up_deg=aoa_deg,
            flip_start_frac=-0.05, flip_duration_frac=0.16,
            deviation_amp_down_deg=deviation_amp_deg,
            deviation_amp_up_deg=-deviation_amp_deg,
            deviation_shape=deviation_shape,
            phase_offset_frac=phase_frac,
        )
        aero = WingAeroModel(span=0.045, mean_chord=0.008, rho=1.225,
                              rotation_axis_frac=0.25, switching_gain=6.0)
        return Wing(name, root_pos, kin, aero)

    wings = {
        "fore_left": make_wing("fore_left", (d_fore, -b_half, 0.0), 90.0, 45.0, 0.0),
        "fore_right": make_wing("fore_right", (d_fore, +b_half, 0.0), 90.0, 45.0, 0.0),
        "hind_left": make_wing("hind_left", (d_hind, -b_half, 0.0), 90.0, 45.0, 0.25),
        "hind_right": make_wing("hind_right", (d_hind, +b_half, 0.0), 90.0, 45.0, 0.25),
    }

    base_angle = {"fore_left": -np.pi / 2, "hind_left": -np.pi / 2,
                  "fore_right": np.pi / 2, "hind_right": np.pi / 2}
    colors = {"fore_left": "tab:blue", "fore_right": "tab:blue",
              "hind_left": "tab:orange", "hind_right": "tab:orange"}

    def tip_and_chord(wing, name, t):
        root = wing.root_pos
        phi = float(np.squeeze(wing.kin.phi(np.array([t]))))
        theta = float(np.squeeze(wing.kin.theta(np.array([t]))))
        alpha = float(np.squeeze(wing.kin.alpha(np.array([t]))))
        R = wing.aero.R

        tip_angle = base_angle[name] + phi
        horiz_r = R * np.cos(theta)
        tip = root + np.array([horiz_r * np.cos(tip_angle),
                                horiz_r * np.sin(tip_angle),
                                R * np.sin(theta)])

        span_vec = tip - root
        s_hat = span_vec / (np.linalg.norm(span_vec) + 1e-12)
        up = np.array([0.0, 0.0, 1.0])
        c0 = np.cross(s_hat, up)
        if np.linalg.norm(c0) < 1e-6:
            c0 = np.cross(s_hat, np.array([1.0, 0.0, 0.0]))
        c0 = c0 / np.linalg.norm(c0)
        chord_dir = _rotate_about_axis(c0, s_hat, alpha)
        chord_len = wing.aero.c_bar * 4.0
        chord_p1 = tip - 0.5 * chord_len * chord_dir
        chord_p2 = tip + 0.5 * chord_len * chord_dir

        f = wing.forces(np.array([t]))
        lift = abs(float(np.squeeze(f["lift"])))
        return root, tip, chord_p1, chord_p2, lift

    t_path = np.linspace(0, T, 120)
    ghost_paths = {}
    for name, wing in wings.items():
        pts = np.array([tip_and_chord(wing, name, tt)[1] for tt in t_path])
        ghost_paths[name] = pts

    n_frames = int(n_cycles * frames_per_cycle)
    t_frames = np.linspace(0, n_cycles * T, n_frames, endpoint=False)
    lift_scale = 4.0

    fig = plt.figure(figsize=(14, 7))
    axes = [fig.add_subplot(121, projection="3d"), fig.add_subplot(122, projection="3d")]
    
    axes[0].view_init(elev=22, azim=-60)
    axes[1].view_init(elev=15, azim=-90)

    axes[0].set_title("Dragonfly wingbeat (Perspective View)")
    axes[1].set_title("Dragonfly wingbeat (Side View)")

    for ax in axes:
        lim = 0.07
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x -- fore(+)/aft(-), m")
        ax.set_ylabel("y -- right(+)/left(-), m")
        ax.set_zlabel("z -- up, m")

        body_t = np.linspace(0, 2 * np.pi, 60)
        ax.plot(0.02 * np.cos(body_t), 0.006 * np.sin(body_t), np.zeros_like(body_t),
                color="gray", linewidth=1)

        for name, pts in ghost_paths.items():
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=colors[name], alpha=0.2, linewidth=1)

    span_lines = {ax: {name: ax.plot([], [], [], linewidth=2.0, color=colors[name])[0] for name in wings} for ax in axes}
    chord_lines = {ax: {name: ax.plot([], [], [], linewidth=1.2, color="black")[0] for name in wings} for ax in axes}
    time_text = axes[0].text2D(0.02, 0.95, "", transform=axes[0].transAxes)

    def update(i):
        t = t_frames[i]
        artists = []
        for name, wing in wings.items():
            root, tip, c1, c2, lift = tip_and_chord(wing, name, t)
            for ax in axes:
                span_lines[ax][name].set_data_3d([root[0], tip[0]], [root[1], tip[1]], [root[2], tip[2]])
                span_lines[ax][name].set_linewidth(1.5 + lift_scale * lift * 1e3)
                chord_lines[ax][name].set_data_3d([c1[0], c2[0]], [c1[1], c2[1]], [c1[2], c2[2]])
                artists.extend([span_lines[ax][name], chord_lines[ax][name]])
        
        time_text.set_text(f"t = {t*1e3:.1f} ms  (t/T = {t/T:.2f})")
        artists.append(time_text)
        return artists

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False, interval=1000 / fps)
    out_path = f"{OUT}/dragonfly_wingbeat_3d_{deviation_shape}.gif"
    anim.save(out_path, writer=PillowWriter(fps=fps))
    
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    demo_single_wing()
    demo_hover()
    demo_animation()