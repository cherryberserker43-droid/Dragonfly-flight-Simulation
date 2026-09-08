# Dragonfly-flight-Simulation
Description: Simulates four-wing dragonfly hover mechanics, calculating lift, drag, and added mass with synchronized 3D Matplotlib animations.

A blade-element quasi-steady aerodynamic simulation of a four-winged flapping insect flyer (Dragonfly) in Python. Implements key physical mechanisms for hovering flight, including translational forces, rotational circulation, added mass, and wake capture, along with 3D synchronized visualization.

 Theoretical Foundations
 This simulation implements mechanisms from foundational insect aerodynamics literature:Sane & Dickinson (2001, J. Exp. Biol. 204): Translational force coefficients, kinematic parameterization ($\phi$, $\alpha$, $\theta$), and added mass forces.Sane & Dickinson (2002, J. Exp. Biol. 205): Rotational circulation using theoretical $C_r = \pi (0.75 - x_0)$ form for a quarter-chord rotation axis.Nabawy (2023, J. R. Soc. Interface 20): Wake capture modeled as a transient boost to the 3D lift-curve slope via a switching function.

Code Architecture (dragonfly_sim.py
)All dependencies, physics, and demonstrations live inside a single executable script structured in dependency order:kinematics (WingKinematics): Generates sweep ($\phi$), pitch ($\alpha$), and stroke deviation ($\theta$) from 6 physical parameters. Uses analytic derivatives (_SmoothedTriangleGenerator) for $\dot{\phi}$ and $\ddot{\phi}$ to prevent numerical artifacts at stroke reversal.aerodynamics (WakeCaptureLiftModel, WingAeroModel): Per-wing blade-element force calculations covering translation, wake capture, rotational circulation, and added mass.dragonfly (Wing, Dragonfly, simulate_hover): Assembles a 4-wing vehicle (fore-left/right, hind-left/right) with independent phase lags and a toy 2-DOF rigid-body integrator.demos:demo_single_wing(): Validates single-wing aerodynamic force profiles.demo_hover(): Evaluates 4-wing vertical force balance and roll moment asymmetries.demo_animation(): Generates synchronized perspective & side-view 3D animations displaying wing sweep, pitch flips ($\alpha$), figure-8 tip trajectories ($\theta$), and dynamic lift-based line weighting.

Outputs:
single_wing_force_traces.png — Lift & drag breakdown over 2 stroke cycles.
hover_force_balance.png — Net vertical force vs. dragonfly body weight.
dragonfly_wingbeat_3d_figure8.gif — Animated 3D side-by-side visualization of wing kinematics.

Simplifications
Geometry: Uniform rectangular wings (no taper or spanwise twist).
Force Projections: Rotational circulation and added mass are modeled normal to the stroke plane rather than fully resolved along chord-normal coordinates.
Flight Regime: Pure hover (no forward flight / free-stream velocity vector).
Rigid-Body Model: Toy decoupled Euler integrator for $z$ (altitude) and $\phi$ (roll).
Interference: No fore/hind wing aerodynamical downwash interaction modeled (wings are treated independently)

Validations
Force Traces: Single-wing force output matches qualitative trends in Sane & Dickinson (2001) — steady mid-stroke lift, reversal dips, and wake capture spikes.
Geometric Moments: Radial second moment of area numerical integration verified against analytic $\frac{1}{3}$.

Next Steps
Fruit-Fly Scale Validation: Re-run at 25 cm model scale in mineral oil to verify absolute force magnitudes (C 
L ∼1.8–2.0).
Morphometrics: Swap placeholder geometry with measured dragonfly data (e.g., Wakeling & Ellington).
Aerodynamic Interference: Add Hu & Deng-style fore/hind wing downwash interaction.
Flight Dynamics: Upgrade toy integrator to full 6-DOF quaternion rigid-body dynamics.

