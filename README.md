# Dragonfly-flight-Simulation
Description: Simulates four-wing dragonfly hover mechanics, calculating lift, drag, and added mass with synchronized 3D Matplotlib animations.


Dragonfly Quasi-Steady Flight Simulation -- v1 

A basic (non-RL) blade-element quasi-steady aerodynamic simulation of a four-winged flapping flyer, implementing the mechanisms discussed across:

Sane & Dickinson (2001), J. Exp. Biol. 204 -- translational force coefficients, kinematic parameterization (phi, alpha, theta), added mass.
Sane & Dickinson (2002), J. Exp. Biol. 205 -- rotational circulation (theoretical Cr = pi*(0.75 - x0) form used here, consistent with their measured ~1.5-1.6 at quarter-chord rotation axis).
Nabawy (2023), J. R. Soc. Interface 20 -- wake capture, modelled as a transient boost to the 3D lift-curve slope via a switching function.

This version merges the original four modules (kinematics.py, aerodynamics.py, dragonfly.py, simulate.py) into a single file, dragonfly_sim.py, and adds a 3D wingbeat animation that didn't exist before. See "What changed from v0" below for the full list.

File
dragonfly_sim.py -- everything: kinematics, aerodynamics, the 4-wing assembly, the toy rigid-body integrator, and three demos. Organized as four clearly-labeled sections in dependency order (each section is the original module's content, unchanged except where noted below):
kinematics -- WingKinematics: generates phi(t), alpha(t), theta(t) for one wing from six physically-meaningful parameters. phi_dot/phi_ddot are analytic (_SmoothedTriangleGenerator) to avoid numerical-differentiation artifacts near stroke reversal.
aerodynamics -- WakeCaptureLiftModel, WingAeroModel: per-wing blade-element force model (translation + wake capture, rotational circulation, added mass).
dragonfly -- Wing, Dragonfly, simulate_hover: assembles 4 wings (fore-left/right, hind-left/right) with independent phase lag, plus a simple rigid-body integrator.
demos -- demo_single_wing(), demo_hover(), demo_animation(): the two original validation demos, plus a new 3D wingbeat animation.

Run: python3 dragonfly_sim.py Outputs go to /mnt/user-data/outputs/:

single_wing_force_traces.png
hover_force_balance.png
dragonfly_wingbeat_3d_figure8.gif (or ..._oval.gif, depending on deviation_shape)
What changed from v0 (the original four-file version)

Nothing about the physics, signs, or default parameters was changed -- this was a straight merge plus two small cleanups plus one new feature:

Merged into one file. Removed the from kinematics import ... / from aerodynamics import ... / from dragonfly import ... lines; everything lives in dragonfly_sim.py now, in the same dependency order as before.
Dead code removed. WingKinematics.alpha() had a placeholder line that was computed and immediately overwritten by the next line -- removed. No behavior change.
A stale doc claim was actually verified. The original aerodynamics.py docstring claimed the blade-element sum "will converge to r_hat_2^2 = 1/3 as a sanity check" for a rectangular wing, but nothing in the code ever checked that. Added _sanity_check_rhat2_moment(), run and printed at the start of demo_single_wing() -- it does converge to 1/3 (confirmed: printed value is 0.333333).
New: demo_animation(). Did not exist before -- the original code only ever plotted the resulting forces, never the wing motion. This demo draws the actual kinematics in 3D:
phi (sweep) -- the tip arcs side to side about its root.
theta (deviation) -- lifts the tip out of that arc; with deviation_shape="figure8" (the default) this produces a genuine self-crossing figure-8 tip path -- each half-stroke dips one way then swings back the other way and crosses itself at the half-stroke midpoint (verified directly by tracing theta(t): e.g. 0° → -19° → 0° → +19° → 0° within a single half-stroke). Pass deviation_shape="oval" instead for the simpler single-loop-per-cycle path (one hump per half-stroke, no self-crossing).
alpha (pitch) -- drawn as a short chord segment at the tip, rotated about the spanwise axis via Rodrigues' rotation formula. Watch it snap through its flip at each stroke reversal.
Wing linewidth scales with instantaneous |lift| (from the same wing.forces(t) used elsewhere), so you can see which wing is loaded hardest at each instant.
Important: this demo builds its own set of wings with a nonzero stroke-deviation amplitude (deviation_amp_deg, default 20°) purely so the flapping and tip path are visible. demo_single_wing() and demo_hover() are untouched and still use deviation_amp = 0 exactly as in the original -- the force/validation numbers printed by those two functions are unaffected.
Known simplifications (read before trusting numbers) -- unchanged from v0
Rectangular wing planform (uniform chord along span) -- real dragonfly wings taper. Swap c_hat(r_hat) in WingAeroModel.__init__ for a real chord distribution when you have one.
Spanwise-uniform angle of attack -- no wing twist modeled.
Rotational circulation and added-mass forces are treated as purely "lift-like" (normal to stroke plane) for simplicity. Only the translation+wake-capture forces are properly split into lift/drag using the Sane & Dickinson geometric convention. Refining this means resolving each force along the wing-chord-normal direction and re-projecting into lift/drag/radial per their "Data analysis" section.
No forward-flight / free-stream velocity term -- hovering only. To extend: vector-add a free-stream velocity to the local blade-element velocity before computing angle of attack and dynamic pressure.
Rigid-body dynamics are a toy: decoupled scalar vertical position
roll angle, explicit Euler integration, a placeholder moment of inertia. Not a real 6-DOF treatment (no proper rotation matrices, no coupling between axes). Good enough to sanity-check "does it hover" and "does an asymmetry produce a turning moment in the right direction" -- not good enough to trust an actual trajectory from.
Wing geometry (span, chord, mass, frequency) are representative guesses, not measured dragonfly values. Swap in real morphometric data when you have it (Wakeling & Ellington's dragonfly papers are the natural source).
Four-wing aerodynamic interference (downwash between fore/hind pairs) is NOT modeled -- each wing's force is computed as if the other three don't exist. This is the Hu & Deng-type correction you'll want to add next; right now all 4 wings are aerodynamically independent, just summed.
The 3D animation is a drawing convenience, not new physics. The stroke-plane-horizontal / theta-is-vertical convention used to place the wing tip in 3D space for plotting is a simplification for visualization only -- it does not feed back into or change any force calculation.
Validated against the papers so far
Single-wing force trace (single_wing_force_traces.png) reproduces the qualitative pattern described in Sane & Dickinson (2001) Results: near-constant lift through mid-stroke, a dip approaching reversal, a positive "wake capture" overshoot immediately after reversal, settling back to steady state -- and drag cleanly flipping sign at each reversal with a transient overshoot in the new direction.
The rectangular-wing second moment of area (r_hat_2^2) numerically converges to the analytic value of 1/3, confirmed at runtime (new in this version -- see point 3 above).
NOT yet validated: absolute force magnitudes / CL, CD values against the papers' measured numbers (different species, different scale -- would need to rerun with fruit-fly-scale wing dimensions and Reynolds-matched fluid properties for a real apples-to-apples check).
Suggested next steps, in order
Re-run demo_single_wing() at fruit-fly scale (25 cm model wing, mineral oil density/viscosity per Sane & Dickinson's actual setup) and check whether computed CL, CD land near their reported peak values (CL ~1.8-2.0, L/D max ~0.8 at alpha=30°/Phi=180°). This is the single most valuable validation step before trusting dragonfly-scale numbers.
Add real dragonfly wing morphometrics (span, chord distribution, mass, wingbeat frequency) from Wakeling & Ellington instead of the current placeholder values.
Add fore/hind wake interference correction (Hu & Deng-style downwash coupling) -- currently the biggest missing "dragonfly-specific" piece.
Replace the toy rigid-body integrator with a proper 6-DOF treatment (quaternion or rotation-matrix based) if you want trustworthy flight trajectories rather than just force/moment sanity checks.
Only after 1-4: consider forward-flight extension, RL action space, etc.
If you extend demo_animation() further: overlaying the actual force vector (lift direction/magnitude) at each wingtip, or animating the simulate_hover() body trajectory itself (currently only the wings are animated, not the body's vertical/roll response), would be natural next visual additions.
