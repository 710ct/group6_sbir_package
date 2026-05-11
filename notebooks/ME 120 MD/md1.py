import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import os
import io
import numpy as np
import warnings
import MDAnalysis as mda
import py3Dmol

def simulate_apple_fall(total_time: float = 10, mass: float = 0.3, initial_velocity: float = -0.1, height: float = 553, timestep: float = 0.05):
    """
    Parameters
    ----------
    total_time: float

    mass: float

    initial_velocity: float

    height: float

    timestep: float
    """
    # Setup the figure and axes...
    fig, ax = plt.subplots(figsize=(8,8))

    ## Adjust axes limits according to the problem. Here we don't need more than a couple of meters left or right, and 600 meters up
    ax.set(xlim=(-2, 2), ylim=(0, height + 50), xlabel='Position, meters', ylabel='Height, meters', title='Apple falling from CN tower')

    gravity = 9.8

    # setting a timestep to be 50 ms
    # dt = 0.05 #s
    num_timesteps = int(total_time // timestep)

    # Allocating arrays for 2D problem
    v = np.zeros((num_timesteps + 1, 2))
    r = np.zeros((num_timesteps + 1, 2))
    f = np.zeros((num_timesteps + 1, 2))

    # initial conditions:
    r[0] = np.array([0., height])
    v[0] = np.array([initial_velocity, 0.])

    # the only force is gravity
    f[:] = np.array([0., mass * -gravity])

    ## Run dynamics:
    for n in range(num_timesteps):
        v[n+1] = v[n] + f[n] / mass * timestep
        r[n+1] = r[n] + v[n+1] * timestep

    ## drawing the first data point  
    scat = ax.scatter(r[0,0], r[0,1], marker='o', c='g', s=200)

    ## animating
    def animate(i):
        scat.set_offsets(r[i])

    ani = animation.FuncAnimation(fig, func=animate, frames=num_timesteps)
    plt.close()

    # Check if the results directory exists
    if not os.path.exists('results'):
        # Create the results directory
        os.makedirs('results')
    ani.save('results/apple_fall.mp4', fps=1//timestep, writer='ffmpeg')

    # # Use if ffmpeg is not installed
    # ani.save('results/apple_fall.gif', fps=1//timestep)

def simulate_three_particles(total_time: float = 10, mass: float = 1.0, ks: int = 5, r0: float = 1.0, timestep: float = 0.05):
    # Setup the figure and axes...
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set(xlim=(-3.5, 3.5), ylim=(-3.5, 3.5), ylabel='meters', xlabel='meters', title='3-Body problem')

    # parameters of the problem
    # T = 10. #s
    # m = 1.0 #kg
    # ks = 5 #N/m
    # r0 = 1. #m

    # setting a timestep to be 50 ms
    # dt = 0.05 #s
    num_timesteps = int(total_time / timestep)

    # Allocating arrays for 2D problem: first axis - time. second axis - particle's number. third - coordinate
    v = np.zeros((num_timesteps + 1, 3, 2))
    r = np.zeros((num_timesteps + 1, 3, 2))
    f = np.zeros((num_timesteps + 1, 3, 2))

    # initial conditions for 3 particles:
    r[0,0] = np.array([0., 2.])
    r[0,1] = np.array([2., 0.])
    r[0,2] = np.array([-1., 0.])

    def compute_forces(n):
        '''The function computes forces on each particle at time step n'''
        for i in range(3):
            for j in range(3):
                if i != j:
                    rij = r[n, i] - r[n, j]
                    rij_abs = np.linalg.norm(rij)
                    f[n, i] -= ks * (rij_abs - r0) * rij / rij_abs 
    
    ## Run dynamics:
    for n in range(num_timesteps):
        compute_forces(n)
        v[n+1] = v[n] + f[n] / mass * timestep
        r[n+1] = r[n] + v[n+1] * timestep

    ## drawing and animating 
    scat = ax.scatter(r[0,:,0], r[0,:,1], marker='o', c=['b', 'k', 'r'], s=1000)

    def animate(i):
        scat.set_offsets(r[i])

    ani = animation.FuncAnimation(fig, animate, frames=num_timesteps)
    plt.close()
    # Check if the results directory exists
    if not os.path.exists('results'):
        # Create the results directory
        os.makedirs('results')
    ## this function will create a lot of *.png files in a folder '3Body_frames'
    ani.save('results/3particles.mp4', fps=1 // timestep)

    # # Use if ffmpeg is not installed
    # ani.save('results/3particles.gif', fps=1 // timestep)
    
    


class KeepOpenStringIO(io.StringIO):
    def close(self):
        pass

def visualize(
    u,
    step=10,
    max_frames=200,
    # protein look
    protein_color="spectrum",
    cartoon_style="trace",      # "trace" is the most robust
    cartoon_thickness=0.35,
    # water hint (fast + obvious)
    water_hint=True,
    water_atoms=1500,           # sampled once; constant across frames
    water_opacity=0.55,
    water_scale=0.65,
    water_color="deepskyblue",
    seed=0,
    # view options
    width=900,
    height=550,
    interval=10
):
    """
    Protein as cartoon; water as a sampled cloud of water oxygens.
    Works for trajectories (subsampled by step and capped by max_frames).
    """
    view = py3Dmol.view(width=width, height=height)

    # --- protein selection ---
    ag_prot = u.select_atoms("protein")
    if len(ag_prot) == 0:
        ag_prot = u.atoms  # fallback

    # --- water hint: oxygen atoms only, sampled once ---
    ag_water = None
    if water_hint:
        ag_water_all = u.select_atoms("resname HOH WAT SOL and name O*")
        if len(ag_water_all) > 0:
            if len(ag_water_all) > water_atoms:
                rng = np.random.default_rng(seed)
                idx = rng.choice(len(ag_water_all), size=water_atoms, replace=False)
                ag_water = ag_water_all[idx]
            else:
                ag_water = ag_water_all

    ag_show = ag_prot if ag_water is None else (ag_prot + ag_water)

    pdb_frames = []
    frame_idx = 0

    # Hide the super-noisy MDAnalysis PDB warnings (harmless here)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module=r"MDAnalysis\.coordinates\.PDB")

        for ts_i, ts in enumerate(u.trajectory[::step]):
            if ts_i >= max_frames:
                break

            sio = KeepOpenStringIO()
            with mda.Writer(sio, format="PDB", multiframe=False) as W:
                W.write(ag_show)

            frame_text = sio.getvalue()
            atom_lines = [l for l in frame_text.splitlines() if l.startswith(("ATOM", "HETATM"))]

            frame_idx += 1
            pdb_frames.append(f"MODEL     {frame_idx:4d}")
            pdb_frames.extend(atom_lines)
            pdb_frames.append("ENDMDL")

    pdb_data = "\n".join(pdb_frames) + "\n"
    view.addModelsAsFrames(pdb_data, "pdb")

    # --- styles ---
    # Key fix: don't rely on {"protein": True}; instead cartoon everything that's NOT water.
    water_sel = {"or": [{"resn": "HOH"}, {"resn": "WAT"}, {"resn": "SOL"}]}
    nonwater_sel = {"not": water_sel}

    view.setStyle(
        nonwater_sel,
        {"cartoon": {"color": protein_color, "style": cartoon_style, "thickness": cartoon_thickness}}
    )

    if water_hint and ag_water is not None:
        view.addStyle(
            water_sel,
            {"sphere": {"scale": water_scale, "opacity": water_opacity, "color": water_color}}
        )

    view.zoomTo()
    if frame_idx > 1:
        view.animate({"loop": "forward", "interval": interval})

    return view