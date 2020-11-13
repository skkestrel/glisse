# GLISSE (GPU Long-term Integrator for Solar System Evolution)

This package provides an integrator designed to simulate large numbers of massive particles over large timescales.
This code works best with systems with small numbers of massive bodies (e.g. in solar systems with <10 planets),
and currently does not support close encounter handling between planet-planet or planet-particle interactions.

## Compilation

Simply type make in the project root directory to compile the integrator, or make cpu to compile a CPU-only version.
The makefile may need to be edited to target your specfic GPU architecture. The current makefile has been tested with the 1080ti and the Tesla K20.

## Usage

There are two required files by GLISSE: one, a configuration file; and two, an initial state file.
The configuration file should contain one entry per line. Lines beginning with a pound sign (#) are ignored.
Below shows a list of configurable options of the integrator.

The units used by the integrator are in natural units. By default, G is set to 1. Thus, to use a system with
days and au, the solar mass should be set to 4 * pi^2 / (365.25)^2 = 2.959 x 10^-4. Days and au are the default unit system.

| Name | Description | Default |
| --- | --- | --- |
| Initial-Time | The starting time of the integration. | 0 |
| Time-Step | The timestep to use in the integration. We recommend this should be at most 180 days (converted to the appropriate units) for outer solar system dynamics. | 122 |
| Final-Time | The time to stop the integration. | |
| Time-Block-Size | The timeblock size; the planetary chunk size. The number of timesteps that the GPU will advance in one kernel launch. | 1024 |
| Cull-Radius | Particles are deactivated if they come within this radius of any planet, in natural units. | 0.5 | BG: ALTER TO A MULTIPLE OF PLANETARY HILL SPHERES; COMPUTE ON CPU AT START AND PASS IN.
| CPU-Thread-Count | The number of threads to use in CPU-only mode. | 4 |
| Log-Interval | The integrator will print the current progress every Log-Interval number of timeblocks. 0 to disable. | 10 |
| Status-Interval | The integrator will write the integration status to the file named `status` in the project output directory every Status-Interval number of timeblocks. 0 to disable. See below. | 1 |
| Track-Interval | The integrator will write orbital elements to the integration track every Track-Interval number of timeblocks. See below. 0 to disable. | 0 |
| Resync-Interval | The integrator will sort ("defragment") the GPU particle array every Resync-Interval. This parameter should be increased when Time-Block-Size is small for performance. | 1 |
| Write-Barycentric-Track | The integrator will write barycentric instead of heliocentric orbital elements to the particle tracks if enabled. | 0 |
| Split-Track-File | If zero, the integrator will write particle tracks into a single file named `track' in the output directory. If nonzero, the integrator will write particle tracks to files with a maximum size of Split-Track-File in bytes, named sequentially in a folder named `tracks' in the output directory. | 0 |
| Dump-Interval | The integrator will dump particle and planet states to a folder named `dumps' in the output directory every Dump-Interval number of timeblocks. 0 to disable. | 1000 |
| Keep-All-Dumps | Whether to keep all the dump files. If set to 0, new dump file will cover the old one. | 0 |
| Write-Binary-Output | Whether to write the output state file in binary format. | 0 |
| Read-Binary-Input | Whether to write the input state file in binary format. | 0 |
| Input-File | The absolute path of the input state file to read. | |
| Output-Folder | The absolute path of the output folder. | |
| Read-Input-Momenta | Whether to interpret momenta instead of velocities in the input state file. | 0 |
| Write-Output-Momenta | Whether to write momenta instead of velocities in the output state file. | 0 |

File formats
Input and output states
Planet count
For each planet: 4 lines
	mass
	x y z
	vx vy vz
	id
Particle count
For each particle: 3 lines
	x y z
	vx vy vz
	id deathflags deathtime

### Particle tracks
The particle track is always in binary format and contains a history of particle and planet orbital elements in single-precision.
TODO

### Utility executables
bin/make-state Generate an initial state file from a template planet data file and uniformly sampling orbital elements for particles
bin/convert-state Convert states from different formats, or between different coordinate systems.
For example: bin/convert-state read state.in to-bary write state.bary.in
bin/filter-state Find particles in a state file that satisfy certain criteria, for example, to find all particles with semimajor axis greater than 20 au
bin/track-info Display information about a particle track

### Utility scripts
scripts/plot_history.py provides utilities to plot data form a particle track
