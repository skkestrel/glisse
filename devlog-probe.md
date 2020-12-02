## 2020-11-12
* Added a new parameter "Keep-All-Dumps" to choose whether or not to keep all dump files;
* Set the final time in config.out to be a constant;
* Removed "Output-File".

## 2020-11-13
* Eliminated "CPU-Thread-Count", "Read-Input-Momenta", and "Write-Input-Momenta";
* Sorted configurations by a logical order.
* Output "Particle-Input-File" and "Planet-Input-File" only if they are not empty.

## 2020-11-16
* Replaced "Cull-Raidius" with "Hill-Radius-Factor" (constants that won't update);
* Added "Particle-Inner-Boundary" and "Particle-Outer-Boundary" that defines the inner and outer boundaries for ejecting test particles. (won't affect planets, planets won't be discarded at all in GLISSE?)

Note: It seems that the particle discarded time has been implemented into GPU already.

## TO DO
* Terminate the integrator whenever planets cross.


## LONG-TERM TO DO
* Change discard logic.

## DONE
* Change Cull-Radius documentation to indicate this is in units of Hill radii
* Introduce Cull radius code to compute this at code start and pass to GPU
* Eliminate CPU-Thread-Count ?  (from docs and make sure appears nowhere in code). DONE.
* Eliminate Limit-Particle-Count ?  (Find out what this does first !)
This parameter appears in the main calculation routine (sia4), and can be used to cap the particle number. NO CHANGE WAS MADE.
* What is the Energy-Interval ?
It's the energy output frequency. Output stored in time.out. NO CHANGE WAS MADE.
* Eliminate option of input and output of momenta (rather than velocities) DONE.
* rename Energy-Interval. DONE!
* Terminate if run out of test particles.