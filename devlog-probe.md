## 2020-11-12
* Added a new parameter "Keep-All-Dumps" to choose whether or not to keep all dump files.
* Set the final time in config.out to be a constant.

## TO DO
* Change Cull-Radius documentation to indicate this is in units of Hill radii
* Introduce Cull radius code to compute this at code start and pass to GPU
* Eliminate CPU-Thread-Count ?  (from docs and make sure appears nowhere in code).
* Eliminate Limit-Particle-Count ?  (Find out what this does first !)
* What is the Status-Interval ?
* Eliminate option of input and output of momenta (rather than velocities)
