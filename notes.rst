Experiments
===========

1. Basic chemotaxis, recreating something comparable to Martin's self-propelled droplets.




TODOs
=====

Incorporate D\theta into diffusion and convection in individual.py so that when
the number of segments change, this has limited influence on what the simulation
is modelling.


Model design concerns...
========================

1. Currently chemistry (reactions and present molecules) is synchronized within
/but not between/ individuals. Arguably, if it is a single world I am simulating, I may
want to synchronize chemistry across all individuals. This would allow me to
compare, for instance, how /utilized/ chemistry is more diverse (I expect) when
behaviour is controlled by agents and not just a random walk.

But then, when do avalanches occur? Individuals would only be different because
of chemical concentrations.

   I have decided to go with a single uniform chemistry throughout the
   simulation. I can simulate a single agent if I want to, and then do a
   comparison of that agent to a random-motion agent. I can do that many times
   (sample over many AC). It just doesn't make any logical sense to have
   multiple concurrent laws-of-chemistry going on... It's cleaner this way, but
   may make me use a different AC than Fernando and Rowe's (perhaps).

2. Could use recorded behavioural agent's path to control non-behavioural agent
(it will have precisely the same statistical properties, unlike a random walk)
but will not be driven by feedback of the agent, so will not be smart.

