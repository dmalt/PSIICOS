Sun Nov  5 19:32:19 MSK 2017
============================

I've reproduced Alex's results with components.
Two problem appeared.
1) I've ran the analyses for the real part of cross-spectrum
   but the phase visualization showed that phase lag for the second
   component is very close to pi/2. Then I checked the code for phase computation
   and it appeared that I was projecting only on the first dipole of each node of a
   network instead of computing orientation.

   This strange phase computation was being done with the projected CT.
   When I run it with the unprojected CT, phases seem normal.
2) The pictures I get with the real part are very nice but they don't get worse
   when I reduce the number of trials. For instance I've ran calculations with only 
   5 trials and got beautiful cross-lateral networks, very similar to what I've seen
   for the complete set of trials.

