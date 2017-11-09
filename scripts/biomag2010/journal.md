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

Thu Nov  9 19:08:19 MSK 2017
============================
I was  bothered by the fact that reducing the number of trials still showed very
reasonable networks.
In order to study how the  number of trials and the level of noise affect performance
I first looked at performance on completely random CT. Both real and imag parts gave
very reasonably looking bundles of connections on 15 000 vertices cortex and a bit worse on
2000. 

Then I ran bootstrap to see if the networks I got are real. First I computed 100 bootstrap
cross-spectrum timeseries on low resolution brain (2000v) with GainSVDTh = 0.001 resulting
into 33 artificial sensors.
is_induced_only parameter was set to true.
Time range was from 0 to 1 s and frequency band was from 16 to 24 Hz.
On each bootstrap iteration I projected CT from signal leakage with pwr_rnk=150,
computed svd and took first left singular vector. Then I did PSIICOS_ScanFast on it and
took 30 most strong connections which I put into Connections container.
Then I averaged them whithin each bootstrap iteration and plotted the result.

Here's what I got:
![image](./pics/bootstrap_LR_all_tr.png)


Then I reduced the number of trials to 30:
![image](./pics/bootstrap_LR_30_tr.png)

And then to 3:
![image](./pics/bootstrap_LR_3_tr.png)
