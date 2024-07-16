This Repository contains code and notes concerning an implementation of the Discrete Dipole Approximation in Julia. 
So far, it has three methods, neither of them fully optimized yet:
1) Direct Inversion (with enormous memory usage),
2) Fixed Point Iteration (using fields to calculate polarizations and vice-versa), following Purcell and Pennypacker, which has yielded the majority of the results so far.
3) Conjugated Gradients, which takes longer but solves cases for which the fixed point method exploded.

Some notes with teoretical derivations of the DDA are provided.

