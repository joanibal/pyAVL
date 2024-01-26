# Optimization 

See [`examples/aircraft/run_opt.py`](https://github.com/joanibal/pyAVL/blob/master/examples/aircraft/run_opt.py) for an example. 

More documentation to follow. 
For now if you have a question post an issue on github.  


## Debugging
-  Some variables (like chord, dihedral, x and z leading edge position) can lead to local minimum. 
   To help fix this add a constraint that keeps the variable monotonic 
