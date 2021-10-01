## POSDEF
This script makes a matrix positive definite. To run you will need to install all requirements in requirements.txt:

```
pip install -r requirements.txt
```

## Notes
The original code had a bug so it always ran 1000 iterations no mater the input for maxiters.
This version does not have that bug and so when incoporating into a new codebase make sure
to update all calls to this function with the appropriate max iters. 
