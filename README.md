## POSDEF
This script makes a matrix positive definite. To run you will need to first install the requirements:

```
pip install -r requirements.txt
```

## Notes
The original code had a bug so it always ran 1000 iterations no matter the input for max iterations.
This version does not have that bug and so when incoporating into a new codebase make sure
to update all calls to this function with the appropriate max iterations. 
