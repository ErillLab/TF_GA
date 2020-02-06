## Branches related
- master: You should NOT develop on this branch.
- develop: This is where the most recent code should be. This is no necessary tested and could have bugs.

## Code related

Please, use explivative variable names.

I like to use prefixes to specify complex variable types, so a new programmer knows the data-type of a variable.

Example:
Use dName for a dictionary or aName for an array 
```python
dMyDictionary = {}
aMyArray = []
```

And do not use magic numbers so new programmers have a better understanding of the code:
```python
#WRONG
for i in range(10):
  print("Loop: {}".format(i))


#RIGHT
numberOfLoops = 10
for i in range(numberOfLoops):
  print("Loop: {}".format(i))

```


