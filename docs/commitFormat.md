To establish and "easy to read" commit message when searching for logs I would like to standarize the format of the commits.
```
[CHANGE-TYPE] scope: Little description of the change
```
## CHANGE-TYPE: This is about the reason of the commit.

- FEAT for new features on the code.
- FIX for for fixes on the code.


## scope: Name of the file without extension.

If I edited the file pObject.py the scope would be pObject.
To commit a change done on the project structure, use the scope `project`.
To commit a change on different files, add all the files to the scope.

## Little description of the change.

Please use the imperative form and ensure that it's not too long to read.

Examples:
```
[FEAT] pObject: Implement the mutation function
[FIX] pObject, cObject: Remove the variable that stores the fit
```
