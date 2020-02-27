# Node object
import random
DEFAULT_LENGTH = 5

class Node:
    
    def __init__(self):
        self.ID = random.randint(0,100)

    def doNothing(self):
        print("Doing nothing... from Node Object.")

