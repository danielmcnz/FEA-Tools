# this is used so a class can be used as a type within its own class (i.e. using the Vec2 type within the Vec2 class)
from __future__ import annotations
import math

class Vec2:
    def __init__(self, x, y):
        self.x : float = float(x)
        self.y : float = float(y)

        self.length : float = math.sqrt(self.x ** 2 + self.y ** 2)


    def __repr__(self):
        print("Vector 2 -> " + str(self.x) + " : " + str(self.y))

    
    def __str__(self):
        return "(" + str(self.x) + " : " + str(self.y) + ")"


    def __eq__(self, other : Vec2):
        return self.x == other.x and self.y == other.y


    def __bool__(self, other : Vec2):
        return self.__eq__(self, other)


    def __add__(self, other : Vec2):
        return Vec2(self.x + other.x, self.y + other.y)
    

    def __sub__(self, other : Vec2):
        return Vec2(self.x - other.x, self.y - other.y)


    def __mul__(self, other):

        if type(other) == Vec2:
            return Vec2(self.x * other.x, self.y * other.y)
        elif type(other) == int or type(other) == float:
            return Vec2(other * self.x, other * self.y)
        

    def __ge__(self, other : Vec2):
        return Vec2(self.x >= other.x, self.y >= other.y)
    
    
    def __gt__(self, other : Vec2):
        return Vec2(self.x > other.x, self.y > other.y)
    

    def __le__(self, other : Vec2):
        return Vec2(self.x <= other.x, self.y <= other.y)
    

    def __lt__(self, other : Vec2):
        return Vec2(self.x < other.x, self.y < other.y)


    def __truediv__(self, other):
        if type(other) == Vec2:
            if (self.x != 0 and other.x == 0) or (self.y != 0 and other.y == 0):
                raise ZeroDivisionError("Cannot divide by zero")
            
            if other.x == 0:
                x = 0
            else:
                x = self.x / other.x
            if other.y == 0:
                y = 0
            else:
                y = self.y / other.y
            return Vec2(x, y)
        
        elif type(other) == int or type(other) == float:
            if other == 0 and self.x != 0:
                raise ZeroDivisionError("Cannot divide by zero")
            return Vec2(self.x / other, self.y / other)
    

    def abs(self):
        return Vec2(abs(self.x), abs(self.y))
    

    def normalize(self) -> Vec2:

        if(self.x == 0 and self.y == 0):
            length = 0
            x = 0
            y = 0
        else:
            length = math.sqrt(self.x ** 2 + self.y ** 2)
            x = self.x / length
            y = self.y / length

        vec = Vec2(x, y)
        vec.length = length

        return vec


    def perpendicular(self):
        return Vec2(-self.y, self.x)