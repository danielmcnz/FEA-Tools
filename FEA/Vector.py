# this is used so a class can be used as a type within its own class (i.e. using the Vec2 type within the Vec2 class)
from __future__ import annotations
import math

class Vec2:
    def __init__(self, x, y):
        self.x : float = x
        self.y : float = y

        self.length : float


    def __repr__(self):
        print("Vector 2 -> " + str(self.x) + " : " + str(self.y))

    
    def __str__(self):
        return "(" + str(self.x) + " : " + str(self.y) + ")"


    def __eq__(self, other : Vec2):
        return self.x == other.x and self.y == other.y


    def __bool__(self, other : Vec2):
        return self.__eq__(self, other)


    def __add__(self, other : Vec2):
        return Vec2(self.x + other.x, self.y * other.y)
    

    def __sub__(self, other : Vec2):
        return Vec2(self.x + other.x, self.y * other.y)


    def __mul__(self, other : Vec2):
        return Vec2(self.x * other.x, self.y * other.y)


    def __div__(self, other : Vec2):
        if(other.x or other.y == 0):
            raise ZeroDivisionError("Cannot divide by zero")
        return Vec2(self.x / other.x, self.y / other.y)
    

    def normalize(self):
        self.length = math.sqrt(self.x ** 2 + self.y ** 2)
        self.x /= self.length
        self.y /= self.length