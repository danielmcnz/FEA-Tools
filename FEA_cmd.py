from FEA.FEA import *

import re
from typing import List

PROMPT = ">> "

HELP = """    
    -----------------------
    get info from an element
    ------------------------
    syntax -> element([index]).[variable]
    e.g. -> element(1).local_stiffness
    
    To check available info from an element, enter element.help.

    ------------------------
    get info from a structure
    ------------------------
    syntax -> structure.[variable]
    e.g. -> structure.total_stiffness
    e.g. -> structure.elements

    To check available info from a structure, enter structure.help.

    ------------------------
    display
    ------------------------
    syntax -> display.[element or structure]
    e.g. display.element(0)
    e.g. display.structure

    Will display the deflected and original structure / element.
"""


def get_n_elements() -> int:

    n_elements = input("number of elements? ")

    try:
        n_elements = int(n_elements)
    except ValueError:
        print("please enter a number")

    return n_elements

def process_matrix(assembly_matrix : str) -> np.ndarray:
    elements = assembly_matrix.split(',')

    np_mat : np.ndarray

    mat : list = []
    cur_row : list = []

    depth : int = 0
    new_row = False

    for element in elements:
        while(element.startswith('[')):
            depth += 1
            element = element.removeprefix('[')
        while(element.endswith(']')):
            depth -= 1
            new_row = True
            element = element.removesuffix(']')
            
        cur_row.append(int(element))

        if(new_row):
            mat.append(cur_row.copy())
            cur_row.clear()
            new_row = False
    
    np_mat = np.array(mat)

    return np_mat


def check_integer_input(prompt : str):

    out = input(prompt)
    try:
        out = int(out)

        return out
    except ValueError:
        print("input must be an integer")
    

def create_elements(n_elements : int):

    elements = []

    for i in range(n_elements):
        assembly_matrix_str = input("assembly matrix of element {i} (in format: [[1,2,3],[4,5,6]]): ".format(i=i))
        E = check_integer_input("youngs modulus of element {i}: ".format(i=i))
        L = check_integer_input("length of element {i}: ".format(i=i))
        A = check_integer_input("area of element {i}: ".format(i=i))
        angle = check_integer_input("angle of element {i}: ".format(i=i))

        assembly_matrix  = process_matrix(assembly_matrix_str)

        I = check_integer_input("Intertia of element {i}: ".format(i=i))
        UDL = check_integer_input("UDL of element {i}: ".format(i=i))
        point_load = check_integer_input("point load of element {i}: ".format(i=i))

        element = Element(assembly_matrix, E, I, L, A, angle, UDL, point_load)
        elements.append(element)

    return elements


def create_structure(elements) -> Structure:
    external_force_vector_str = input('External force vector? ')

    external_force_vector = process_matrix(external_force_vector_str)

    structure = Structure(elements, external_force_vector)

    return structure


def get_element_info(elements):
    pass

def get_structure_info(structure : Structure):
    pass

def process_input_info(structure : Structure):
    action = input(PROMPT)

    action_elems = action.split('.')
    

    if(action == "help"):
        print(HELP)
    elif("element" in action_elems[0]):
        pass
    elif(action_elems[0] == "structure"):
        pass
    elif(action_elems[0] == "display"):
        nodes_str = input("Enter nodal matrix or element/struture: ")
        nodes = process_matrix(nodes_str)

        display_mag = check_integer_input("display magnitude (how large do you want to scale the deflections by?): ")
        resolution = check_integer_input("resolution of plot: ")

        if(action_elems[1] == "structure"):
            structure.plot_structure(nodes, display_mag, resolution)
        elif("element" in action_elems[1]):
            start = action_elems[1].index("(")
            end = action_elems[1].index(")")

            index = int(action_elems[1][start+1:end])

            structure.elements[index].plot_element(nodes, int(display_mag), int(resolution))


def main():
    n_elements = get_n_elements()

    elements = create_elements(n_elements)

    structure = create_structure(elements)

    process_input_info(structure)


main()