import random
import sys
import matplotlib
import json
import numpy as np

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

from .Element import Element
from .Structure import Structure

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class ProcessInput():
    def __init__(self):
        pass


    @staticmethod
    def str_to_matrix(string : str) -> np.ndarray:
        elements = string.split(',')

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
                
            cur_row.append(float(element))

            if(new_row):
                mat.append(cur_row.copy())
                cur_row.clear()
                new_row = False
        
        np_mat = np.array(mat)

        return np_mat
    

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.main_layout = QGridLayout(self)
        self.widget = QWidget()

        # ----------------- #
        #  Menu Dock        #
        self.title = QLabel("Finite Element Analysis 2D")

        # ----------------- #
        #  Right Dock       #
        self.right_layout = QGridLayout()
        self.right_widget = QWidget()
        self.right_widget_scroll = QScrollArea()

        self.n_elements : int = 1
        self.external_force_vector = QLineEdit("[[0],[140e3],[0]]")
        self.elements = []

        self.add_element_btn = QPushButton("Add Element")
        self.calculate_fea_btn = QPushButton("Calculate FEA")
        

        # ----------------- #
        #  Central Widget   #
        self.central_layout = QGridLayout()

        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.nodes = QLineEdit("[[[0,10],[0,0]],[[0,0],[10,0]]]")
        self.plot_displacement_magnitude = QLineEdit("100")
        self.plot_resolution = QLineEdit("20")
        self.update_plot_btn = QPushButton("Update Plot") 
        
        self.update_ui()


    def update_right_pane(self):

        for i in range(self.n_elements):

            if(i == len(self.elements)):
                element_info = {}

                element_info['name'] = QLabel("Element {}".format(i + 1))
                element_info['asm mat text'] = QLabel("Assembly Matrix")
                element_info['asm mat input'] = QLineEdit("[[0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,0,0]]")
                element_info['E text'] = QLabel("Young's Modulus")
                element_info['E input'] = QLineEdit("200e9")
                element_info['I text'] = QLabel("Interia")
                element_info['I input'] = QLineEdit("5e-4")
                element_info['L text'] = QLabel("Length")
                element_info['L input'] = QLineEdit("10")
                element_info['A text'] = QLabel("Area")
                element_info['A input'] = QLineEdit("1e-5")
                element_info['angle text'] = QLabel("Angle")
                element_info['angle input'] = QLineEdit("-90")
                element_info['UDL text'] = QLabel("UDL")
                element_info['UDL input'] = QLineEdit("0")
                element_info['LVL text'] = QLabel("LVL")
                element_info['LVL input'] = QLineEdit("0")
                element_info['point load len text'] = QLabel("Point Load Length")
                element_info['point load len input'] = QLineEdit("0")
                element_info['point load angle text'] = QLabel("Point Load Angle")
                element_info['point load angle input'] = QLineEdit("0")
                element_info['point load mag text'] = QLabel("Point Load Magnitude")
                element_info['point load mag input'] = QLineEdit("0")
                element_info['point load x sign text'] = QLabel("Point Load x sign")
                element_info['point load x sign input'] = QLineEdit("0")
                element_info['point load y sign text'] = QLabel("Point Load y sign")
                element_info['point load y sign input'] = QLineEdit("0")
                
                self.elements.append(element_info)

            height = i * 28

            self.right_layout.addWidget(self.elements[i]['name'],                   height + 2, 0, 1, 2)
            self.right_layout.addWidget(self.elements[i]['asm mat text'],           height + 3, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['asm mat input'],          height + 3, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['E text'],                 height + 4, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['E input'],                height + 4, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['I text'],                 height + 5, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['I input'],                height + 5, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['L text'],                 height + 6, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['L input'],                height + 6, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['A text'],                 height + 7, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['A input'],                height + 7, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['UDL text'],               height + 8, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['UDL input'],              height + 8, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['LVL text'],               height + 9, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['LVL input'],              height + 9, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['angle text'],             height + 10, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['angle input'],            height + 10, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load len text'],    height + 11, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load len input'],   height + 11, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load angle text'],  height + 12, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load angle input'], height + 12, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load mag text'],    height + 13, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load mag input'],   height + 13, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load x sign text'], height + 14, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load x sign input'],height + 14, 1, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load y sign text'], height + 15, 0, 1, 1)
            self.right_layout.addWidget(self.elements[i]['point load y sign input'],height + 15, 1, 1, 1)

        self.right_layout.addWidget(QLabel("External Force Vector"), 1, 0, 1, 1)
        self.right_layout.addWidget(self.external_force_vector, 1, 1, 1, 1)

        self.right_layout.addWidget(self.calculate_fea_btn, 0, 0, 1, 1)
        self.right_layout.addWidget(self.add_element_btn, 0, 1, 1, 1)


    def update_central_pane(self):
        self.central_layout.addWidget(self.toolbar, 0, 0)
        self.central_layout.addWidget(self.canvas, 1, 0)

        self.central_layout.addWidget(QLabel("Nodes"), 2, 0)
        self.central_layout.addWidget(self.nodes, 2, 1)
        self.central_layout.addWidget(QLabel("Displacement Magnitude"), 3, 0)
        self.central_layout.addWidget(self.plot_displacement_magnitude, 3, 1)
        self.central_layout.addWidget(QLabel("Resolution"), 4, 0)
        self.central_layout.addWidget(self.plot_resolution, 4, 1)

        self.central_layout.addWidget(self.update_plot_btn, 5, 0)


    def update_ui(self):
        self.setWindowTitle("FEA 2D")
        self.setStyleSheet("background-color: white;")
        
        self.title.setStyleSheet("font-size: 30px; font-weight: bold;")
        self.title.setAlignment(Qt.AlignCenter)

        self.update_central_pane()
        self.update_right_pane()

        self.update_plot_btn.clicked.connect(self.update_plot)
        self.add_element_btn.clicked.connect(self.add_element)
        self.calculate_fea_btn.clicked.connect(self.calculate_fea)

        self.right_widget.setLayout(self.right_layout)

        self.right_widget_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.right_widget_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.right_widget_scroll.setWidgetResizable(True)
        self.right_widget_scroll.setWidget(self.right_widget)
        
        self.main_layout.addWidget(self.right_widget_scroll, 0, 1, 1, 1)
        self.main_layout.addLayout(self.central_layout, 0, 0, 1, 1)
        self.widget.setLayout(self.main_layout)

        self.setMenuWidget(self.title)
        self.setCentralWidget(self.widget)
    

    @Slot()
    def update_plot(self, x_data, y_data):
        self.canvas.axes.cla()
        self.canvas.axes.plot(x_data, y_data)
        self.canvas.draw()


    @Slot()
    def add_element(self):
        self.n_elements += 1

        self.update_right_pane()

    
    @Slot()
    def calculate_fea(self):

        elems = []

        for element in self.elements:
            elem = Element(
                np.array(json.loads(element['asm mat input'].text())),
                (float(element['E input'].text())), 
                (float(element['I input'].text())), 
                (float(element['L input'].text())), 
                (float(element['A input'].text())), 
                (float(element['angle input'].text())),
                (float(element['UDL input'].text())),
                (float(element['LVL input'].text())),
                (
                    (float(element['point load len input'].text())),
                    (float(element['point load angle input'].text())),
                    (float(element['point load mag input'].text())),
                    (float(element['point load x sign input'].text())),
                    (float(element['point load y sign input'].text())),
                )
            )
            elems.append(elem)
            
        external_force_vector = np.array(json.loads(self.external_force_vector.text()))

        structure = Structure(elems, external_force_vector)

        structure.solve()

        nodes = np.array(json.loads(self.nodes.text()))
        displacement_magnitude = int(float(self.plot_displacement_magnitude.text()))
        resolution = int(float(self.plot_resolution.text()))

        for element in structure.elements:
            print(element.global_force)

        self.canvas.axes.cla()
        structure.plot_structure(nodes, displacement_magnitude, resolution, self.canvas.axes)
        self.canvas.draw()
        



if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = MainWindow()

    w.resize(800, 600)
    w.show()

    sys.exit(app.exec())