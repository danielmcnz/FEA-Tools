import random
import sys
import matplotlib

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import *

from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.x_data = [x for x in range(0, 10)]
        self.y_data = [random.randint(0, 10) for x in range(0, 10)]

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
        self.elements = []

        self.add_element_btn = QPushButton("Add Element")
        self.calculate_fea_btn = QPushButton("Calculate FEA")
        

        # ----------------- #
        #  Central Widget   #
        self.central_layout = QGridLayout()

        self.canvas = MplCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.update_plot_btn = QPushButton("Update Plot") 
        
        self.update_ui()


    def update_right_pane(self):

        for i in range(self.n_elements):
            element_info = {}

            element_info['name'] = QLabel("Element {}".format(i + 1))
            element_info['E_text'] = QLabel("Young's Modulus")
            element_info['E_input'] = QLineEdit("200e9")
            element_info['I text'] = QLabel("Interia")
            element_info['I input'] = QLineEdit("1e-5")
            element_info['L text'] = QLabel("Length")
            element_info['L input'] = QLineEdit("10")
            element_info['A text'] = QLabel("Area")
            element_info['A input'] = QLineEdit("4e-5")

            height = i * len(element_info)

            self.right_layout.addWidget(element_info['name'],     height + 1, 0, 1, 2)
            self.right_layout.addWidget(element_info['E_text'],   height + 2, 0, 1, 1)
            self.right_layout.addWidget(element_info['E_input'],  height + 2, 1, 1, 1)
            self.right_layout.addWidget(element_info['I text'],   height + 3, 0, 1, 1)
            self.right_layout.addWidget(element_info['I input'],  height + 3, 1, 1, 1)
            self.right_layout.addWidget(element_info['L text'],   height + 4, 0, 1, 1)
            self.right_layout.addWidget(element_info['L input'],  height + 4, 1, 1, 1)
            self.right_layout.addWidget(element_info['A text'],   height + 5, 0, 1, 1)
            self.right_layout.addWidget(element_info['A input'],  height + 5, 1, 1, 1)

        self.right_layout.addWidget(self.calculate_fea_btn, 0, 0, 1, 1)
        self.right_layout.addWidget(self.add_element_btn, 0, 1, 1, 1)


    def update_ui(self):
        self.setWindowTitle("FEA 2D")
        self.setStyleSheet("background-color: white;")
        
        self.title.setStyleSheet("font-size: 30px; font-weight: bold;")
        self.title.setAlignment(Qt.AlignCenter)

        self.central_layout.addWidget(self.toolbar, 0, 0)
        self.central_layout.addWidget(self.canvas, 1, 0)
        self.central_layout.addWidget(self.update_plot_btn, 2, 0)

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
    def update_plot(self):
        self.y_data = [random.randint(0, 10) for x in range(0, 10)]

        self.canvas.axes.cla()
        self.canvas.axes.plot(self.x_data, self.y_data)
        self.canvas.draw()


    @Slot()
    def add_element(self):
        self.n_elements += 1

        self.update_right_pane()

    
    @Slot()
    def calculate_fea(self):
        pass


if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = MainWindow()

    w.resize(800, 600)
    w.show()

    sys.exit(app.exec())