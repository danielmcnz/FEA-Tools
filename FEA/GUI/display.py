import random
import sys
import matplotlib

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import (QApplication, QLabel, QPushButton,
                               QVBoxLayout, QWidget, QMainWindow)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# from __feature__ import snake_case, true_property

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        sc = MplCanvas(self, width=5, height=4, dpi=100)
        sc.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])
        self.setCentralWidget(sc)
        
        self.show()

class MyWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        self.hello = [
            "Hello world",
            "Hallo Welt"
        ]

        self.button = QPushButton("Click me!")
        self.message = QLabel("Hello World")
        self.message.setAlignment(Qt.AlignCenter)

        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.message)
        self.layout.addWidget(self.button)
        self.button.clicked.connect(self.magic)

    @Slot()
    def magic(self):
        self.message.setText(random.choice(self.hello))

if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = MainWindow()

    widget = MyWidget()
    widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())