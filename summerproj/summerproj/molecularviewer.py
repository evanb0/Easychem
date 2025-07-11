import sys
import os
from PySide6.QtWidgets import QFrame, QVBoxLayout, QLabel
from PySide6.QtCore import Qt


# NOTE: THIS WILL BE SKIPPED FOR NOW IDFK HOW TO DO THIS WILL COME BACK TO (VERY) LATER (EMPHASIS ON VERY)
# As of version 2 this has still not been changed at all.

class MolecularViewer(QFrame):
    """Molecular viewer placeholder widget
    As of now, simply just creates a box that will allow for
    A view of the molecule, but currently there is no way 
    To actually implement a molecule, so kinda useless rn"""
    def __init__(self):
        super().__init__()
        self.setFrameStyle(QFrame.Box) 
        self.setMinimumSize(400, 300) 
        self.setStyleSheet("background-color: #008080; border: 2px solid gray;") 
        
        layout = QVBoxLayout()
        label = QLabel("Molecular Viewer")
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("color: white; font-size: 14px;")
        layout.addWidget(label)
        self.setLayout(layout)