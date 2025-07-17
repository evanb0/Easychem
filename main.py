import sys
import os

from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QGridLayout, QLabel, QComboBox, 
                               QLineEdit, QPushButton, QTextEdit, QFrame,
                               QGroupBox, QCheckBox, QSpinBox, QTabWidget)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from molecularviewer import MoleculeViewer
from parameterpanel import ParameterPanel 
from filepanel import FilePanel   

# The main 'Screen' we see if you will. Basically standard from the book, I just put it into classes rather than running it directly as they do in the book

class EasyChemMainWindow(QMainWindow):
    """Main application window"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EasyChem")
        self.setMinimumSize(1000, 700) 
        self.setup_ui()
        
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout()
        
        # Right panel - Molecular viewer
        right_panel = QVBoxLayout()
        viewer_label = QLabel("Molecular Analysis")
        viewer_label.setAlignment(Qt.AlignCenter)
        viewer_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 5px; font-color;")
        
        self.molecular_viewer = MoleculeViewer()
        
        right_panel.addWidget(viewer_label)
        right_panel.addWidget(self.molecular_viewer)
        
        # Left Panel - Parameters
        left_panel = QVBoxLayout()
        param_label = QLabel("Calculation Parameters")
        param_label.setAlignment(Qt.AlignCenter)
        param_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 5px;")
        
        from inputmodel import InputModel
        self.input_model = InputModel()
        self.parameter_panel = ParameterPanel(model=self.input_model)
        
        left_panel.addWidget(param_label)
        left_panel.addWidget(self.parameter_panel)
        
        # Middel Panel - File operations
        middle_panel = QVBoxLayout()
        file_label = QLabel("File Operations")
        file_label.setAlignment(Qt.AlignCenter)
        file_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 5px; font-color;")
        
        self.file_panel = FilePanel(model=self.input_model, parameter_panel=self.parameter_panel)
        
        middle_panel.addWidget(file_label)
        middle_panel.addWidget(self.file_panel)
        
        # Add panels to main layout
        main_layout.addLayout(left_panel, 2)  
        main_layout.addLayout(middle_panel, 1)  
        main_layout.addLayout(right_panel, 2)  
        
        central_widget.setLayout(main_layout)





class EasyChemApp(QApplication):
    """Main application class that inherits from QApplication"""
    
    def __init__(self, argv):
        """Constructor that initializes the application"""
        super().__init__(argv)  # Call parent constructor with command line arguments
        self.setApplicationName("EasyChem")  # Set the application name
        
        # Set application style 
        self.setStyle('Fusion')
        # Set global background color to orange
        self.setStyleSheet("QWidget { background-color: #7785AC; }")
        # Create main window instance
        self.main_window = EasyChemMainWindow()
        
        
    def run(self):
        """Method to run the application and start the event loop"""
        self.main_window.show()  # Show the main window
        return self.exec()       # Start the event loop and return exit code
    

def main():
    """Main entry point function that creates and runs the application"""
    app = EasyChemApp(sys.argv)  # Create application instance with command line arguments
    sys.exit(app.run())          # Run the application and exit with its return code 

if __name__ == "__main__":
    main()  

    