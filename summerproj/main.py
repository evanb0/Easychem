
#IMPORT OF CLI/FILE MOVEMENT MODULES
import sys
import os

# Import of QtWidget Modules, added as I needed. Same for Other things added below.
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                               QHBoxLayout, QGridLayout, QLabel, QComboBox, 
                               QLineEdit, QPushButton, QTextEdit, QFrame,
                               QGroupBox, QCheckBox, QSpinBox, QTabWidget)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from molecularviewer import MolecularViewer # CUSTOM SUBCLASS MADE BY ME
from parameterpanel import ParameterPanel # CUSTOM SUBCLASS MADE BY ME
from filepanel import FilePanel # CUSTOM SUB CLASS MADE BY ME

#ALL OTHER SUB CLASSES HAVE BEEN IMPORTED, BUT I DONT ACTUALLY KNOW
#HOW THIS WILL RUN ON OTHER SYSTEMS, WILL PROBABLY NEED TO SPEAK TO JULIA ABOUT HOW THIS SHIT ACUTALLY WORKS
# CUZ IM NOT SURE

class EasyChemMainWindow(QMainWindow):
    """Main application window"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EasyChem")
        self.setMinimumSize(1000, 700) # Arbituary size, May need reusing
        self.setup_ui()
        
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout()
        
        # Right panel - Molecular viewer
        right_panel = QVBoxLayout()
        viewer_label = QLabel("Molecular Structure")
        viewer_label.setAlignment(Qt.AlignCenter)
        viewer_label.setStyleSheet("font-size: 16px; font-weight: bold; margin: 5px;")
        
        self.molecular_viewer = MolecularViewer()
        
        right_panel.addWidget(viewer_label)
        right_panel.addWidget(self.molecular_viewer)
        
        # Left Panel - Parameters
        left_panel = QVBoxLayout()
        param_label = QLabel("Calculation Parameters")
        param_label.setAlignment(Qt.AlignCenter)
        param_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 5px;")
        
        self.parameter_panel = ParameterPanel()
        
        left_panel.addWidget(param_label)
        left_panel.addWidget(self.parameter_panel)
        
        # Middel Panel - File operations
        middle_panel = QVBoxLayout()
        file_label = QLabel("File Operations")
        file_label.setAlignment(Qt.AlignCenter)
        file_label.setStyleSheet("font-size: 20px; font-weight: bold; margin: 5px; font-color;")
        
        self.file_panel = FilePanel()
        
        middle_panel.addWidget(file_label)
        middle_panel.addWidget(self.file_panel)
        
        # Add panels to main layout
        main_layout.addLayout(left_panel, 2)  # 2/5 of width
        main_layout.addLayout(middle_panel, 2)  # 2/5 of width
        main_layout.addLayout(right_panel, 1)  # 1/5 of width
        
        central_widget.setLayout(main_layout)



# The book normally does not get an app running in this specific way,
# But the process is basically the same.

class EasyChemApp(QApplication):
    """Main application class that inherits from QApplication"""
    
    def __init__(self, argv):
        """Constructor that initializes the application"""
        super().__init__(argv)  # Call parent constructor with command line arguments
        self.setApplicationName("EasyChem")  # Set the application name
        
        # Set application style 
        # I like Fusion and there was also a good bit on it in the book so went with it
        self.setStyle('Fusion')
        
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

# I didn't see this in the book anywhere, But I think this is a safer way to run GUI apps.
# 
if __name__ == "__main__":
    main()  # Call the main function to start the application

    