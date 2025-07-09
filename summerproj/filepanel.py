import sys
import os
from PySide6.QtWidgets import (
    QWidget, 
    QVBoxLayout, 
    QGroupBox, 
    QHBoxLayout, 
    QLabel, 
    QLineEdit, 
    QPushButton, 
    QTextEdit
)

from PySide6.QtCore import QThread, Signal

class FilePanel(QWidget):
    """File input/output panel"""
    def __init__(self):
        super().__init__()
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Input file section
        input_group = QGroupBox("Input Files")
        input_layout = QVBoxLayout()
        
        input_row = QHBoxLayout()
        input_row.addWidget(QLabel("File"))
        self.input_file_edit = QLineEdit()
        self.input_file_edit.setPlaceholderText("Select File")
        input_row.addWidget(self.input_file_edit)
        self.browse_input_btn = QPushButton("Browse")
        input_row.addWidget(self.browse_input_btn)
        input_layout.addLayout(input_row)
        
        input_group.setLayout(input_layout)
        input_group.setStyleSheet("font-weight: bold;")
        
        # Output directory section
        output_group = QGroupBox("Output")
        output_layout = QVBoxLayout()
        output_group.setStyleSheet("font-weight: bold;")
        
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Directory"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory")
        output_row.addWidget(self.output_dir_edit)
        self.browse_output_btn = QPushButton("Browse")
        output_row.addWidget(self.browse_output_btn)
        output_layout.addLayout(output_row)
        
        
        output_group.setLayout(output_layout)
        
        # Make button
        self.make_btn = QPushButton("Make Input Files")
        self.make_btn.setMinimumHeight(40)
        self.make_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        
        # Progress/Status area
        status_group = QGroupBox("Status")
        status_layout = QVBoxLayout()
        
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(100)
        self.status_text.setPlaceholderText("Ready to create input files...")
        self.status_text.setReadOnly(True)
        status_layout.addWidget(self.status_text)
        
        status_group.setLayout(status_layout)
        
        # Add all components to layout
        layout.addWidget(input_group)
        layout.addWidget(output_group)
        layout.addWidget(self.make_btn)
        layout.addWidget(status_group)
        layout.addStretch()
        
        self.setLayout(layout)