import sys
import os
from PySide6.QtWidgets import QWidget, QComboBox, QLabel, QGroupBox, QVBoxLayout, QHBoxLayout, QSpinBox, QSlider
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont



class ParameterPanel(QWidget): 
    """Parameter selection panel
    This class deals with our two calculation softwares, Gaussian & ORCA
    Included are some of the most common Functionals/basis sets for each I could find 
    (NOTE: WILL PROBABLY HAVE TO CHECK WITH JULIA TO SEE IF THIS IS GOOD.)
    """
    def __init__(self):
        super().__init__()
        self.setup_software_data()
        self.setup_ui()
    
    def setup_software_data(self):
        """Define software-specific functionals and basis sets"""
        # Orca functional and basis set definitions
        self.ORCA_FUNCTIONALS = ["HF", "MP2", "CCSD", "CCSD(T)", "BLYP", "PBE",
            "revPBE", "B3LYP", "MO6L", "MO62X", "B97-3C"
   
        ]
        
        self.ORCA_BASIS_SETS = [
            '6-31G(d)', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ',
            'aug-cc-pVTZ', 'aug-cc-pVQZ', 'def2-SVP', 'def2-TZVP', 'def2-QZVP',
            'def2-TVZPP', 'def2-QZVPP', 'def2-TZVPPD', 'def2-QZVPPD', 'ma-def2-SVP',
            'ma-def2-TZVP', 'ma-def2-QZVP',
        ]
        
        # Gaussian functional and basis set definitions
        self.GAUSSIAN_FUNCTIONALS =['AM1', 'PM3', 'RHF', 'B3LYP', 'WB97XD','MP2', 'CCSD']
        
        self.GAUSSIAN_BASIS_SETS = [
            'STO-3G', '3-21G', '6-31G(d)', '6-31G(d,p)', 'LANL2DZ',
            'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-p5VZ', 'cc-p6VZ',
            'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-p5VZ', 'aug-cc-p6VZ',
            'Def2SV', 'Def2TZV', 'Def2QZV', 'Def2SVP', 'Def2TZVP',
            'Def2QZVP', 'Def2SVPP', 'Def2TZVPP', 
        ]
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        # Software selection
        software_group = QGroupBox("Software")
        software_group.setStyleSheet("font-weight: bold")
        software_layout = QHBoxLayout()
        
        self.software_combo = QComboBox()
        self.software_combo.addItems(["Orca", "Gaussian"]) # A combo box containing Our softwares
        # These names should not be modified from how they currently are. Errors will occur if so.
        #If for any reason, e.g to add more calculation parameters, Specification is also needed
        #within on_software_changed function to prevent errors.
        self.software_combo.currentTextChanged.connect(self.on_software_changed)
        software_layout.addWidget(QLabel("Software:",))
        software_layout.addWidget(self.software_combo)
        software_group.setLayout(software_layout)
        
        # Functional selection - Selection for our Functions/Basis sets
        functional_group = QGroupBox("Method")
        functional_group.setStyleSheet("font-weight: bold")
        functional_layout = QVBoxLayout()
        
        func_row = QHBoxLayout()
        func_row.addWidget(QLabel("Functional"))
        self.functional_combo = QComboBox()
        func_row.addWidget(self.functional_combo)
        functional_layout.addLayout(func_row)
        
        basis_row = QHBoxLayout()
        basis_row.addWidget(QLabel("Basis Set"))
        self.basis_combo = QComboBox()
        basis_row.addWidget(self.basis_combo)
        functional_layout.addLayout(basis_row)
        
        functional_group.setLayout(functional_layout)
        
        # Initialize with default software (Orca)
        self.on_software_changed("Orca") # Inital state when app is launched will be ORCA.
        
        # Processor count selection
        processor_group = QGroupBox("Computational Resources")
        processor_group.setStyleSheet('font-weight: bold')
        processor_layout = QVBoxLayout()
        
        proc_row = QHBoxLayout()
        proc_row.addWidget(QLabel("Number of Processors:"))
        self.nprocs_label = QLabel("1")
        proc_row.addWidget(self.nprocs_label)
        processor_layout.addLayout(proc_row)
        
        self.nprocs_slider = QSlider(Qt.Horizontal)
        self.nprocs_slider.setRange(1, 60)
        self.nprocs_slider.setValue(1)
        self.nprocs_slider.valueChanged.connect(self.on_nprocs_changed)
        processor_layout.addWidget(self.nprocs_slider)
        
        processor_group.setLayout(processor_layout)

        #Memory 
        memory_group = QGroupBox("Memory")
        memory_group.setStyleSheet('font-weight: bold')
        memory_layout = QVBoxLayout()

        mem_row = QHBoxLayout()
        mem_row.addWidget(QLabel('Memory:'))
        self.mem_label = QLabel("1")
        mem_row.addWidget(self.mem_label)
        memory_layout.addLayout(mem_row)

        self.mem_slider = QSlider(Qt.Horizontal)
        self.mem_slider.setRange(1, 16)
        self.mem_slider.setValue(1)
        self.mem_slider.valueChanged.connect(self.on_mem_changed)
        memory_layout.addWidget(self.mem_slider)

        memory_group.setLayout(memory_layout)
        
        # Solvation model
        solvation_group = QGroupBox("Solvation")
        solvation_group.setStyleSheet('font-weight: bold')
        solvation_layout = QVBoxLayout()
        
        solv_row = QHBoxLayout()
        solv_row.addWidget(QLabel("Solvent Model"))
        self.solvation_combo = QComboBox()
        self.solvation_combo.addItems(["None", "PCM", "SMD", "COSMO"])
        solv_row.addWidget(self.solvation_combo)
        solvation_layout.addLayout(solv_row)
        
        solvent_row = QHBoxLayout()
        solvent_row.addWidget(QLabel("Solvent"))
        self.solvent_combo = QComboBox()
        self.solvent_combo.addItems(["None", "Water", "Methanol", "Ethanol", "Acetone", "DMSO"])
        solvent_row.addWidget(self.solvent_combo)
        solvation_layout.addLayout(solvent_row)
        
        solvation_group.setLayout(solvation_layout)
        
        # Calculation type -> What do you want to calculate? OPT/FRQS ETC.
        calc_group = QGroupBox("Calculation Type")
        calc_group.setStyleSheet('font-weight: bold')
        calc_layout = QVBoxLayout()
        
        calc_row = QHBoxLayout()
        calc_row.addWidget(QLabel("Job Type:"))
        self.calc_combo = QComboBox()
        self.calc_combo.addItems(["Single Point", "Geometry Optimization", "Frequency", "TD-DFT"])
        calc_row.addWidget(self.calc_combo)
        calc_layout.addLayout(calc_row)
        
        # Additional options -> IDK WHAT OTHER PARAMATERS WOULD BE USED, NEED TO DISCUSS WITH JULIA
        
        charge_spin_row = QHBoxLayout()
        charge_spin_row.addWidget(QLabel("Charge:"))
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(0)
        charge_spin_row.addWidget(self.charge_spin)
        
        charge_spin_row.addWidget(QLabel("Multiplicity:"))
        self.multiplicity_spin = QSpinBox()
        self.multiplicity_spin.setRange(1, 10)
        self.multiplicity_spin.setValue(1)
        charge_spin_row.addWidget(self.multiplicity_spin)
        calc_layout.addLayout(charge_spin_row)
        
        calc_group.setLayout(calc_layout)
        
        # Add all groups to main layout
        layout.addWidget(software_group)
        layout.addWidget(functional_group)
        layout.addWidget(processor_group)
        layout.addWidget(memory_group)
        layout.addWidget(solvation_group)
        layout.addWidget(calc_group)
        layout.addStretch()
        
        self.setLayout(layout)
    
    def on_software_changed(self, software_name):
        """Update functional and basis set options based on selected software"""
        # Clear existing items
        self.functional_combo.clear()
        self.basis_combo.clear()
        
        # Purpose of this is to ensure you are not trying to run an ORCA basis set within
        # A gaussian input file, we want our inputs to actually work when complete!
        if software_name == "Orca":
            # Add Orca functionals and basis sets
            self.functional_combo.addItems(self.ORCA_FUNCTIONALS)
            self.basis_combo.addItems(self.ORCA_BASIS_SETS)
        elif software_name == "Gaussian":
            # Add Gaussian functionals and basis sets
            self.functional_combo.addItems(self.GAUSSIAN_FUNCTIONALS)
            self.basis_combo.addItems(self.GAUSSIAN_BASIS_SETS)

    def on_nprocs_changed(self, value):
        """Update the processor count label when slider changes"""
        self.nprocs_label.setText(str(value))
    
    def get_current_functional(self):
        """Get the current functional selection"""
        return self.functional_combo.currentText() # Returns current selection 
    
    def get_current_basis_set(self):
        """Get the current basis set selection"""
        return self.basis_combo.currentText() # Returns current selection of basis set.
    
    def get_current_nprocs(self):
        """Get the current processor count selection"""
        return self.nprocs_slider.value() # Returns current processor count
    
    def on_mem_changed(self, value):
        """Update Memory when slider changes"""
        self.mem_label.setText(str(value))

    def get_current_mem(self):
        """Get Current Memory allocation"""
        return self.mem_slider.value()