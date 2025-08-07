import os
import sys
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QDialog,
    QMessageBox, QFileDialog, QGroupBox, QTextEdit, QApplication, QSlider, QFormLayout
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import QEventLoop, QUrl, Signal, Qt
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import py3Dmol
import math


class MoleculeViewer(QWidget):
    """
    Main widget for Molecule Viewer application.
    Displays molecules from SDF files, with the 3D view integrated directly.
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecule Viewer")
        self.molecules = []
        self.setup_ui()
        # Initially hide the 3D viewer group
        self.viewer_group.setVisible(False) 

    def setup_ui(self):
        # This container holds all the options, but NOT the 3D viewer itself
        self.options_container = QWidget()
        options_layout = QVBoxLayout(self.options_container)
        options_layout.setContentsMargins(0, 0, 0, 0)
        options_layout.setSpacing(10)
        
        # File Information Section
        file_group = QGroupBox("File Information")
        file_layout = QVBoxLayout()
        self.file_info_label = QLabel("No file selected")
        self.file_info_label.setWordWrap(True)
        file_layout.addWidget(self.file_info_label)
        file_group.setLayout(file_layout)
        options_layout.addWidget(file_group)

        # Display Options Section
        display_group = QGroupBox("Display Options")
        display_layout = QVBoxLayout()

        combo_style = """
            QComboBox { 
                background-color: #edf2f4; 
                color: #2b2d42; 
                border: 1px solid #8d99ae; 
                border-radius: 3px; 
                padding: 2px; 
                combobox-popup: 0;
            }
            QComboBox QAbstractItemView {
                background-color: #edf2f4;
                color: #2b2d42;
                selection-background-color: #8d99ae;
                selection-color: #edf2f4;
            }
        """

        style_row = QHBoxLayout()
        style_row.addWidget(QLabel("Display Style:"))
        self.style_selector = QComboBox()
        self.style_selector.setStyleSheet(combo_style)
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'Surface'])
        self.style_selector.currentIndexChanged.connect(self.render_selected_molecule)
        style_row.addWidget(self.style_selector)
        display_layout.addLayout(style_row)

        molecule_row = QHBoxLayout()
        molecule_row.addWidget(QLabel("Molecule:"))
        self.molecule_selector = QComboBox()
        self.molecule_selector.setStyleSheet(combo_style)
        self.molecule_selector.currentIndexChanged.connect(self.on_molecule_change)
        self.molecule_selector.setMaxVisibleItems(10)
        molecule_row.addWidget(self.molecule_selector)
        display_layout.addLayout(molecule_row)
        
        #Zoom slider (this actually has to be dragged, otherwise it doesnt work idk why)
        # zoom_row = QHBoxLayout()
        # zoom_row.addWidget(QLabel("Zoom:"))
        # self.zoom_slider = QSlider(Qt.Horizontal)
        # self.zoom_slider.setMinimum(100)
        # self.zoom_slider.setMaximum(200)
        # self.zoom_slider.setValue(120)
        # self.zoom_slider.setTickInterval(10)
        # self.zoom_slider.setSingleStep(10)
        # self.zoom_label = QLabel("1.2x")
        # self.zoom_slider.valueChanged.connect(self.update_zoom_label)
        # self.zoom_slider.sliderReleased.connect(self.render_selected_molecule)
        
        # zoom_row.addWidget(self.zoom_slider)
        # zoom_row.addWidget(self.zoom_label)
        # display_layout.addLayout(zoom_row)

        # // List slider  \\
        
        display_group.setLayout(display_layout)
        options_layout.addWidget(display_group)

        # Buttons Section
        button_layout = QHBoxLayout()
        
        self.export_sdf_button = QPushButton("Export Selected")
        self.export_sdf_button.setMinimumHeight(30)
        self.export_sdf_button.clicked.connect(self.on_export_sdf)
        self.export_sdf_button.setEnabled(False)

        #Button for exporting lowest energy molecules
        self.export_lowest_energy_molecules_btn = QPushButton("Export 10 Lowest Energy Molecules")
        self.export_lowest_energy_molecules_btn.setMinimumHeight(30)
        self.export_lowest_energy_molecules_btn.clicked.connect(self.export_lowest_energy_molecules)
        self.export_lowest_energy_molecules_btn.setEnabled(False)

        button_layout.addWidget(self.export_sdf_button)
        button_layout.addWidget(self.export_lowest_energy_molecules_btn)
        options_layout.addLayout(button_layout)

        # Molecular Information Section
        self.info_group = QGroupBox("Molecular Information")
        self.info_layout = QFormLayout()
        
        self.formula_label = QLabel("N/A") # all states N/A before file imported.
        self.weight_label = QLabel("N/A")
        self.atoms_label = QLabel("N/A")
        self.bonds_label = QLabel("N/A")
        self.charge_label = QLabel("N/A")
        self.mult_label = QLabel("N/A")
        self.energy_label = QLabel("N/A")

        self.info_layout.addRow("Formula:", self.formula_label)
        self.info_layout.addRow("Weight:", self.weight_label)
        self.info_layout.addRow("Atoms:", self.atoms_label)
        self.info_layout.addRow("Bonds:", self.bonds_label)
        self.info_layout.addRow("Charge:", self.charge_label)
        self.info_layout.addRow("Multiplicity:", self.mult_label)
        self.info_layout.addRow("(MMFF) Energy:", self.energy_label)
        
        self.info_group.setLayout(self.info_layout)
        options_layout.addWidget(self.info_group)
        
        options_layout.addStretch()

        # 3D Viewer Section - This is the part that will be hidden/shown
        self.viewer_group = QGroupBox("3D Viewer")
        viewer_layout = QVBoxLayout()
        self.web_view = QWebEngineView()
        viewer_layout.addWidget(self.web_view)
        self.viewer_group.setLayout(viewer_layout)

    # def update_zoom_label(self, value):
    #     zoom_factor = value / 100.0
    #     self.zoom_label.setText(f"{zoom_factor:.1f}x")

    def load_molecules_from_file(self, filename):
        if not filename or not filename.lower().endswith('.sdf'):
            self.clear_molecules()
            return
        
        try:
            supplier = Chem.SDMolSupplier(filename)
            self.molecules = []
            
            for i, mol in enumerate(supplier):
                if mol is not None:
                    self.molecules.append((mol, i))

            if not self.molecules:
                QMessageBox.warning(self, "Error", "No valid molecules found in the SDF file.")
                self.clear_molecules()
                return

            self.viewer_group.setVisible(True)
            file_basename = os.path.basename(filename)
            self.file_info_label.setText(f"File: {file_basename}\nMolecules found: {len(self.molecules)}")

            self.molecule_selector.currentIndexChanged.disconnect()
            self.molecule_selector.clear()
            
            for i, (mol, original_idx) in enumerate(self.molecules):
                mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Molecule {original_idx + 1}"
                if not mol_name.strip():
                    mol_name = f"Molecule {original_idx + 1}"
                
                try:
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    label = f"{mol_name} ({formula})"
                except:
                    label = mol_name
                
                self.molecule_selector.addItem(label)

            self.molecule_selector.currentIndexChanged.connect(self.on_molecule_change)

            self.export_sdf_button.setEnabled(True)
            self.export_lowest_energy_molecules_btn.setEnabled(True)
            
            if self.molecules:
                self.molecule_selector.setCurrentIndex(0)
                self.on_molecule_change(0)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load molecules from SDF file.\n\n{str(e)}")
            self.clear_molecules()
            self.viewer_group.setVisible(False)

    def clear_molecules(self):
        self.molecules = []
        self.file_info_label.setText("No file selected")
        self.molecule_selector.clear()
        self.export_sdf_button.setEnabled(False)
        self.export_lowest_energy_molecules_btn.setEnabled(False)
        self.web_view.setHtml("")
        self.viewer_group.setVisible(False)
        self.formula_label.setText("N/A")
        self.weight_label.setText("N/A")
        self.atoms_label.setText("N/A")
        self.bonds_label.setText("N/A")
        self.energy_label.setText("N/A")

    def on_molecule_change(self, index):
        if index < 0 or not self.molecules or index >= len(self.molecules):
            self.formula_label.setText("N/A")
            self.weight_label.setText("N/A")
            self.atoms_label.setText("N/A")
            self.bonds_label.setText("N/A")
            self.charge_label.setText("N/A")
            self.mult_label.setText("N/A")
            self.energy_label.setText("N/A")
            return
        
        mol, original_idx = self.molecules[index]
        
        try:
            # Display basic molecule info
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            num_atoms = mol.GetNumAtoms()
            num_bonds = mol.GetNumBonds()
            charge = Chem.GetFormalCharge(mol)
            mult = self.getMult(mol)
            
            mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Molecule {original_idx + 1}"
            if not mol_name.strip():
                mol_name = f"Molecule {original_idx + 1}"
            
            self.formula_label.setText(formula)
            self.weight_label.setText(f"{mol_weight:.2f} g/mol")
            self.atoms_label.setText(str(num_atoms))
            self.bonds_label.setText(str(num_bonds))
            self.charge_label.setText(str(charge))
            self.mult_label.setText(str(mult))
            self.energy_label.setText("Calculating...")

            QApplication.processEvents()

            # Calculate and display the energy
            try:
                # First check for 3D coordinates, if not present, embed them
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    
                # Perform a quick energy minimization
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
                ff.Minimize()
                energy = ff.CalcEnergy()
                
                # Update the displayed info with the energy
                self.energy_label.setText(f"{energy:.4f} kcal/mol")

            except Exception as energy_e:
                self.energy_label.setText("Error")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error getting molecule details: {str(e)}")
            self.formula_label.setText("Error")
            self.weight_label.setText("Error")
            self.atoms_label.setText("Error")
            self.bonds_label.setText("Error")
            self.charge_label.setText("Error")
            self.mult_label.setText("Error")
            self.energy_label.setText("Error")
            
        self.render_selected_molecule()

    def getMult(self, mol):
        """2S + 1, where S is the total spin, 1 unpaired electron = 1/2"""
        unpaired_electrons = 0
        for atom in mol.GetAtoms():
            unpaired_electrons += atom.GetNumRadicalElectrons()
        return int(unpaired_electrons / 2 + 1)

    def render_selected_molecule(self):
        if not self.molecules:
            return

        index = self.molecule_selector.currentIndex()
        if index < 0 or index >= len(self.molecules):
            return

        mol, original_idx = self.molecules[index]
        style = self.style_selector.currentText()
        
        if not mol.GetNumConformers():
            try:
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                QMessageBox.critical(self, "Error", "Could not generate 3D coordinates for this molecule.")
                return
        
        mol_block = Chem.MolToMolBlock(mol)
        
        viewer = py3Dmol.view(width='100%', height='100%')
        viewer.addModel(mol_block, 'mol')

        if style == 'Stick':
            viewer.setStyle({'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({'stick': {}, 'sphere': {}})
        elif style == 'Surface':
            viewer.setStyle({'stick': {}})
            viewer.addSurface(py3Dmol.VDW, {'opacity': 0.8})
        
        # zoom_factor = self.zoom_slider.value() / 100.0
        viewer.zoomTo()
        # viewer.zoom(zoom_factor)
        html = viewer._make_html()
        self.web_view.setHtml(html)

    def on_export_sdf(self):
        if not self.molecules:
            QMessageBox.warning(self, "Error", "No molecules available. Please select an SDF file in File Operations.")
            return

        current_index = self.molecule_selector.currentIndex()
        if current_index == -1:
            QMessageBox.warning(self, "Error", "No molecule selected. Please select a molecule to export.")
            return

        mol, original_idx = self.molecules[current_index]

        filename, _ = QFileDialog.getSaveFileName(self, "Save SDF File", "molecule.sdf", "SDF Files (*.sdf)")

        if filename:
            try:
                writer = SDWriter(filename)
                writer.write(mol)
                writer.close()
                
                export_message = f"SUCCESS: SDF file exported!\n\n"
                export_message += f"File saved to: {filename}\n"
                export_message += f"Molecule: {self.molecule_selector.currentText()}"
                
                QMessageBox.information(self, "Export Successful", export_message)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export SDF file.\n\n{str(e)}")

    def export_lowest_energy_molecules(self):
        if not self.molecules:
            QMessageBox.warning(self, "Error", "No molecules available. Please select an SDF file.")
            return

        QMessageBox.information(self, "Calculation", "Calculating energies for all molecules... This may take a moment.")
        QApplication.processEvents()

        molecules_with_energy = []
        for mol, original_idx in self.molecules:
            try:
                # Ensure the molecule has 3D coordinates
                if mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

                # Perform energy minimization and get the energy
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
                ff.Minimize()
                energy = ff.CalcEnergy()
                
                # Store the molecule object and its calculated energy
                molecules_with_energy.append((mol, energy))
            except Exception:
                # If energy calculation fails, skip this molecule
                continue

        if not molecules_with_energy:
            QMessageBox.critical(self, "Error", "No molecules could be processed for energy calculation.")
            return

        # Sort molecules by energy
        molecules_with_energy.sort(key=lambda x: x[1])

        # Get the top 10 lowest energy molecules, or fewer if not available
        lowest_10_molecules = [item[0] for item in molecules_with_energy[:10]]

        filename, _ = QFileDialog.getSaveFileName(self, "Save Lowest Energy Molecules", "lowest_energy_molecules.sdf", "SDF Files (*.sdf)")
        
        if filename:
            try:
                writer = SDWriter(filename)
                for mol in lowest_10_molecules:
                    writer.write(mol)
                writer.close()
                
                export_message = f"SUCCESS: Exported {len(lowest_10_molecules)} lowest energy molecules!\n\n"
                export_message += f"File saved to: {filename}"
                QMessageBox.information(self, "Export Successful", export_message)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export SDF file.\n\n{str(e)}")
        else:
            QMessageBox.information(self, "Export Cancelled", "Export cancelled by user.")