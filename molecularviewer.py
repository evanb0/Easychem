"""
File now actually does something. This does not use QPREP, And it may not be entirely what you are after either, But it is something

I couldn't find a way to import .sdf files and get an accurate representation, so basically main process of this is as followed:
EITHER BY USER INPUT OR STRIPPED FROM A SPECIFIC .sdf file (more can be added), the smiles of the molecules is obtained, then
10 Energy conformers are given. Idk how RDKit orders these (prob on doc page and i just missed it) but the most stable is usually the first one to view.
Number can also be updated from 10, but ik u said 10 in ur example of what to add so I just went with that for now.

It works pretty well, I only had it crash on me once, and that was a real bad molecule (like 10 8C rings stacked on top of eachother)

For some molecules (if they're real bad) i don't think it will generate conformers, idk how to change this but if it looks reasonable
i think we good 

Like I said This does not use QPREP, I think it could be implemented but if you're happy with this then that is good.

Def is not perfect, zooming in is a lil wonky (usually better to view the molecule in full screen) and sometimes it opens the 
molecule window twice idk why, prob an easy fix im just blind rn 

UI also needs to be sorted, but that can be done whenever.

"""
from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QDialog,
    QMessageBox, QFileDialog
)
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtCore import QEventLoop, QUrl
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import py3Dmol
import sys
import tempfile
import io
from PIL import Image
from PySide6.QtGui import QFont, QPixmap


class Molecule3DWindow(QDialog):
    """
    Initial styling of the window, displaying it etc.
    """
    def __init__(self, mol_block, style, conformer_label=""): 
        super().__init__()
        # Set the window title dynamically based on the conformer label (e.g u are viewing conf 1/ or conf 9 etc.)
        self.setWindowTitle(f"3D Molecule Viewer - {conformer_label}")
        self.resize(800, 600)
        layout = QVBoxLayout()
        self.web_view = QWebEngineView()
        layout.addWidget(self.web_view)
        self.setLayout(layout)

        # Initialize viewer with specified dimensions and add to viewer
        viewer = py3Dmol.view(width=1000, height=1000) # random dimensions i chose idk
        viewer.addModel(mol_block, 'mol')

        # Apply the selected display style
        if style == 'Stick':
            viewer.setStyle({'stick': {}})
        elif style == 'Ball and Stick':
            viewer.setStyle({'stick': {}, 'sphere': {}})
        elif style == 'Surface':
            viewer.setStyle({'stick': {}})
            viewer.addSurface(py3Dmol.VDW, {'opacity': 0.85}) # Pretty sure this shows ED regions within the molecule
        else:
            viewer.setStyle({'stick': {}})

        # Zoom to fit the molecule in the view
        viewer.zoomTo()
        html = viewer._make_html()
        self.web_view.setHtml(html)


class MoleculeViewer(QWidget):
    """
    Main widget for RDKit Molecule Viewer application.
    Allows users to input SMILES, generate conformers, view them in 3D,
    and export them as SDF files.
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RDKit Molecule Viewer")
        self.resize(500, 300)


        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(15, 15, 15, 15)
        main_layout.setSpacing(12)

        form_layout = QVBoxLayout()
        form_layout.setSpacing(10)

        button_layout = QHBoxLayout()
        button_layout.setSpacing(10)

        export_layout = QHBoxLayout()
        export_layout.setSpacing(10)

        file_load_layout = QHBoxLayout()

        font_label = QFont("Arial", 10, QFont.Bold)

        self.smiles_label = QLabel("SMILES Input:")
        self.smiles_label.setFont(font_label)
        self.smiles_input = QLineEdit()
        self.smiles_input.setPlaceholderText("Enter or load SMILES string...")

        self.load_smiles_button = QPushButton("Load SMILES from File")
        self.load_smiles_button.clicked.connect(self.load_smiles_from_file)

        self.style_label = QLabel("Display Style:")
        self.style_label.setFont(font_label)
        self.style_selector = QComboBox()
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'Surface'])

        self.conformer_label = QLabel("Conformer:")
        self.conformer_label.setFont(font_label)
        self.conformer_selector = QComboBox()
        self.conformer_selector.currentIndexChanged.connect(self.on_conformer_change)

        self.draw_button = QPushButton("Show 3D Molecule")
        self.draw_button.setStyleSheet("padding: 8px; font-weight: bold;")
        self.draw_button.clicked.connect(self.on_draw_click)

        self.export_sdf_button = QPushButton("Export as .sdf")
        self.export_sdf_button.clicked.connect(self.on_export_sdf)

        # Assemble layouts
        form_layout.addWidget(self.smiles_label)
        form_layout.addWidget(self.smiles_input)

        file_load_layout.addWidget(self.load_smiles_button)
        form_layout.addLayout(file_load_layout)

        form_layout.addWidget(self.style_label)
        form_layout.addWidget(self.style_selector)
        form_layout.addWidget(self.conformer_label)
        form_layout.addWidget(self.conformer_selector)

        button_layout.addWidget(self.draw_button)
        export_layout.addWidget(self.export_sdf_button)

        main_layout.addLayout(form_layout)
        main_layout.addLayout(button_layout)
        main_layout.addLayout(export_layout)

        self.setLayout(main_layout)

    def generate_conformers_with_energy(self, smiles, max_confs=10):
        """
        Generates 3D conformers for a given SMILES string and calculates their energies.
        Uses RDKit's EmbedMultipleConfs for conformer generation and MMFF/UFF for energy calculation.
        IDK How accurate this is tbh, Probably worth you testing to see if it is good
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None

        mol = Chem.AddHs(mol)  # Add hydrogens to the molecule
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=max_confs)

        results = []
        for conf_id in conf_ids:
            try:
                # Try to use MMFF force field for energy calculation -> May be more accurate way than this, idk
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
                energy = ff.CalcEnergy()
            except:
                try:
                    # Fallback to UFF force field if MMFF fails -> Worse than above way (i think, but prob better for weirder molec)
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                    energy = ff.CalcEnergy()
                except:
                    energy = None  # If both fail, energy is N/A (molecule is trash)
            results.append((mol, conf_id, energy))

        return results

    def on_draw_click(self):
        """
        Handles the 'Show 3D Molecule' button click.
        Generates conformers from the SMILES input, adds them to be selected,
        and displays the first conformer in a new 3D viewer window.
        """
        smiles = self.smiles_input.text().strip()
        style = self.style_selector.currentText()

        results = self.generate_conformers_with_energy(smiles)
        if not results:
            QMessageBox.warning(self, "Invalid SMILES", "Could not generate conformers for the given SMILES.")
            return

        self.mol_conformers = results

        # Populate conformer selector with energies
        self.conformer_selector.clear()
        for i, (_, _, energy) in enumerate(self.mol_conformers):
            label = f"Conformer {i+1}"
            if energy is not None:
                label += f" (energy: {energy:.2f} kcal/mol)"
            else:
                label += " (energy: N/A)"
            self.conformer_selector.addItem(label)

        # Show first conformer in new window
        self.show_selected_conformer()

    def on_conformer_change(self, index):
        """
        Handles the change in the conformer selector.
        Displays the newly selected conformer in the 3D viewer window.
        """
        if index < 0 or not self.mol_conformers:
            return
        self.show_selected_conformer()

    def show_selected_conformer(self):
        """
        Displays the currently selected conformer in a new 3D viewer window.
        """
        index = self.conformer_selector.currentIndex()
        if index < 0 or index >= len(self.mol_conformers):
            return
        mol, conf_id, _ = self.mol_conformers[index]
        style = self.style_selector.currentText()
        mol_block = Chem.MolToMolBlock(mol, confId=conf_id)
        
        # Get the conformer label to pass to the 3D viewer window
        conformer_label = self.conformer_selector.currentText()

        viewer_window = Molecule3DWindow(mol_block, style, conformer_label) # Pass conformer_label
        viewer_window.exec_()

    def on_export_sdf(self):
        """
        Handles the 'Export as .sdf' button click.
        Exports the currently selected conformer as an SDF file to a user-specified location.
        """
        if not self.mol_conformers:
            QMessageBox.information(self, "No Molecule", "Generate a molecule first.")
            return

        # Get the currently selected conformer
        current_index = self.conformer_selector.currentIndex()
        if current_index == -1: # No conformer selected
            QMessageBox.warning(self, "No Conformer Selected", "Please select a conformer to export.")
            return

        mol, conf_id, _ = self.mol_conformers[current_index]

        # Open file dialog to get save location from user
        filename, _ = QFileDialog.getSaveFileName(self, "Save SDF File", "molecule.sdf", "SDF Files (*.sdf)")

        if filename:
            try:
                # Write the selected conformer to the chosen file
                writer = SDWriter(filename)
                writer.write(mol, confId=conf_id)
                writer.close()
                QMessageBox.information(self, "Export Successful", f"SDF file saved to: {filename}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to save SDF file: {e}")


    def load_smiles_from_file(self):
        """
        Handles the 'Load SMILES from File' button click.
        Opens a file dialog, reads an SDF or XYZ file, extracts the SMILES string,
        and adds the SMILES extracted to smiles input box. 

        WHEN I tested this, .XYZ seemed to be a lil buggy, may be worth testing but i might've just done smthing dumb idk, 
        .sdf def works though
        """
        file_filters = "Chemical Files (*.sdf *.xyz);;SDF Files (*.sdf);;XYZ Files (*.xyz)"
        filename, _ = QFileDialog.getOpenFileName(self, "Select Molecule File", "", file_filters)

        if filename:
            try:
                mol = None
                # Load molecule based on file extension
                if filename.lower().endswith('.sdf'):
                    mol = Chem.MolFromMolFile(filename)
                elif filename.lower().endswith('.xyz'):
                    mol = Chem.MolFromXYZFile(filename)

                if mol:
                    # Convert loaded molecule to SMILES and set it in the input field
                    smiles = Chem.MolToSmiles(mol)
                    self.smiles_input.setText(smiles)
                    QMessageBox.information(self, "SMILES Loaded", f"SMILES extracted and loaded from {filename}")
                else:
                    QMessageBox.warning(self, "Error Loading File", "Could not load molecule from the selected file. Ensure it's a valid SDF or XYZ file.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred: {e}")