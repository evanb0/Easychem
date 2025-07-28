from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QDialog,
    QMessageBox, QFileDialog, QGroupBox
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
    def __init__(self, mol_block, style, coord_system, distance_labels, conformer_label=""): 
        super().__init__()
        # Set the window title dynamically based on the conformer label (e.g u are viewing conf 1/ or conf 9 etc.)
        self.setWindowTitle(f"3D Molecule Viewer - {conformer_label}")
        self.showMaximized()  # Open maximized (windowed fullscreen)
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

        # Add coordinate system if requested
        if coord_system in ['XYZ Axes', 'Both']:
            # Add XYZ coordinate axes
            # X-axis (red)
            viewer.addLine({'start': {'x': -5, 'y': 0, 'z': 0}, 'end': {'x': 5, 'y': 0, 'z': 0}, 'color': 'red', 'linewidth': 3})
            viewer.addLabel('X', {'position': {'x': 5.5, 'y': 0, 'z': 0}, 'fontColor': 'red', 'fontSize': 14, 'fontOpacity': 0.8})
            
            # Y-axis (green)
            viewer.addLine({'start': {'x': 0, 'y': -5, 'z': 0}, 'end': {'x': 0, 'y': 5, 'z': 0}, 'color': 'green', 'linewidth': 3})
            viewer.addLabel('Y', {'position': {'x': 0, 'y': 5.5, 'z': 0}, 'fontColor': 'green', 'fontSize': 14, 'fontOpacity': 0.8})
            
            # Z-axis (blue)
            viewer.addLine({'start': {'x': 0, 'y': 0, 'z': -5}, 'end': {'x': 0, 'y': 0, 'z': 5}, 'color': 'blue', 'linewidth': 3})
            viewer.addLabel('Z', {'position': {'x': 0, 'y': 0, 'z': 5.5}, 'fontColor': 'blue', 'fontSize': 14, 'fontOpacity': 0.8})

        if coord_system in ['Grid Lines', 'Both']:
            # Add grid lines for better spatial reference
            grid_range = 10
            grid_step = 2
            grid_color = 'gray'
            grid_opacity = 0.3
            
            # XY plane grid
            for i in range(-grid_range, grid_range + 1, grid_step):
                if i != 0:  # Don't overlap with main axes
                    # Lines parallel to X-axis
                    viewer.addLine({
                        'start': {'x': -grid_range, 'y': i, 'z': 0}, 
                        'end': {'x': grid_range, 'y': i, 'z': 0}, 
                        'color': grid_color, 'linewidth': 1, 'opacity': grid_opacity
                    })
                    # Lines parallel to Y-axis
                    viewer.addLine({
                        'start': {'x': i, 'y': -grid_range, 'z': 0}, 
                        'end': {'x': i, 'y': grid_range, 'z': 0}, 
                        'color': grid_color, 'linewidth': 1, 'opacity': grid_opacity
                    })

        # Add distance labels if requested
        if distance_labels != 'None':
            self.add_distance_labels(viewer, mol_block, distance_labels)

        # Zoom to fit the molecule in the view
        viewer.zoomTo()
        html = viewer._make_html()
        self.web_view.setHtml(html)

    def add_distance_labels(self, viewer, mol_block, distance_type):
        """Add distance labels between atoms"""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            import math
            
            # Parse the mol block to get atom positions
            mol = Chem.MolFromMolBlock(mol_block)
            if not mol or not mol.GetNumConformers():
                return
                
            conf = mol.GetConformer()
            atoms = []
            
            # Get atom positions and symbols
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                atoms.append({
                    'idx': i,
                    'symbol': atom.GetSymbol(),
                    'pos': {'x': pos.x, 'y': pos.y, 'z': pos.z}
                })
            
            if distance_type == 'Bond Lengths':
                # Only show distances for bonded atoms
                for bond in mol.GetBonds():
                    idx1 = bond.GetBeginAtomIdx()
                    idx2 = bond.GetEndAtomIdx()
                    
                    pos1 = atoms[idx1]['pos']
                    pos2 = atoms[idx2]['pos']
                    
                    # Calculate distance
                    dx = pos2['x'] - pos1['x']
                    dy = pos2['y'] - pos1['y']
                    dz = pos2['z'] - pos1['z']
                    distance = math.sqrt(dx*dx + dy*dy + dz*dz)
                    
                    # Midpoint for label
                    mid_x = (pos1['x'] + pos2['x']) / 2
                    mid_y = (pos1['y'] + pos2['y']) / 2
                    mid_z = (pos1['z'] + pos2['z']) / 2
                    
                    # Add distance label
                    label_text = f"{distance:.2f}Å"
                    viewer.addLabel(label_text, {
                        'position': {'x': mid_x, 'y': mid_y, 'z': mid_z},
                        'fontColor': 'black',
                        'fontSize': 10,
                        'fontOpacity': 0.8,
                        'backgroundColor': 'white',
                        'backgroundOpacity': 0.7
                    })
                    
            elif distance_type == 'All Distances': # Not recommended, but I added it in case it may be useful.
                max_distance = 5.0  # Only show distances up to 5 Angstroms to avoid clutter
                
                for i in range(len(atoms)):
                    for j in range(i + 1, len(atoms)):
                        pos1 = atoms[i]['pos']
                        pos2 = atoms[j]['pos']
                        
                        # Calculate distance
                        dx = pos2['x'] - pos1['x']
                        dy = pos2['y'] - pos1['y']
                        dz = pos2['z'] - pos1['z']
                        distance = math.sqrt(dx*dx + dy*dy + dz*dz)
                        
                        # Only show distances within reasonable range
                        if distance <= max_distance:
                            # Midpoint for label
                            mid_x = (pos1['x'] + pos2['x']) / 2
                            mid_y = (pos1['y'] + pos2['y']) / 2
                            mid_z = (pos1['z'] + pos2['z']) / 2
                            
                            # Add distance label
                            label_text = f"{distance:.2f}Å"
                            viewer.addLabel(label_text, {
                                'position': {'x': mid_x, 'y': mid_y, 'z': mid_z},
                                'fontColor': 'purple',
                                'fontSize': 8,
                                'fontOpacity': 0.6,
                                'backgroundColor': 'lightyellow',
                                'backgroundOpacity': 0.5
                            })
                            
        except Exception as e:
            print(f"Error adding distance labels: {e}")
            # Continue without distance labels if there's an error


class MoleculeViewer(QWidget):
    """
    Main widget for RDKit Molecule Viewer application.
    Allows users to input SMILES, generate conformers, view them in 3D,
    and export them as SDF files.
    Updated to match the UI style of ParameterPanel and FilePanel.
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RDKit Molecule Viewer")
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        # SMILES Input Section
        smiles_group = QGroupBox("SMILES Input")
        smiles_group.setStyleSheet("font-weight: bold;")
        smiles_layout = QVBoxLayout()
        
        smiles_row = QHBoxLayout()
        smiles_row.addWidget(QLabel("SMILES:"))
        self.smiles_input = QLineEdit()
        self.smiles_input.setPlaceholderText("Enter or load SMILES string...")
        smiles_row.addWidget(self.smiles_input)
        self.load_smiles_button = QPushButton("Browse")
        self.load_smiles_button.clicked.connect(self.load_smiles_from_file)
        smiles_row.addWidget(self.load_smiles_button)
        smiles_layout.addLayout(smiles_row)
        smiles_group.setLayout(smiles_layout)

        # Display Options Section
        display_group = QGroupBox("Display Options")
        display_group.setStyleSheet("font-weight: bold;")
        display_layout = QVBoxLayout()

        style_row = QHBoxLayout()
        style_row.addWidget(QLabel("Display Style:"))
        self.style_selector = QComboBox()
        self.style_selector.addItems(['Stick', 'Ball and Stick', 'Surface'])
        style_row.addWidget(self.style_selector)
        display_layout.addLayout(style_row)

        conformer_row = QHBoxLayout()
        conformer_row.addWidget(QLabel("Conformer:"))
        self.conformer_selector = QComboBox()
        self.conformer_selector.currentIndexChanged.connect(self.on_conformer_change)
        conformer_row.addWidget(self.conformer_selector)
        display_layout.addLayout(conformer_row)

        display_group.setLayout(display_layout)

        # Buttons Section
        button_layout = QHBoxLayout()
        
        self.draw_button = QPushButton("Generate Conformers")
        self.draw_button.setMinimumHeight(40)
        self.draw_button.setStyleSheet("font-size: 12px; font-weight: bold; background-color: #4287f5; color: white;")
        self.draw_button.clicked.connect(self.on_draw_click)
        
        self.view_button = QPushButton("Show 3D Molecule")
        self.view_button.setMinimumHeight(30)
        self.view_button.setStyleSheet("font-size: 12px; font-weight: bold; background-color: #4287f5; color: white;")
        self.view_button.clicked.connect(self.show_selected_conformer)
        self.view_button.setEnabled(False)
        
        self.export_sdf_button = QPushButton("Export SDF")
        self.export_sdf_button.setMinimumHeight(30)
        self.export_sdf_button.setStyleSheet("font-size: 10px; font-weight: bold; background-color: #f5b942; color: black;")
        self.export_sdf_button.clicked.connect(self.on_export_sdf)
        self.export_sdf_button.setEnabled(False)

        button_layout.addWidget(self.view_button)
        button_layout.addWidget(self.draw_button)
        button_layout.addWidget(self.export_sdf_button)

        # Status Display Section (matching FilePanel style)
        status_group = QGroupBox("Status / Information")
        status_layout = QVBoxLayout()
        from PySide6.QtWidgets import QTextEdit
        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(200)
        self.status_text.setPlaceholderText("Ready to generate molecular conformers")
        self.status_text.setReadOnly(True)
        status_layout.addWidget(self.status_text)
        status_group.setLayout(status_layout)

        # Final Layout
        layout.addWidget(smiles_group)
        layout.addWidget(display_group)
        layout.addLayout(button_layout)
        layout.addWidget(status_group)
        layout.addStretch()
        self.setLayout(layout)

    def generate_conformers_with_energy(self, smiles, max_confs=10):
        """
        Generates 3D conformers for a given SMILES string and calculates their energies.
        Uses RDKit's EmbedMultipleConfs for conformer generation and MMFF/UFF for energy calculation.
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
        Handles the 'Generate Conformers' button click.
        Generates conformers from the SMILES input, adds them to be selected,
        and displays information in the status area.
        """
        smiles = self.smiles_input.text().strip()
        
        if not smiles:
            self.status_text.setPlainText("Error: No SMILES string provided. Please enter a SMILES string.")
            return

        self.status_text.setPlainText("Processing... Generating conformers. \n\nThis may take a moment depending on the complexity of your molecule.")
        self.draw_button.setEnabled(False)
        
        try:
            results = self.generate_conformers_with_energy(smiles)
            if not results:
                self.status_text.setPlainText("Error: Could not generate conformers for the given SMILES.\n\nPlease check that the SMILES string is valid.")
                self.draw_button.setEnabled(True)
                return

            # Sort conformers by energy (lowest energy first)
            # Handle cases where energy might be None by putting them at the end
            results_sorted = sorted(results, key=lambda x: x[2] if x[2] is not None else float('inf'))
            self.mol_conformers = results_sorted

            # Temporarily disconnect the signal to prevent double triggering
            self.conformer_selector.currentIndexChanged.disconnect()
            
            # Populate conformer selector with energies (now sorted by energy)
            self.conformer_selector.clear()
            for i, (_, _, energy) in enumerate(self.mol_conformers):
                if i == 0 and energy is not None:
                    label = f"Conformer {i+1} (MOST STABLE)"
                else:
                    label = f"Conformer {i+1}"
                
                if energy is not None:
                    label += f" (energy: {energy:.2f} kcal/mol)"
                else:
                    label += " (energy: N/A)"
                self.conformer_selector.addItem(label)

            # Set the first (most stable) conformer as selected
            self.conformer_selector.setCurrentIndex(0)

            # Reconnect the signal
            self.conformer_selector.currentIndexChanged.connect(self.on_conformer_change)

            # Update status with success information, highlighting most stable conformer
            lowest_energy = self.mol_conformers[0][2] if self.mol_conformers[0][2] is not None else "N/A"
            success_message = f"SUCCESS: Generated {len(results)} conformers!\n\n"
            success_message += f"SMILES: {smiles}\n"
            success_message += f"Display Style: {self.style_selector.currentText()}\n"
            success_message += f"Number of conformers: {len(results)}\n"
            success_message += f"📌 Most stable conformer energy: {lowest_energy}"
            if lowest_energy != "N/A":
                success_message += " kcal/mol"
            success_message += "\n\n"
            success_message += "Conformers are sorted by energy (most stable first).\n"
            success_message += "Click 'Show 3D Molecule' to visualize the selected conformer,\n"
            success_message += "or click 'Export SDF' to save it."
            
            self.status_text.setPlainText(success_message)
            
            # Trigger the conformer change event to update status with first conformer info
            self.on_conformer_change(0)
            
            # Enable other buttons
            self.view_button.setEnabled(True)
            self.export_sdf_button.setEnabled(True)
            
        except Exception as e:
            self.status_text.setPlainText(f"Error: Failed to generate conformers.\n\n{str(e)}")
        
        finally:
            self.draw_button.setEnabled(True)

    def on_conformer_change(self, index):
        """
        Handles the change in the conformer selector.
        Updates the status display with information about the selected conformer.
        """
        if index < 0 or not hasattr(self, 'mol_conformers') or not self.mol_conformers:
            return
        
        # Update status to show selected conformer info
        _, _, energy = self.mol_conformers[index]
        conformer_info = f"Selected: {self.conformer_selector.currentText()}\n\n"
        conformer_info += f"SMILES: {self.smiles_input.text().strip()}\n"
        conformer_info += f"Display Style: {self.style_selector.currentText()}\n"
        if energy is not None:
            conformer_info += f"Energy: {energy:.2f} kcal/mol\n\n"
        else:
            conformer_info += "Energy: N/A\n\n"
        conformer_info += "Click 'Show 3D Molecule' to visualize this conformer,\n"
        conformer_info += "or click 'Export SDF' to save it as an SDF file."
        
        self.status_text.setPlainText(conformer_info)

    def show_selected_conformer(self):
        """
        Displays the currently selected conformer in a new 3D viewer window.
        """
        if not hasattr(self, 'mol_conformers') or not self.mol_conformers:
            self.status_text.setPlainText("Error: No conformers available. Please generate conformers first.")
            return

        index = self.conformer_selector.currentIndex()
        if index < 0 or index >= len(self.mol_conformers):
            self.status_text.setPlainText("Error: No conformer selected. Please select a conformer first.")
            return

        mol, conf_id, _ = self.mol_conformers[index]
        style = self.style_selector.currentText()
        mol_block = Chem.MolToMolBlock(mol, confId=conf_id)
        
        conformer_label = self.conformer_selector.currentText()
        
        coord_system = 'XYZ Axes'  # Could also use 'XYZ Axes', 'Grid Lines', or 'Both'
        distance_labels = 'Bond Lengths'  # Could also use 'Bond Lengths' or 'All Distances'

        viewer_window = Molecule3DWindow(mol_block, style, coord_system, distance_labels, conformer_label)
        viewer_window.exec_()


    def on_export_sdf(self):
        """
        Handles the 'Export SDF' button click.
        Exports the currently selected conformer as an SDF file to a user-specified location.
        """
        if not hasattr(self, 'mol_conformers') or not self.mol_conformers:
            self.status_text.setPlainText("Error: No conformers available. Please generate conformers first.")
            return

        # Get the currently selected conformer
        current_index = self.conformer_selector.currentIndex()
        if current_index == -1: # No conformer selected
            self.status_text.setPlainText("Error: No conformer selected. Please select a conformer to export.")
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
                
                # Update status with export success
                export_message = f"SUCCESS: SDF file exported!\n\n"
                export_message += f"File saved to: {filename}\n"
                export_message += f"Conformer: {self.conformer_selector.currentText()}\n"
                export_message += f"SMILES: {self.smiles_input.text().strip()}"
                
                self.status_text.setPlainText(export_message)
            except Exception as e:
                self.status_text.setPlainText(f"Error: Failed to export SDF file.\n\n{str(e)}")

    def load_smiles_from_file(self):
        """
        Handles the 'Browse' button click.
        Opens a file dialog, reads an SDF or XYZ file, extracts the SMILES string,
        and adds the SMILES extracted to smiles input box. 
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
                    
                    # Update status display
                    file_basename = filename.split('/')[-1]  # Get just the filename
                    load_message = f"SMILES loaded from file: {file_basename}\n\n"
                    load_message += f"SMILES: {smiles}\n"
                    load_message += f"File path: {filename}\n\n"
                    load_message += "Click 'Generate Conformers' to create 3D conformers for this molecule."
                    
                    self.status_text.setPlainText(load_message)
                else:
                    self.status_text.setPlainText("Error: Could not load molecule from the selected file.\n\nEnsure it's a valid SDF or XYZ file.")
            except Exception as e:
                self.status_text.setPlainText(f"Error: Failed to load file.\n\n{str(e)}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec_())
