"""
This file is our main area, and where most issues occured. Most things that went wrong happened here, and I imagine that will continue to be the case with further
Updates of the code. So Far, I have only ran codes with ORCA in mind. If you do end up reading the code over the weekend, and want something to do, I would suggest
Attempting to run a GAUSSIAN script. It will probably return an error, and I will probably have to work on this, but idrk gaussian that well so on monday I think
I'm gonna look at our github page to view it anyway, because I kinda just tried following the documentation and hoped for the best.

Another thing, I have NOT USED RDKIT at all, Python can read files pretty well, however I do realise this is kind of a 'brute-force' way to get the code to
work properly for specific file formats, so if u think this is bad I will probably have to look up how stuff in RDKIT works.
"""


# Import system and OS libraries
import sys
import os

# Import necessary PySide6 Qt widgets and core modules for GUI construction and threading
from PySide6.QtWidgets import (
    QWidget, 
    QVBoxLayout, 
    QGroupBox, 
    QHBoxLayout, 
    QLabel, 
    QLineEdit, 
    QPushButton, 
    QTextEdit,
    QFileDialog
)
from PySide6.QtCore import QThread, Signal

# Supported input file formats for selection in the file dialog.
FILE_FILTERS = [
    "Structured Data Files (*.sdf)",
    "XYZ (*.xyz)",
    "Gaussian Input Files (*.com *.gjf)",
]

class FilePanel(QWidget):
    """A panel to manage input file selection, output directory, and generate input scripts."""

    def __init__(self, model=None, parameter_panel=None):
        super().__init__()
        self.model = model  # External model containing computational parameters
        self.parameter_panel = parameter_panel  # External widget for additional input parameters
        self.setup_ui()  # Set up the GUI layout

    def setup_ui(self):
        """Set up the layout and widgets in the panel."""
        layout = QVBoxLayout()

        # Input File Section
        input_group = QGroupBox("Input Files")
        input_layout = QVBoxLayout()
        
        input_row = QHBoxLayout()
        input_row.addWidget(QLabel("File"))

        self.input_file_edit = QLineEdit()
        self.input_file_edit.setPlaceholderText("Select File")
        input_row.addWidget(self.input_file_edit)

        self.browse_input_btn = QPushButton("Browse")
        self.browse_input_btn.clicked.connect(self.get_filename)
        input_row.addWidget(self.browse_input_btn)

        input_layout.addLayout(input_row)
        input_group.setLayout(input_layout)
        input_group.setStyleSheet("font-weight: bold;")

        # Output Directory Section
        output_group = QGroupBox("Output")
        output_layout = QVBoxLayout()
        output_group.setStyleSheet("font-weight: bold;")

        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Directory"))

        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory")
        output_row.addWidget(self.output_dir_edit)

        self.browse_output_btn = QPushButton("Browse")
        self.browse_output_btn.clicked.connect(self.get_output_directory)
        output_row.addWidget(self.browse_output_btn)

        output_layout.addLayout(output_row)
        output_group.setLayout(output_layout)

        # Generate Button
        self.make_btn = QPushButton("Make Input Files")
        self.make_btn.setMinimumHeight(40)
        self.make_btn.setStyleSheet("font-size: 14px; font-weight: bold;")
        self.make_btn.clicked.connect(self.generate_input_file)

        # Status Display Section
        status_group = QGroupBox("Status")
        status_layout = QVBoxLayout()

        self.status_text = QTextEdit()
        self.status_text.setMaximumHeight(200)
        self.status_text.setPlaceholderText("Ready to create input files")
        self.status_text.setReadOnly(True)  # Ensure status box is always read-only
        status_layout.addWidget(self.status_text)

        status_group.setLayout(status_layout)

        # Final Layout of page
        layout.addWidget(input_group)
        layout.addWidget(output_group)
        layout.addWidget(self.make_btn)
        layout.addWidget(status_group)
        layout.addStretch()
        self.setLayout(layout)

    def generate_input_file(self):
        """Generate the input script file based on selected parameters."""
        try:
            # Gather parameters from UI and model
            params = self.collect_input_params()
            script, filename = self.create_input_script(params)
            output_path = os.path.join(params['output_dir'], filename)

            # Write the generated script to the specified output file
            with open(output_path, 'w') as f:
                f.write(script)

            # Display success message
            self.status_text.setPlainText(f"Input file generated successfully:\n{output_path}")
        except Exception as e:
            # Show detailed traceback if an error occurs
            import traceback
            tb = traceback.format_exc()
            self.status_text.setPlainText(f"Error generating input file: {str(e)}\n\nTraceback:\n{tb}")

    def collect_input_params(self):
        """Collect all input values from the UI and model."""
        params = {}

        # Extract values from the model if provided
        if self.model:
            params['software'] = self.model.software().lower()
            params['functional'] = self.model.functional()
            params['basis'] = self.model.basisSet()
            params['nprocs'] = self.model.nprocs()
            params['mem'] = self.model.mem()

        # Get values from the UI
        input_file = self.input_file_edit.text().strip()
        output_dir = self.output_dir_edit.text().strip()
        params['input_file'] = input_file
        params['output_dir'] = output_dir
        params['charge'] = 0
        params['multiplicity'] = 1

        # If parameter panel is available, get job type and solvent info
        if self.parameter_panel:
            params['job_type'] = self.parameter_panel.get_current_job_type()
            params['solvent_model'] = self.parameter_panel.get_current_solvent_model()
            params['solvent'] = self.parameter_panel.get_current_solvent()
        else:
            # Use default values
            params['job_type'] = 'Geometry Optimization'
            params['solvent_model'] = 'None'
            params['solvent'] = 'None'

        # Infer molecule name from file name
        params['molecule_name'] = os.path.splitext(os.path.basename(input_file))[0] if input_file else 'input'
        return params

    def create_input_script(self, params):
        """Generate the actual input script and determine filename based on software type."""
        coords = self.extract_coordinates(params['input_file'])

        # ORCA or Gaussian script generation
        if params['software'] == 'orca':
            script = self.make_orca_input(params, coords)
            filename = params['molecule_name'] + '.inp'  # ORCA uses .inp
        else:
            script = self.make_gaussian_input(params, coords)
            filename = params['molecule_name'] + '.com'  # Gaussian uses .com

        return script, filename

    def extract_coordinates(self, input_file):
        """Extract atomic coordinates from .xyz or .sdf file."""
        if not input_file or not os.path.isfile(input_file):
            return ''
        ext = os.path.splitext(input_file)[1].lower()
        coords = ''
        try:
            if ext == '.xyz':
                with open(input_file) as f:
                    lines = f.readlines()
                coords = ''.join(lines[2:])  # Skip the first two lines (atom count and comment)
            elif ext == '.sdf':
                with open(input_file) as f:
                    lines = f.readlines()
                counts_line = lines[3]
                try:
                    num_atoms = int(counts_line[:3])
                except Exception:
                    num_atoms = 0
                atom_lines = lines[4:4+num_atoms]
                # SDF atom line: x y z symbol ...
                coord_lines = []
                for line in atom_lines:
                    parts = line.split()
                    if len(parts) >= 4:
                        x, y, z, symbol = parts[0], parts[1], parts[2], parts[3]
                        coord_lines.append(f"{symbol} {x} {y} {z}\n")
                coords = ''.join(coord_lines)
        except Exception:
            coords = ''  # If any error occurs, return an empty string
        return coords

    def make_orca_input(self, params, coords):
        """Generate an ORCA-compatible input script with correct job and solvent keywords."""
        # Always include Opt for Frequency and TD-DFT jobs -> Currently this is forced, I just added this to test if multi-processing jobs work, however
        # I think it's also pretty common to run a geometry optimsation with your calculations anyway, so may be worth keeping.
        full_job_script = {
            'Single Point': [],
            'Geometry Optimization': ['Opt'],
            'Frequency': ['Opt', 'Freq'],
            'TD-DFT': ['Opt', 'TDDFT']
        }
        job_types = full_job_script.get(params.get('job_type', 'Geometry Optimization'), ['Opt'])
        job_kw_str = ''.join(f" {kw}" for kw in job_types)

        # Construct solvent block depending on the solvent model
        solvent_block = ''
        model = params.get('solvent_model', 'None')
        solvent = params.get('solvent', 'None')

        if model in ('PCM', 'SMD') and solvent != 'None':
            solvent_block = f"\n%cpcm\n  smd true\n  SMDsolvent \"{solvent}\"\nend"
        elif model == 'COSMO' and solvent != 'None':
            solvent_block = f"\n%cosmo\n  epsilon 999\nend" 
            # For now, I think this is incorrect. I think I need to add a parameter to specify dielectric constant of
            # Solvent being used depending on the model, but this can be done later

        # Final ORCA input file content -> Right now this Forces the user into TightSCF, Apparently this is the default for geometry optimisations so rn its just
        # Kinda a lock here to test running scripts actually works, ofcourse however, it can also be manually updated.
        return f"""! {params['functional']} {params['basis']} TightSCF{job_kw_str}
%pal nprocs {params['nprocs']} end
%maxcore {params['mem']}000{solvent_block}
* xyz {params['charge']} {params['multiplicity']}
{coords}*
"""
    # THEORETICALLY this code should work the same way as the ORCA script above, however I've never actually ran GAUSSIAN code, so it may just fail.
    def make_gaussian_input(self, params, coords):
        """Generate a Gaussian-compatible input script with correct job and solvent keywords."""
        # Always include Opt for Frequency and TD-DFT jobs
        full_job_script = {
            'Single Point': [],
            'Geometry Optimization': ['Opt'],
            'Frequency': ['Opt', 'Freq'],
            'TD-DFT': ['Opt', 'TD(NStates=10)']
        }
        job_types = full_job_script.get(params.get('job_type', 'Geometry Optimization'), ['Opt'])
        job_kw_str = ''.join(f" {kw}" for kw in job_types)

        # Construct solvent keyword block
        model = params.get('solvent_model', 'None')
        solvent = params.get('solvent', 'None')
        scrf = ''
        if model == 'PCM' and solvent != 'None':
            scrf = f" SCRF=(PCM,Solvent={solvent})"
        elif model == 'SMD' and solvent != 'None':
            scrf = f" SCRF=(SMD,Solvent={solvent})"
        elif model == 'COSMO' and solvent != 'None':
            scrf = f" SCRF=(COSMO)"

        # Final Gaussian input file content
        return f"""%nprocshared={params['nprocs']}
%mem={params['mem']}GB
#p {params['functional']}/{params['basis']}{job_kw_str}{scrf}

{params['molecule_name']}

{params['charge']} {params['multiplicity']}
{coords}
"""
    # All of this stuff is just pretty standard from the book.
    def get_filename(self):
        """Open file dialog to select a structure file."""
        initial_filter = FILE_FILTERS[0]
        filters = ';;'.join(FILE_FILTERS)

        filename, selected_filter = QFileDialog.getOpenFileName(
            self,
            filter=filters,
            selectedFilter=initial_filter
        )

        if filename:
            self.input_file_edit.setText(filename)
            self.display_file_info(filename)

    def display_file_info(self, filename):
        """Update the status box with selected file details."""
        try:
            file_basename = os.path.basename(filename)
            status_text = f"Selected file: {file_basename}\n"
            status_text += f"Full path: {filename}\n"
            status_text += "-" * 50 + "\n"
            self.status_text.setPlainText(status_text)
        except Exception as e:
            error_msg = f"Error reading file: {os.path.basename(filename)}\nError: {str(e)}"
            self.status_text.setPlainText(error_msg)

    def get_output_directory(self):
        """Open dialog to select output folder."""
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory",
            "",
            QFileDialog.ShowDirsOnly
        )

        if directory:
            self.output_dir_edit.setText(directory)
            self.update_output_status(directory)

    def update_output_status(self, directory):
        """Append output directory info to status box."""
        try:
            dir_basename = os.path.basename(directory) or directory
            current_text = self.status_text.toPlainText()
            output_info = f"\nOutput directory: {dir_basename}\n"
            output_info += f"Full path: {directory}\n"
            output_info += "-" * 50 + "\n"

            if current_text and not current_text.startswith("Ready to create"):
                updated_text = current_text + "\n" + output_info
            else:
                updated_text = output_info

            self.status_text.setPlainText(updated_text)
        except Exception as e:
            error_msg = f"Error with output directory: {directory}\nError: {str(e)}"
            self.status_text.append(error_msg)
