"""
This file is our main area, and where most issues occured. Most things that went wrong happened here, and I imagine that will continue to be the case with further
Updates of the code. So Far, I have only ran codes with ORCA in mind. If you do end up reading the code over the weekend, and want something to do, I would suggest
Attempting to run a GAUSSIAN script. It will probably return an error, and I will probably have to work on this, but idrk gaussian that well so on monday I think
I'm gonna look at our github page to view it anyway, because I kinda just tried following the documentation and hoped for the best.

Another thing, I have NOT USED RDKIT at all, Python can read files pretty well, however I do realise this is kind of a 'brute-force' way to get the code to
work properly for specific file formats, so if u think this is bad I will probably have to look up how stuff in RDKIT works.


// I suggested RDKit because I thought extracting the SDF details to mol objects might be easier to work with later on when we want to have more interactive visualisation \\
// AQME's QPREP takes in SDF files anyhow so before we do have more to do we don't really need to worry - we just need to pass the file name as a parameter to the aqme thread \\
// I also added some test sdf files to the repo, one has a few molecules and the other has a few thousand, the latter took 15.5 seconds when i ran it so not bad\\
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
from qprepworker import QPrepWorker

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
        # self.make_btn.clicked.connect(self.generate_input_file)
        self.make_btn.clicked.connect(self.run_qprep)

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

    # def test_qprep(self):
    #     from aqme.qprep import qprep
    #     print("TEST: checking if qprep works in the main thread (spoiler: it does)")
    #     params = self.collect_input_params()
    #     qprep(
    #         files=params['input_file'],
    #         program=params['software'],
    #         qm_input=params['functional'] + " " + params['basis'],
    #         mem=f"{params['mem']}GB",
    #         nprocs=params['nprocs'],
    #     )

    def run_qprep(self):
        """Run the QPrep worker thread to generate input files."""
        # since we have to pass the params to the worker:
        params = self.collect_input_params()
        if not params['input_file']:
            self.status_text.setPlainText("Error: No input file selected. Please select a file.")
            return
        self.thread = QThread()
        self.worker = QPrepWorker(params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        # print("Debug: Starting QPrep worker thread with parameters:", params)
        self.thread.start()
        self.make_btn.setEnabled(False)
        self.thread.finished.connect(lambda: self.make_btn.setEnabled(True))
        

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
        params['charge'] = 0 # //unsure about this one, why is it being hardcoded? (i might be blind here)\\
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

        # // whoever wrote qprep decided to be lazy so for the whole qm_input we have to figure out the solvent block ourselves — thankfully you already did some of the work so we can just make a solvent_block string in params \\
        # // If you look in the orca 6.0 manual, table 7.29 on page 987 (1005 in the pdf) has all the implicit solvation methods that we can implement. I imagine gaussian will have something similar. The below will have to be changed eventually to make sure it's compatible w the software\\
        if params['solvent_model'] != 'None' and params['solvent'] != 'None':
            solvent_blocks = {
            'orca': {
                'PCM': f"CPCM({params['solvent']})",
                'SMD': f"SMD({params['solvent']})",
                'COSMO': f"ddCOSMO({params['solvent']})"
            },
            'gaussian': {
                'PCM': f"SCRF=(PCM,Solvent={params['solvent']})",
                'SMD': f"SCRF=(SMD,Solvent={params['solvent']})",
            }
            }
            params["solvent_block"] = solvent_blocks.get(params['software'], {}).get(params['solvent_model'], '')
        else:
            params["solvent_block"] = ''

        # // same as what you already did below\\
        job_types = {
            'Single Point': '',
            'Geometry Optimization': 'opt',
            'Frequency': 'freq',
            'Opt+Freq': 'opt freq', # // need to add this as an option \\
        }
        # now just adding the keyword to params to pass it to qm_input (i know this is giving brute force and we could absolutely make it better)
        params["keywords"] = job_types.get(params['job_type'], "")

        # Infer molecule name from file name
        params['molecule_name'] = os.path.splitext(os.path.basename(input_file))[0] if input_file else 'input'
        return params

    # // For now I have commented this out because I didn't need it to test qprep \\
    # def extract_coordinates(self, input_file):
    #     """Extract atomic coordinates from .xyz or .sdf file."""
    #     if not input_file or not os.path.isfile(input_file):
    #         return ''
    #     ext = os.path.splitext(input_file)[1].lower()
    #     coords = ''
    #     try:
    #         if ext == '.xyz':
    #             with open(input_file) as f:
    #                 lines = f.readlines()
    #             coords = ''.join(lines[2:])  # Skip the first two lines (atom count and comment)
    #         elif ext == '.sdf':
    #             with open(input_file) as f:
    #                 lines = f.readlines()
    #             counts_line = lines[3]
    #             try:
    #                 num_atoms = int(counts_line[:3])
    #             except Exception:
    #                 num_atoms = 0
    #             atom_lines = lines[4:4+num_atoms]
    #             # SDF atom line: x y z symbol ...
    #             coord_lines = []
    #             for line in atom_lines:
    #                 parts = line.split()
    #                 if len(parts) >= 4:
    #                     x, y, z, symbol = parts[0], parts[1], parts[2], parts[3]
    #                     coord_lines.append(f"{symbol} {x} {y} {z}\n")
    #             coords = ''.join(coord_lines)
    #     except Exception:
    #         coords = ''  # If any error occurs, return an empty string
    #     return coords

#     def make_orca_input(self, params, coords):
#         """Generate an ORCA-compatible input script with correct job and solvent keywords."""
#         # Always include Opt for Frequency and TD-DFT jobs -> Currently this is forced, I just added this to test if multi-processing jobs work, however
#         # I think it's also pretty common to run a geometry optimsation with your calculations anyway, so may be worth keeping.
#         full_job_script = {
#             'Single Point': [],
#             'Geometry Optimization': ['Opt'],
#             'Frequency': ['Opt', 'Freq'],
#             'TD-DFT': ['Opt', 'TDDFT']
#         }
#         job_types = full_job_script.get(params.get('job_type', 'Geometry Optimization'), ['Opt'])
#         job_kw_str = ''.join(f" {kw}" for kw in job_types)

#         # Construct solvent block depending on the solvent model
#         solvent_block = ''
#         model = params.get('solvent_model', 'None')
#         solvent = params.get('solvent', 'None')

#         if model in ('PCM', 'SMD') and solvent != 'None':
#             solvent_block = f"\n%cpcm\n  smd true\n  SMDsolvent \"{solvent}\"\nend"
#         elif model == 'COSMO' and solvent != 'None':
#             solvent_block = f"\n%cosmo\n  epsilon 999\nend" 
#             # For now, I think this is incorrect. I think I need to add a parameter to specify dielectric constant of
#             # Solvent being used depending on the model, but this can be done later

#         # Final ORCA input file content -> Right now this Forces the user into TightSCF, Apparently this is the default for geometry optimisations so rn its just
#         # Kinda a lock here to test running scripts actually works, ofcourse however, it can also be manually updated.
#         return """
#         {params['functional']} {params['basis']} TightSCF{job_kw_str} {solvent_block}
#         """
# #     f"""! {params['functional']} {params['basis']} TightSCF{job_kw_str}
# # %pal nprocs {params['nprocs']} end
# # %maxcore {params['mem']}000{solvent_block}
# # * xyz {params['charge']} {params['multiplicity']}
# # {coords}*

    # THEORETICALLY this code should work the same way as the ORCA script above, however I've never actually ran GAUSSIAN code, so it may just fail.
#     def make_gaussian_input(self, params, coords):
#         """Generate a Gaussian-compatible input script with correct job and solvent keywords."""
#         # Always include Opt for Frequency and TD-DFT jobs
#         full_job_script = {
#             'Single Point': [],
#             'Geometry Optimization': ['Opt'],
#             'Frequency': ['Opt', 'Freq'],
#             'TD-DFT': ['Opt', 'TD(NStates=10)']
#         }
#         job_types = full_job_script.get(params.get('job_type', 'Geometry Optimization'), ['Opt'])
#         job_kw_str = ''.join(f" {kw}" for kw in job_types)

#         # Construct solvent keyword block
#         model = params.get('solvent_model', 'None')
#         solvent = params.get('solvent', 'None')
#         scrf = ''
#         if model == 'PCM' and solvent != 'None':
#             scrf = f" SCRF=(PCM,Solvent={solvent})"
#         elif model == 'SMD' and solvent != 'None':
#             scrf = f" SCRF=(SMD,Solvent={solvent})"
#         elif model == 'COSMO' and solvent != 'None':
#             scrf = f" SCRF=(COSMO)"

#         # Final Gaussian input file content
#         return f"""%nprocshared={params['nprocs']}
# %mem={params['mem']}GB
# #p {params['functional']}/{params['basis']}{job_kw_str}{scrf}

# {params['molecule_name']}

# {params['charge']} {params['multiplicity']}
# {coords}
# """
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

# // Will probably have to change this slightly to work better with the worker thread. \\
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
