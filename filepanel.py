import os
import sys
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog
)
from PySide6.QtCore import QThread, Signal
from qprepworker import QPrepWorker

FILE_FILTERS = [
    "Structured Data Files (*.sdf)",
    "XYZ (*.xyz)",
    "Gaussian Input Files (*.com *.gjf)",
]

class FilePanel(QWidget):
    """Panel to manage input file selection, output directory, and generate input scripts."""
    
    # Signal to notify when a new file is selected
    fileSelected = Signal(str)

    def __init__(self, model=None, parameter_panel=None, molecular_viewer=None):
        super().__init__()
        self.model = model
        self.parameter_panel = parameter_panel
        self.molecular_viewer = molecular_viewer
        self.setup_ui()
        self.connect_model_signals()
        
        if self.molecular_viewer:
            self.fileSelected.connect(self.molecular_viewer.load_molecules_from_file)

    def setup_ui(self):
        layout = QVBoxLayout()

        # Input File Section
        input_group = QGroupBox("Input Files")
        input_layout = QVBoxLayout()
        input_row = QHBoxLayout()
        input_row.addWidget(QLabel("File"))
        self.input_file_edit = QLineEdit()
        self.input_file_edit.setPlaceholderText("Select File")
        self.input_file_edit.setStyleSheet("QLineEdit { background-color: #edf2f4; color: #2b2d42; border: 1px solid #8d99ae; border-radius: 3px; padding: 5px; }")
        input_row.addWidget(self.input_file_edit)
        self.browse_input_btn = QPushButton("Browse")
        self.browse_input_btn.clicked.connect(self.get_filename)
        input_row.addWidget(self.browse_input_btn)
        input_layout.addLayout(input_row)
        input_group.setLayout(input_layout)

        # Output Directory Section
        output_group = QGroupBox("Output")
        output_layout = QVBoxLayout()
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Directory"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory")
        self.output_dir_edit.setStyleSheet("QLineEdit { background-color: #edf2f4; color: #2b2d42; border: 1px solid #8d99ae; border-radius: 3px; padding: 5px; }")
        output_row.addWidget(self.output_dir_edit)
        self.browse_output_btn = QPushButton("Browse")
        self.browse_output_btn.clicked.connect(self.get_output_directory)
        output_row.addWidget(self.browse_output_btn)
        output_layout.addLayout(output_row)
        output_group.setLayout(output_layout)

        # Buttons
        self.make_btn = QPushButton("Make Input Files")
        self.make_btn.setMinimumHeight(40)
        self.make_btn.clicked.connect(self.run_qprep)
        self.preview_btn = QPushButton("Preview Input")
        self.preview_btn.setMinimumHeight(30)
        self.preview_btn.clicked.connect(self.preview_input)
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.preview_btn)
        button_layout.addWidget(self.make_btn)
        self.slurm_btn = QPushButton("SLURM")
        self.slurm_btn.setMinimumHeight(30)
        self.slurm_btn.clicked.connect(self.generate_slurm_script)
        button_layout.addWidget(self.slurm_btn)

        # Status Display
        status_group = QGroupBox("Status / Preview")
        status_layout = QVBoxLayout()
        self.status_text = QTextEdit()
        self.status_text.setStyleSheet("QTextEdit { background-color: #edf2f4; color: #2b2d42; border: 1px solid #8d99ae; border-radius: 3px; padding: 5px; }")
        self.status_text.setMaximumHeight(200)
        self.status_text.setPlaceholderText("Ready to create input files")
        self.status_text.setReadOnly(True)
        status_layout.addWidget(self.status_text)
        status_group.setLayout(status_layout)

        # Final Layout
        layout.addWidget(input_group)
        layout.addWidget(output_group)
        layout.addLayout(button_layout)
        layout.addWidget(status_group)
        layout.addStretch()
        self.setLayout(layout)
        self.update_button_color()

    def connect_model_signals(self):
        if self.model:
            try:
                if hasattr(self.model, 'softwareChanged'):
                    self.model.softwareChanged.connect(self.update_button_color)
                elif hasattr(self.model, 'software_changed'):
                    self.model.software_changed.connect(self.update_button_color)
                elif hasattr(self.model, 'dataChanged'):
                    self.model.dataChanged.connect(self.update_button_color)
                elif hasattr(self.model, 'modelChanged'):
                    self.model.modelChanged.connect(self.update_button_color)
            except Exception:
                self.setup_fallback_timer()
        else:
            self.setup_fallback_timer()

    def setup_fallback_timer(self):
        from PySide6.QtCore import QTimer
        self.last_software = None
        self.timer = QTimer()
        self.timer.timeout.connect(self.check_software_change)
        self.timer.start(500)

    def check_software_change(self):
        if self.model:
            current_software = self.model.software().lower()
            if current_software != self.last_software:
                self.last_software = current_software
                self.update_button_color()

    def update_button_color(self):
        style = ""
        slurm_style = ""
        
        self.preview_btn.setStyleSheet(style)
        self.make_btn.setStyleSheet(style)
        self.slurm_btn.setStyleSheet(slurm_style)

    def preview_input(self):
        params = self.collect_input_params()
        if not params['input_file']:
            self.status_text.setPlainText("Error: No input file selected. Please select a file first.")
            self.preview_btn.setEnabled(True)
            return
        try:
            if params['software'].lower() == 'orca':
                preview_content = self.generate_orca_preview(params)
            elif params['software'].lower() == 'gaussian':
                preview_content = self.generate_gaussian_preview(params)
            else:
                preview_content = "Error: Unsupported software type."
            self.status_text.setPlainText(f"=== INPUT FILE PREVIEW ===\n\n{preview_content}\n\n{'='*50}\n\nThis is a preview of the input file that will be generated.\nClick 'Make Input Files' to create the actual files.")
        except Exception as e:
            self.status_text.setPlainText(f"Error generating preview: {str(e)}")
        self.preview_btn.setEnabled(True)

    def generate_orca_preview(self, params):
        calc_line = f"! {params['functional']} {params['basis']}"
        if params['keywords']:
            calc_line += f" {params['keywords']}"
        if params['solvent_block']:
            calc_line += f" {params['solvent_block']}"
        preview = f"{calc_line}\n\n"
        preview += f"%pal nprocs {params['nprocs']} end\n"
        preview += f"%maxcore {params['mem']}000\n\n"
        preview += f"* xyz {params['charge']} {params['multiplicity']}\n"
        preview += "[Molecular coordinates will be extracted from input file]\n*\n"
        return preview

    def generate_gaussian_preview(self, params):
        preview = f"%nprocshared={params['nprocs']}\n"
        preview += f"%mem={params['mem']}GB\n"
        route = f"#p {params['functional']}/{params['basis']}"
        if params['keywords']:
            route += f" {params['keywords']}"
        if params['solvent_block']:
            route += f" {params['solvent_block']}"
        preview += f"{route}\n\n"
        preview += f"{params['molecule_name']}\n\n"
        preview += f"{params['charge']} {params['multiplicity']}\n"
        preview += "[Molecular coordinates will be extracted from input file]\n\n"
        return preview

    def run_qprep(self):
        params = self.collect_input_params()
        if not params['input_file']:
            self.status_text.setPlainText("Error: No input file selected. Please select a file.")
            return
        self.status_text.setPlainText("Processing... Generating input files. \n\nThis may take a moment depending on the size of your input file.")
        self.thread = QThread()
        self.worker = QPrepWorker(params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.worker.resultReady.connect(self.on_qprep_success)
        self.worker.errorReady.connect(self.on_qprep_error)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()
        self.make_btn.setEnabled(False)
        self.preview_btn.setEnabled(False)
        self.thread.finished.connect(lambda: self.make_btn.setEnabled(True))
        self.thread.finished.connect(lambda: self.preview_btn.setEnabled(True))

    def on_qprep_success(self, message):
        params = self.collect_input_params()
        output_dir = params.get('output_dir', 'Current directory')
        success_message = f"SUCCESS: Input files generated successfully!\n\n"
        success_message += f"Output location: {output_dir}\n"
        success_message += f"Software: {params['software']}\n"
        success_message += f"Method: {params['functional']}/{params['basis']}\n"
        success_message += f"Memory: {params['mem']}GB\n"
        success_message += f"Processors: {params['nprocs']}\n"
        if params['solvent_model'] != 'None':
            success_message += f"🧪 Solvent: {params['solvent']} ({params['solvent_model']})\n"
        success_message += f"\n{message}"
        self.status_text.setPlainText(success_message)

    def on_qprep_error(self, error_message):
        self.status_text.setPlainText(f"ERROR: Failed to generate input files.\n\n{error_message}")

    def collect_input_params(self):
        params = {}
        if self.model:
            params['software'] = self.model.software().lower()
            params['functional'] = self.model.functional()
            params['basis'] = self.model.basisSet()
            params['nprocs'] = self.model.nprocs()
            params['mem'] = self.model.mem()
        input_file = self.input_file_edit.text().strip()
        output_dir = self.output_dir_edit.text().strip()
        params['input_file'] = input_file
        params['output_dir'] = output_dir
        
        if self.parameter_panel:
            try:
                params['charge'] = int(self.parameter_panel.get_current_charge())
            except Exception:
                params['charge'] = 0
            try:
                params['multiplicity'] = int(self.parameter_panel.get_current_multiplicity())
            except Exception:
                params['multiplicity'] = 1
        else:
            params['charge'] = 0
            params['multiplicity'] = 1
            
        if self.parameter_panel:
            params['job_type'] = self.parameter_panel.get_current_job_type()
            params['solvent_model'] = self.parameter_panel.get_current_solvent_model()
            params['solvent'] = self.parameter_panel.get_current_solvent()
        else:
            params['job_type'] = 'Geometry Optimization'
            params['solvent_model'] = 'None'
            params['solvent'] = 'None'
        if params['solvent_model'] != 'None' and params['solvent'] != 'None':
            solvent_blocks = {
                'orca': {
                    'CPCM': f"CPCM({params['solvent']})",
                    'SMD': f"SMD({params['solvent']})",
                    'COSMORS': f"COSMORS({params['solvent']})",
                    'DRACO': f"CPCM({params['solvent']}) DRACO",
                },
                'gaussian': {
                    'PCM': f"SCRF=(PCM,Solvent={params['solvent']})",
                    'SMD': f"SCRF=(SMD,Solvent={params['solvent']})",
                    'IEFPCM': f"SCRF=(IEFPCM,Solvent={params['solvent']})",
                    'CPCM': f"SCRF=(CPCM,Solvent={params['solvent']})",
                }
            }
            params["solvent_block"] = solvent_blocks.get(params['software'], {}).get(params['solvent_model'], '')
        else:
            params["solvent_block"] = ''
        job_types = {
            'Single Point': '',
            'Geometry Optimization': 'opt',
            'Frequency': 'freq',
            'Opt+Freq': 'opt freq',
        }
        params["keywords"] = job_types.get(params['job_type'], "")
        params['molecule_name'] = os.path.splitext(os.path.basename(input_file))[0] if input_file else 'input'
        return params

    def get_filename(self):
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
            
            self.fileSelected.emit(filename)

    def display_file_info(self, filename):
        try:
            file_basename = os.path.basename(filename)
            file_size = os.path.getsize(filename)
            file_ext = os.path.splitext(filename)[1].lower()
            status_text = f"Selected file: {file_basename}\n"
            status_text += f"Full path: {filename}\n"
            status_text += f"File size: {file_size:,} bytes\n"
            status_text += f"File type: {file_ext}\n"
            
            if file_ext == '.sdf':
                try:
                    from rdkit import Chem
                    supplier = Chem.SDMolSupplier(filename)
                    mol_count = len([mol for mol in supplier if mol is not None])
                    status_text += f"Molecules in file: {mol_count}\n"
                    status_text += "→ Molecules loaded in Molecular Analysis panel\n"
                except Exception:
                    status_text += "→ Could not read molecule count\n"
            
            status_text += "-" * 50 + "\n"
            status_text += "Click 'Preview Input' to see what the generated input file will look like,\n"
            status_text += "or click 'Make Input Files' to create the actual files."
            
            if file_ext == '.sdf':
                status_text += "\n\nFor SDF files: View molecules in the Molecular Analysis panel on the right."
            
            self.status_text.setPlainText(status_text)
        except Exception as e:
            error_msg = f"Error reading file: {os.path.basename(filename)}\nError: {str(e)}"
            self.status_text.setPlainText(error_msg)

    def get_output_directory(self):
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
        try:
            dir_basename = os.path.basename(directory) or directory
            current_text = self.status_text.toPlainText()
            output_info = f"\nOutput directory: {dir_basename}\n"
            output_info += f"Full path: {directory}\n"
            output_info += "-" * 50 + "\n"
            output_info += "Ready to generate input files!"
            if current_text and not current_text.startswith("Ready to create"):
                updated_text = current_text + "\n" + output_info
            else:
                updated_text = output_info
            self.status_text.setPlainText(updated_text)
        except Exception as e:
            error_msg = f"Error with output directory: {directory}\nError: {str(e)}"
            self.status_text.append(error_msg)

    def generate_slurm_script(self):
        params = self.collect_input_params()
        input_file = params.get('input_file', '').strip()
        nprocs = params.get('nprocs', 1)
        software = params.get('software', 'orca').lower()
        if not input_file or not os.path.isfile(input_file):
            self.status_text.setPlainText("Error: No valid input file selected for SLURM script generation.")
            return
        input_basename = os.path.basename(input_file)
        base_name = os.path.splitext(input_basename)[0]
        output_dir = params.get('output_dir', '').strip()
        if output_dir:
            script_dir = output_dir
        else:
            script_dir = os.path.dirname(input_file)
        if software == 'gaussian':
            script_filename = "submitgaussian.txt"
        else:
            script_filename = "submitorca.txt"
        script_path = os.path.join(script_dir, script_filename)
        if software == 'gaussian':
            inp_basename = base_name + '.inp'
            out_basename = base_name + '.inp'
            script_content = f'''#!/bin/bash --login
#SBATCH -p multicore # (or --partition) Single-node multicore
#SBATCH -n {nprocs} # (or --ntasks=) Number of cores (2--40)
#SBATCH -t 4-0
# Load g16 for the CPU type our job is running on
module load gaussian/g16c01_em64t_detectcpu
## Set up scratch dir (please do this!)
export GAUSS_SCRDIR=/scratch/$USER/gau_temp_$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR
## Say how much memory to use (4GB per core)
export GAUSS_MDEF=$((SLURM_NTASKS*4))GB
## Inform Gaussian how many cores to use
export GAUSS_PDEF=$SLURM_NTASKS
$g16root/g16/g16 < {inp_basename} > {out_basename}
'''
        else:
            inp_basename = base_name + '.inp'
            script_content = f'''#!/bin/bash --login

#SBATCH -p multicore    # (or --partition=) Submit to the AMD Genoa nodes
#SBATCH -n {nprocs}            # (or --ntasks=) Number of cores (2--168). Must match the number in your ORCA input file!
#SBATCH -t 4-0          # Wallclock timelimit (4-0 is 4 days, max permitted is 7-0)

module purge
module load apps/binapps/orca/6.0.1-avx2

$ORCA_HOME/orca {inp_basename}  > results.${{SLURM_JOB_ID}}.txt
'''
        try:
            with open(script_path, 'w') as f:
                f.write(script_content)
            self.status_text.setPlainText(f"SLURM script generated: {script_path}")
        except Exception as e:
            self.status_text.setPlainText(f"Error writing SLURM script: {str(e)}")