from PySide6.QtCore import  Signal, QObject, Slot
from aqme.qprep import qprep # //i think this was giving some grief when it was imported in run so i moved it here and the memory allocation error ceased for now\\


class QPrepWorker(QObject):
    resultReady = Signal(str)
    errorReady = Signal(str)
    finished = Signal()
    def __init__(self, params):
        super().__init__()
        self.params = params


    # Changed your code slightly, it now works for import/export anywhere but if QREP isnt meant to work like this then can prob just revert.
    @Slot()
    def run(self):
        try:
            import shutil, os
            input_path = self.params['input_file']
            testfiles_dir = os.path.join(os.path.dirname(__file__), 'testfiles')
            if not os.path.exists(testfiles_dir):
                os.makedirs(testfiles_dir)
            if not os.path.abspath(input_path).startswith(os.path.abspath(testfiles_dir)):
                basename = os.path.basename(input_path)
                dest_path = os.path.join(testfiles_dir, basename)
                shutil.copy2(input_path, dest_path)
                self.params['input_file'] = dest_path

            print("***********************")
            print("Debug: Running qprep with parameters:", self.params)
            print("***********************")

            qprep(files=self.params['input_file'],
                program=self.params['software'],
                qm_input=self.params["keywords"] + ' ' + self.params['functional'] + " " + self.params['basis'] + " " + self.params['solvent_block'],
                mem = f"{self.params['mem']}GB",
                nprocs=self.params['nprocs'],
            )

            # Move generated files to user-selected output directory if specified
            output_dir = self.params.get('output_dir', '').strip()
            qcalc_dir = os.path.join(os.path.dirname(__file__), 'QCALC')
            if output_dir and os.path.isdir(output_dir) and os.path.abspath(output_dir) != os.path.abspath(qcalc_dir):
                for fname in os.listdir(qcalc_dir):
                    src = os.path.join(qcalc_dir, fname)
                    dst = os.path.join(output_dir, fname)
                    if os.path.isfile(src):
                        import shutil
                        shutil.move(src, dst)

            self.finished.emit()
            self.resultReady.emit("Input files generated successfully.")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.errorReady.emit(f"Error: {str(e)}\n\nTraceback:\n{tb}")