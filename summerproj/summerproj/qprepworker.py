# Displays a message if the file is submitted successfully, and also displays an error if not.

from PySide6.QtCore import  Signal, QObject, Slot
from aqme.qprep import qprep # //i think this was giving some grief when it was imported in run so i moved it here and the memory allocation error ceased for now\\


class QPrepWorker(QObject):
    resultReady = Signal(str)
    errorReady = Signal(str)
    finished = Signal()
    def __init__(self, params):
        super().__init__()
        self.params = params

    @Slot()
    def run(self):
        try:
            # qprep(**self.params)
            # //QPREP will not work with just a dictionary, it needs to be in the form of keyword arguments
            # see the documentation for more info on this)\\

            print("***********************")
            print("Debug: Running qprep with parameters:", self.params)
            print("***********************")

            qprep(files=self.params['input_file'],
                program=self.params['software'],
                qm_input=self.params["keywords"] + ' ' + self.params['functional'] + " " + self.params['basis'] + " " + self.params['solvent_block'],
                mem = f"{self.params['mem']}GB",
                nprocs=self.params['nprocs'],
            )
            self.finished.emit()
            self.resultReady.emit("Input files generated successfully.")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.errorReady.emit(f"Error: {str(e)}\n\nTraceback:\n{tb}")
