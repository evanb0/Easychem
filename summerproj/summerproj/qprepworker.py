# Displays a message if the file is submitted successfully, and also displays an error if not.

from PySide6.QtCore import QThread, Signal

class QPrepWorker(QThread):
    resultReady = Signal(str)
    errorReady = Signal(str)
    def __init__(self, params):
        super().__init__()
        self.params = params
    def run(self):
        try:
            from aqme.qprep import qprep
            qprep(**self.params)
            self.resultReady.emit("Input files generated successfully.")
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.errorReady.emit(f"Error: {str(e)}\n\nTraceback:\n{tb}")
