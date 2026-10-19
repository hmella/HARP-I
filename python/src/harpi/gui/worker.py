"""Running long computations off the UI thread.

A QThread is used rather than a separate process, since the arrays are
hundreds of megabytes and pickling them would cost more than the compute.
The FFTs, the BLAS solves and the unwrapper all release the GIL.

Cancellation is checked at frame boundaries. A single dense solve cannot be
interrupted, and in 3-D one frame can take a minute or more.
"""

from __future__ import annotations

import threading
import traceback

from PySide6.QtCore import QObject, QThread, Signal

__all__ = ["Worker", "run_in_background"]


class Worker(QObject):
    """Runs ``fn(progress=..., cancelled=...)`` and reports back."""

    progress = Signal(int, int, str)     # done, total, stage
    finished = Signal(object)
    failed = Signal(str, str)            # message, traceback

    def __init__(self, fn):
        super().__init__()
        self._fn = fn
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    @property
    def cancelled(self) -> bool:
        return self._cancel.is_set()

    def run(self):
        try:
            result = self._fn(
                progress=lambda i, n, s: self.progress.emit(int(i), int(n), str(s)),
                cancelled=lambda: self._cancel.is_set(),
            )
        except Exception as exc:  # noqa: BLE001 - reported to the user
            self.failed.emit(str(exc), traceback.format_exc())
        else:
            self.finished.emit(result)


def run_in_background(parent, fn, on_done, on_fail=None, on_progress=None):
    """Start ``fn`` on a worker thread. Returns ``(thread, worker)``.

    The caller must keep both alive until ``finished`` fires.
    """
    thread = QThread(parent)
    worker = Worker(fn)
    worker.moveToThread(thread)

    thread.started.connect(worker.run)
    if on_progress is not None:
        worker.progress.connect(on_progress)

    def _done(result):
        thread.quit()
        on_done(result)

    def _fail(msg, tb):
        thread.quit()
        if on_fail is not None:
            on_fail(msg, tb)

    worker.finished.connect(_done)
    worker.failed.connect(_fail)
    thread.finished.connect(worker.deleteLater)
    thread.start()
    return thread, worker
