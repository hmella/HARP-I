"""PySide6 interface. Requires the ``gui`` extra: ``pip install 'harpi[gui]'``."""

__all__ = ["main", "MainWindow"]


def __getattr__(name):
    if name in __all__:
        from .app import MainWindow, main

        return {"main": main, "MainWindow": MainWindow}[name]
    raise AttributeError(name)
