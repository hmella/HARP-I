# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""PySide6 interface. Requires the ``gui`` extra: ``pip install 'harpi[gui]'``."""

__all__ = ["main", "MainWindow"]


def __getattr__(name):
    if name in __all__:
        from .app import MainWindow, main

        return {"main": main, "MainWindow": MainWindow}[name]
    raise AttributeError(name)
