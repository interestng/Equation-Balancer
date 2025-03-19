# Chemical Equation Balancer

Chemical Equation Balancer is a Python GUI application that automatically balances chemical equations and visually displays the process. It uses Tkinter for the user interface and PyCairo for drawing a custom table that shows element counts. When the balanced coefficient differs from the original value, the updated number is displayed in sky blue with an arrow indicating the change.

## Features

- **Input Interface:** Enter unbalanced chemical equations (reactants and products) via a simple GUI.
- **Automatic Balancing:** Uses Sympy's linear algebra capabilities to compute the appropriate coefficients.
- **Visual Display:** A custom-drawn table displays the element counts for each compound. Updated numbers are highlighted in sky blue.
- **Standalone Executable:** Easily build a standalone `.exe` using PyInstaller.

## Dependencies

- **Python 3.10+**
- **Tkinter** (bundled with Python)
- **Sympy** for balancing equations
- **PyCairo** for advanced drawing
- **Pillow** (PIL) for image processing and integration with Tkinter
- **PyInstaller** (optional, for creating a standalone executable)

You can install the required libraries using pip:

```bash
pip install sympy pycairo Pillow
