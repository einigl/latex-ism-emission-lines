# latex-ism-emission-lines

This code implements a conversion of emission lines formatted according to the convention of the Meudon PDR code (<http://ism.obspm.fr/pdr.html>).

## Getting started

This code is designed in a very simple way, in a single Python file `ism_lines_helpers.py`, so it can be easily copied into your project. Only built-in libraries are used so the only dependency is `Python >= 3.0`.

## Lines formatted according to the Meudon PDR standard

The emission line are supposed to be formatted according to the Meudon PDR standard `species_highlevels__lowlevels`.

For instance, `h2_v0_j3__v0_j1` represents the transition $J=3-2$ of $H_2$ at $\nu=0$.

## Usage

The main function is `line_to_latex` which returns a LaTeX printable version of a formatted emission line.

Other available functions are the following
- `molecule_and_transition`: segment the formatted emission line in order to return the formatted chemical species and energy transition,
- `molecule`: returns only the formatted chemical species,
- `transition`: returns only the formatted energy transition,
- `is_line_of`: returns whether the line belongs to a given molecule,
- `filter_molecules`: returns the sublist of a formatted emission lines list containing only lines of species contained in a given list of chemical species,
- `molecules_among_lines`: returns a list of chemical species present in a list of formatted emission lines (without duplicates),
- `molecule_to_latex`: returns a LaTeX printable version of a formatted chemical species,
- `transition_to_latex`: returns a LaTeX printable version of a formatted energy transition.
- `remove_hyperfine`: return a shorter version of the formatter line without hyperfine levels.
- `is_hyperfine`: returns whether a line corresponds to a degenerate energy level or whether two lines belong to the same hyperfine structure.

After moving the file `ism_lines_helpers.py` in a directory that is in your Python path (if it is not the case, consider using the command `sys.path.append(path_of_containing_folder)`), you can import any of the above functions. For instance
```python
from ism_lines_helpers import line_to_latex
```

## Examples

```python
line_to_latex("h2_v0_j3__v0_j1")
>>> '$H_2$ ($\\nu=0$, $J=2$ $\\to$ $\\nu=0$, $J=0$)'
```

$H_2$ ($\nu=0$, $J=2$ $\to$ $\nu=0$, $J=0$)

```python
line_to_latex("h2o_j1_ka1_kc1__j0_ka0_kc0")
>>> '$H_2O$ ($J=1$, $k_a=1$, $k_c=1$ $\\to$ $J=0$, $k_a=0$, $k_c=0$)'
```

$H_2O$ ($J=1$, $k_a=1$, $k_c=1$ $\to$ $J=0$, $k_a=0$, $k_c=0$)

```python
line_to_latex("h_el3d_j5_2__el1s_j1_2")
>>> '$H$ ($3d$, $J=\\frac{5}{2}$ $\\to$ $1s$, $J=\\frac{1}{2}$)'
```

$H$ ($3d$, $J=\frac{5}{2}$ $\to$ $1s$, $J=\frac{1}{2}$)

More examples are available in `examples.ipynb`, especially for the other functions. An HTML version of this notebook is also available as documentation.

## Customization

This code is designed so that it can be enriched with new features. You can add or modify

- the LaTeX writing of chemical species,
- the different eligible aliases of some chemical species,
- the LaTeX writing of energies.

To modify these variables, simply edit the Python file and modify the following dictionaries.

```python
# Molecular names to LaTeX
_molecules_to_latex = {
    "h": "H",
    "h2": "H_2",
    # ...
}

# Molecular names aliases
_molecules_aliases = {
    "13co": "13c_o",
    "c18o": "c_18o",
    # ...
}

# Energy to LaTeX
_energy_to_latex = {
    "j": "J={}",
    "v": "\\nu={}",
    # ...
}

# Electronic state to LaTeX
_elstate_to_latex = {
    "s": "{}s",
    "p": "{}p",
    # ...
}

# Literal to LaTeX
_literal_to_latex = {
    "pp": "p=+",
    "pm": "p=-"
}
```