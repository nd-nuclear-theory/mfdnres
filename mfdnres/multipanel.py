"""multipanel.py

    Provide tools for generating multipanel arrays of axes.

    Recovers functionality from SciDraw Figure and Multipanel combination.

    Reference:

        M. A. Caprio, Comput. Phys. Commun. 171, 107 (2005).

    Mark A. Caprio
    University of Notre Dame

    - 05/23/21 (mac): Created.

"""

import itertools

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

################################################################
# helper functions
################################################################

# TODO provide panel lettering helper function

def panel_index(dimensions,panel_indices,direction="x",shift=0):
    """ Recover counting index for panel.

    Arguments:

        dimensions (tuple of int): array dimensions (num_rows,num_cols) of multipanel grid

        panel_indices (tuple of int): array indices (row,col) in multipanel grid

        direction (str, optional): most-rapidly varying index ("x" or "y") for iteration

        shift (int, optional): shift to panel indexing

    Returns:

       panel_index (int): panel index

    """
    row, col = panel_indices
    if direction=="x":
        panel_index = row*dimensions[1]+col
    elif direction=="y":
        panel_index = col*dimensions[0]+row
    else:
        raise(ValueError("unrecognized direction {}".format(direction)))

    panel_index += shift
    return panel_index

def panel_letter(
        dimensions,panel_indices,
        direction="x",
        base = "a",
        shift = 0,
        delimiters = ("(",")"),
):
    """ Generate letter for panel.

    Underlying calculation is carried out by panel_index.

    Recovers functionality from SciDraw Multipanel and PanelLetter.

    TODO: Add analogs to PanelLetterOrigin and PanelLetterDimensions, via panel_index.

    Arguments:

        dimensions (tuple of int): array dimensions (num_rows,num_cols) of multipanel grid

        panel_indices (tuple of int): array indices (row,col) in multipanel grid

        base (str, optional): base panel letter (must be single character)

        direction (str, optional): most-rapidly varying index ("x" or "y") for iteration

        shift (int, optional): shift to panel indexing

 (tuple of str, optional): (left,right) delimeter strings for panel letter

    Returns:

       panel_letter (str): panel letter

    """

    panel_letter = chr(ord(base)+panel_index(dimensions,panel_indices,direction=direction,shift=shift))
    decorated_panel_letter = "{delimiters[0]}{panel_letter}{delimiters[1]}".format(panel_letter=panel_letter,delimiters=delimiters)
    return decorated_panel_letter

def grid_iterator(dimensions):
    """Provide iterator over rows and columns (or n-dimensional Cartesian grid).

    Example:

        >>> g=grid_iterator((2,3))
        >>> list(g)
        [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    Arguments:

        dimensions (tuple): (num_rows, num_cols, ...)

    Returns:

        (iterable): iterator over (row, col, ...) in range 0..num_rows, etc.

    """
    return itertools.product(*map(range,dimensions))

def multipanel_fig_gs(
        dimensions=(1,1),
        canvas_size=None, panel_size=(3.,2.),
        canvas_margin=((1.,1.),(1.,1.)),
        x_panel_sizes=None, y_panel_sizes=None,
        x_panel_gaps=0., y_panel_gaps=0.
):
    """Generate figure and gridspec such as to provide user-controlled panel sizes.

    Recovers functionality from SciDraw Figure and Multipanel combination.

    Note: Recovering support for different relative gaps between different
    panels would require abandoning mpl's Subplots/GridSpec paradigm.  Would
    have to directly construct Axes object on appropriate rectangle in figure
    coordinates (then perhaps register the axes with fig.add_subplot(ax)?).

    Arguments:
        dimensions (tuple of int, optional): array dimensions (num_rows,num_cols) of multipanel grid
        panel_size (tuple of float, optional): physical dimensions (px,py) of single axes
        canvas_size (tuple of float, optional): physical dimensions (lx,ly) of the main canvas area, deduced from panel_size if None
        canvas_margin (tuple of tuple of float, optional): additional margin ((dx0,dx1),(dy0,dy1)) around the main canvas area
        x_panel_sizes, y_panel_sizes (list of float, optional): relative panel sizes
        x_panel_gaps, y_panel_gaps (float, optional): relative gap sizes

    Returns:
        fig (mpl.figure.Figure): figure
        gs (mpl.gridspec.GridSpec): gridspec

    """

    if canvas_size is None:
        canvas_size = np.array(panel_size)*np.array(tuple(reversed(dimensions)))
        canvas_size += np.array(panel_size)*np.array((x_panel_gaps,y_panel_gaps))*(np.array(tuple(reversed(dimensions)))-1)
    x_main_length, y_main_length = canvas_size
    x_margin, y_margin = canvas_margin
    x_margin_0, x_margin_1 = x_margin
    y_margin_0, y_margin_1 = y_margin
    x_total_length = x_main_length+x_margin_0+x_margin_1
    y_total_length = y_main_length+y_margin_0+y_margin_1
    
    padded_canvas_size = (x_total_length, y_total_length)
    fig = plt.figure(figsize=padded_canvas_size)
    gs = fig.add_gridspec(
        *dimensions,
        left=x_margin_0/x_total_length, right=1-x_margin_1/x_total_length,
        bottom=y_margin_0/y_total_length, top=1-y_margin_1/y_total_length,
        width_ratios=x_panel_sizes, height_ratios=y_panel_sizes,
        wspace=x_panel_gaps, hspace=y_panel_gaps
    )

    return fig, gs

def main():    
    pass

if __name__ == "__main__":
    main()
