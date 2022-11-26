"""multipanel.py

    Provide tools for generating multipanel arrays of axes.

    Recovers functionality from SciDraw Figure and Multipanel combination.

    Reference:

        M. A. Caprio, Comput. Phys. Commun. 171, 107 (2005).

    Mark A. Caprio
    University of Notre Dame

    - 05/23/21 (mac): Created.
    - 05/25/21 (mac): Move in suppress_interior_labels from data.py.
    - 12/26/21 (mac): Add support for secondary axes in suppress_interior_labels.

"""

import itertools

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

################################################################
# helper functions
################################################################

def panel_index(dimensions,panel_indices,direction="horizontal",shift=0):
    """ Recover counting index for panel.

    Arguments:

        dimensions (tuple of int): array dimensions (num_rows,num_cols) of multipanel grid

        panel_indices (tuple of int): array indices (row,col) in multipanel grid

        direction (str, optional): most-rapidly varying index ("horizontal" or "vertical") for iteration

        shift (int, optional): shift to panel indexing

    Returns:

       panel_index (int): panel index

    """
    row, col = panel_indices
    if direction=="horizontal":
        panel_index = row*dimensions[1]+col
    elif direction=="vertical":
        panel_index = col*dimensions[0]+row
    else:
        raise(ValueError("unrecognized direction {}".format(direction)))

    panel_index += shift
    return panel_index

def panel_letter(
        dimensions,panel_indices,
        direction="horizontal",
        base = "a",
        shift = 0,
        delimiters = ("(",")"),
):
    """ Generate letter for panel.

    Underlying calculation is carried out by panel_index.

    Recovers functionality from SciDraw Multipanel and PanelLetter.

    TODO: Add analogs to SciDraw's PanelLetterOrigin and PanelLetterDimensions, via panel_index.

    Arguments:

        dimensions (tuple of int): array dimensions (num_rows,num_cols) of multipanel grid

        panel_indices (tuple of int): array indices (row,col) in multipanel grid

        base (str, optional): base panel letter (must be single character)

        direction (str, optional): most-rapidly varying index ("horizontal" or "vertical") for iteration

        shift (int, optional): shift to panel indexing

        delimiters (tuple of str, optional): (left,right) delimeter strings for panel letter

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

    if x_panel_sizes is None:
        x_panel_sizes = dimensions[1]*[1.]
    if y_panel_sizes is None:
        y_panel_sizes = dimensions[0]*[1.]
    total_relative_sizes = np.array([sum(x_panel_sizes),sum(y_panel_sizes)])
    if canvas_size is None:
        ## canvas_size = np.array(panel_size)*np.array(tuple(reversed(dimensions)))
        canvas_size = np.array(panel_size)*total_relative_sizes
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

################################################################
# multipanel utilities
################################################################

def suppress_interior_labels(
        ax,
        axis="both", show_axis_label=False, show_tick_labels=False,
        secondary_x_axis=None, secondary_y_axis=None
):
    """ Suppress axis and tick labels on interior axes.

    Arguments:

        ax (mpl.axes.Axes): axes object

        axis (str, optional): which axis to act on ("x", "y", or "both")

        show_axis_label (bool, optional): whether or not to still permit axis label

        show_tick_labels (bool, optional): whether or not to still permit tick labels

        secondary_x_axis, secondary_y_axis (mpl.axes._secondary_axes.SecondaryAxis, optional): secondary axes

    """

    # From ax.is_last_row()...
    #
    # In matplotlib 3.4.3:
    #
    # MatplotlibDeprecationWarning: The is_last_row
    # function was deprecated in Matplotlib 3.4 and will be removed two minor
    # releases later. Use ax.get_subplotspec().is_last_row() instead.
    #
    # But, in matplotlib 3.3.0, this new interface is not yet available:
    #
    # AttributeError: 'SubplotSpec' object has no attribute 'is_last_row'

    # TODO: Rather than annihilating the ticks and labels, simply control their
    # visibility, using ax.tick_params(right|left|...=bool,
    # labelright|labelleft|...=bool).
    
    if (axis in {"x","both"}) and (not ax.get_subplotspec().is_last_row()):
        if not show_axis_label:
            ax.set_xlabel(None)
        if not show_tick_labels:
            ax.set_xticklabels([])
    if (axis in {"x","both"}) and (not ax.get_subplotspec().is_first_row()) and (secondary_x_axis is not None):
        if not show_axis_label:
            secondary_x_axis.set_xlabel(None)
        if not show_tick_labels:
            secondary_x_axis.set_xticklabels([])
    if (axis in {"y","both"}) and (not ax.get_subplotspec().is_first_col()):
        if not show_axis_label:
            ax.set_ylabel(None)
        if not show_tick_labels:
            ax.set_yticklabels([])
    if (axis in {"y","both"}) and (not ax.get_subplotspec().is_last_col()) and (secondary_y_axis is not None):
        if not show_axis_label:
            secondary_y_axis.set_ylabel(None)
        if not show_tick_labels:
            secondary_y_axis.set_yticklabels([])

################################################################
# main
################################################################

def main():    
    pass

if __name__ == "__main__":
    main()
