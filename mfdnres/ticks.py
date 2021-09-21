"""ticks.py

    Custom tick mark generation for linear, log, and general nonlinear axes.

    Port from CustomTicks V2.1.1 (August 26, 2020) for Mathematica (originally
    MCAxes package, January 10, 2003).

    Reference:

        M. A. Caprio, Comput. Phys. Commun. 171, 107 (2005).

    Mark A. Caprio
    University of Notre Dame

    - 04/08/21 (mac): Created.
    - 09/07/21 (mac): Add set_ticks().

"""

import warnings

import matplotlib as mpl
##import matplotlib.pyplot as plt

################################################################
# helper functions
################################################################

def in_range(x,interval):
    """ Range test."""
    (a,b) = interval
    return a<=x<=b

def approx_in_range(x,interval,tolerance=1e-10):
    """ Approximate range test."""
    (a,b) = interval
    return (a-tolerance)<=x<=(b+tolerance)


# Options[LinTicks]={
#   (* tick options *)
#   ExtraTicks->{},
#   TickPreTransformation->Identity,TickPostTransformation->Identity,
#   ShowFirst->True,ShowLast->True,ShowTickLabels->True,ShowMinorTickLabels->False,ShowMinorTicks->True,
#   TickLabelStart->0,TickLabelStep->1,TickRange->{-Infinity,Infinity},TickLabelRange->{-Infinity,Infinity},TickLabelFunction->Automatic,
#   DecimalDigits->Automatic,
#   MajorTickLength->0.010,MinorTickLength->0.005,MinorTickLabelStyle->None,TickLengthScale->1,TickDirection->In,TickReverse->False,
#   MajorTickStyle->{},MinorTickStyle->{},MinorTickIndexRange->{1,Infinity},MinorTickIndexTransformation->Identity,
#   TickTest->(True&),TickLabelTest->(True&),
# 
#   (* pass through to FixedPointForm *)
#   NumberSigns->Automatic,  (* suppress space before positive numbers *)
#   NumberPoint->".",(* SignPadding is pretty irrelevant to ticks, so omit... *)
# 
#   (* debugging *)
#   Debug->False
# };

# LinTicks[RawMajorCoordList_List,RawMinorCoordList_List,Opts___]:=
#   Module[
#     {
#       FullOpts= Flatten[{Opts,Options[LinTicks]}],
#       MajorCoordList,
#       LabeledCoordList,
#       MinorCoordList,
#       MajorTickList,
#       MinorTickList,
#       UsedTickLabelFunction,
#       DefaultTickLabelFunction,
#       TickValue,
#       TickPosition,
#       TickLabel,
#       TickLength,
#       TickStyle,
#       TickLabelStyleWrapper,
#       i
#     },
# 
#     (* make major ticks *)
#     MajorCoordList=Select[
#       (TickPreTransformation/.FullOpts)/@Union[RawMajorCoordList,ExtraTicks/.FullOpts],
#       (ApproxInRange[(TickRange/.FullOpts),#]&&(TickTest/.FullOpts)[#])&
# 		   ];
#     LabeledCoordList=Flatten[Table[
#       TickValue=MajorCoordList[[i]];
#       If[
# 	(ShowTickLabels/.FullOpts)
# 	&&ApproxInRange[(TickLabelRange/.FullOpts),TickValue]
# 	&&(Mod[i-1,(TickLabelStep/.FullOpts)]==(TickLabelStart/.FullOpts))
# 	&&((i!=1)||(ShowFirst/.FullOpts))
# 	&&((i!=Length[MajorCoordList])||(ShowLast/.FullOpts))
# 	&&(TickLabelTest/.FullOpts)[TickValue],
# 	TickValue,
# 	{}
#       ],
#       {i,1,Length[MajorCoordList]}
# 		     ]];
#     DefaultTickLabelFunction=FixedTickFunction[
#       LabeledCoordList,
#       DecimalDigits/.FullOpts,
#       FilterRules[FullOpts,Options[FixedPointForm]]
# 			     ];
#     UsedTickLabelFunction=Switch[
#       (TickLabelFunction/.FullOpts),
#       Automatic,(#2&),
#       _,(TickLabelFunction/.FullOpts)
# 			  ];
#     TickLength=ResolveTickLength[
#       (MajorTickLength/.FullOpts),(TickLengthScale/.FullOpts),
#       (TickDirection/.FullOpts),(TickReverse/.FullOpts)
# 	       ];
#     TickStyle=(MajorTickStyle/.FullOpts);
#     MajorTickList=Table[
# 
#       (* calculate tick value *)
#       TickValue=MajorCoordList[[i]];
# 
#       (* calculate coordinate for drawing tick *)
#       TickPosition=(TickPostTransformation/.FullOpts)[TickValue];
# 
#       (* construct label, or null string if it should be suppressed -- if: tick labels shown,  tick is in TickLabelRange, in designated modular cycle if only a cycle of major ticks are to be labeled, not explicitly suppressed as first or last label; will only then be used if tick is also in TickRange *)
#       TickLabel=If[
# 	(ShowTickLabels/.FullOpts)
# 	&&ApproxInRange[(TickLabelRange/.FullOpts),TickValue]
# 	&&(Mod[i-1,(TickLabelStep/.FullOpts)]==(TickLabelStart/.FullOpts))
# 	&&((i!=1)||(ShowFirst/.FullOpts))
# 	&&((i!=Length[MajorCoordList])||(ShowLast/.FullOpts)),
# 	UsedTickLabelFunction[TickValue,DefaultTickLabelFunction[TickValue]],
# 	""
# 		];
# 
#       (* make tick *)
#       {TickPosition,TickLabel,TickLength,TickStyle},
# 
#       {i,1,Length[MajorCoordList]}
# 		  ];
# 
#     (* make minor ticks *)
#     MinorCoordList=Select[
#       (TickPreTransformation/.FullOpts)/@RawMinorCoordList,
#       (ApproxInRange[(TickRange/.FullOpts),#]&&(TickTest/.FullOpts)[#])&
# 		   ];
#     TickLength=ResolveTickLength[
#       (MinorTickLength/.FullOpts),(TickLengthScale/.FullOpts),
#       (TickDirection/.FullOpts),(TickReverse/.FullOpts)
# 	       ];
#     TickStyle=(MinorTickStyle/.FullOpts);
#     TickLabelStyleWrapper=With[
#       {UsedMinorTickLabelStyle=(MinorTickLabelStyle/.FullOpts)},
#       If[
# 	UsedMinorTickLabelStyle=!=None,
# 	(Style[#,UsedMinorTickLabelStyle]&),
# 	Identity
#       ]
# 			  ];
#     MinorTickList=If[(ShowMinorTicks/.FullOpts),
# 		     Table[
# 
# 		       (* calculate tick value *)
# 		       TickValue=MinorCoordList[[i]];
# 
# 		       (* calculate coordinate for drawing tick *)
# 		       TickPosition=(TickPostTransformation/.FullOpts)[TickValue];
# 
# 		       (* construct label, or null string if it should be suppressed -- if: tick labels shown, *minor* tick labels shown, tick is in TickLabelRange; will only then be used if tick is also in TickRange *)
# 		       TickLabel=If[
# 			 (ShowTickLabels/.FullOpts)&&(ShowMinorTickLabels/.FullOpts)
# 			 &&ApproxInRange[(TickLabelRange/.FullOpts),TickValue],
# 			 TickLabelStyleWrapper[UsedTickLabelFunction[TickValue,DefaultTickLabelFunction[TickValue]]],
# 			 ""
# 				 ];
# 
# 		       (* make tick *)
# 		       {TickPosition,TickLabel,TickLength,TickStyle},
# 
# 		       {i,1,Length[MinorCoordList]}
# 		     ],
# 		     {}
# 		  ];
# 
#     (* combine tick lists*)
#     If[
#       (Debug/.FullOpts),
#       Print[RawMajorCoordList];
#       Print[RawMinorCoordList];
#       Print[Join[MajorTickList,MinorTickList]]
#     ];
# 
#     Join[MajorTickList,MinorTickList]
# 
#   ];

def linear_ticks(
        x1,x2,major_spacing,minor_subdivisions,
        minor_index_range=None,
        minor_index_transformation=None
):
    """Construct lists of major and minor tick locations.

    This is primarily intended as a "helper" function for set_linear_ticks() and
    should not normally need to be called by the user.

    Arguments:

        x1 (float): starting value

        x2 (float): ending bound

        major_spacing (float): spacing between major ticks

        minor_subdivisions (int): subdivisions for minor ticks between major ticks

        minor_index_range (tuple of int, optional): range of minor tick indices
            to include (used in log tick generation)

        minor_index_transformation (callable, optional): transformation on minor tick
            indices (used in log tick generation)

    Returns:

        major_tick_list (list of float): list of major tick locations

        minor_tick_list (list of float): list of minor tick locations (excludes
            major tick locations)

    """

    # construct preliminary tables of ticks
    
    # These are indexed by:
    #
    #     major_index=0,1,...,max_major_index
    #
    #     minor_index=0,1,...,minor_subdivisions-1
    #
    #  where minor_index=0 gives the major tick (excluded), though there should
    #  actually be no minor ticks after last major tick.

    max_major_index = round((x2-x1)/major_spacing)
    
    major_location_list = [
        x1+major_index*major_spacing
        for major_index in range(0,max_major_index+1)
        ]

    minor_location_list = [
        (
            x1+major_index*major_spacing
            +(minor_index if minor_index_transformation is None else minor_index_transformation(minor_index))
            *major_spacing/minor_subdivisions
        )
        for major_index in range(0,max_major_index+1)
        for minor_index in range(1,(minor_subdivisions-1)+1)
        if minor_index_range is None or in_range(minor_index,minor_index_range)
        ]

    # restrict tick ranges
    
    # There are generally ticks to be suppressed at the upper end, since the
    # major tick index rounds up to the next major tick (for safety in
    # borderline cases where truncation might fail), and the loop minor tick
    # index iterates for a full series of minor ticks even after the last major
    # tick.

    major_location_list = list(filter(lambda x : approx_in_range(x,(x1,x2)), major_location_list))
    minor_location_list = list(filter(lambda x : approx_in_range(x,(x1,x2)), minor_location_list))

    return major_location_list, minor_location_list

def linear_tick_locations(*args,**kwargs):
    """ Legacy name for linear_ticks.
    """
    warnings.warn("linear_tick_locations() has been renamed to linear_ticks()", DeprecationWarning)
    return linear_ticks(*args,**kwargs)
                          
def set_ticks(
        ax,axis,major_minor_ticks
):
    """Apply ticks to axis.

    Beware that applying ticks beyond the axis range (limits) will stretch the
    limits, which may need to then be reset.

    FUTURE: provide support for tick label formatting

    Arguments:

        ax (mpl.Axes.axes): axes

        axis (str): axis on which to operate ("x" or "y")

        major_minor_ticks (tuple):
            major_ticks (list): major tick locations
            minor_ticks (list): minor tick locations

    """

    if axis == "x":
        ticks_setter = ax.set_xticks
    elif axis == "y":
        ticks_setter = ax.set_yticks
    else:
        raise(ValueError("Unrecognized direction ({})".format(direction)))

    (major_ticks,minor_ticks) = major_minor_ticks
    ticks_setter(major_ticks)
    ticks_setter(minor_ticks, minor=True)
    
def main():    

    print(linear_tick_locations(0,10,2,4))

if __name__ == "__main__":
    main()
