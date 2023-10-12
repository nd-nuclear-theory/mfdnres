""" Test code for mfdn_format_7_pm.

    Language: Python 3
    Mark A. Caprio
    Patrick J. Fasano
    University of Notre Dame

"""

import mfdnres

if (__name__ == "__main__"):

    filename = r"MFDn.res.Z4.N5.JISP16.Nmin1.Nm13.hw20.0.La500.St06.tol1e-6"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_pm")
    print(filename)
    print(info)

    filename = r"MFDn.res.A07.Z3.N4.Nm02.chi2bSMSI0B_nocd_srg0400ho40J_HO016"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_pm_lenpic")
    print(filename)
    print(info)
