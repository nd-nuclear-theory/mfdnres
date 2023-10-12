""" Test code for mfdn_format_7_ho.

    Language: Python 3
    Mark A. Caprio
    Patrick J. Fasano
    University of Notre Dame

"""

import mfdnres

if (__name__ == "__main__"):

    filename = r"run0000-mfdn-Z2-N6-Daejeon16-coul1-hw05.000-a_cm20-Nmax02-Mj0.0-lan500-tol1.0e-06-natorb-no0.res"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"run0000-mfdn-Z2-N6-Daejeon16-coul1-hw05.000-a_cm20-Nmax02x-Mj0.0-lan500-tol1.0e-06.res"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"runpjf0015-mfdn15-Z3-N4-JISP16-coul1-hw20.000-a_cm40-Nmax02-Mj0.5-lan1000-tol1.0e-06-natorb-no0.res"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"runpjf0069-mfdn15-Z2-N1-Daejeon16-coul1-hw22.500-a_cm0-Nmax14-Mj0.5-lan400-tol1.0e-06-natorb-J00.5-g0-n01-no0.res"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)

    filename = r"runmac0746-mfdn15-Z6-N8-Daejeon16-coul1-hw15.000-a_cm50-Nmax04-Mj0.0-lan800-tol1.0e-06-J00.0-g0-n01-U3SpSnS-dlan1200.lanczos"
    info = mfdnres.input.parse_filename(filename, filename_format="mfdn_format_7_ho")
    print(filename)
    print(info)
