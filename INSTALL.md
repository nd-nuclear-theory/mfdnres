# mfdnres installation guide #

Mark A. Caprio, Patrick J. Fasano
Department of Physics, University of Notre Dame

+ 05/24/19 (mac): Created.
+ 05/18/22 (mac): Update basic examples.
+ 07/02/26 (mac): Add environment variable MFDNRES_RESULTS_DIR.
+ 09/20/26 (mac): Remove pip install --editable flag.
----------------------------------------------------------------

# 1. Retrieving and installing source

  Change to the directory where you want the repository to be installed,
  e.g.,

  ~~~~~~~~~~~~~~~~
  % cd ~/code
  ~~~~~~~~~~~~~~~~

  Clone the `mfdnres` repository.

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git clone https://github.com/nd-nuclear-theory/mfdnres.git
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Change your working directory to the repository for the following steps:
  
  ~~~~~~~~~~~~~~~~
  % cd mfdnres
  ~~~~~~~~~~~~~~~~

  If you want the bleeding-edge, potentially broken version, check out the
  `develop` branch:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % git checkout -t origin/develop
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Set up the package in your `PYTHONPATH` by running `pip`:

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % python3 -m pip install --user .
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Note that the `.` here means to install the Python package defined by the code
  in the current directory.  If you are actively developing `mfdnres` itself,
  you may want to pass the `--editable` flag to `pip`, so that your edits take
  effect immediately, without your needing to run `pip install` again.

  This basic installation does not check that you have certain dependencies
  (matplotlib or Pandas) installed.  These are only needed if you are using the
  analysis and plotting machinery, not for the basic operation of `mfdnres`.  If
  you want to make sure that these dependencies are also installed, do instead
  (or in addition):
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % python3 -m pip install --user --editable ".[analysis]"
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  a. Subsequently updating source:

  ~~~~~~~~~~~~~~~~
  % git pull
  % python3 -m pip install --user --editable .
  ~~~~~~~~~~~~~~~~

  This subsequent `pip install`, when updating the source code, is a precaution
  in case, e.g., the package dependencies have changed.

# 2. Environment variables

  For most purposes, you do not need to set any environment variables.  However,
  if you wish to use the function `mfdnres.res_file_directory` to construct path
  names to your data files, you will need to set the `MFDNRES_RESULTS_DIR` environment
  variable, to point to the parent directory where various users' results files
  are stored.  This assumes a directory structure such as the following:

      results/

        alice/
          mfdn/
            runalice0001/
            runalice0002/
          spncci/
            runalice0003/
            runalice0004/

        bob/
          mfdn/
            runbob0001/
            runbob0002/

  Then you should set `MFDNRES_RESULTS_DIR` to point to the `results/`
  directory.  E.g., if you have downloaded all the results to to a directory
  named `results` under your home directory, you would just set
  `MFDNRES_RESULTS_DIR` to point to `${HOME}/results`.

  If your default shell is csh, define initialization as follows (adjusting
  directory names to match your own choices as appropriate):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  setenv MFDNRES_RESULTS_DIR ${HOME}/results
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Alternatively, if your default shell is bash:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  export MFDNRES_RESULTS_DIR=${HOME}/results
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3. Basic tests

  Basic test scripts to read and parse results files may be found in
  `doc/examples/read_res_test`.

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cd mfdnres
  % cd doc/examples/read_res_test
  % python3 read_res_mfdn.py
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  See also the tutorial in `doc/examples/tutorial`.
