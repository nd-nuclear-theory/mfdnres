# mfdnres installation guide #

Mark A. Caprio, Patrick J. Fasano  
Department of Physics, University of Notre Dame

+ 05/24/19 (mac): Created.
+ 05/18/22 (mac): Update basic examples.

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
  % python3 -m pip install --user --editable .
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Note that the `.` here means to install the Python package defined by the code
  in the current directory.

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
  if you wish to use the function mfdnres.res_file_directory() to construct path
  names to your data files, you will need to set the GROUP_HOME environment
  variable, to point to the parent directory where various users' results files
  are stored.  This assumes a directory structure such as the following:

    /afs/crc.nd.edu/group/nuclthy/

      results/

        alice/
          runalice0001/
          runalice0002/

        bob/
          runbob0001/
          runbob0002/

  Then you should set GROUP_HOME to point to the directory which contains the
  results/ directory, e.g., here /afs/crc.nd.edu/group/nuclthy.  Or, you have
  downloaded all the results to to a directory named results under your home
  directory, you would just set GROUP_HOME to point to ${HOME}.

  In your csh initialization file, define initialization as follows
  (adjusting directory names to match your own choices as
  appropriate):

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  setenv GROUP_HOME /afs/crc.nd.edu/group/nuclthy
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  or

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  setenv GROUP_HOME ${HOME}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Alternatively, if your default shell is bash:
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  export GROUP_HOME=/afs/crc.nd.edu/group/nuclthy
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  or

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mfdnres
  export GROUP_HOME=${HOME}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3. Basic tests

  Basic test scripts to read and parse results files may be found in
  `docs/examples/read_res_test`.

  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % cd mfdnres
  % cd docs/examples/read_res_test
  % python3 read_res_mfdn.py
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  See also the tutorial in `docs/examples/tutorial`.