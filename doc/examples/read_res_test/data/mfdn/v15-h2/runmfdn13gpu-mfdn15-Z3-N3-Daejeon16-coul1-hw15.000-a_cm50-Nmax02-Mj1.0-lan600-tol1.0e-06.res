 
[MFDn]
Version = 15
Revision = v15b01-205-gdb2400a
Platform = 
Username = 
ndiags   =        1
MPIranks  =        1
OMPthreads =       32
 
[PARAMETERS]
 
[Basis]
Nprotons  =        3
Nneutrons =        3
# Orbitals read in from file
TwoMj  =        2
parity =        1
Nmin   =        0
Nmax   =        2
DeltaN =        2
WTmax  =     4.1000
 
[Interaction]
# reading complete Hamiltonian TBMEs in H2 format
Hrank  =        2
TBMEfile = tbme-H.bin
 
[Many-body matrix]
dimension  =                711
numnonzero =              28238
 
[Observables]
numTBops   =       15
# TBME files for additional operators
TBMEfile(2) = tbme-H.bin
TBMEfile(3) = tbme-Ncm.bin
TBMEfile(4) = tbme-Tintr.bin
TBMEfile(5) = tbme-Tcm.bin
TBMEfile(6) = tbme-VNN.bin
TBMEfile(7) = tbme-VC.bin
TBMEfile(8) = tbme-L2.bin
TBMEfile(9) = tbme-Sp2.bin
# TBME file for relative R2 operator
TBMEfile(1) = tbme-rrel2.bin
 
 
[RESULTS]
# following Condon-Shortley phase conventions
# for conventions of reduced matrix elements and
# for Wigner-Eckart theorem, see Suhonen (2.27)
 
[Energies]
# Seq       J    n      T        Eabs        Eexc        Error    J-full
    1      1.0   1   0.000      -25.578      0.000     0.61E-05    1.0000
    2      3.0   1   0.000      -24.148      1.429     0.55E-05    3.0000
    3      2.0   1   1.000      -19.983      5.594     0.91E-05    2.0000
    4      2.0   2   0.001      -18.975      6.602     0.15E-04    2.0000
 
[Oscillator quanta]
# Seq    J    n      T      Amp(N)^2 
    1   1.0   1   0.000      0.8818      0.1182    
    2   3.0   1   0.000      0.8867      0.1133    
    3   2.0   1   1.000      0.8937      0.1063    
    4   2.0   2   0.001      0.8271      0.1729    
 
[Angular momenta]
# Seq    J    n      T      <L^2>      <S^2>      <Sp^2>     <Sn^2>     <J^2>      <T^2>
    1   1.0   1   0.000     0.2358     1.9653     0.7664     0.7663     2.0000     0.0001
    2   3.0   1   0.000     6.1210     2.0880     0.7827     0.7824    12.0000     0.0001
    3   2.0   1   1.000     4.8333     0.7182     0.7813     0.7814     6.0000     1.9996
    4   2.0   2   0.001     5.9736     2.0309     0.7697     0.7694     6.0000     0.0006
 
[Relative radii]
# Seq    J    n      T      r(p)       r(n)    r(matter)     r_pp       r_pn       r_nn
    1   1.0   1   0.000     2.0904     2.0830     2.0867     0.9681     1.5773     0.9641
    2   3.0   1   0.000     2.0438     2.0364     2.0401     0.9406     1.5492     0.9366
    3   2.0   1   1.000     2.0651     2.0562     2.0606     0.9502     1.5651     0.9454
    4   2.0   2   0.001     2.1809     2.1724     2.1767     1.0056     1.6507     1.0009
 
[Other 2-body observables]
# Seq    J    n      T     see header for TBME file names of observables
    1   1.0   1   0.000     -25.5777         0.902046E-08      75.1700          11.2500         -102.499          1.75138         0.235849         0.766417    
    2   3.0   1   0.000      6.12105         0.782653         -19.9833         0.236560E-08      76.5757          11.2500         -98.3325          1.77357    
    3   2.0   1   1.000     -91.2318          1.68472          5.97360         0.769692          0.00000          0.00000          0.00000          0.00000    
    4   2.0   2   0.001      0.00000          0.00000          0.00000          0.00000          0.00000          0.00000          0.00000          0.00000    
 
[Occupation probabilities]
# orb      n_rad   l   2*j
    1   pro    0    0    1     1.936        1.935        1.940        1.940    
    2   pro    0    1    1    0.2451       0.3324E-01   0.9336E-01   0.4515    
    3   pro    0    1    3    0.7704       0.9891       0.9261       0.5264    
    4   pro    1    0    1    0.1409E-01   0.1913E-01   0.1864E-01   0.1263E-01
    5   pro    0    2    3    0.1600E-01   0.1178E-01   0.9370E-02   0.2172E-01
    6   pro    0    2    5    0.1260E-01   0.7383E-02   0.7137E-02   0.1664E-01
    7   pro    1    1    1    0.3003E-02   0.1193E-14   0.1609E-02   0.1895E-01
    8   pro    1    1    3    0.6625E-03   0.3446E-03   0.1220E-02   0.1722E-02
    9   pro    0    3    5    0.2060E-02   0.9655E-04   0.1621E-02   0.8197E-02
   10   pro    0    3    7    0.1279E-14   0.4015E-02   0.1025E-02   0.2746E-02
# sum of proton occ.prob.     3.000        3.000        3.000        3.000    
 
    1   neu    0    0    1     1.935        1.933        1.939        1.939    
    2   neu    0    1    1    0.2436       0.3357E-01   0.8485E-01   0.4463    
    3   neu    0    1    3    0.7730       0.9896       0.9360       0.5338    
    4   neu    1    0    1    0.1480E-01   0.2025E-01   0.1950E-01   0.1267E-01
    5   neu    0    2    3    0.1577E-01   0.1177E-01   0.9168E-02   0.2127E-01
    6   neu    0    2    5    0.1244E-01   0.7210E-02   0.7027E-02   0.1652E-01
    7   neu    1    1    1    0.2708E-02   0.1246E-14   0.1228E-02   0.1755E-01
    8   neu    1    1    3    0.4611E-03   0.1682E-03   0.8896E-03   0.1474E-02
    9   neu    0    3    5    0.2060E-02   0.9847E-04   0.1495E-02   0.8285E-02
   10   neu    0    3    7    0.1551E-14   0.3970E-02   0.1109E-02   0.2716E-02
# sum of neutron occ.prob.    3.000        3.000        3.000        3.000    
 
