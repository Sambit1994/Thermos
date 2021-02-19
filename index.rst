.. Automated NFS-G4(MP2) documentation.

==================
   NFS-G4(MP2) 
==================

The Original G4(MP2) protocol is inherent to Gaussian software package. This remodelled code runs the normal G4(MP2) protocol using ORCA suite of package along with the devised NFS-G4(MP2) protocols for high speed-accuracy trade-off. The script also provide the opportunity to the user to explore the affects of selection of Theory/Basis-set combination in the composite model and its overall affect on the prediction of the Thermochemical properties.

Current version is |release|.

.. _GitHub: https://github.com/Sambit1994/FAcT.github.io/


Key Features
============

- Prediction of user defined Thermochemical properties including, Enthalpy of Formation, Ionization Potential, Atomization Energy, Electron Affinity, Proton Affinity and Binding Energy.

- User defined Theory-Basis set combination at various levels.

- Geometry optimization free Energy calculation.

- HLC free Energy calculation.

Script/Dependent-libraries Installations
========================================

Installing running Script
-------------------------

.. code-block:: bash

   $ pip install ----

----statement----

.. code-block:: bash

   $ pip install -----

----statement----

Installing dependent libraries
------------------------------

.. code-block:: bash

   $ pip install ----

----statement----

Requirements and Outcome
========================

Input Files
-----------

The script requires the following files to run and generate the required thermochemical properties. 

- ``geom.xyz``: The cartesian coordinates (Angstrom unit) of the molecule along with information about the total no. of atoms, charge and multiplicity. One such example for methane molecule is given below.

.. code-block:: text

             5 [Total no. of atoms]
             0           1 [charge and multiplicity]
    C        0.00000000     0.00000000     0.00000000 [Cartesian Coordinates] 
    H        1.09336000     0.00000000     0.00000000
    H       -0.36445000     0.00000000    -1.03083000
    H       -0.36445000    -0.97188000     0.34361000
    H       -0.36445000     0.97188000     0.34361000

- ``orca_g4mp2.inp``: This input file contains the technical details for the various steps involved in the composite procedure. One such example for estimating `Enthalpy of Formation` using the normal G4(MP2) protocol is :ref:`G4MP2-inp`. It also allows the user to manipulate the settings as per their choice. For detailed information, please go to the tutorial section.

- ``read_geom_freq.dat``: The user can choose to determine the G4(MP2) energy on a pre-optimized geometry. In that case, they have to provide this input file along with the above two, contating information regarding the optimzed geometry and vibrational frequencies. To be noted, all vibrational frequencies must be included to run the script, i.e., 3N-5 for linear systems and 3N-6 for non-linear systems. One such example for pre-optimized methane molecule is given below:

.. code-block:: text

             5 [Total no. of atoms]
   Optimized atomic coordinates (Angstrom) [Title line]
         C     -0.00000123    -0.00000000    -0.00000079 [Optimized geometry coordinates]
         H      0.68957985     0.00051101    -0.84657788
         H     -1.02803745     0.00063322    -0.36792717
         H      0.16896006    -0.89208896     0.60649423
         H      0.16951220     0.89094474     0.60802020
  1343.36 [Vibrational frequencies at the optimized geometry]
  1343.36
  1343.38
  1564.11
  1564.17
  3042.96
  3156.57
  3156.61
  3156.67

Output Files
------------

The script produces three output files, ``ORCA_G4MP2.com``, ``ORCA_G4MP2.out``, and ``Thermochemistry.out``. The ``ORCA_G4MP2.com`` file compiles all the ORCA input files created by the script on-the-fly --- for the various steps involved in the composite procedure. Similarly, the output files generated from these corresponding input files are collected in ``ORCA_G4MP2.out``. The ``Thermochemistry.out`` file collects all the required energies and energy corrections involved in the G4(MP2) protocol along with the predicted value of the thermochemical property asked by the user. One such output files genrerated for methane molecule is given here. :ref:`G4MP2-com`, :ref:`G4MP2-out`, and :ref:`G4MP2-Thermo`. 

Tutorial
========

Technical Details
-----------------

The detailed information involved in the ``orca_g4mp2.inp`` are provided here: :ref:`G4MP2-Technical-Details`.


The G4(MP2) Variants
---------------------

The various devised G4(MP2) variants on the basis of the approximations used in the different steps of calculations are presented below:
	
+------------+-----------------+-----------+------------+------------+-----------+-----------+  
|   Method   |    Opt          |         CCSD(T)        |           MP2          |     HF    |
+============+=================+===========+============+============+===========+===========+
|            | RIJCOSX         |  DLPNO    |     RIJK   |   RI-MP2   |     RIJK  |     RIJK  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| G4(MP2)opt |       |tk|      |   |cr|    |     |cr|   |   |cr|     |     |cr|  |     |cr|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x3.16      |       |tk|      |   |tk|    |     |tk|   |   |cr|     |     |cr|  |     |cr|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x1.07      |       |tk|      |   |cr|    |     |cr|   |   |tk|     |     |tk|  |     |cr|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x1.20      |       |tk|      |   |cr|    |     |cr|   |   |cr|     |     |cr|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x1.29      |       |tk|      |   |cr|    |     |tk|   |   |cr|     |     |tk|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x1.32      |       |tk|      |   |cr|    |     |tk|   |   |tk|     |     |tk|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x7.20      |       |tk|      |   |tk|    |     |tk|   |   |cr|     |     |cr|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x1.32 |r8| |       |tk|      |   |cr|    |     |cr|   |   |tk|     |     |tk|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| x3.91      |       |tk|      |   |tk|    |     |tk|   |   |tk|     |     |tk|  |     |cr|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+
| G4(MP2)XP  |       |tk|      |   |tk|    |     |tk|   |   |tk|     |     |tk|  |     |tk|  |
+------------+-----------------+-----------+------------+------------+-----------+-----------+

The ``orca_g4mp2.inp`` files containing the technical details of the involved calculations are: 

- G4(MP2)opt: :ref:`G4MP2-inp-opt` 
- x3.16: :ref:`G4MP2-inp-x3.16` 
- x1.07: :ref:`G4MP2-inp-x1.07`
- x1.20: :ref:`G4MP2-inp-x1.20`
- x1.29: :ref:`G4MP2-inp-x1.29`
- x1.32: :ref:`G4MP2-inp-x1.32`
- x7.20: :ref:`G4MP2-inp-x7.20`
- x1.32 |r8|: :ref:`G4MP2-inp-x1.32-RI8`
- x3.91: :ref:`G4MP2-inp-x3.91`
- G4(MP2)XP: :ref:`G4MP2-inp-XP`

Demo
----

After collecting all the required files (see `Input Files`_), and selecting the required thermochemical property (see `Thermochemical Property Determination <file:///S:/TIFR/Sphinx/SPHINX/html/Technical-Details.html#thermochemical-property-determination>`_)

, the user just has to run the binary from the corresponding directory. The following shows one such example when HOF estimation is done for `CH4` using normal G4(MP2) method. For extended information of the running input files, please refer to :ref:`G4MP2-out`. 

.. code-block:: bash

   $ /apps/orca_g4mp2/g4mp2.x &

Geometry Optimization & Frequency calculation
.............................................

The calculation will start from the geomtery optimization and subsequent frequency calculation (at **GTBAS3** basis set level) and an input file (`input.com`) will be created with the necessary information to run the calculation:

.. code-block:: text 

    !   B3LYP/G Grid7 Printbasis  TightOpt  Freq
    * xyz   0  1
    C     0.00000000     0.00000000     0.00000000
    H     1.09336000     0.00000000     0.00000000
    H    -0.36445000     0.00000000    -1.03083000
    H    -0.36445000    -0.97188000     0.34361000
    H    -0.36445000     0.97188000     0.34361000
    *

Once the geometry optimization and frequency calculation is over, the script will extract the optimized geometry from the corresponding output file (`input.out`). The coordinates of the optimized geometry goes through an in-built scheme to verify the whether the molecule is linear or not. Depending on that the required *vibrational frequencies* are collected from the `input.out` file and are saved in ``freq.txt``.

The optimized coordinates and scaled frequencies are also subsequently get printed on the :ref:`G4MP2-Thermo` file 

.. code-block:: text

     Optimized atomic coordinates (Angstrom)
           C      0.00000024    -0.00000002     0.00000000
           H      1.07698196    -0.18012947    -0.00002129
           H     -0.52882274    -0.95534548     0.00013221
           H     -0.27411608     0.56763359    -0.89161936
           H     -0.27404603     0.56784158     0.89150844
     
     Unscaled harmonic wavenumbers (cm^-1)
            1323.51
            1323.52
            1323.52
            1541.17
            1541.17
            2998.34
            3110.48
            3110.49
            3110.50


After geometry optimization, the optimized coordinates are used for rest of the single point energy calculations.

CCSD(T) Energy calculation
..........................

The earlier `input.com` and `input.out` files are discarded and new files are created for each of the single point calculations.

.. code-block:: text

    ! CCSD(T) Printbasis
    * xyz   0  1
    C     0.00000024    -0.00000002     0.00000000
    H     1.07698196    -0.18012947    -0.00002129
    H    -0.52882274    -0.95534548     0.00013221
    H    -0.27411608     0.56763359    -0.89161936
    H    -0.27404603     0.56784158     0.89150844
    *

From the CCSD(T) step, both the CCSD(T) energy and inherent MP2 energy (at small basis **GTBAS1**) are collected

.. code-block:: text

            -------------------------   --------------------
            FINAL SINGLE POINT ENERGY       -40.354687160506 [CCSD(T) Energy] 
            -------------------------   --------------------

            Initial E(tot)                             ...    -40.331248094 [Inherent MP2 energy in the CCSD(T) step]


MP2 Energy calculation
......................

CCSD(T) calculation is followed by MP2 energy calculation at a large basis-set (**GTMP2LargeXP**)

.. code-block:: text

    ! MP2 Printbasis
    * xyz   0  1
    C     0.00000024    -0.00000002     0.00000000
    H     1.07698196    -0.18012947    -0.00002129
    H    -0.52882274    -0.95534548     0.00013221
    H    -0.27411608     0.56763359    -0.89161936
    H    -0.27404603     0.56784158     0.89150844
    *

The MP2 energy is collected from the corresponding `input.out` file and is substracted to the onde from the CCSD(T) step to calculate DE(MP2). Along with the MP2 energy the inherent HF energy is also collected, required to determine DE(HF). 

.. code-block:: text

            -------------------------   --------------------
            FINAL SINGLE POINT ENERGY       -40.405646162405 [MP2 energy at large basis]
            -------------------------   --------------------

            ----------------
            TOTAL SCF ENERGY
            ----------------

            Total Energy       :          -40.21213373 Eh           -1094.22779 eV [Inherent HF energy]

HF Energy calculation
.....................

Two HF energy calculations are carried out, first with a small basis **GFHFB3** and the next with a large basis-set **GFHFB4**

.. code-block:: text

    ! HF Printbasis
    * xyz   0  1
    C     0.00000024    -0.00000002     0.00000000
    H     1.07698196    -0.18012947    -0.00002129
    H    -0.52882274    -0.95534548     0.00013221
    H    -0.27411608     0.56763359    -0.89161936
    H    -0.27404603     0.56784158     0.89150844
    *

Using the extrapolation scheme, the HF energy at basis-set limit is determined (see `HF Energies Determination <file:///S:/TIFR/Sphinx/SPHINX/html/Technical-Details.html#hf-energies-determination>`_) and is substracted from the HF energy obtained from **GTMP2LargeXP** basis set level to determine DE(HF).  

After each sucessive step, the ``Thermochemistry.out`` file keeps on updating and prints the computation time involved for each step

.. code-block:: text

     * Geometry optimization/frequencies done in             151.62 s
     ** Elapsed time =             152.26 s
     
     * CCSD(T) done in               1.91 s
     ** Elapsed time =             154.18 s
     
     * MP2/L done in                 6.91 s
     ** Elapsed time =             161.11 s
     
     * HF/VTZ done in                3.18 s
     ** Elapsed time =             164.32 s
     
     * HF/VQZ done in               15.46 s
     ** Elapsed time =             179.87 s


Once all the calculation steps are over for a system, the various energies and energy corrections as well as the calculated G4(MP2) energies are printed on the ``Thermochemistry.out`` file. The ``Elapsed time`` shows the time involved in the execution of each of the individual steps.

.. code-block:: text

               Temperature=     298.150000               Pressure=       1.000000
                    E(ZPE)=       0.044157             E(Thermal)=       0.047028
                       HLC=      -0.037888                     SO=       0.000000
                E(CCSD(T))=     -40.354687
                   DE(MP2)=      -0.074398                 DE(HF)=      -0.004870
                G4MP2(0 K)=     -40.427686           G4MP2 Energy=     -40.424814
            G4MP2 Enthalpy=     -40.423870      G4MP2 Free Energy=       0.000000


- ``E(ZPE)``: Scaled Zero Point Vibrational Energy

- ``E(Thermal)``: Molecular Thermal energy corrections, determined from translational, rotational and vibrational degrees of freedom

- ``HLC``: Higher Level Correction (The semi-empirical correction)

- ``SO``: The Spin-Orbit correction 

- ``E(CCSD(T))`` The CCSD(T) energy

- ``DE(MP2)`` The MP2 energy correction

- ``DE(HF)`` The HF energy correction

- ``G4MP2(0 K)`` The total G4(MP2) energy (:math:`U_{0}`). The energy is used for estimating various thermochemical properties.

- ``G4MP2 Energy`` The G4(MP2) energy after thermal correction (:math:`U_{T}`) 

- ``G4MP2 Enthalpy`` The G4(MP2) enthalpy (:math:`H_{T}`)  

- ``G4MP2 Free Energy`` Free energy. **not getting estimated**

For estimating HOF, after G4(MP2) energy is estimated for the molecule is done, the same for the constituent atoms are determined. Similar to the molecule, the dfferent single point energies for each of the constituent atom are carried (except geometry optimization and frequency calculation). 

The ``Thermochemistry.out`` file prints the energies of the individual atoms once they are done.

.. code-block:: text

     * G4MP2 for reference atoms*
     
     * ATOM:   C
     
     * Geometry optimization skipped for atom
     
     * CCSD(T) done in               1.99 s
     ** Elapsed time =             182.08 s
     
     * MP2/L done in                 2.27 s
     ** Elapsed time =             184.46 s
     
     * HF/VTZ done in                2.63 s
     ** Elapsed time =             187.17 s
     
     * HF/VQZ done in                3.93 s
     ** Elapsed time =             192.22 s
     
               Temperature=     298.150000               Pressure=       1.000000
                    E(ZPE)=       0.000000             E(Thermal)=       0.001416
                       HLC=      -0.013971                     SO=      -0.000140
                E(CCSD(T))=     -37.751279
                   DE(MP2)=      -0.025410                 DE(HF)=      -0.003476
                G4MP2(0 K)=     -37.794277           G4MP2 Energy=     -37.792861
            G4MP2 Enthalpy=     -37.791917      G4MP2 Free Energy=       0.000000

     * G4(MP2) done
     ** Elapsed time =             192.27 s
     
     * ATOM:   H
     
     * Geometry optimization skipped for atom
     
     * CCSD(T) done in               2.12 s
     ** Elapsed time =             194.44 s
     
     * MP2/L done in                 1.49 s
     ** Elapsed time =             195.95 s
     
     * HF/VTZ done in                1.65 s
     ** Elapsed time =             197.63 s
     
     * HF/VQZ done in                1.58 s
     ** Elapsed time =             199.24 s
     
               Temperature=     298.150000               Pressure=       1.000000
                    E(ZPE)=       0.000000             E(Thermal)=       0.001416
                       HLC=      -0.002115                     SO=       0.000000
                E(CCSD(T))=      -0.498233
                   DE(MP2)=      -0.001585                 DE(HF)=      -0.000161
                G4MP2(0 K)=      -0.502094           G4MP2 Energy=      -0.500677
            G4MP2 Enthalpy=      -0.499733      G4MP2 Free Energy=       0.000000
     
     * G4(MP2) done
     ** Elapsed time =             199.26 s

After the required calculation, the ``SUMMARY`` section prints the estimated data of the thermochemical property chosen by the user (Here HOF)

.. code-block:: text

     *---------* 
     * SUMMARY * 
     *---------* 
     
     * Atomization energy (Atoms - mol) =                0.62503407 Hartree
     * Atomization energy (Atoms - mol) =               17.00804473 eV
     * Atomization energy (Atoms - mol) =              392.21486014 kcal/mol
     * Atomization energy (Atoms - mol) =             1641.02697481 kj/mol
     * Atomization energy (Atoms - mol), NO HLC =              382.51544470 kcal/mol
     
     ** Elapsed time =             199.26 s
     
     *---------* 
     * SUMMARY * 
     *---------* 
     
     * Heat of formation =               -0.02806423 Hartree
     * Heat of formation =               -0.76366657 eV
     * Heat of formation =              -17.61057079 kcal/mol
     * Heat of formation =              -73.68262820 kj/mol
     * Heat of formation, NO HLC =               -7.91115535 kcal/mol

For HOF, both `Atomization Energy` and `Enthapy of formation` are printed in different units, i.e., *Hartree*, *eV*, *kcal/mol*, and *kj/mol*. 


.. |tk| unicode:: U+2714
.. |cr| unicode:: U+2718
.. |r8| unicode:: U+01C2


Source code
===========

The project is hosted on GitHub_

Please feel free to file an issue on the `bug tracker
<https://github.com/>`_ if you have found a bug
or have some suggestion in order to improve the library.

Dependencies
============

Communication channels
======================

Contributing
============

