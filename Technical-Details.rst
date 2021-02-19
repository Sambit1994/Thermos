.. _G4MP2-Technical-Details:

Technical-Details
=================

.. currentmodule:: G4MP2

.. code-block:: text

Libraries and Machine Details
-----------------------------

.. code-block:: text

   openmpi_dir ="/apps/openmpi-3.0.0_install/"
   orca_dir ="/home/project_g4mp2_2019/ORCA_4_2/orca_4_2_0_linux_x86-64_openmpi314/"
   install_dir ="/apps/orca_g4mp2/"
  
   maxcore_mb = 32000
   nproc = 1

This segment contains information regarding the absolute path to different libraries necessary for running the script. The ``openmpi_dir`` contains the path to the directory containing Open MPI (Message Passing Interface) library for parallel computation. The variable ``orca_dir`` defines the absolute path to the directory containing ORCA libraries. The ``install_dir`` specifies the directory containing binary files of the script as well as the custom Basis-sets used in the calculation.

The option ``maxcore_mb`` defines the assigned memory for the ORCA calculations and ``nproc`` determines the number of cores assigned to each of the calculation.

Geometry Optimization
---------------------

.. code-block:: text

   method_opt_freq = "B3LYP/G "
   basis_opt_freq = "GTBAS3"
   custombasis_opt_freq = .true.
   scalfac = 0.9854
  
   String_Opt = "TightOpt"
  
   MGGA = .false.
  
   FROZEN_GEOM = .false.

The details involving geometry optimization of the required molecule as well as subsequent frequency calculation of the optimized structure are defined here. The ``method_opt_freq`` option specifies the method of choice for geometry optimization and **B3LYP/G** (Gaussian implementaion of B3LYP) DFT functional is the recommended method for G4(MP2) study. 

The ``basis_opt_freq`` option selects the custom basis from the library directory (``install_dir``) to use along with **B3LYP/G** functional for geometry optimization, i.e., **GTBAS3** [modified form of **6-31G(2df,p)**], while ``custombasis_opt_freq =`` **.true.** uses the selected custom basis-set in geometry optimization and frequency calculation. The parameter ``scalfac`` decides the *scaling factor* to determine the Zero Point Vibrational Energy (ZPVE) from the calculated vibarational frequency, which is specific for each *Method/Basis-set* combination. For **B3LYP/6-31G(2df,p)** the recommended scaling factor for ZPVE is 0.9854.

The user has the libery to run the geometry optimization with *Method/Basis-set* of his choice as well as other additional specifications by changing the predefined options. One such example with **wB97X-D3** functional and **def2-TZVP** basis set is presented below.  

.. code-block:: text

   method_opt_freq = "wB97X-D3 def2-TZVP Grid7 Printbasis "
   basis_opt_freq = ""
   custombasis_opt_freq = .false.
   scalfac = 0.9791

   "The script by default runs frequency calculation after geometry optimization"

The ``String_Opt`` selects the convergence criteria for the geometry optmization, i.e., **LooseOpt**, **NormalOpt**, **TightOpt**, or **VeryTightOpt**. The ``MGGA`` option is utilized for methods having no atalytical hessian like meta-GGA functionals, **M06**, **M062X**, **TPSS0**, correlated methods like **MP2**. For those methods, the user can decide to find Hessian and frequency numerically by keeping ``MGGA =`` **.true.**. Switching the option to true will add **NumHess** and **NumFreq** options to the geometry optimization step to carry out the ORCA implemented numerical hessian and frequency calculation.

If the user decides to skip the geometry optimization and frequency calculation and use a pre-determined geomtry to estimate the G4(MP2) energy, they have to keep the option ``FROZEN_GEOM`` to **.true.**. Along with this they have to provide the ``read_geom_freq.dat`` file containing the details about geometry and vibrational frequencies, please see `Input Files <https://moldis-group.github.io/Pople/#input-files>`_ 


The devised G4(MP2)XP protocol utilizes *RIJCOSX* approximation to fasten the optimization calculation with **GTBAS3** as orbital basis and **def2/J** as auxiliary basis-set and is presented below.

.. code-block:: text

   method_opt_freq = "B3LYP/G Grid7 RIJCOSX def2/J GridX9 Printbasis "
   basis_opt_freq = "GTBAS3"
   custombasis_opt_freq = .true.
  
   String_Opt = "TightOpt"


CCSD(T) Energy Determination
----------------------------

.. code-block:: text
  
   method_ccsdt = "CCSD(T) "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .false.

   switch_read_RIMP2_small = .false.  

   switch_RIMP2_small=.false.
   method_mp2_s = " "
   basis_mp2_s = ""
   custombasis_mp2_s = .true.

The CCSD(T) calculation determines the reference energy in G4(MP2) theory, to which the different energy corrections are added. The method of choice for reference energy along with other additional informations are provided in ``method_ccsdt``. The ``basis_ccsdt`` option selects the custom basis for CCSD(T) step, i.e., **GTBAS1** [modified form of **6-31G(d)**]. To create the ORCA input file for CCSD(T) calcuation with the **GTBAS1** custom basis, the user has to keep ``custombasis_ccsdt`` to **.true.**. 

The user can choose other basis-sets along with additional specifications, by making changes in the above options. For an instance, the example below uses **def2-TZVP** basis set to carry out the CCSD(T) calculation.

.. code-block:: text
  
   method_ccsdt = "CCSD(T) def2-TZVP "
   basis_ccsdt = ""
   custombasis_ccsdt = .false.

The devised highly in-expensive G4(MP2)-XP protocol utilizes *DLPNO* approximation in the CCSD(T) and the same with **TightPNO** threshold settings can be invoked as follows:

.. code-block:: text
  
   method_ccsdt = "DLPNO-CCSD(T) def2-SVP/C TightPNO Printbasis "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .true.

where **GTBAS1** is used as orbital basis and **def2-SVP/C** as the auxiliary basis.

The next set of options is used to manipulate the MP2 energy at the **GTBAS1** basis-set level.

The :math:`\Delta E(MP2)` correction energy in G4(MP2) theory utilized MP2 energy from two basis-sets levels, **GTMP2LargeXP** (Large, see `MP2 Energy Determination`_) and **GTBAS1** (small). The latter is determined from the inherent MP2 energy determined in the CCSD(T) calculation itself. However, the MP2 energy obtained in the DLPNO-CCSD(T) calculation considers local approimation and can provide slightly deviating energy from that of Normal MP2 or RI-MP2. However, we made sure irrespective of the approximation, for MP2 energy calculation step (see `MP2 Energy Determination`_) , the level of theory remains intact for MP2 energy at either basis set level. 

The user can decide to do so through two ways. 

1. The user can choose to implement an additional MP2 calculation in the CCSD(T) step along with the inherent MP2 and then would need to select the ``switch_read_RIMP2_small`` option to **true**. For and instance, both G4(MP2)-x3.16 and G4(MP2)XP uses DLPNO approximation for CCSD(T) with **GTBAS1** basis-set, while MP2 energy calculation at **GTMP2LargeXP** basis-set uses different approximations. Thus, for G4(MP2)-x3.16 which calculates Normal MP2 energy with **GTMP2LargeXP** basis-set, the MP2 energy with **GTBAS1** basis-set should also be Normal MP2 energy, rather than the local MP2 energy inherent in DLPNO-CCSD(T) calculation. So to achieve the same, the following can be done:

.. code-block:: text
  
   method_ccsdt = "DLPNO-CCSD(T) MP2 def2-SVP/C RIJK def2/JK TightPNO Printbasis "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .true.

   switch_read_RIMP2_small = .true.  

   switch_RIMP2_small=.false.
   method_mp2_s = " "
   basis_mp2_s = ""
   custombasis_mp2_s = .true.

`RIJK` implement the `RI` approximation in the SCF step as well while using **def2/JK** as auxiliary basis-set. Similarly, in case of G4(MP2)XP, that uses RI-MP2 approximation for determining MP2 energy with **GTMP2LargeXP** basis-set, an additional RI-MP2 calculation is done in the CCSD(T) step as follows:

.. code-block:: text
  
   method_ccsdt = "DLPNO-CCSD(T) RI-MP2 def2-SVP/C RIJK def2/JK TightPNO Printbasis "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .true.

   switch_read_RIMP2_small = .true.  

   switch_RIMP2_small=.false.
   method_mp2_s = " "
   basis_mp2_s = ""
   custombasis_mp2_s = .true.

2. The user can also decide to run an independent MP2 calculation with **GTBAS1** basis-set and would need to invoke the ``switch_RIMP2_small`` option to **true**. For G4(MP2)-x3.16, it is done as:

.. code-block:: text
  
   method_ccsdt = "DLPNO-CCSD(T) def2-SVP/C RIJK def2/JK TightPNO Printbasis "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .true.

   switch_read_RIMP2_small = .false.  

   switch_RIMP2_small=.true.
   method_mp2_s = "MP2 RIJK def2/JK Printbasis"
   basis_mp2_s = "GTBAS1"
   custombasis_mp2_s = .true.

The same for G4(MP2)XP is done as follows:

.. code-block:: text
  
   method_ccsdt = "DLPNO-CCSD(T) def2-SVP/C RIJK def2/JK TightPNO Printbasis "
   basis_ccsdt = "GTBAS1"
   custombasis_ccsdt = .true.
   DLPNO_CCSDT = .true.

   switch_read_RIMP2_small = .false.  

   switch_RIMP2_small=.true.
   method_mp2_s = "RI-MP2 def2-SVP/C RIJK def2/JK Printbasis"
   basis_mp2_s = "GTBAS1"
   custombasis_mp2_s = .true.


Although both method-1 and 2 will provide the same accuracy, with method-2 an additional HF/SCF calculation is added, which although uses negligible computation time for small molecules, but for big molecular systems, it can add to the total time to a large extent. Thus, we perferred to rely on method-1 to manipulate the MP2 energy determined with **GTBAS1** basis-set.

MP2 Energy Determination
------------------------

.. code-block:: text

   method_mp2 = "MP2 "
   basis_mp2 = "GTMP2largeXP"
   custombasis_mp2 = .true.
   flag_RIMP2 = .false.
   flag_DLPNOMP2 = .false.

The above block of variables modulate the MP2 energy determined with the large basis-set, **GTMP2largeXP** [modified form of **6-311G(3df,2p)** basis-set]. Similar to the earlier steps, the user can decide to rely on other basis set than the custom basis **GTMP2LargeXP** and need to change the variables ``method_mp2``, ``basis_mp2``, and ``custombasis_mp2``. Similar example with **def2-QZVP** basis set is given below:

.. code-block:: text

   method_mp2 = "MP2 def2-QZVP "
   basis_mp2 = ""
   custombasis_mp2 = .false.

The variables ``flag_RIMP2`` and ``flag_DLPNOMP2`` are used to activate the `RI` and `DLPNO` approximation in the MP2 calculation step. Subsequent changes are also needed to be done in the ``method_mp2`` variable as well. The examples below show how RI-MP2 and DLPNO-MP2 calculations are triggered in the script.

RI-MP2 with **def2-QZVP** orbital basis and **def2-QZVPPD/C** auxiliary basis. 

.. code-block:: text

   method_mp2 = "RI-MP2 def2-QZVP def2-QZVPPD/C "
   basis_mp2 = ""
   custombasis_mp2 = .false.
   flag_RIMP2 = .true.
   flag_DLPNOMP2 = .false.

DLPNO-MP2 with **def2-QZVP** orbital basis and **def2-QZVPPD/C** auxiliary basis. 

.. code-block:: text

   method_mp2 = "DLPNO-MP2 def2-QZVP def2-QZVPPD/C "
   basis_mp2 = ""
   custombasis_mp2 = .false.
   flag_RIMP2 = .false.
   flag_DLPNOMP2 = .true.

The devised G4(MP2)XP protocol utilizes RI approximations in this step at both SCF iteration and MP2 correlation calculation while using the custom basis **GTMP2LargeXP** as orbital basis-set. `RIJK` invokes the approximation in the SCF step and uses **def2/JK** as auxiliary basis-set, while **def2-TZVP/C** is used as auxillary basis set for correlation calculation.

.. code-block:: text

   method_mp2 = "RI-MP2 def2-TZVP/C RIJK def2/JK Printbasis "
   basis_mp2 = "GTMP2largeXP"
   custombasis_mp2 = .true.
   flag_RIMP2 = .true.
   flag_DLPNOMP2 = .false.

HF Energies Determination
--------------------------

.. code-block:: text

   method_hf3 = "HF "
   basis_hf3 = "GFHFB3"
   custombasis_hf3 = .true.
  
   method_hf4 = "HF "
   basis_hf4 = "GFHFB4"
   custombasis_hf4 = .true.

   HF_CBS_default      = .true.
   HF_CBS_orca_23_def2 = .false.
   HF_CBS_orca_34_def2 = .false.
   HF_CBS_orca_23_cc   = .false.
   HF_CBS_orca_34_cc   = .false.

The modulation of HF calculations to determine HF energy at basis set limit (:math:`HF_{limit}`) is done here. The G4(MP2) theory uses a two-point extrapolation scheme to evaluate (:math:`HF_{limit}`) and uses HF energies at two levels of basis-set size. ``method_hf3`` and ``method_hf4`` define the two levels of HF calculations; Triple :math:`\zeta` and Quadruple :math:`\zeta` level.

``basis_hf3`` and ``basis-hf4`` selects the custom basis-sets for the corresponding HF calculations, which are **GFHFB3** (modified form of **aug-cc-pV(T+d)Z**) and **GFHFB4** (modified form of **aug-cc-pV(Q+d)Z**) basis-sets. keeping ``custombasis_hf3`` and ``custombasis_hf4`` varibales to **.true.**, selects those custom basis sets from the library directory (``install_dir``) to run the HF calculations.

The keyword ``HF_CBS_default`` is kept **.true.** to use the extrapolation formula (given below) to calculate :math:`HF_{limit}` from the individual energies.

HF_CBS_default:

.. math::
       
    HF_{limit} =  \frac{ HF_{Q} - [HF_{T} * \exp{(-1.63)}] } {1 - \exp{(-1.63)}}

In case of using different basis-sets, the corresponding changes has to be done in both the methods section and the choice of extrapolation equation. The ``HF_CBS_orca_23_def2`` and ``HF_CBS_orca_34_def2`` variables select the following extrapolation equations to determine :math:`HF_{limit}` when *Karlsruhe* basis sets of size Double :math:`\zeta`::Triple :math:`\zeta` or Triple :math:`\zeta`::Quadruple :math:`\zeta` levels are used.

HF_CBS_orca_23_def2

.. math::
  
    HF_{limit} =  \frac{ HF_{T} - [HF_{D} * (\exp{(-10.39)*(\sqrt{3} - \sqrt{2})})]}  {1 - (\exp{(-10.39)*(\sqrt{3} - \sqrt{2})})} 

HF_CBS_orca_34_def2

.. math::
  
    HF_{limit} =  \frac{ HF_{Q} - [HF_{T} * (\exp{(-7.88)*(\sqrt{4} - \sqrt{3})})]}  {1 - (\exp{(-7.88)*(\sqrt{4} - \sqrt{3})})} 

Similarly, the user can select to run the HF calculation using *correlation-consistent* basis sets of size Double :math:`\zeta`::Triple :math:`\zeta` or Triple :math:`\zeta`::Quadruple :math:`\zeta` levels, which uses the extrapolation equations respectively as given below. 
      
HF_CBS_orca_23_cc

.. math::
  
    HF_{limit} =  \frac{ HF_{T} - [HF_{D} * (\exp{(-4.42)*(\sqrt{3} - \sqrt{2})})]}  {1 - (\exp{(-4.42)*(\sqrt{3} - \sqrt{2})})} 

HF_CBS_orca_34_cc

.. math::
  
    HF_{limit} =  \frac{ HF_{Q} - [HF_{T} * (\exp{(-5.46)*(\sqrt{4} - \sqrt{3})})]}  {1 - (\exp{(-5.46)*(\sqrt{4} - \sqrt{3})})} 

An example presenting HF equations using **def2-SVPD** (Double :math:`\zeta`) and **def2-TZVPPD** (Triple :math:`\zeta`) basis sets is given below.

.. code-block:: text

   method_hf3 = "HF def2-SVPD "
   basis_hf3 = ""
   custombasis_hf3 = .false.
  
   method_hf4 = "HF def2-TZVPPD "
   basis_hf4 = ""
   custombasis_hf4 = .false.

   HF_CBS_default      = .false.
   HF_CBS_orca_23_def2 = .true.
   HF_CBS_orca_34_def2 = .false.
   HF_CBS_orca_23_cc   = .false.
   HF_CBS_orca_34_cc   = .false.


The G4(MP2)XP protocol uses `RIJK` approximation to faciliate the HF calculations while using the custom basis-sets, as shown below

.. code-block:: text

   method_hf3 = "HF RIJK def2/JK Printbasis "
   basis_hf3 = "GFHFB3"
   custombasis_hf3 = .true.
  
   method_hf4 = "HF RIJK def2/JK Printbasis "
   basis_hf4 = "GFHFB4"
   custombasis_hf4 = .true.

   HF_CBS_default      = .true.
   HF_CBS_orca_23_def2 = .false.
   HF_CBS_orca_34_def2 = .false.
   HF_CBS_orca_23_cc   = .false.
   HF_CBS_orca_34_cc   = .false.
 

Thermochemical Property Determination
-------------------------------------

.. code-block:: text

   calc_HF = .true.

   calc_AE = .false.

   calc_IP = .false.
   verticalIP = .false.
  
   calc_EA = .false.
   verticalEA = .false.
  
   calc_PA = .false.
  
   calc_BE = .false.
  
The user can decide to estimate one of the given thermochemical property using G4(MP2) energies, i.e., `Enthalpy of Formation` (``calc_HF``), `Atomization Energy` (``calc_AE``), `Ionization Potential` (``calc_IP``), `Electron Affinity` (``calc_EA``), `Proton Affinity` (``calc_PA``), and `Binding Energy` (``calc_BE``) and need to make only that variable to **.true.**, while keeping the rest as **.false.**.

Although, the `Ionization Energy` and `Electron Affinity` of various systems estimated in the *G3/05* dataset are `adiabatic`, the script allows the user to predict the respective `vertical` property if interested.

In case of none of the above thermochemical properties are selected, the output file (`Thermochemistry.out <https://moldis-group.github.io/Pople/Thermochemistry.html>`_ ) will only have information regarding the G4(MP2) energies of the corresponding systems.


Additional Details
------------------

.. code-block:: text

   HLCeqZERO = .false.

   SO_3rdrow_mols = .true.

   conv_scf = "VeryTight"

   SCFDIIS  = .false.
   SOSCF    = .false.
  
   LSHIFT  = .false.
  
   optdiis = .false.
  
   restart_cc  = .false.
   restart_mp2 = .false.
   restart_hf3 = .false.
   restart_hf4 = .false.

The G4(MP2) protocol uses additional semi-empirical parameters---*Higher Level Correction* (**HLC**) to estimate the energy of the system. However, the user can choose to compute the energy as well as the desirable thermochemical property without the inclusion of the empirical correction. For that, they need to keep the variabe ``HLCeqZERO`` to **.true.**.

``SO_3rdrow_mols`` is used for specific 3rd-row atomic molecules in *doublet* state and include the assigned *spin-orbit* correction term to the total G4(MP2) energy. The molecules are **BrO**, **SeH**, **KBr+**, **AsH+**, **HBr+**, **BrF+**, **NaBr+**, and **Br2+**.

The rest of the variables are used to manipulate the SCF calculations in the various steps involved in the composite procedure. The variable ``conv_scf`` is used to assign the convergence tolerance for the SCF determination and in-case of failure of SCF converegence, the tolerance can be reduced.

ORCA depends upon the robust `Direct Inversion in Iterative Subspace` (*DIIS*) technique for finding convergence for SCF equation which in turn considers a set no. (usually 5-7) of Fock matrices. In case convergence is not achieved, the user can decide to increase the no. of Fock matrices (to 15) by keeping ``SCFDISS`` to **.true.**. The `Second Order SCF` (*SOSCF*) leads to faster SCF convergence and is useful when *DIIS* is stuck. Switching ``SOSCF`` to **.true.** will set the maximum iteration in SOSCF to 12. ``LSHIFT`` is used to modify the `levelshift` for faster convergence.

The `Coupled-Perturbed SCF` (*CPSCF*) equations are solved during second derivative calculations (required for vibrational frequencies) and the default solver is **Pople**. However, in case of troubled convergence, the user can swith the option ``optdiis`` to **.true.**, for using **DIIS** solver for sloving *CPSCF* equations.

**Question: are these options only for optimization or all the single points?**

The rest set of variables, i.e., ``restart_cc``, ``restart_mp2``, ``restart_hf3``, and ``restart_hf4`` are used when CCSD(T), MP2, HF/Triple :math:`\zeta` and(or) HF/Quadruple :math:`\zeta` calculation(s) end up with error. Switching the respective variables to **.true.** will launch an additional HF calculation with the respective basis-set using the positively charged molecular species prior to the original single point calculation(s). The optimized SCF is then used for relaunching the previously failed single point calculation(s).
