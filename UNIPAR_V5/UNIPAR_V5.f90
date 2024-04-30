!  UNIPAR_V5.f90 
!**************************************************************************************************
! UNIPAR_V5.f90 copyrightable work product created by faculty at the University of Florida. 
!
! © 2024 the University of Florida Research Foundation, Inc. All commercial rights reserved. All publications and presentations referencing the UF works should include
!   the University of Florida Research Foundation, Inc. as the copyright holder and reference the UF creators.
!
!  PROGRAM: UNIfied Partitioning and Aerosol Phase Reaction (UNIPAR) 
!
!	UNIPAR copyrightable work product created by faculty at the University of Florida. © 2024 the University of Florida Research Foundation, Inc. All commercial rights reserved. 
!	All publications and presentations referencing the UF works should include the University of Florida Research Foundation, Inc. as the copyright holder and reference the UF creators.
!
!  PURPOSE: To simulate SOA mass based on multiphase reactions of organic species 
!  Main Program for SOA model to call Submodules for SOA_H5 (aerosol phase reactions) and SOA_AE3 (partitioning)
!  Contact Information: 
!	Dr. Myoseon Jang, University of Florida at mjang@ufl.edu
!
!	This UNIPAR code is to simulate alkane SOA. If paramters for other SOA precursors are needed,
!	please contact Dr. Jang.
!  
!  PHRCSOA is the previous version of UNIPAR, coded by Yunseok Im and Jang
!  		Im et al., "Simulation of aromatic SOA formation using the lumping model 
!  		integrated with explicit gas-phase kinetic mechanisms and aerosol-phase reactions" ACP, doi:10.5194/acp-14-4013-2014
! 
!  Extended to 51 lumping on August 1, 2016 (Isoprene SOA)
!  	 	Beardsley and Jang, "Simulating the SOA formation of isoprene from partitioning and aerosol phase reactions 
!    	in the presence of inorganics” DOI:10.5194/acp-16-5993-2016, ACP, 2016
!  Updated to dynamic features (alpha, MW, O:C, and H-bonding )
!		“Simulation of SOA Formation from the Photooxidation of Monoalkylbenzenes 
! 		in the Presence of Aqueous Aerosols Containing Electrolytes under Various NOx Levels” ACP, 19, 5719-5735, 2019.
!  
!	Demonstration of UNIPAR performance for various SOA precusors: 11 aromatics (toluene, benzene, xylenes, trimethylbenzenes, and naphthalene),
!  	linear and branched alkanes from C9-C24, biogneics (isoprene, terpene, sesquiterpene for reactions with OH, O3, and NO3), phenol, and cresol
!
!  		Choi and Jang "Dual roles of inorganic aqueous phase on SOA growth from benzene and phenol" egusphere-2023-2461, 2024
!
!		Madhu et al., “Modeling the influence of carbon branching structure on SOA formation via multiphase reactions of alkanes” 
!  		ACP, https://doi.org/10.5194/egusphere-2023-1500, 2024
! 		
!		Madhu, et al., “Modeling the influence of a chain length on SOA formation via multiphase reactions of alkanes” 
!		ACP, 23, 1661–1675, 2023.
!		
!		Han and Jang “Modeling Diurnal Variation of Biogenic SOA Formation via Multiphase Reaction of Biogenic Hydrocarbons” 
!		acp-2022-327, ACP, 23, 1209–1226, 2023.
!
!		Choi and Jang "Suppression of the phenolic SOA formation in the presence of electrolytic inorganic seed." 
!		STOTEN (doi.org/10.1016/j.scitotenv.2022.158082), 851, 158082, 2022
!
!		Han and Jang, “Prediction of Secondary Organic Aerosol from the Multiphase Reaction of Gasoline Vapor by Using Volatility–Reactivity Base Lumping” 
!		ACP, 22, 625–639, DOI: https://doi.org/10.5194/acp-22-1-2022, 2022
!
!		Yu et al. “Simulation of monoterpene SOA formation by multiphase reactions using explicit mechanisms” 
!		ACS Earth and Space Chemistry, 5, 1455–1467, 2021 
!		
!	Demonstration of UNIPAR in regional scale
!		Yu,et al. "Secondary Organic Aerosol Formation via Multiphase Reaction of Hydrocarbons in Urban Atmosphere Using the CAMx Model Integrated with the UNIPAR model", 
!		doi.org/10.5194/acp-2021-1002, 22, 9083-9098, 2022
!	
!		Jo, et al., "Multiphase reactions of hydrocarbons into an air quality model with CAMx-UNIPAR: Impacts of humidity and NOx on secondary organic aerosol formation in the Southern USA."
! 		Journal of Advances in Modeling Earth System, DOI: 10.22541/au.170629218.86079671/v1, 2024
!	
!		Jo, et al., “CAMx-UNIPAR Simulation of SOA Mass Formed from Multiphase Reactions of Hydrocarbons under Urban Atmospheres in Northern California”, 
!		https://doi.org/10.5194/acp-24-487-2024, 24, 487–508, 2024
! 
!***************************************************************************************************

    program UNIPAR_V5
	use SOA_AE3
	use SOA_H5
	implicit none
	
	Character (len = 4) precursor			  ! index of target aromatics for input parameter
								  ! AR01: Toluene   
                                  ! AR02: Ethylbenzene  
                                  ! AR03: Propylbeneze
							      ! AR04: o-Xylene
                                  ! AR05: m-Xylene      
                                  ! AR06: p-Xylene, 
								  ! AR07: 135TMB    
                                  ! AR08: 123TMB 	    
                                  ! AR09: 124TMB
								  ! AR10: Benzene 
								  ! AA01: updated Benzene with phenol
								  ! AA02: updated toluene with cresol
                                  ! TA01: a_pinene with autoxidation (PRAM)
                                  ! TA02: b_pinene with autoxidation (PRAM)     
                                  ! TA03: d_limonene with autoxidation (PRAM)
                                  ! TP01: a_pinene   
                                  ! TP02: b_pinene       
                                  ! TP03: d_limonene
                                  ! IS01: isoprene
                                  ! SP01: Sesquiterpene
								  ! AU##: Linar alkane unified parameters (## represents carbon length)
								  ! BU##: Branched alkane unified parameters (## represents carbon length. Specific structure labeled inside file)
                                  ! BB01: phenol
                                  ! BB02: o-cresol
                                  ! BA01: updated phenol (aqueous reaction included)
								  ! BA02: updated o-cresol (aqueous reaction included)
								  ! BAO3: catechol
								  ! PA01: Napthalene


    Real C_number			      ! carbon # to calculate HCs concentration [ppbC] in aging scale

!********************************************************************************
!*                       Declaration of all variables							*
!********************************************************************************
	Integer Max_array			  ! array size
		parameter(Max_array = 500)! 500*DT=1500 min = 25 hour	
	Integer N_cyc				  ! Number of computational cycle (load from the input file)
	Integer DT					  ! Time step[min] 
		parameter(DT = 6)		  ! Previous time DT=3 based on SMPS. Since Feb/2018 DT=6
	Integer i,j,k                   ! j: computional cycle index (time variable); k: output time index
	Integer NCVAP				  ! number of SVOC species that partition
		parameter(NCVAP = 51)
!	CHARACTER (len = 5) space     ! 12
!        parameter(space = "  ")   ! predefine space to delimit output data 
	CHARACTER (len = 99) input_file 
								  ! 15
	CHARACTER (len = 99) output_file
								  ! 22
    CHARACTER (len = 99) log_file                
	REAL Test(0:Max_array, NCVAP)
	REAL frac_OMHi
	
! *** Environmental variables ***				
	Integer in_phase_status       ! inorganic phases: 0-->liquid phase; 1-->solid phase (default status: 0)
	REAL RH(0:Max_array)	      ! relative humidity [0 ~ 1]
	REAL ERH(0:Max_array)	      ! Efflorscence relative humidity [0 ~ 1]
	REAL DRH(0:Max_array)	      ! Deliquescence relative humidity [0 ~ 1]   ! assuming DRH = ERH + 0.50  
	REAL temp(0:Max_array)	      ! temperature [K]
 	REAL P_amb			          ! air pressure [Pa]
		parameter (P_amb = 101325)
    REAL dROG(0:Max_array)	      ! (Delta ROG) concentration of reacted VOCs during time step[ppb] 
	REAL RO2(0:Max_array)		  ! RO2 concentration of during time step [ppb]
	REAL HO2(0:Max_array)		  ! HO2 concentration of during time step [ppb]
	REAL HC_initial		          ! Initial hydrocarbon (HC) concentration [ppb]
	REAL aging_Scale_fresh		  ! lower boundary of Log((HO2+RO2)/ROG) determined at 9:30AM using the sunlight on June/20/2012 (i.e., -10.24 for ethylbenzene) 
	REAL aging_Scale_aged		  ! upper bondary of Log((HO2+RO2)/ROG) determined at 4:00 pm using the sunlight on June/20/2012 (i.e., -3.012 for ethylbenzene) 
	REAL aging_Scale_current(0:Max_array)	
							      ! current(HO2+RO2)/ROG value in log scale
	REAL aging_Scale(0:Max_array) ! [(current - fresh)/(aged - fresh)], scale: 0~1
	REAL R_V1N_LNOx			  
		parameter (R_V1N_LNOx = 15)
								  ! need to be updated, since different HC has different upper and lower limits
								  ! low NOx condtion: upper bound of the VOC/NOx ratio (25) for the calculation of the NOx parameter at urban areas !originally: 25
	REAL R_V1N_HNOx 			  
		parameter (R_V1N_HNOx = 2)
								  ! high NOx condition: lower bound of the VOC/NOx ratio (2) for the calculation of the NOx parameter 
	REAL R_V1N				   	  ! VOC:NOx ratio 								  
	REAL NOx_Scale  			  ! [(LNOx - current)/(LNOx - HNOx)], HC/NOx scale: 0~1, applied to the calculation of the dynamic MW and O:C array 	

! ** HC ppb/NOx levels for precursor parameters
    Integer NOx_Levels			  ! Number of NOx levels
		parameter(NOx_Levels = 9) 	
    REAL HC_NOx_Level(NOx_Levels)
    REAL Aging_Scale_Upper(NOx_Levels)  ! aging scale for the aged condition at a given NOx level
    REAL Aging_Scale_Lower(NOx_Levels)  ! aging scale for the fresh condition at a given NOx level
    
! ** Molecular weight **
	REAL MW_LNOx_fresh(NCVAP)     ! MW array at the low NOx under the fresh condition 
	REAL MW_HNOx_fresh(NCVAP)     ! MW array at the high NOx under the fresh condition 
	REAL MW_fresh(NCVAP)     	  ! MW array at a given NOx under the fresh condition	
	REAL MW_LNOx_aged(NCVAP)	  ! MW array at the low NOx under the aged condition
	REAL MW_HNOx_aged(NCVAP) 	  ! MW array at the high NOx under the aged condition
	REAL MW_aged(NCVAP)     	  ! MW array at a given NOx under the aged condition		
	REAL MW(NCVAP)                ! dynamic MW of each lumped groups [g/mol]
	REAL MW_save(0:Max_array,NCVAP)
								  ! Saving dynamic MW of each lumped groups at the end of main loop (this sort of set up prevents the excess change on the other two modules) 
	REAL MW_m1(NCVAP)   		  ! Inverse SVOC MW's [ mol / g ]        
	REAL MW_or(0:Max_array)       ! Average Molecular weight of organic aerosol (OMP+OMH)
	REAL MW_omp(0:Max_array)      ! Average Molecular weight of OMP (monomer)
	REAL MW_omh(0:Max_array)      ! Average Molecular weight of OMH (oligomer)
        real mw_rog               ! Molecular weight of reactive organic gases
!	REAL mwoldOMP                 ! Molecular weight of pre-exisiting SOA [ g / mol ]
!       parameter (mwoldOMP = 200)! *************** fix this value
!	REAL mwOMH                    ! molecular weight of pre-existing hetero mass [g/mol]
!	    parameter (mwOMH = 300)

! ** Oxygen to carbon (O:C) ratio **
	REAL O1C_LNOx_fresh(NCVAP)    ! O:C array at the low NOx under the fresh condition 
	REAL O1C_HNOx_fresh(NCVAP)    ! O:C array at the high NOx under the fresh condition 
	REAL O1C_fresh(NCVAP)     	  ! O:C array at a given  NOx under the fresh condition 	
	REAL O1C_LNOx_aged(NCVAP)	  ! O:C array at the low NOx under the aged condition 
	REAL O1C_HNOx_aged(NCVAP) 	  ! O:C array at the high NOx under the aged condition 
	REAL O1C_aged(NCVAP)     	  ! O:C array at a given  NOx under the aged condition 
	REAL O1C(NCVAP)	 			  ! dynamic O:C ratio of each lumped group 
	REAL O1C_save(0:Max_array,NCVAP)
								  ! Saving dynamic O:C ratio of each lumped groups at the end of main loop	
	REAL O1Ci(NCVAP)			  ! sum up molar fraction weighted O:C ratio for all the lumping species   
    REAL O1CHi(NCVAP)             ! sum up molar fraction weighted O:C ratio for OMH 
	REAL O1C_SOA(0:Max_array)	  ! O:C ratio of aerosol ( sum up(O1C(i))/MW of SOA )
    REAL O1C_OMH(0:Max_array)	  ! O:C ratio of OMH ( sum up(O1CH(i))/MW of OMH )
    REAL MW_SOA(0:Max_array)
	
! ** Hydrogen Bonding for the calculation of the acitivity coeff of organic in the inorganic phase **  PS: -OH = 1; -C(=O)OH = 1.6
	REAL HBonding_LNOx_fresh(NCVAP)   ! H-bonding array at the low NOx under the fresh condition 
	REAL HBonding_HNOx_fresh(NCVAP)   ! H-bonding at the high NOx under the fresh condition 
	REAL HBonding_fresh(NCVAP)     	  ! H-bonding at a given  NOx under the fresh condition 	
	REAL HBonding_LNOx_aged(NCVAP)	  ! H-bonding at the low NOx under the aged condition 
	REAL HBonding_HNOx_aged(NCVAP) 	  ! H-bonding at the high NOx under the aged condition 
	REAL HBonding_aged(NCVAP)     	  ! H-bonding at a given  NOx under the aged condition 
	REAL HBonding(NCVAP)	 		  ! dynamic H-bonding of each lumped group 
	REAL HBonding_save(0:Max_array,NCVAP)
								      ! Saving dynamic H-bonding of each lumped groups at the end of main loop	
! ** coefficents in polynomial equations to produce stoichimetric arrays at different NOx and againg conditions
    REAL coeff1_LNOx_Fresh(NCVAP)
    REAL coeff2_LNOx_Fresh(NCVAP)
    REAL coeff3_LNOx_Fresh(NCVAP)
    REAL coeff4_LNOx_Fresh(NCVAP)
    REAL coeff1_LNOx_Aged(NCVAP)
    REAL coeff2_LNOx_Aged(NCVAP)
    REAL coeff3_LNOx_Aged(NCVAP)
    REAL coeff4_LNOx_Aged(NCVAP)
    
    REAL coeff1_HNOx_Fresh(NCVAP)
    REAL coeff2_HNOx_Fresh(NCVAP)
    REAL coeff3_HNOx_Fresh(NCVAP)
    REAL coeff4_HNOx_Fresh(NCVAP)
    REAL coeff1_HNOx_Aged(NCVAP)
    REAL coeff2_HNOx_Aged(NCVAP)
    REAL coeff3_HNOx_Aged(NCVAP)
    REAL coeff4_HNOx_Aged(NCVAP)
    Character(LEN = 6) Header_NCVAP

	REAL prod_mole_conc(NCVAP) 	  ! mole concentration of the species in aerosol phase	
    REAL OMH_mole_conc(NCVAP)     ! mole concnetration of the species in OMH
	REAL act_or_in(NCVAP)		  ! activity coefficients of lumping species in inorganic aerosol 	(using semi-emperical model compound from AIOMFAC)
	REAL act_or_in_save(0:Max_array,NCVAP)
								  ! Saving dynamic activity coeff of organics in inorganic phase of each lumped groups at the end of main loop	

!******************************************************************************
!    		    Parameters for Heterogeneous Mass Estimation				  *
!******************************************************************************
	REAL alpha_fresh(NCVAP)
								  ! stoiciometric coeff. of fresh gas pahse at each step as a function of VOC/NOx ratio
	REAL alpha_aged(NCVAP)	      
								  ! stoiciometric coeff. of aged system at each step as a function of VOC/NOx ratio
	REAL alpha_dyn(0:Max_array, NCVAP)	          
								  ! dynamic stoiciometric coeff. calculated from combining fresh and aged alpha fractioned with aging_Scale  
	REAL alpha_auto(0:Max_array, NCVAP)	! dynamic stoiciometric coeff. calculated from autoxidation and MCM
    REAL Rho_om			          ! ensity of organic matter (OMT)
		parameter (Rho_om = 1.38)		
	REAL VP(NCVAP)		          ! saturated vapor pressure in [mmHg]
	REAL Cin(NCVAP)		          ! solubilized SVOC concentration in inorganic phase[ug/m3]
	REAL Cin_sum(0:Max_array)	  ! Sum-up of solubilized SVOC concentration in inorganic phase[ug/m3]
	REAL OMo					  ! Preexisting organic particular matter [ug/(m**3)]  			
	REAL dOMHi(NCVAP)	          ! individual OMH formed during time interval (DT) [ug/( m**3 )]
	REAL dOMH(0:Max_array)        ! OMH formed during DT: sum of dOMHi(i) at a given time step[ug/( m**3 )]
	REAL OMHi(0:Max_array, NCVAP) ! Accumulated individual OMH components after time step [ug/(m**3), OMHi(j,i)]
	REAL OMH(0:Max_array)	      ! Accumulted OMH after time step [ug/( m**3 ), OMH(j)]: sum of OMHi(i)		
	REAL M_in(0:Max_array)	      ! Mass of inorganic aerosols [ug/m3]
	REAL in_mol2ug(NCVAP)         ! Unit conversion factor for Cin [ug/m3 of air per mol/L of organic in inorg aerosol] 
	REAL FS(0:Max_array)	      ! acidity descriptor for inorganic aerosol, FS[H+]=(SO4/(SO4+NH4))
	REAL C_proton(0:Max_array)    ! proton concentration in inorganic aerosol [mol/L]   
	REAL SO4(0:Max_array)		  ! Total sulfate concentration [umol/m3 of air]
	REAL SO4_OS(0:Max_array)      ! Consumed Sulfate concentration by OS formation	[umol/m3 of air]
	REAL OS(0:Max_array)		  ! Accumulated consumed Sulfate concentration by OS formation for each timestep [ug/m3 of air]
	REAL NH4(0:Max_array)		  ! Total NH4+ concentraion [umol/m3 of air] in aerosol
	REAL dSO4(0:Max_array)		  ! SO4 formed during DT [umol/m3 of air]    ! use dSO4(0) and dNH4(0) to represent the initial concentration
	REAL dNH4(0:Max_array)		  ! NH4+ formed during DT [umol/m3 of air] in aerosol
	REAL SO4_free(0:Max_array)	  ! Sulfate ion which can participate in Organosulfate formation (umol/m3)
	REAL N_OS					  ! Sum of # of Functional group associated with organosulfate formation
								  ! Weighing facgtor:
								  ! Medium(2),Fast(4), MultiAlcohol(4([least volatile],4[2nd least vol],3[3rd least vol]) 
                                  ! weighing factor for fast should be reduced due to polymerization.    
    REAL dN_OS
    REAL f_OS					  ! Factor for acidity drop by Organosulfate formation
		parameter (f_OS = 0.04)   ! 0.017, 0.071 (CF and Yun). 0.04 (04222021 updated by testing various SOA by Jang )
	REAL OS_Convfct				  ! Conversion factor for OS formation.
								  ! remained SO4 ratio after OS formation.
	REAL X_water_in(0:Max_array)  ! Mass fraction of water in inorganic aerosol; water/(water+salt)
	REAL M_salt(0:Max_array)	  ! Mass of salt (ug/m3)
	REAL M_water(0:Max_array)	  ! Mass of water (ug/m3)	
	REAL MW_in					  ! Molecular weight of inorganic phases [g/mol]
	REAL V_in 					  ! Volume of inorganic aerosol [cm3]
	REAL Vin			    	  ! Molar volume of inorganic phase.[cm3/mol]
	REAL Rho_in					  ! density of inorganic phase.[g/cm3]
	REAL beta1(0:Max_array, NCVAP)! coefficient1 for OMH calculation
	REAL beta2(0:Max_array, NCVAP)! coefficient2 for OMH calculation 
	REAL kac(0:Max_array, NCVAP)  ! kH*M_in/1000+ko
	REAL Kor(0:Max_array, NCVAP)  ! partioning coeff. between organic aerosol and the gas phase
	REAL Kin(NCVAP)		          ! partioning coeff. between inorganic aerosol and gas	
	REAL Kin_save(0:Max_array,NCVAP)		          
								  ! partioning coeff. between inorganic aerosol and gas
	REAL org1sulf            	  ! Bulk organic (in inorganic phase) to inorganic sulfur mass ratio (org:sulf) is equal to OMH:SO4_free    
	
!******************************************************************************
!    			 Parameters for Partitioning Mass Estimation				  *
!******************************************************************************
	INTEGER LOGDEV
	INTEGER LAYER                 ! model layer number
    INTEGER NPSPCS                ! number of aerosol precursor species
    REAL OMP(0:Max_array)         ! new value of OMP after time step [ug/( m**3 )]
	REAL OMPi(0:Max_array, NCVAP) ! new values of OMP components after time step [ug/( m**3 )]
	REAL OMT(0:Max_array)         ! Total OMT = OMH + OMP
    
! *** Define the descriptive numbers ***
	PARAMETER (Layer=1, NPSPCS=1, LOGDEV=3)

! *** module containing all information and subroutines used in SUBROUTINE NEWT ***
    logical first_time
    data first_time / .true. /
    save first_time

! *** Gas/aerosol partitioning parameters ***
    real cstar(NCVAP)	          ! Effective saturation concentrations of SVOC's [ ug / m**3 ] at 298 K
						          ! [ ug / m **3 ] / [ ug / m**3 ]
    real rgas1                    ! reciprocal of universal gas constant = 1/R
        parameter( rgas1 = 1.0/8.314510 )
    real hvap(NCVAP)	          ! h_vap: enthalpy of vaporization [ J / mol ]156. before correction		
	real hvap_c(NCVAP)	          ! h_vap: enthalpy of vaporization [ J / mol ]156  after correction
	real fac_a		          	  ! enthalphy correctionfactor
		parameter(fac_a = 1.5)
	real fac_b                 	  ! enthalphy correctionfactor
		parameter(fac_b = 0.5)
                                  ! fixed value of h_vap * rgas1  (enthlpy of vap/R)
                                  ! parameter( hfac = 40.0e3 * rgas1 )						        

! *** Reference temperature and pressure *** 
    real T_ref                    ! reference temperature
		parameter( T_ref = 298.0) ! [ K ]
    real T_refm1                  ! inverse of reference temperature
		parameter( T_refm1 = 1.0 / T_ref )
    real P_ref                    ! reference pressure [ Pa ]
        parameter( P_ref = 101325.0 ) !760/101325 = 133.32 
    real P_refm1                  ! inverse of reference pressure
		parameter( P_refm1 = 1.0 / P_ref )

! *** Unit conversion factors at reference temperature and pressure ***
    real dROG_ppb2ug              ! [ ppb per ug/m3 ] for ROG's
        save dROG_ppb2ug
	
! *** Variables used to adjust cstar as a function of temperature ***
    real convfac_298              ! P/RT = mole/V at 1 mmHg and 298 K [mole/m**3]
		parameter( convfac_298 = P_ref * rgas1 * T_refm1 )
    real convfac, convfacm1       ! convfac:(T1/T2)*(P in ambient/P og reference pressure) 
                                  ! difference of conversion factor from standard state
                                  ! convfacm1 was in the previous SOA model
    real tt1,tt2                  ! temperature-related factors (tt1 = T1/T2, tt2=1/T1 -1/T2)
    real tempcorr(NCVAP)          ! temperature correction factor for cstar
	
! *** Variables used in organic equilibrium calculations ***
    real totrog(NCVAP)            ! dROG concentrations mapped to each SVOC [ ug / m**3 ]
    real c0(NCVAP)                ! cstar at AIRTEMP [ ug / m**3 ]
	real Ctb(0:Max_array, NCVAP)  ! SVOC conc. for each lumping group before the current time step [ug/m**3]
	real Ctf(NCVAP)               ! SVOC conc. for each lumping group after the current time step [ug/m**3]
	real Ctf_sum				  ! total SVOC conc. after the current time step [ug/m**3] (for redistribution of all SVOC)
	real ctotf(NCVAP)	          ! Ctf - Cin [ug/m**3] (SVOC only for gas/org partitioning)
    real caer(NCVAP)              ! SVOC conc in aerosol phase [ ug / m**3 ]
    real totorgsw                 ! molar concentration of POA [ u-mole / m**3 ]
    real totorg                   ! SOA + POA before time step [ u-mole / m**3 ]
    real threshold                ! criterion for establishing Gas/Part equil.
    real threshmin                ! small positive number
		parameter( threshmin = 1.0E-19)
    real conmin                   ! concentration lower limit for SVOC in aerosol
		parameter ( conmin = 1.0e-30 )
    real CTOLMIN
		parameter ( CTOLMIN = 1.E-06 )
		
! *** Variable internal to NEWT subroutine ***
    logical check                 ! flag to indicate if NEWT subroutine
                                  ! converged to a spurious root
								  
! *** Initialization of the varibles for the first cycle calculation ***							  
	N_OS = 0.00                   ! define N_OS initial value
    dN_OS = 0.00 
	OS(0) = 0.00				  ! define initial OS mass [ug/m**3]
	Ctf(NCVAP) = 0.00		      ! define Ctf initial value
	in_phase_status = 0   		  ! intial inorganic phase: liquid phase, since most exp start at high RH condition
    O1C_SOA(0) = 0.3              ! Initial O:C ratio of organic aerosol  
    O1C_OMH(0) = 0.8              ! Initial O:C ratio of OMH
    SO4_OS(0) = 0

!##############################################################################	
!     Lumping structure for all precursors and model paramters for each precusor 
!##############################################################################	
	
! Lumping Structure
	!  |  i1  i2  i3  i4  i5  i6  i7  i8  
	!-----------------------------------	
	!VF|   1   2   3   4   5   6   7   8
	! F|   9  10  11  12  13  14  15  16
	! M|  17  18  19  20  21  22  23  24
	! S|  25  26  27  28  29  30  31  32
	! P|  33  34  35  36  37  38  39  40  
	!MA|  41  42  43  44  45  46  47  48
    ! R|                      49  50  51  
    !  |                    EPOX MGLY GLY
	

! *** Vapor pressure (mmHg) ***  
	VP	= (/1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &
			1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &
			1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &
			1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &
            1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &  
            1E-08,	1E-06,	1E-05,	1E-04,	1E-03,	1E-02,	1E-01,  1.0,  &
				                                    1E-02,	 28.0,  43.0  /) !43 !MGLY 28.0 !4E-08 for 1P original
! Explicit vapor pressures for group 49-51
! Group 49 = 0.01, Group 50 = 27.6 mmHg, Group 51 = 42.7 mmHg at 298K

! *** Enthalpy of Vaporization of lumping speceis (J/mol) ***  
   hvap	= (/140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, & 
				140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, & 
				140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, &
                140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, & 
				140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, &
                140E3,	106E3,	96E3,	89E3,	82E3,	58E3,  58E3,  58E3, & 
				                                        58E3,  58E3,  58E3 /)

!	hvap	= (/140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, & 
!				140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, &
!				140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, &
 !               140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, &
!				140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, &
!                140E3,	121E3,	109E3,	98E3,	88E3,	57E3,	57E3,   57E3, &
!				                                        57E3, 	57E3,   57E3 /)

!"Zhao, L.; Ni, N.; Yalkowsky, S.H. A modification of Trouton¡¯s rule by simple molecular parameters  for hydrocarbon compounds. Ind. Eng. Chem. Res. 1999, 38, 324-327	

! Enthalpy of vaporization (hvap) correction      
	hvap_c = hvap/((hvap/hvap(38)*fac_a)**fac_b)			

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!     READ the input file 													   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! read input the file from the input_parameter file to choose precusor parameter and input data	
!    open(1, File='Input_parameter.dat')									  
!        READ (1,*)
!	    READ (1,*) precursor
!	    READ (1,*)
!	    READ (1,*) input_file
!	    READ (1,*)
!	    READ (1,*) output_file
!    CLOSE(1)

  !read input parameters fron *.job file
    READ (*,'(20x,a)') precursor
    READ (*,'(20x,a)') input_file
    READ (*,'(20x,a)') output_file
    READ (*,'(20x,a)') log_file   
    
	open(2, FILE= input_file)	  ! read data files
    ! *** Read HO2, RO2, dROG, Temp, RH, dSO4, dNH4 ***
	    READ (2,*)
	    READ (2,*) N_cyc              ! Number of iteration 
	    READ (2,*)
	    READ (2,*) OMo                ! preexsiting OM 
	    READ (2,*)
	    READ (2,*) R_V1N              ! Initial HC/NOx ratio
	    READ (2,*)		
	    READ (2,*) HC_initial         ! initial ROG concentration [ppb]
	    READ (2,*) 
	    DO j = 0, N_cyc               ! (e.g., 101 data lines and N_cyc = 100, thus j needs one more line: j = 0, can be treated as initail condition)
	        READ (2,*) HO2(j), RO2(j), dROG(j), Temp(j), RH(j), dSO4(j), dNH4(j)
								      ! aging parameter was introduced 02/17/18
								      ! ROG(j) = dROG(j)*20
	    End Do	
    CLOSE(2)

    open(3, FILE= log_file)	  ! log files
!    open(3, FILE= 'log.dat')	  ! log files
	open(4, file= output_file)	  ! write results files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Model parameters for each precursor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	open(5, File= 'SOA_PAR/SOAPAR_'//precursor//'.dat')	
								  ! read input files
                                  
        ! *** Read carbon number and MW of precusor; aging scale, MW, O:C, and H-bonding of lumped spwceis; soichiometric arrays of lumping species ***
        DO j = 0, 6               ! read 7 commend lines
	        READ (5,*) 
        end do	
    
        READ (5,*) C_Number             ! carbon number of precursor  
        READ (5,*)
	    READ (5,*) MW_rog               ! MW of precursor 
        READ (5,*)                      ! read aging factor for aged and fresh conditions of the oxidation with OH 
        READ (5,*)                      ! read VOC/NOx ratio
        
        DO j = 1, NOx_Levels
            READ (5,*) HC_NOx_Level(j), Aging_Scale_Lower(j), Aging_Scale_Upper(j)
        End DO
	   
    ! Read pysicochemical parameters of lumping speceis (MW, O:C, and H-bonding) and coefficients (a1, a2, a3, and a4) for polynomial equation to produce dynamic alpha arrays)
        READ (5,*)                      ! read comments for lumping species (LNOx_Fresh) 
        DO j = 1, NCVAP                 !! read model parameters for low NOx fresh
            READ (5,*) Header_NCVAP, MW_LNOx_Fresh(j), O1C_LNOx_Fresh(j), HBonding_LNOx_Fresh(j), &
            coeff1_LNOx_Fresh(j), coeff2_LNOx_Fresh(j), coeff3_LNOx_Fresh(j), coeff4_LNOx_Fresh(j)  
        ENDDO
    
        READ (5,*)                      ! read comments for lumping species    (LNOx_aged) 
        DO j = 1, NCVAP                 !! read model parameters for low NOx aged
            READ (5,*) Header_NCVAP, MW_LNOx_Aged(j), O1C_LNOx_Aged(j), HBonding_LNOx_Aged(j), &
            coeff1_LNOx_Aged(j), coeff2_LNOx_Aged(j), coeff3_LNOx_Aged(j), coeff4_LNOx_Aged(j)  
        ENDDO
    
        READ (5,*)                      ! read comments for lumping species (HNOx_Fresh) 
        DO j = 1, NCVAP                 !! read model parameters for high NOx fresh
            READ (5,*) Header_NCVAP, MW_HNOx_Fresh(j), O1C_HNOx_Fresh(j), HBonding_HNOx_Fresh(j), &
            coeff1_HNOx_Fresh(j), coeff2_HNOx_Fresh(j), coeff3_HNOx_Fresh(j), coeff4_HNOx_Fresh(j)  
        ENDDO
    
        READ (5,*)                      ! read comments for lumping species    (HNOx_aged) 
        DO j = 1, NCVAP                 !! read model parameters for high NOx aged
            READ (5,*) Header_NCVAP, MW_HNOx_Aged(j), O1C_HNOx_Aged(j), HBonding_HNOx_Aged(j), &
            coeff1_HNOx_Aged(j), coeff2_HNOx_Aged(j), coeff3_HNOx_Aged(j), coeff4_HNOx_Aged(j)  
        ENDDO
                              
    CLOSE(5)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!     Determine aging scales for upper and lower caps												   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!		
     IF (R_V1N .gt. (HC_NOx_Level(1)*0.5+HC_NOx_Level(2)*0.5)) then
        aging_Scale_fresh = Aging_Scale_Lower(1)
        aging_Scale_aged = Aging_Scale_Upper(1)
    ELSEIF (R_V1N .le. (HC_NOx_Level(1)*0.5+HC_NOx_Level(2)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(2)*0.5+HC_NOx_Level(3)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(2)
        aging_Scale_aged = Aging_Scale_Upper(2)
    ELSEIF (R_V1N .le. (HC_NOx_Level(2)*0.5+HC_NOx_Level(3)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(3)*0.5+HC_NOx_Level(4)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(3)
        aging_Scale_aged = Aging_Scale_Upper(3)
    ELSEIF (R_V1N .le. (HC_NOx_Level(3)*0.5+HC_NOx_Level(4)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(4)*0.5+HC_NOx_Level(5)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(4)
        aging_Scale_aged = Aging_Scale_Upper(4)
    ELSEIF (R_V1N .le. (HC_NOx_Level(4)*0.5+HC_NOx_Level(5)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(5)*0.5+HC_NOx_Level(6)*0.5)) then
        aging_Scale_fresh = Aging_Scale_Lower(5)
        aging_Scale_aged = Aging_Scale_Upper(5)
    ELSEIF (R_V1N .le. (HC_NOx_Level(5)*0.5+HC_NOx_Level(6)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(6)*0.5+HC_NOx_Level(7)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(6)
        aging_Scale_aged = Aging_Scale_Upper(6)
    ELSEIF (R_V1N .le. (HC_NOx_Level(6)*0.5+HC_NOx_Level(7)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(7)*0.5+HC_NOx_Level(8)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(7)
        aging_Scale_aged = Aging_Scale_Upper(7)
    ELSEIF (R_V1N .le. (HC_NOx_Level(7)*0.5+HC_NOx_Level(8)*0.5) .AND. R_V1N .gt. (HC_NOx_Level(8)*0.5+HC_NOx_Level(9)*0.5 )) then
        aging_Scale_fresh = Aging_Scale_Lower(8)
        aging_Scale_aged = Aging_Scale_Upper(8)
    ELSE
        aging_Scale_fresh = Aging_Scale_Lower(9)
        aging_Scale_aged = Aging_Scale_Upper(9)   
    ENDIF

! *** Initial guess of Molecular weight of organic aerosol ***
	
    MW_omp(0) = 200 
	MW_omh(0) = 250
	MW_or(0) = 250	
    
! ***  upper and lower  alpha arrays at two aging levels as a functional of NOx levels ***
    DO i = 1, NCVAP  
        
        IF (R_V1N .ge. 14) then     ! HC/NOx > 14
            !calculate alpha_fresh array
            alpha_fresh(i)= coeff1_LNOX_fresh(i)*R_V1N**3+coeff2_LNOX_fresh(i)*R_V1N**2 &
			+coeff3_LNOX_fresh(i)*R_V1N+coeff4_LNOX_fresh(i)
            !calculate alpha_aged array
            alpha_aged(i)= coeff1_LNOX_aged(i)*R_V1N**3+coeff2_LNOX_aged(i)*R_V1N**2 &
			+coeff3_LNOX_aged(i)*R_V1N+coeff4_LNOX_aged(i)
        ELSE        ! HC/NOx < 14
            !calculate alpha_fresh array
            alpha_fresh(i)= coeff1_HNOX_fresh(i)*R_V1N**3+coeff2_HNOX_fresh(i)*R_V1N**2 &
			+coeff3_HNOX_fresh(i)*R_V1N + coeff4_HNOX_fresh(i)
            !calculate alpha_aged array
            alpha_aged(i)= coeff1_HNOX_aged(i)*R_V1N**3+coeff2_HNOX_aged(i)*R_V1N**2 &
			+coeff3_HNOX_aged(i)*R_V1N + coeff4_HNOX_aged(i)
        ENDIF
    ENDDO
 write (3,*) alpha_fresh(i), alpha_aged(i)
!#####################################################################################
!  Estimation of dynamic aging scales
!#####################################################################################

DO j=1, N_cyc					  ! j index represents the time step, i index represents 51 lumping group, j starts from 1 

! *** calculate aging parameter (aging_Scale) ***
	
    if ((HO2(j) + RO2(j)) .le. 0.000) then
        aging_Scale_current(j) = 0.000        
    else 
        aging_Scale_current(j) = log10((HO2(j) + RO2(j)) / (HC_initial * C_number))
                                      ! the accumulated [HO2] and [RO2] to represent the aging status of the reaction system 			
	end if
    							  
	if (aging_Scale_current(j) .lt. aging_Scale_fresh) then
		aging_Scale_current(j) = aging_Scale_fresh
									  ! set up the low boundary
	else if (aging_Scale_current(j) .gt. aging_Scale_aged) then
		aging_Scale_current(j) = aging_Scale_aged
	                                  ! set up the high boundary
	end if
	
	aging_Scale(j) = (aging_Scale_current(j) - aging_Scale_fresh)/(aging_Scale_aged - aging_Scale_fresh) 
                                      ! calulate the aging parameter, scale: 0~1

! sensitivity test of the aging
!	aging_Scale(j) = 0 				  ! all the time fresh
!	aging_Scale(j) = 1 				  ! all the time aged

!###################################################################################################
!  Estimation of dyanmic alpha(N_Cyc,NCVAP), MW(N_Cyc,NCVAP), O:C(N_Cyc,NCVAP), HBonding(N_Cyc,NCVAP)
!###################################################################################################

    do i = 1, NCVAP          
        alpha_dyn(j,i) = alpha_fresh(i)*(1-aging_Scale(j)) + alpha_aged(i)*aging_Scale(j)
            if (alpha_dyn(j,i).lt.0.00) then 
            alpha_dyn(j,i) = 0.00
            end if    
    end do
	
! *** calculate NOx scale (NOx_Scale) ***
	if (R_V1N .lt. R_V1N_HNOx) then
		R_V1N = R_V1N_LNOx
									  ! set up the low boundary
	else if (R_V1N .gt. R_V1N_LNOx) then
		R_V1N = R_V1N_LNOx
	                                  ! set up the high boundary
	end if	

	NOx_Scale = (R_V1N_LNOx - R_V1N)/(R_V1N_LNOx - R_V1N_HNOx)
									  ! 1-(current-HNOx)/(LNOx-HNOx), scale: 0~1
							  
! *** calculate dyanmic MW and dynamic O:C ratio for each time step ***
	!dynamic MW calculation 
	MW_fresh(:) =  MW_LNOx_fresh(:)*(1-NOx_Scale) + MW_HNOx_fresh(:)*NOx_Scale
	MW_aged(:) = MW_LNOx_aged(:)*(1-NOx_Scale) + MW_HNOx_aged(:)*NOx_Scale	
	MW(:) = MW_fresh(:)*(1-aging_Scale(j)) + MW_aged(:)*aging_Scale(j)
	
    !dynamic O:C calculation 
	O1C_fresh(:) =  O1C_LNOx_fresh(:)*(1-NOx_Scale) + O1C_HNOx_fresh(:)*NOx_Scale
	O1C_aged(:) = O1C_LNOx_aged(:)*(1-NOx_Scale) + O1C_HNOx_aged(:)*NOx_Scale
	O1C(:) = O1C_fresh(:)*(1-aging_Scale(j)) + O1C_aged(:)*aging_Scale(j)

    !dynamic H-bonding calculation 
	HBonding_fresh(:) =  HBonding_LNOx_fresh(:)*(1-NOx_Scale) + HBonding_HNOx_fresh(:)*NOx_Scale
	HBonding_aged(:) = HBonding_LNOx_aged(:)*(1-NOx_Scale) + HBonding_HNOx_aged(:)*NOx_Scale
	HBonding(:) = HBonding_fresh(:)*(1-aging_Scale(j)) + HBonding_aged(:)*aging_Scale(j)
	
	do i = 1, NCVAP                	  ! Force the MW to be 1 when it smaller than 1 and O:C ratio to be 1.2 when it greater than 1.2.
		if (MW(i) .lt. 1.00) then 
			MW(i) = 1.00
		end if
		
		if (O1C(i) .gt. 1.2) then	  ! The upper boundary of O:C -- 1.2, determined from the AIOMFAC regression for 26 organic compunds
			O1C(i) = 1.2
		else if (O1C(i) .le. 0.0) then
			O1C(i) = 0.001   		  ! The lower boundary of O:C -- if group is not exist, set up its O:C to be 0.001 
		end if						  ! (very nonpolar, large a_coeff, small Kp --> negeligible contribution to SOA)
		
		if (HBonding(i) .lt. 0.000) then 
			HBonding(i) = 0.000
		end if		
    end do

!#####################################################################################
!  Estimation of total organic concentrations : Ctb (j), Tempcorr(j), and c0(j)
!#####################################################################################    
    
! *** Set unit conversion and inverse MW
    dROG_ppb2ug =  mw_rog * convfac_298/1000  !from ppb to (ug/m3)   10**-9 MW*10**6 (P/RT)
    MW_m1(:) = 1.0 / MW(:)
	
! ***	The calculation of saturated concentration (cstar, ug/m3) of each lumping group 
	cstar = VP(:) * MW(:) * rgas1 * T_refm1 * 10**6*133.32 !ALN P Pa

! *** set temperature factors ***
	tt1 = T_ref / temp(j-1)           ! T1/T2
	tt2 = T_refm1 - 1.0 / temp(j-1)	  ! 1/T1 -1/T2
	convfac = tt1 * P_amb * P_refm1   ! (T1/T2)*(P in ambient/P og reference pressure) Conc correction based on stand.state
      
! *** set SVOC concentrations  [ ppb ] -> [ ug / m**3 ]
!     and initialize Ctb from SVOCS
 
   !Saturated gas concentration of SVOC
	tempcorr = tt1 * exp( hvap_c * rgas1 * tt2 )   ! K1/K2 = T1/T2 * exp(delH/RT1)/(exp(delH/RT2))  
	c0 = cstar * tempcorr						   ! c0 = 1/Kp only for SOA growth when activity coeff =1    

!******************************************************************************** 
! estimation of the initial SVOC available for the SOA formation.
! Ctb is estimated in the beggining as shown below (to produce OMH) 
! Ctf is to estimate in the previous step after OMH (call subroutine HET)
!********************************************************************************

	Ctb(j,:) = Ctf(:) + alpha_dyn(j,:)*dROG(j)*dROG_ppb2ug*convfac !Total SVOC    !Ctf ---> Ctf(:)   
!*** test *** summing up all of the Ctf and new dROG, and then redistributed by alpha_dyn      
	!Ctf_sum = SUM(Ctf(:))
	!Ctb(j,:) = alpha_dyn(j,:)* ( dROG(j)*dROG_ppb2ug*convfac + Ctf_sum )
	
!#####################################################################
!                    Heterogeneous SOA (OMH)                         #
!#####################################################################

!**********************************************************************
!                    FS, ERH(j), and DRH(j) calculation 
!**********************************************************************
! FS is a numerical descritor and used to semiempirically describe other thermodynamic parameters together with RH.
! For example, aerosol water content, water activity in aerosol or ERH as function of RH and FS

	NH4(j) = NH4(j-1) + dNH4(j) !umol/m3
	SO4(j) = SO4(j-1) + dSO4(j) !umol/m3

	if (SO4(j) .le. 0) then	! NO inorganic aerosol case
    
		FS(j) = 0.33
		M_in(j) = 0
		Rho_in = 1
		MW_in = 18
		Vin = 18
	
	else					! In the presence of inorganic aerosol
	    FS(j) = SO4(j) / (SO4(j)+NH4(j))	! Update FS for describing the composition of the inorganic phase.
!print *, FS(2)
!pause			
        ERH(j) = -4.8536*FS(j)**2 + 2.4022*FS(j) + 0.0509
		    ! ERH regression using FS for the SO4-NH4-H2O aerosol. 
		    !references
		    !(1) Atmos. Chem. Phys. Discuss., 2, 2449-2487, 2002
		    !    A novel model to predict the physical state of atmospheric H2SO4/NH4/H2O aerosol particles
		    !(2) Spann and Richardson(1985)
		    !    http://www.sciencedirect.com/science/article/pii/0004698185900721

!******************************************************************************
!*          Aerosol Acidity, aerosol water content  and OS formation          *
!******************************************************************************    		  
		if (RH(j) .gt. ERH(j).AND. in_phase_status .eq. 0) then		! modified for inorganci phase status
			
			SO4_free(j) = SO4(j) - 0.5*NH4(j)
			
			if (SO4_free(j) .le. 0 ) then !(NH4)2SO4
				SO4_free(j) = 0
				OS_Convfct = 0         ! Conversion factor for OS formation
    
			else
				!OS_Convfct = 1 - 1/(1 + f_OS*(N_OS/SO4_free(j)))  ! OS_Convfct = OS ratio (0 - 1). = OS(j)/SO4(j-1)
                OS_Convfct = 1 - 1/(1 + f_OS*(dN_OS/SO4_free(j)))  ! OS_Convfct = OS ratio (0 - 1). = OS(j)/SO4(j-1)
																  ! 0 = NO OS formation; All SO4 remained
																  ! 1 = 100 % OS formation; NO SO4 remained
			!! f_OS is the coefficient for the OS formation related to the ratio of the total reactive functional group to free SO4                                                   
            end if		
			
			SO4_OS(j) = SO4_free(j)*OS_Convfct  ! Concentration of the consumed SO4 by OS formation [umol/m3 of air]
			
			if (SO4_OS(j) .gt. SO4_free(j)) then!restriction: if (SO4_OS > SO4_free)
				SO4_OS(j) = SO4_free(j) 			
			end if
			
			SO4(j) = SO4(j) - SO4_OS(j)         !umol/m3, the sulfate that is not related to OS
			
			FS(j) = SO4(j) / (SO4(j)+NH4(j))	!updating inorganic composition factor.

			!***set the FS upper and lower boundarys
			if (FS(j) .le. 0.3333) then
				FS(j) = 0.3333		 		! constrain water mass not less than 0.3333		
			else if (FS(j) .ge. 1) then
				FS(j) = 1 					! constrain water mass not greater than 1		
			end if	
			
! *** Mass fraction of water in inorganic aerosol ***        
!			if (FS(j) .gt. 0.55) then	! regression of water content in the E-AIM II thermodynamic model using RH and FS (RH 0.1 ~ 0.9, temp = 298.15K)
!				X_water_in(j) = 0.3431*FS(j)**1.1*EXP(0.9736*RH(j)/FS(j)**0.95) !Mass fraction of water
!			else
!				X_water_in(j) = FS(j)**0.151*RH(j)**1.7872
!			end if 

! water mass fraction equation 
			X_water_in(j) = 0.30073566 - 2.6161978*FS(j) + 0.7194458*RH(j) + 5.26241655*FS(j)**2 &
				- 0.36628446*RH(j)**2 -2.679446*FS(j)**3 + 0.76097502*RH(j)**3 &
				+ 0.404578295*FS(j)*RH(j) - 0.17713611*FS(j)**2*RH(j) - 0.727581490539455*FS(j)*RH(j)**2
						   
		else	!RH <= ERH, no inorganic phase reactions

			X_water_in(j) = 0		 		! otherwise, keep the solid phase, kac = 0
			OS_Convfct = 0			
			in_phase_status = 1	
			
        end if		
		
		! Accumulated OS mass [ug/m**3]
		OS(j) = OS(j-1) + SO4_OS(j)*96		
		
		! Mass of Salt (SO4, NH4) Unit: ug/m3
		M_salt(j) = SO4(j)*96 + NH4(j)*18

		! Mass of water, Unit: ug/m3 (Known: inorganic mass (M_salt) and X_water_in)				   
		if (X_water_in(j) .le. 0) then
			M_water(j) = 0.00		 		! constrain water mass not less than 0		
		else if (X_water_in(j) .ge. 1) then
			X_water_in(j) = 1 		! constrain water mass not greater than 1		
		end if		
		
		M_water(j) = X_water_in(j) / (1 - X_water_in(j)) * M_salt(j)
		
		! Mass of inorganic aerosol, Unit: ug/m3
		M_in(j) = M_water(j) + M_salt(j)
		
		! Molecular weight of inorganic phase(g/mol)
		MW_in = M_in(j) /(M_water(j)/18 + SO4(j) + NH4(j))
				
		! *******  density of inorganic aerosol [g/cm3] *********
		Rho_in = -(0.494-(-0.937*(1-FS(j))**2 + 0.1534*(1-FS(j))))*RH(j) + 1.5757+0.7907*(1-FS(j))**2 - 0.0161*(1-FS(j))
		! regression from the E-AIM thermodynamic model
		! Vin = MW_in / Rho_in      ! Molar volume of inorganic phase (cm3/mol)
		! Volume of inorganic aerosol 
		! V_in = M_in/Rho_in
		
    end if

!NOTE: to include nitrate in future
!AC = EXP(-8.0437*RH-2.9984*ln(O:C)+0.0368*MW-0.9704*HB+12.1897*fanion-1.7258*fN)
!Where fanion = ([SO4]+[NO3])/([SO4]+[NO3]+[NH4]), mol/mol
!                fN = [NO3]/([NO3]+[SO4]), mol/mol
!                RH = 0-1

! *** calculate the organics to sulfate ratio in inorganic-phase for each lumping group
	Cin_sum(j) = sum(Cin(1:NCVAP))                           ! sum-up concentration of organic in inorganic phase
	
	if (SO4(j) .le. 0.00 ) then 			 ! restrict the value to prevent error
		org1sulf = 0.00 
	else
		org1sulf = Cin_sum(j)/(SO4(j)*96)    
	end if
		
! sensitivity test of the activity coefficent of organics in inorganic phase	
!	act_or_in (:) = act_or_in (:) / 2

!*********************************************************************
	call HET (j, DT, NCVAP, Max_array, N_cyc, RH, ERH, DRH, in_phase_status,MW, MW_in, Rho_in, beta1, beta2, kac, Kor, Kin, Vin, &
     		Cin, c0, FS, Rho_om, Ctb, Ctf, OMHi, OMT, OMo, M_in, C_proton, dOMH, dN_OS, TEST, MW_or, act_or_in, &
               O1C_SOA, SO4_OS, MW_omp, aging_Scale_current, O1C, HBonding, O1C_OMH, C_number, precursor)

	do i = 1, NCVAP
		OMH(j) = OMH(j) + OMHi(j,i)	! New Total OMH, initial OMH = 0. OMHi(j,i) is from SOA_H module
	end do

!##############################################################################	
!     Partitioning SOA (OMP)
!##############################################################################	

! *** Begin OMP solution code ***
!**************************************************************************
!                   Threshold to operate OMP
!  When either enough amount of POA is present or when the gas is saturated
!  or when inorganic phase exist, OMH begins 
!**************************************************************************
    threshold = 0.0   ! initialization 
    do i = 1, NCVAP
		threshold  = threshold  +  Ctb(j,i) / c0(i)   !test whether the gas phase is saturated
	end do
	
! *** initialize totorgsw and totorg (mole conc. of POA and pre-existing SOA) ***
    totorgsw = OMo/MW_omp(j-1) + OMH(j)/MW_omh(j-1) ! molar concentration of POA
! *** totorgsw = OMo + OMH(j)/MW_omh ***
	totorg = totorgsw + OMP(j-1)/MW_omp(j-1)   ! molar concentration of (SOA + POA) before time step

! *** check if gas/particle equilibrium can be established ***
    if ((threshold .gt. 1.0) .or. (totorgsw .gt. threshmin)) then
	
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! *** Initial guess of caer is computed as follows:
!     From eqn (8) of Schell et al. (2001)
!       caer = ctotf - c0 * (caer/MW) / totorg
!                                       ! totorg: SOA + POA before time step [ u-mole / m**3 ]            
!     Assuming totorg doesn't change during the timestep,
!       caer * (1 + c0/MW / totorg) = ctotf     (ctotf is the conc. of the total SVOC (gas+paticle))
!       caer = ctotf / ( 1 + c0/MW / totorg)
!
!     REVISISON ALN 1/07
!	Assuming no change in totorg during a time step and OMH is included as
!             caer = ctotf - OLDPOAi - c0*(caer/MW) / totorg    
! where OLDPOA is included in totorg and OLDPOAi is the OMH attributable to product i
!             caer = (ctotf - oldpoa) / (1 + c0/MW/totorg)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      	
    ctotf = Ctf - Cin					 ! Ctf - Cin [ug/m**3] (SVOC only for gas/org partitioning)

	caer = ctotf/(1+ c0*MW_m1/totorg)    ! initial guess for the Newtonian model 
    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c *** calculate new SOA by partitioning. This method uses a globally *** c
!c     convergent Newton-Raphson method coded by Dr Benedikt Schell		  c
!c     to solve the nonlinear quadratic system shown in eqn 8 of 		  c
!c     Schell et al:													  c
!c        A(i)  * caer(i) ** 2 + B(i) * caer(i) + C(i) = 0.0,			  c
!c        where B(i) contains the sum of all caer(j),					  c
!c       for j not equal to i.											  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   print *, "caer before", caer(:)
        
         call NEWT( LAYER, caer, NCVAP, check, &
                 ctotf, c0, MW_m1, totorgsw ) 
         
         if( check ) then
!  Try again with initial guess of 50/50 gas/aerosol split.  If NEWT solution has problem, change ratio from 0.5 to others. 
            do i = 1, NCVAP
             caer(i) = 0.1 * ctotf(i)  
            end do
            
            call NEWT( LAYER, caer, NCVAP, check,  &
                 ctotf, c0, MW_m1, totorgsw ) 
					    
            if (check) then                           
                write(3,*) ' problem in NEWT at LAYER = ', LAYER
            end if            
        end if
	
! *** Constrain caer to values between CONMIN and ctotf
        do i = 1, NCVAP
            if (caer(i) .lt. CONMIN ) then
                If(caer(i) .lt. 0.0 ) then
                    write(logdev,*) ' caer negative '
                end if
                caer(i) = CONMIN
            end if

           if( ( caer(i) - ctotf(i) ) .gt. TOLMIN ) then
!            if( ( caer(i) - ctotf(i) ) .gt. CTOLMIN ) then
!	                  

                write(LOGDEV,*) ' i = ', i
                write(LOGDEV,*) '  caer(i) = ', caer(i)
                write(LOGDEV,*) ' ctotf(i) = ', ctotf(i)
!		        
                caer(i) = ctotf(i)	 
                write(3,*) '  caer(i) = ', caer(i)
                write(3,*) ' ctotf(i) = ', ctotf(i)
		 
		    end if
        end do

    else
! *** threshold not exceeded; no material transferred to aerosol phase
        do i = 1, NCVAP
           caer(i) = 0.0
        end do

    end if                     ! check on equilibrium threshold
    
!##########################################################################################
!   Calculate total SOA concentrations (OMT) at the end of time step
!##########################################################################################
! *** calculation of partitioning organic matter (OMP) ***
        do i= 1, NCVAP
		    OMPi(j,i) = caer(i) + Cin(i)
		    OMP(j) = OMP(j) + OMPi(j,i)
        end do

! *** calculation of average molecular weight of OMP ************	
		!** initialization **
		MW_omp = 0
		MW_omh = 0		
		
		do i= 1, NCVAP	
			frac_OMHi = OMhi(j,i)/OMh(j)
			MW_omp(j) =  MW_omp(j) + OMPi(j,i)/OMP(j) * MW(i)  !sum of mass fraction(i) * MW(i)
			MW_omh(j) =  MW_omh(j) +  OMhi(j,i)/OMH(j) * MW(i)*2 !assume dimerization ==> MW were doubled
			
		end do
			
		if (OMP(j) .eq. 0) then
			MW_omp(j) = 200
		end if
		
		if (OMH(j) .eq. 0) then
			MW_omh(j) = 300      
		end if
		
        OMT(j) = OMo + OMP(j) + OMH(j)

		MW_or(j) = ( (OMP(j) + OMo) * MW_omp(j) + OMH(j) * MW_omh(j) ) / OMT(j)
		
		if ((OMH(j) .eq. 0.0) .or. (OMP(j) .eq. 0.0) .or. (OMT(j) .eq. 0.0)) then 
			MW_or(j) = 250.0
		end if 
		
! *** calculation of SOA O:C ratio ***   		 
        do i= 1, NCVAP
			prod_mole_conc(i) = OMpi(j,i)/MW(i) + OMHi(j,i)/(MW(i)*2)
            OMH_mole_conc(i) = OMHi(j,i)/(MW(i)*2)
			O1Ci(i) = O1C(i) * prod_mole_conc(i)
            O1CHi(i) = O1C(i) * OMHi(j,i)/(MW(i)*2)
        end do	
		
			O1C_SOA(j) = sum(O1Ci(1:NCVAP))/ sum(prod_mole_conc(1:NCVAP))
            O1C_OMH(j) = sum(O1CHi(1:NCVAP))/ sum(OMH_mole_conc(1:NCVAP) )
        
        if (O1C_SOA(j) .le. 0.3) then
            O1C_SOA(j) = 0.3 !0.45 for TMBs  0.35 for o-xyl
        end if
        
        if (O1C_OMH(j) .le. 0.35) then
          O1C_OMH(j) = 0.35 
        end if
		
! *** Section to save the dynamic MW and O:C for printing out ***		
	MW_save(j,:) = MW(:)
	O1C_save(j,:) = O1C(:)
	HBonding_save(j,:) = HBonding(:)
	act_or_in_save(j,:) = act_or_in(:)
	Kin_save(j,:) = Kin(:)


End do

!##############################################################################	
!#             Writing simulation results into the file (summary) 			  #
!##############################################################################	

! *** Formating the output file *** 
	100 Format(T1,'Time',T7,'OMT',T15,'OMAR',T25,'OMP',T35,'GLY',T43,'M_in',T50,'FS',T56,'ERH',T63,'Proton',&
        T71,'M_water',T79,'Cin_sum',T90,'OS'/)
	102 Format(T1,'=====',T7,'=======',T15,'=========',T25,'=========',T35,'=======',T43,'======',T50,'=====',&
        T56,'======',T63,'=======',T71,'=======',T79,'=========',T90,'======='/)
    105 Format(T1,I5,T7,F7.3,T15,ES9.2,T25,ES9.2,T35,F7.3,T43,F6.1,T50,F5.2,T56,F6.3,T63,F7.3,T71,F7.2,T79,ES9.2,T90,F7.4)

! For sensitivity test
!    100 Format(T1,'OMT',T9,'OMAR',T17,'OMP'/)
!	 102 Format(T1,'=======',T9,'=======',T17,'======='/)
!    105 Format(T1,F7.3,T9,F7.3,T17,F7.3)
! End

	200 Format(T1,8ES12.3)
	201	Format(T1,8F9.5)				!For the input of alpha
	202	Format(T1,8F9.3)				!For the input of MW and O:C
	205 Format(T1,8ES12.3)	

! *** Main SOA output ***
	write(4,*) "Simulation results"
	write(4,100,advance='no') 			! Titles
	write(4,102,advance='no')			! Deliminations
	do j = 0, N_cyc
		write(4,105) (j*DT),OMT(j),OMH(j),OMP(j),OMHi(j,51),M_in(j),FS(j),ERH(j),C_proton(j),M_water(j),Cin_sum(j),OS(j)
!        write(4,105) OMT(j),OMH(j),OMP(j)   !For sensitivity test
	end do 

! *** Output for test ***
	write(4,*) "free SO4 (ug/m3)"
	j=j-1
	write(4,*) j, SO4_free(j)*96
 
! *** out put structure ***
	! Lumping Structure
	!  |  i1  i2  i3  i4  i5  i6  i7  i8  
	!-----------------------------------	
	!VF|   1   2   3   4   5   6   7   8
	! F|   9  10  11  12  13  14  15  16
	! M|  17  18  19  20  21  22  23  24
	! S|  25  26  27  28  29  30  31  32
	! P|  33  34  35  36  37  38  39  40  
	!MA|  41  42  43  44  45  46  47  48
    ! R|  49  50  51  
    !  |EPOX MGLY GLY

	
		j = 50			!7:00 EST for TMB, since the relatively fast decay
		k = 110			!6 hours from the beginning
	
		!Fresh condition
		write(4,*) 
		write(4,*) (j*DT/60)," hours from the beginning"
		write(4,*) "OMAR"
		write(4,200) OMHi(j,:)   
		write(4,*) "OMp"
		write(4,200) OMpi(j,:)
		write(4,*) "Ctb" 	! SVOC conc. before the current time step [ug/m**3]
		write(4,200) Ctb(j,:)
		write(4,*) "alpha_dyn"
		write(4,201) alpha_dyn(j,:)
		write(4,*) "MW"
		write(4,202) MW_save(j,:)
		write(4,*) "O1C"
		write(4,202) O1C_save(j,:)	
		write(4,*) "activity_coeff_organics_in_inorganics"
		write(4,205) act_or_in_save(j,:)
		write(4,*) "Kin"
		write(4,205) Kin_save(j,:)
		write(4,*) "kac"
		write(4,205) kac(j,:)	
		write(4,*) "Hydrogen bonding term"
		write(4,202) HBonding_save(j,:)
		write(4,*) "SOA O:C ratio = ", O1C_SOA(j)
        write(4,*) "OMH O:C ratio = ", O1C_OMH(j)
        write(4,*) "aging_scale", aging_Scale(j)
        write(4,*) "aging_scale_current", aging_Scale_current(j)
        write(4,*) "MW_SOA = ", MW_or(j)
		
		!Highly aged condition
		write(4,*)		
		write(4,*) (k*DT/60)," hours from the beginning"
		write(4,*) "OMAR"
		write(4,200) OMHi(k,:)   
		write(4,*) "OMp"
		write(4,200) OMpi(k,:)
		write(4,*) "Ctb" 	! SVOC conc. before the current time step [ug/m**3]
		write(4,200) Ctb(k,:)
		write(4,*) "alpha_dyn"
		write(4,201) alpha_dyn(k,:)
		write(4,*) "MW"
		write(4,202) MW_save(k,:)
		write(4,*) "O1C"
		write(4,202) O1C_save(k,:)
		write(4,*) "activity_coeff_organics_in_inorganics"
		write(4,205) act_or_in_save(k,:)
		write(4,*) "Kin"
		write(4,205) Kin_save(k,:)
		write(4,*) "kac"
		write(4,205) kac(k,:)
		write(4,*) "Hydrogen bonding term"
		write(4,202) HBonding_save(k,:)
		write(4,*) "SOA O:C ratio = ", O1C_SOA(k)
        write(4,*) "OMH O:C ratio = ", O1C_OMH(k)
        write(4,*) "aging_scale", aging_Scale(k)
        write(4,*) "aging_scale_current", aging_Scale_current(k)
        write(4,*) "MW_SOA = ", MW_or(k)

 CLOSE (4)
 CLOSE (3)
end program UNIPAR_V5
