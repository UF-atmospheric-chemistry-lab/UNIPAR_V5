!############################################################################################
!  PROGRAM: UNIfied Partitioning and Aerosol Phase Reaction (UNIPAR) OMH module SOA_H5
!
!  PURPOSE: To simulate SOA mass based on multiphase reactions of organic species 
!  Submodule of UNIPAR_V5 for aerosol phase reactions
!
!  Contact Information: 
!	Dr. Myoseon Jang, University of Florida at mjang@ufl.edu
!
!	This UNIPAR code is to simulate alkane SOA. If paramters for other SOA precursors are needed,
!	please contact Dr. Jang.
!  
!	Heterogeneous reactions are based on a flow reactor data collected using various model compounds at different humidity
!	and seed composition (acidity). 
!
!	Jang et al., "Semiempirical Model for Organic Aerosol Growth by Acid-Catalyzed Heterogeneous Reactions of Organic Carbonyls" 
!	Env. Sci. Tech., doi: 10.1021/es048977h
!
!	Rate constants for heterogeneous reactions were updated due to solubility calculation. Updated solubility parameters were
!	applied to flow reactor data from Jang et al. (2005) and used to regenerate rate constants. Updated rate constants were 
!	validated with data from various SOA precursor hydrocarbons.
!  
!  PHRCSOA is the previous version of UNIPAR, coded by Yunseok Im
!  		Im et al., "Simulation of aromatic SOA formation using the lumping model 
!  		integrated with explicit gas-phase kinetic mechanisms and aerosol-phase reactions" ACP, doi:10.5194/acp-14-4013-2014
! 
!  Extended to 51 lumping on August 1, 2016 (Isoprene SOA)
!  	 	Beardsley and Jang, "Simulating the SOA formation of isoprene from partitioning and aerosol phase reactions 
!    	in the presence of inorganics” DOI:10.5194/acp-16-5993-2016, ACP, 2016
!  Updated to dynamic features (alpha, MW, O:C, and H-bonding )
!		Zhou et al., “Simulation of SOA Formation from the Photooxidation of Monoalkylbenzenes 
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

 Module SOA_H5

	implicit none  
	
	contains

	SUBROUTINE HET( &
                   j, DT, NCVAP, Max_array, N_cyc, RH, ERH, DRH, in_phase_status, MW, MW_in, Rho_in, &
                    beta1, beta2, kac, Kor, Kin, Vin, Cin, c0, Fs, Rho_om, Ctb, Ctf, OMHi, OMT, OMo, &
                    M_in, C_proton, dOMH, dN_OS, TEST, MW_or, act_or_in, &
                    O1C_SOA, SO4_OS, MW_omp, aging_scale_current, O1C, HBonding, O1C_OMH, C_number, precursor)

	logical seed
          parameter(seed = .FALSE.)
								  !with inorganic seed ->TRUE, SA from SO2 ->.FALSE.
	Integer j, NCVAP, i, Max_array, DT, N_cyc
	REAL test(0:Max_array, NCVAP)
	Character (len = 4) precursor
    REAL RH(0:Max_array) 		  ! [0 ~ 1]
	REAL ERH(0:Max_array)	      ! Efflorscence relative humidity [0 ~ 1]
	REAL DRH(0:Max_array)	      ! Deliquescence relative humidity [0 ~ 1]  
	Integer in_phase_status       ! inorganic phases: 0-->liquid phase; 1-->solid phase (default status: 0)
    REAL O1C(NCVAP)	 			  ! dynamic O:C ratio of each lumped group 
    REAL HBonding(NCVAP)
    REAL MW(NCVAP), Cin(NCVAP), c0(NCVAP), Csum(NCVAP)
	REAL Ctb(0:Max_array, NCVAP), Ctf(NCVAP)
	REAL OMT(0:Max_array)
		
	REAL OMHi(0:Max_array, NCVAP), dOMHi(NCVAP), dOMH(0:Max_array), N_OS(0:Max_array)

	REAL MW_or(0:Max_array)     ! Average Molecular weight of organic aerosol (OMP+OMH)
	REAL Rho_om, Rho_in		      ! density of inorganics(M_in)
	REAL Fs(0:Max_array)	      ! inorganic composition factor 
	REAL M_in(0:Max_array)	      ! Mass of inorganic aerosol [ug/m3]
	REAL Cor(NCVAP)		          ! SVOC conc. in organic aerosol phase [mol/L of OM]	
	REAL koo(NCVAP)		          ! oligomerization rate constants of oxidized products in organic phase [L/mol/s]
	REAL kac_o(NCVAP)	      	  ! Intercept of acid-catalyzed reaction rate constants of oxidized products in inorganic phase [L/mol/s]
	REAL pkBH(NCVAP)	          ! inverse log of KBH+(carbonyl protonation eq. constant)
	REAL II(NCVAP)                ! reactivity based on the structure lumping constant  !
	REAL X				          ! excess acidity based on RH
    REAL XACID			          ! excess acidity based on RH and H2SO4
	REAL act_W_in		          ! activity of water in inorganic aerosol
	REAL C_W_in                   ! concentration of water in inorganic aerosol
	REAL C_proton(0:Max_array)    ! concentration of proton in inorganic aerosol (constant Mseed version)
	REAL O1C_SOA(0:Max_array)	  ! O:C ratio of aerosol ( sum up(O1C(i))/MW of SOA )
    REAL O1C_OMH(0:Max_array)	  ! O:C ratio of OMH ( sum up(O1CH(i))/MW of OMH )
    REAL O1C_SOA_vis(0:Max_array) ! to control viscosity
    REAL kH(NCVAP)		          ! Heterogeneous rate constant for each lumped product group
	REAL kac(0:Max_array, NCVAP)  ! kH*M_in/1000+ko
	REAL term_inor				  ! log of water activity x proton conc for inorganic aerosol
	REAL act_or_in(NCVAP)		  ! activity coefficients of lumping species in inorganic aerosol (diluted condition/Under saturation)
	REAL act_or_in_sat(NCVAP)	  ! activity coefficients of lumping species in inorganic aerosol (when saturated)
	REAL act_or_or(NCVAP)		  ! activity coefficients of lumping species in organic aerosol
	REAL Kor(0:Max_array, NCVAP)  ! partioning coeff. between organic aerosol and gas
	REAL Kin(NCVAP)		          ! partioning coeff. between inorganic aerosol and gas
	REAL Koi(NCVAP)		          ! partioning coeff. between orgnic and inorganic aerosol
	REAL Vor			          ! morlar volume of organic medium (average)
		parameter(Vor = 95)
	REAL Vin			          ! morlar volume or inorganic medium.
	REAL t_crt(NCVAP)	          !(critical time) the time which .
	REAL beta1(0:Max_array,NCVAP) ! coefficient1 for OMH calculation (related with Oligomerizations (Organic phase))
	REAL beta2(0:Max_array,NCVAP) ! coefficient2 for OMH calculation (related with Acid-catalyzed reactions (inorganic phase))
	REAL f_in (NCVAP)			  ! beta2 / (beta1 + beta2) 
								  ! fraction of OM_AC to the Total OMH (OM_AC + OM_olg)
	REAL beta3(NCVAP)	          ! coefficient3 for OMH calculation
    REAL beta4(NCVAP)	          ! coefficient4 for OMH calculation
	REAL in_mol2ug(NCVAP)         ! [ ug/m3 of air per umol/L of inorganic ] for Cin
	REAL or_mol2ug(NCVAP)         ! [ ug/m3 of air per umol/L of inorganic ] for Cor
	REAL MW_in                    ! moleculr weight of inorganic medium [ g/mol ]   
	REAL OMo					  ! Preexisting organic particular matter [ug/(m**3)] 
    Real C_number
!New addtition
    REAL ko(NCVAP)				  ! oligomerization rate constants of oxygenated products in organic phase [L/mol/s]
    REAL dN_OS
    REAL SO4_OS(0:Max_array)      ! Consumed Sulfate concentration by OS formation	[umol/m3 of air]
    REAL MW_omp(0:Max_array)                   ! Average Molecular weight of OMH (oligomer)
    REAL aging_Scale_current(0:Max_array) ! [(current - fresh)/(aged - fresh)], scale: 0~1
    
!******************************************************************
!                       Lumping Matrix 							  *
!******************************************************************

	! Lumping Structure
	!  | i1	i2	i3	i4	i5	i6  i7  i8  
	!----------------------------------	
	!VF|  1	 2	 3	 4	 5	 6   7   8
	! F|  9	10	11	12  13  14  15  16
	! M| 17	18  19  20  21  22  23  24
	! S| 25	26	27	28	29	30  31  32
	! P| 33 34  35  36  37  38  39  40  
	!MA| 41 42  43  44  45  46  47  48
    ! R|                    49  50  51
    !49 is for IEPOX; 50 is for Mglyoxal; and 51 is for glyoxal

!################################################################
!           Rate constants for aerosol phase reactions
!#################################################################   
!*************** Basicity of organic species *********

       pkBH	= (/-1.5,	-1.5,	-1.5,	-1.5,	-1.5,	-1.5, -1.5, -1.5, &
                -1.5,	-1.5,	-1.5,	-1.5,	-1.5,   -1.5, -1.5, -1.5, &
				-1.5,	-1.5,	-1.5,	-1.5,	-1.5,	-1.5, -1.5, -1.5, &
				-2.0,	-2.0,	-2.0,	-2.0,	-2.0,	-2.0, -2.0, -2.0, &
				 0.0,	 0.0,	 0.0,	 0.0,	 0.0,	 0.0,  0.0,  0.0, &
				 0.0,	 0.0,	 0.0,	 0.0,	 0.0,	 0.0,  0.0,  0.0, &	
						                                -1.5,  -1.5, -1.5	/)  
 
! Reactivity scale of organic species in aerosol phase. IEPOX is for isoprene.  MA is high in isoprene (middle reactivity)	   

    II = (/	8.5,	8.5,	8.5,	8.5,	8.5,	8.5,  8.5,  8.5, &           	
            8.2,	8.2,	8.2,	8.2,	8.2,	8.2,  8.2,  8.2, &	
            4.0,	4.0,	4.0,	4.0,	4.0,	4.0,  4.0,  4.0, &
        	1.0,	1.0,	1.0,	1.0,	1.0,	1.0,  1.0,  1.0, &
        	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,  0.0,  0.0, &
        	0.0,    0.0,    0.0,    0.0,    0.0,    0.0,  0.0,  0.0, &
		                    		               8.2,  9.0,  12.5  /)
! IEPOX is treated as F 
    

!*************** activity coeffienct of organic compounds in the organic phase  *********
	act_or_or 	=  (/1,	1,	1,	1,	1,	1,	1,  1, &
                    1,	1,	1,	1,	1,	1,	1,  1, &
					1,	1,	1,	1,	1,	1,	1,  1, &
      				1,	1,	1,	1,	1,	1,	1,  1, &
      				1,	1,	1,	1,	1,	1,	1,  1, &
      				1,	1,	1,	1,	1,	1,	1,  1, &
                                        1,	1,  1  /)
!*************** activity coeffienct of organic compounds in the inorganic phase  *********
	do i = 1, NCVAP
		
        if (O1C(i) .ge. 0.00 .AND. O1C(i) .le. 0.8) then  
            act_or_in(i)=exp(-1.342*HBonding(i)+0.0471*(MW(i))-3.435*log(O1C(i))-1.440*(FS(j))-0.0306*(100*RH(j))) ! activity regression 10-10000
         
        else if (O1C(i) .gt. 0.8) then     
           act_or_in(i)=exp(-1.1209*HBonding(i)+0.0354*(MW(i))-2.7036*log(O1C(i))-0.3295*(FS(j))-0.0216*(100*RH(j))) ! Zhou et al., (2019), ACP
        end if
        
        If (act_or_in(i) .ge. 1E+8) then
            act_or_in(i) = 1E+8    
        else if(act_or_in(i) .lt. 1)then
			act_or_in(i) = 1
		end if
	end do

!******** Proton concnetration in the inorganic phase ******
    If (FS(j) .gt. 0.55) then	
        C_proton(j) = FS(j)*(-16.051*RH(j)**3 + 21.096*RH(j)**2 - 16.782*RH(j) + 12.399)   !original
        !regression from E-AIM II (at 298.15K) with the correction at higher acidic conditions
    else
        C_proton(j) = 15**(-1/FS(j)) * (-759.264*RH(j)**3 + 105.556*RH(j)**2 - 18.524*RH(j) + 743.116)
        !regression from E-AIM II (at 298.15K) with correction at lower acidic conditions
    End if

!*********** Calculating Excess Acidity (X) ***********************************
	X=(-4*RH(j)**3 + 7*RH(j)**2 -7.41*RH(j) +4.2343)
	Xacid= (0.90689+Fs(j)*8.94391+RH(j)*(-5.98593))*0.322*(1-Fs(j))+X*Fs(j)		
	
!*** Reference for X function  
! 1. Cox, R.A.; Yates, K., Kinetic equations for reactions in concentrated aqueous acids based on the concept of "excess acidity". 
!    Canadian Journal of Chemistry 1979, 57, 2944-2951
! 2. Rochester C. H. Acidity Function, Academic Press, New York, 1970
! 3. Li, J. and Jang, M., Aerosol Acidity Measurement Using Colorimetry Coupled with a Reflectance UV-Visible Spectrometer 
!    Aerosol Sci. Technol., in print, 2012

!*** reference for the sulfuric acid composition vs humidity 
! 1. Liu, B.Y.H.; Levi, J. Generation of aerosols and facilities for exposure experiments 
!    Ann Arbor Science Publishers Inc., Ann Arbor, 1980, pp. 317-336
! 2. Perry, R.H. Perry's Chemical Engineering Handbook, 6th, Mcgraw-Hill Cor.  New York, 1984, pp. 364-366

!***********	Calculating the rate constant (kH) of organic compounds in inorganic phase *********************
	term_inor=LOG10(RH(j)*C_proton(j))      !log of water activity x proton conc for inorganic aerosol 
    O1C_SOA_vis = O1C_SOA - 0.5   
     if (O1C_SOA_vis(j-1) .le. 0) then 
		O1C_SOA_vis = 0
     end if
          
        kac(j,:)=10**(0.25*pkbH + 1.0*Xacid + 0.95*II + term_inor  - 2.58)      ! LLPS-AIOMFAC 1-1000 solubility
        ko(:)= 10**(0.25*pkbH + 0.95*II - 2.58 - 3.8 + 1.2*(1-1/(1+exp((300-MW_or(j-1))*0.005))) &
        + 2.2/(1+exp((0.75-O1C_SOA(j-1))*6)))  ! considering of viscosity of organic species by usibng MW_or and O:C ratios.


!*** for sensitivity test ***  
!	kac(j,:) = kac(j,:) / 2
	!for sensitivity test 
	!ko = ko / 150
    
!#####################################################################
!                          OMH calculation 
!#####################################################################   
			
!**** Partitioning coefficient between organic and inorganic aerosol phases

	! ******* Calculation of Kor and Kin***********************  		
	if (MW_or(j-1) .le. 0.0) then 
		MW_or(j) = 250
	end if
	
	Kor(j,:)= 1/c0 * MW/MW_or(j-1) / act_or_or
	Kin	= 1/c0 * MW/MW_in / act_or_in    !**** partioning coeff of organics in the inorganic phase**** 

	if (RH(j) .gt. ERH(j).AND. in_phase_status .eq. 0) then
									! inorganic has liquid phase
		kac = kac
	else					
		kac(j,:) = 0.000		! Under ERH, no acid-catalyzed reactions
		in_phase_status = 1
		Kin = 0.000  			! Under ERH, Organic compounds are not soluble in the inorganic phase								
    end if
		

! *** Concenrations of organic compounds in inorganic and organic phases (Cin and Cor) ***
    OMT(0) = OMo   
	in_mol2ug = MW*M_in(j)/Rho_in/1000  !Conversion factor for Cin from [mol/L of M_in] to .[ug/m3 of air]
	or_mol2ug = MW*OMT(j-1)/Rho_om/1000 !Conversion factor for Cor from [mol/L of M_in] to .[ug/m3 of air]
	Cor = Kor(j,:)*OMT(j-1)/(1+Kor(j,:)*OMT(j-1)+Kin*M_in(j))*Ctb(j,:) ! before oligomerization
	Cin = Kin*M_in(j)/(1+Kor(j,:)*OMT(j-1)+Kin*M_in(j))*Ctb(j,:)       ! before oligomerization

! *** Constant terms for the analytical solution of OMH ******    
	beta1(j,:)=ko*1000*Rho_om*Kor(j,:)**2*OMT(j-1)/(MW*(1+Kor(j,:)*OMT(j-1)+Kin*M_in(j))**2) ! if OMT = 0, beta1 = 0
    !beta1(j,:)=koo*1000*Rho_om*Kor(j,:)**2*OMT(j-1)/(MW*(1+Kor(j,:)*OMT(j-1)+Kin*M_in(j))**2) ! if OMT = 0, beta1 = 0
	beta2(j,:)=kac(j,:)*1000*Rho_in*Kin**2*M_in(j)/(MW*(1+Kor(j,:)*OMT(j-1)+Kin*M_in(j))**2) ! if M_in = 0, beta2 = 0
	
! for the compounds produced OMH	
	do i=1, 32  !Acid-Catalyzed Heterogeneous reaction(Inorganic Phase) + Oligomerization (Organic Phase)
				!Fast(1-8), Medium(9-16), Slow(17-24) reactivity lumping group				
		dOMHi(i) = (beta1(j,i)+beta2(j,i))*Ctb(j,i)**2*DT*60 &					
      				  /(1+(beta1(j,i)+beta2(j,i))*Ctb(j,i)*DT*60)
		Ctf(i) = Ctb(j,i) - dOMHi(i)  ! SVOC for calculating partiotning mass (to NEWT)          
		Cor(i) = Kor(j,i)*OMT(j-1)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i) ! after oligomerization  
		Cin(i) = Kin(i)*M_in(j)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i)    ! after oligomerization (just check)
    end do

! for the medium reactivity with multifunctinal alcohol (this is not important in aromatic but important in isoprene SOA     
!   do i=41, 48  !Acid-Catalyzed Heterogeneous reaction(Inorganic Phase) + Oligomerization (Organic Phase)
!				!Fast(1-8), Medium(9-16), Slow(17-24) reactivity lumping group	
!		dOMHi(i) = (beta1(j,i)+beta2(j,i))*Ctb(j,i)**2*DT*60 &
 !     				  /(1+(beta1(j,i)+beta2(j,i))*Ctb(j,i)*DT*60)
!		Ctf(i) = Ctb(j,i) - dOMHi(i)  ! SVOC for calculating partiotning mass (to NEWT)          
!		Cor(i) = Kor(j,i)*OMT(j-1)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i) ! after oligomerization  
!		Cin(i) = Kin(i)*M_in(j)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i)    ! after oligomerization (just check)
 !   end do

! for the compounds produced OMH (gly, mgly and epoxydiol)	
	do i=49, 51  !Acid-Catalyzed Heterogeneous reaction(Inorganic Phase) + Oligomerization (Organic Phase)
				!Fast(1-8), Medium(9-16), Slow(17-24) reactivity lumping group	
		dOMHi(i) = (beta1(j,i)+beta2(j,i))*Ctb(j,i)**2*DT*60 &
      				  /(1+(beta1(j,i)+beta2(j,i))*Ctb(j,i)*DT*60)
		Ctf(i) = Ctb(j,i) - dOMHi(i)  ! SVOC for calculating partiotning mass (to NEWT)          
		Cor(i) = Kor(j,i)*OMT(j-1)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i) ! after oligomerization  
		Cin(i) = Kin(i)*M_in(j)/(1+Kor(j,i)*OMT(j-1)+Kin(i)*M_in(j))*Ctf(i)    ! after oligomerization (just check)
    end do    
    
! No Heterogeneous Reaction: Partitioning only groups(33-40)
	do i=33,48       !multialcohol is only for partitioning 11/29/2021 by MJ and JC
		Ctf(i) = Ctb(j,i)
		dOMHi(i) = 0
    end do

! Assume Multi-Alchol(25-30) groups become OMH. (remove on July/12/2016)   
!	do i=25,30
!		Ctf(i) = 0
!		dOMHi(i) = Ctb(j,i)
!	end do
        
!**** Accumulated OMH (previous + new)  ******      
    OMHi(j,:) = OMHi(j-1,:) + dOMHi(:)   
    
    do i=1, NCVAP 
        dOMH(j) = dOMH(j) + dOMHi(i)
    end do

!#################################################################
!            Acidity reduction due to OS formaiton 
!#################################################################   
! Free sulfates can react with organics and form OS. 
! THe formation OS reduces the acidity of aerosol.  In turn, the reduced acidity affect aerosol phase reaction.
! The reactivity weighting factors for OS formation was applied to lumping groups
! for example
! Highly reactive group: 4
! mid-reactivity group: 2
! multifunctional alcohol: 3 and 4
! Weighting factor is applyed to OS_convfct in PHRCSOA (fos)
! The rate constants were adjusted together with OS formation (kac_o)

!******** Weighting factor of Functional group associated with OS formation ********

	!**** Fraction of OMH originating from the inorganic phase reaction 
	do i=1, NCVAP
	if (beta1(j,i) + beta2(j,i) .eq. 0) then
		f_in(i) = 0.0
	else
		f_in(i) = beta2(j,i)/(beta1(j,i)+beta2(j,i))
	!f_in = [ OM_AC / (OM_Olg + OM_AC) ]
	!	OM_Olg: OMH originating from the organic phase reaction (oligomerization)
	!	OM_AC : OMH originating from the inorganic phase reaction (Acid-catalyzed Rxn)
	end if
    end do
	
    If (precursor .eq. 'IS01') then
	        dN_OS =									&
            !Very Fast reactivity lumping group (defaul Factor 4)
				+ dOMHi(1)*f_in(1)/MW(1)*4	&
				+ dOMHi(2)*f_in(2)/MW(2)*4	&
				+ dOMHi(3)*f_in(3)/MW(3)*4	&
				+ dOMHi(4)*f_in(4)/MW(4)*4	&
				+ dOMHi(5)*f_in(5)/MW(5)*4	&
				+ dOMHi(6)*f_in(6)/MW(6)*3	&
                + dOMHi(7)*f_in(7)/MW(7)*2	&
                + dOMHi(8)*f_in(8)/MW(8)*2	&
			!Fast reactivity lumping group (defaul Factor 2)
                + dOMHi(9)*f_in(9)/MW(9)*3	&
				+ dOMHi(10)*f_in(10)/MW(10)*3	&
				+ dOMHi(11)*f_in(11)/MW(11)*3	&
				+ dOMHi(12)*f_in(12)/MW(12)*3	&
				+ dOMHi(13)*f_in(13)/MW(13)*3	&
				+ dOMHi(14)*f_in(14)/MW(14)*3	&
                + dOMHi(15)*f_in(15)/MW(15)*3	&
                + dOMHi(16)*f_in(16)/MW(16)*2	&
            !medium reactivity lumping group (default Factor 1)
				+ dOMHi(17)*f_in(17)/MW(17)*4   &
				+ dOMHi(18)*f_in(18)/MW(18)*2  	&
				+ dOMHi(19)*f_in(19)/MW(19)*3   &
				+ dOMHi(20)*f_in(20)/MW(20)*2	&
				+ dOMHi(21)*f_in(21)/MW(21)*2	&
				+ dOMHi(22)*f_in(22)/MW(22)*1	&
                + dOMHi(23)*f_in(23)/MW(23)*1	&
                + dOMHi(24)*f_in(24)/MW(24)*2	&   
			!multialcohol S lumping group (alcohols in isoprene products)
				+ dOMHi(25)*f_in(25)/MW(25)*1   &
				+ dOMHi(26)*f_in(26)/MW(26)*1  	&
				+ dOMHi(27)*f_in(27)/MW(27)*1   &
				+ dOMHi(28)*f_in(28)/MW(28)*1	&
				+ dOMHi(29)*f_in(29)/MW(29)*1	&
				+ dOMHi(30)*f_in(30)/MW(30)*1	&
                + dOMHi(31)*f_in(31)/MW(31)*1	&
                + dOMHi(32)*f_in(32)/MW(32)*1	&   
            !multialcohol S lumping group (alcohols in isoprene products)
				+ dOMHi(33)*f_in(33)/MW(33)*2   &
				+ dOMHi(34)*f_in(34)/MW(34)*2  	&
				+ dOMHi(35)*f_in(35)/MW(35)*1   &
				+ dOMHi(36)*f_in(36)/MW(36)*2	&
				+ dOMHi(37)*f_in(37)/MW(37)*1	&
				+ dOMHi(38)*f_in(38)/MW(38)*1	&
                + dOMHi(39)*f_in(39)/MW(39)*1	&
                + dOMHi(40)*f_in(40)/MW(40)*1	&      
            !multialcohol  lumping group (default Factor 3)
				+ dOMHi(41)*f_in(41)/MW(41)*3   &
				+ dOMHi(42)*f_in(42)/MW(42)*3  	&
				+ dOMHi(43)*f_in(43)/MW(43)*3   &
				+ dOMHi(44)*f_in(44)/MW(44)*3	&
				+ dOMHi(45)*f_in(45)/MW(45)*3	&
				+ dOMHi(46)*f_in(46)/MW(46)*3	&
                + dOMHi(47)*f_in(47)/MW(47)*3	&
                + dOMHi(48)*f_in(48)/MW(48)*3	&   
            !Special Reactive lumping group (default factor 4) 
				+ dOMHi(49)*f_in(49)/MW(49)*4	&    
                + dOMHi(50)*f_in(50)/MW(50)*1	&
				+ dOMHi(51)*f_in(51)/MW(51)*2	
                   
    else if ((precursor .eq. 'TA01') .or. (precursor .eq. 'TA02') .or. (precursor .eq. 'TA03')) then
                dN_OS =                         &
               !Very Fast reactivity lumping group by using a-pinene 
                + dOMHi(1)*f_in(1)/MW(1)*3	&
				+ dOMHi(2)*f_in(2)/MW(2)*2	&
				+ dOMHi(3)*f_in(3)/MW(3)*2	&
				+ dOMHi(4)*f_in(4)/MW(4)*2	&
				+ dOMHi(5)*f_in(5)/MW(5)*2	&
				+ dOMHi(6)*f_in(6)/MW(6)*2	&
                + dOMHi(7)*f_in(7)/MW(7)*2	&
                + dOMHi(8)*f_in(8)/MW(8)*2	&
			!Fast reactivity lumping group (defaul Factor 2)
                + dOMHi(9)*f_in(9)/MW(9)*2	&
				+ dOMHi(10)*f_in(10)/MW(10)*2	&
				+ dOMHi(11)*f_in(11)/MW(11)*2	&
				+ dOMHi(12)*f_in(12)/MW(12)*3	&
				+ dOMHi(13)*f_in(13)/MW(13)*2	&
				+ dOMHi(14)*f_in(14)/MW(14)*2	&
                + dOMHi(15)*f_in(15)/MW(15)*2	&
                + dOMHi(16)*f_in(16)/MW(16)*2	&
            !medium reactivity lumping group (default Factor 1)
				+ dOMHi(17)*f_in(17)/MW(17)*1   &
				+ dOMHi(18)*f_in(18)/MW(18)*2  	&
				+ dOMHi(19)*f_in(19)/MW(19)*1   &
				+ dOMHi(20)*f_in(20)/MW(20)*1	&
				+ dOMHi(21)*f_in(21)/MW(21)*1	&
				+ dOMHi(22)*f_in(22)/MW(22)*1	&
                + dOMHi(23)*f_in(23)/MW(23)*1	&
                + dOMHi(24)*f_in(24)/MW(24)*1	&   
			!multialcohol  lumping group (default Factor 1)
				+ dOMHi(41)*f_in(41)/MW(41)*3   &
				+ dOMHi(42)*f_in(42)/MW(42)*3  	&
				+ dOMHi(43)*f_in(43)/MW(43)*3   &
				+ dOMHi(44)*f_in(44)/MW(44)*3	&
				+ dOMHi(45)*f_in(45)/MW(45)*3	&
				+ dOMHi(46)*f_in(46)/MW(46)*3	&
                + dOMHi(47)*f_in(47)/MW(47)*3	&
                + dOMHi(48)*f_in(48)/MW(48)*3	&   
            !Special Reactive lumping group (default factor 4) 
                + dOMHi(49)*f_in(49)/MW(49)*4	&    
                + dOMHi(50)*f_in(50)/MW(50)*1	&
				+ dOMHi(51)*f_in(51)/MW(51)*2
			else if ((precursor .eq. 'BB01') .or. (precursor .eq. 'BB02') .or. (precursor .eq. 'BA01') .or. (precursor .eq. 'BA02') &
			.or. (precursor .eq. 'AA01')) then
                dN_OS =                         &
               !Very Fast reactivity lumping group by using a-pinene 
                + dOMHi(1)*f_in(1)/MW(1)*3	&
				+ dOMHi(2)*f_in(2)/MW(2)*3	&
				+ dOMHi(3)*f_in(3)/MW(3)*2	&
				+ dOMHi(4)*f_in(4)/MW(4)*4	&
				+ dOMHi(5)*f_in(5)/MW(5)*4	&
				+ dOMHi(6)*f_in(6)/MW(6)*4	&
                + dOMHi(7)*f_in(7)/MW(7)*3	&
                + dOMHi(8)*f_in(8)/MW(8)*2	&
			!Fast reactivity lumping group (defaul Factor 2)
                + dOMHi(9)*f_in(9)/MW(9)*2	&
				+ dOMHi(10)*f_in(10)/MW(10)*2	&
				+ dOMHi(11)*f_in(11)/MW(11)*2	&
				+ dOMHi(12)*f_in(12)/MW(12)*2	&
				+ dOMHi(13)*f_in(13)/MW(13)*4	&
				+ dOMHi(14)*f_in(14)/MW(14)*3	&
                + dOMHi(15)*f_in(15)/MW(15)*2	&
                + dOMHi(16)*f_in(16)/MW(16)*2	&
            !medium reactivity lumping group (default Factor 1)
				+ dOMHi(17)*f_in(17)/MW(17)*3   &
				+ dOMHi(18)*f_in(18)/MW(18)*2  	&
				+ dOMHi(19)*f_in(19)/MW(19)*3   &
				+ dOMHi(20)*f_in(20)/MW(20)*1	&
				+ dOMHi(21)*f_in(21)/MW(21)*1	&
				+ dOMHi(22)*f_in(22)/MW(22)*1	&
                + dOMHi(23)*f_in(23)/MW(23)*1	&
                + dOMHi(24)*f_in(24)/MW(24)*2	&   
			!multialcohol  lumping group (default Factor 1)
				+ dOMHi(41)*f_in(41)/MW(41)*5   &
				+ dOMHi(42)*f_in(42)/MW(42)*4  	&
				+ dOMHi(43)*f_in(43)/MW(43)*3   &
				+ dOMHi(44)*f_in(44)/MW(44)*3	&
				+ dOMHi(45)*f_in(45)/MW(45)*3	&
				+ dOMHi(46)*f_in(46)/MW(46)*3	&
                + dOMHi(47)*f_in(47)/MW(47)*3	&
                + dOMHi(48)*f_in(48)/MW(48)*3	&   
            !Special Reactive lumping group (default factor 4) 
                + dOMHi(49)*f_in(49)/MW(49)*4	&    
                + dOMHi(50)*f_in(50)/MW(50)*1	&
				+ dOMHi(51)*f_in(51)/MW(51)*2
			else 
    	!alkanes
			dN_OS =                         &
	        !Very Fast reactivity lumping group by using alkane(default Factor 1) ??Added OS vals for 38*42 lumping groups
              +  dOMHi(1)*f_in(1)/MW(1)*1    &
              +  dOMHi(2)*f_in(2)/MW(2)*1    &
              +  dOMHi(3)*f_in(3)/MW(3)*1    &
              +  dOMHi(4)*f_in(4)/MW(4)*1    &
              +  dOMHi(5)*f_in(5)/MW(5)*2    &
              +  dOMHi(6)*f_in(6)/MW(6)*2    &
              +  dOMHi(7)*f_in(7)/MW(7)*2    &
              +  dOMHi(8)*f_in(8)/MW(8)*2    &
        !Fast reactivity lumping group (defaul Factor 1)"
              +  dOMHi(9)*f_in(9)/MW(9)*1    &
              +  dOMHi(10)*f_in(10)/MW(10)*1    &
              +  dOMHi(11)*f_in(11)/MW(11)*1    &
              +  dOMHi(12)*f_in(12)/MW(12)*1    &
              +  dOMHi(13)*f_in(13)/MW(13)*1    &
              +  dOMHi(14)*f_in(14)/MW(14)*3    &
              +  dOMHi(15)*f_in(15)/MW(15)*1    &
              +  dOMHi(16)*f_in(16)/MW(16)*1    &
        !medium reactivity lumping group (default Factor 1)
              +  dOMHi(17)*f_in(17)/MW(17)*1    &
              +  dOMHi(18)*f_in(18)/MW(18)*1    &
              +  dOMHi(19)*f_in(19)/MW(19)*1    &
              +  dOMHi(20)*f_in(20)/MW(20)*3    &
              +  dOMHi(21)*f_in(21)/MW(21)*1    &
              +  dOMHi(22)*f_in(22)/MW(22)*1    &
              +  dOMHi(23)*f_in(23)/MW(23)*2    &
              +  dOMHi(24)*f_in(24)/MW(24)*1    &
        !Slow reactivity lumping group (default factor 1)"
              +  dOMHi(25)*f_in(25)/MW(25)*1    &
              +  dOMHi(26)*f_in(26)/MW(26)*1    &
              +  dOMHi(27)*f_in(27)/MW(27)*1    &
              +  dOMHi(28)*f_in(28)/MW(28)*1    &
              +  dOMHi(29)*f_in(29)/MW(29)*1    &
              +  dOMHi(30)*f_in(30)/MW(30)*1    &
              +  dOMHi(31)*f_in(31)/MW(31)*1    &
              +  dOMHi(32)*f_in(32)/MW(32)*1    &
        !Partitioning  lumping group (default factor 1)
              +  dOMHi(33)*f_in(33)/MW(33)*1    &
              +  dOMHi(34)*f_in(34)/MW(34)*1    &
              +  dOMHi(35)*f_in(35)/MW(35)*1    &
              +  dOMHi(36)*f_in(36)/MW(36)*1    &
              +  dOMHi(37)*f_in(37)/MW(37)*1    &
              +  dOMHi(38)*f_in(38)/MW(38)*1    &
              +  dOMHi(39)*f_in(39)/MW(39)*1    &
              +  dOMHi(40)*f_in(40)/MW(40)*1    &
        !multialcohol  lumping group (default Factor 1)
              +  dOMHi(41)*f_in(41)/MW(41)*1    &
              +  dOMHi(42)*f_in(42)/MW(42)*1    &
              +  dOMHi(43)*f_in(43)/MW(43)*1    &
              +  dOMHi(44)*f_in(44)/MW(44)*3    &
              +  dOMHi(45)*f_in(45)/MW(45)*1    &
              +  dOMHi(46)*f_in(46)/MW(46)*1    &
              +  dOMHi(47)*f_in(47)/MW(47)*1    &
              +  dOMHi(48)*f_in(48)/MW(48)*1    
    end if
			
    END SUBROUTINE	
 
 END module SOA_H5
      