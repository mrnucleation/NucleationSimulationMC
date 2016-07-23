#SHELL = /bin/sh
CUR_DIR := $(shell pwd)
# ====================================
#        Compiler Options
# ====================================
FC := mpif90
#FC := mpifort
#FC := gfortran
CC := mpicc
#OPTIMIZE_FLAGS := -O3 -xHost -ipo
OPTIMIZE_FLAGS := -O3
#OPTIMIZE_FLAGS += -prof-gen -prof-dir=$(CUR_DIR)/profiling
#OPTIMIZE_FLAGS += -prof-use -prof-dir=$(CUR_DIR)/profiling
#OPEN_MP_FLAGS := -fopenmp
#DEBUGFLAGS := -g -fbacktrace -fcheck=all
#DEBUGFLAGS := -fbounds-check
#DEBUGFLAGS := -check bounds
#DEBUGFLAGS += -heap-arrays 1024
#DEBUGFLAGS += -check bounds -traceback -g
DEBUGFLAGS += -pg 
#DEBUGFLAGS += -ffpe-trap=invalid
#DEBUGFLAGS := -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -fbacktrace
COMPFLAGS := $(OPEN_MP_FLAGS) $(DEBUGFLAGS) $(OPTIMIZE_FLAGS)


# ====================================
#        Directory List
# ====================================

SRC := $(CUR_DIR)/src
MODS := $(CUR_DIR)
OBJ := $(CUR_DIR)/objects
TRIAL := $(CUR_DIR)/Trials
ESUB := $(CUR_DIR)/src/EnergyFunctions
CBMC := $(CUR_DIR)/src/CBMC_Functions
CLUSTSUB := $(CUR_DIR)/src/ClusterCriteriaFunctions
SWAP := $(CUR_DIR)/src/SwapFunctions
INPUTSUB := $(CUR_DIR)/src/InputFunctions
ANALYSIS_SUB := $(CUR_DIR)/src/AnalysisFunctions

#RUN_DIR := $(TRIAL)/Trial1_TransRot_Test
#RUN_DIR := $(TRIAL)/Trial2_Water_Test



# Define file extensions
.SUFFIXES:
.SUFFIXES: .f .f90 .o .mod 
# ====================================
#        Compiler specific commands
# ====================================

#MODFLAGS := -I $(MODS) -J $(MODS)
# ====================================
#        Source Files
# ====================================
MOD_FILES := $(MODS)/acceptrates.mod\
		$(MODS)/bendingfunctions.mod\
		$(MODS)/bondstretchfunctions.mod\
		$(MODS)/cbmc_variables.mod\
		$(MODS)/constants.mod\
		$(MODS)/coordinatetypes.mod\
		$(MODS)/coords.mod\
		$(MODS)/energyCriteria.mod\
		$(MODS)/forcefield.mod\
		$(MODS)/forcefieldfunctions.mod\
		$(MODS)/forcefieldvariabletype.mod\
		$(MODS)/improperanglefunctions.mod\
		$(MODS)/indexingfunctions.mod\
		$(MODS)/interenergy_lj_electro.mod\
		$(MODS)/intraenergy_lj_electro.mod\
		$(MODS)/parallelvar.mod\
		$(MODS)/simparameters.mod\
		$(MODS)/torsionalfunctions.mod\
		$(MODS)/umbrellafunctions.mod\
		$(MODS)/units.mod
MOD_SRC := $(SRC)/Common.f90 \
 		$(SRC)/Units.f90 \
 		$(SRC)/ForceFieldFunctions.f90
SRC_ENERGY := $(ESUB)/Bending_Functions.f90 \
            $(ESUB)/BondStretch_Functions.f90 \
            $(ESUB)/LJ_Electro_Functions.f90 \
            $(ESUB)/Pedone_Functions.f90 \
            $(ESUB)/Intra_LJ_Electro_Functions.f90\
            $(ESUB)/Torsional_Functions.f90 \
            $(ESUB)/Improper_Functions.f90 \
            $(ESUB)/Rosen_Boltz_Functions.f90\
            $(ESUB)/Rosen_Pedone_Functions.f90\
            $(ESUB)/EnergyInterfaceFunctions.f90\
            $(ESUB)/EnergyInterfaceFunctions_Pedone.f90\
            $(ESUB)/EnergyPointers.f90
SRC_CRIT:=  $(SRC)/ClusterCriteria_Energy.f90\
            $(SRC)/ClusterCriteria_Distance.f90
SRC_MAIN := $(SRC)/BasicMovement.f90\
            $(SRC)/AnalysisFunctions.f90\
            $(SRC)/MCMove_Module.f90\
            $(SRC)/DebugFunctions.f90\
            $(SRC)/Main.f90\
 		$(SRC)/WHAM.f90\
 		$(SRC)/RandomNew.f90\
 		$(SRC)/OutputFunctions.f90\
 		$(SRC)/Input_Ultility.f90\
 		$(SRC)/UmbrellaSampling.f\
 		$(SRC)/AngleIntegration.f90\
 		$(SRC)/CoordinateFunctions.f90\
		$(SRC)/ReadInput.f90
SRC_MAIN2:=  $(SRC)/ETableFunctions.f90
SRC_CBMC := $(CBMC)/CBMC.f90\
            $(CBMC)/CBMC_Initialize.f90\
            $(CBMC)/CBMC_ConfigGen.f90\
            $(CBMC)/CBMC_Utility.f90\
            $(CBMC)/CBMC_PartialRegrowth.f90\
            $(CBMC)/CBMC_Rosen_AVBMC_ConfigGen.f90
SRC_SWAP := $(SWAP)/AVBMC_EBias_Rosen.f90\
            $(SWAP)/Exchange.f90
#SRC_SWAP := $(SWAP)/AVBMC_EBias.f90
#SRC_SWAP := $(SWAP)/AVBMC_Uniform.f90
SRC_COMPLETE:= $(SRC_ENERGY) $(SRC_CBMC) $(SRC_MAIN2) $(SRC_SWAP) $(SRC_CRIT) $(MOD_SRC) $(SRC_MAIN) 
# ====================================
#        Object Files
# ====================================
OBJ_TEMP:=$(patsubst $(SRC)/%.f, $(OBJ)/%.o, $(SRC_MAIN))
OBJ_MAIN:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(SRC)/%.f, $(OBJ)/%.o, $(SRC_MAIN2))
OBJ_MAIN2:=$(patsubst $(SRC)/%.f90, $(OBJ)/%.o, $(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(ESUB)/%.f,$(OBJ)/%.o,$(SRC_ENERGY))
OBJ_ENERGY:=$(patsubst $(ESUB)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(SRC)/%.f,$(OBJ)/%.o,$(MOD_SRC))
OBJ_MOD:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(CBMC)/%.f, $(OBJ)/%.o, $(SRC_CBMC))
OBJ_CBMC:=$(patsubst $(CBMC)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(SWAP)/%.f,$(OBJ)/%.o,$(SRC_SWAP))
OBJ_SWAP:=$(patsubst $(SWAP)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))

OBJ_TEMP:=$(patsubst $(SRC)/%.f,$(OBJ)/%.o,$(SRC_CRIT))
OBJ_CRIT:=$(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(OBJ_TEMP))
# ====================================
#        Compile Commands
# ====================================


.f.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
.f90.o :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

.f.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

.f90.mod :     
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<

$(OBJ)/%.o: $(ESUB)/%.f 
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
$(OBJ)/%.o: $(ESUB)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
            
$(OBJ)/%.o: $(CBMC)/%.f 
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<		        
$(OBJ)/%.o: $(CBMC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
            
$(OBJ)/%.o: $(SWAP)/%.f 
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
$(OBJ)/%.o: $(SWAP)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
            
$(OBJ)/%.o: $(SRC)/%.f
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<
$(OBJ)/%.o: $(SRC)/%.f90
		@echo Creating $<
		@$(FC) $(COMPFLAGS) $(MODFLAGS) -c -o $@ $<




# ====================================
#        Compile Commands
# ====================================
default: startUP createMods energyFunctions generalNucleation finale
engOnly: startUP energyFunctions generalNucleation finale
quick: startUP generalNucleation finale
comp_run: startUP createMods generalNucleation run finale
neat: startUP createMods generalNucleation removeObject finale
clean: removeObjects removeExec finale    
    
createMods: $(MOD_SRC)
		@echo =============================================
		@echo            Creating Module Files
		@echo =============================================		
		@echo		
		@echo  -------- Compiling VariablePrecision.f90
		@$(FC) -c $(SRC)/VariablePrecision.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/VariablePrecision.o
		@echo  -------- Compiling Units.f90
		@$(FC) -c $(SRC)/Units.f90  $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/Units.o		
		@echo  -------- Compiling ForceFieldFunctions.f90
		@$(FC) -c $(SRC)/ForceFieldFunctions.f90 $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/ForceFieldFunctions.o				
		@echo  -------- Compiling Common.f90
		@$(FC) -c $(SRC)/Common.f90  $(COMPFLAGS) $(MODFLAGS) -o $(OBJ)/Common.o
		@echo =============================================
		@echo            Creating Object Files
		@echo =============================================	
		@echo  	

            
energyFunctions: $(OBJ_CRIT) $(OBJ_ENERGY) 
		@$(FC) $(COMPFLAGS)  $< -c
      
        
generalNucleation:  $(OBJ_CRIT) $(OBJ_ENERGY) $(OBJ_MAIN2) $(OBJ_MOD) $(OBJ_CBMC) $(OBJ_SWAP) $(OBJ_MAIN) 
		@echo =============================================
		@echo     Compiling and Linking Source Files
		@echo =============================================	
		@$(FC) $(COMPFLAGS) $(MODFLAGS)  $^ -o $@ 		
		
run:
		@echo =============================================
		@echo            Running Code
		@echo =============================================		
		cd $(RUN_DIR)\
		&& ../../generalNucleation.exe	
		
runMPI:
		@echo =============================================
		@echo            Running Code
		@echo =============================================		
		cd $(RUN_DIR)\
		&& mpirun -np 2 ../../generalNucleation.exe		
		
startUP:
		@echo ==================================================================
		@echo ---------------------- Begin ---------------------------------		
		@echo Current Directory:$(CUR_DIR)		
		@echo Compiler and Flags used:	$(FC) $(COMPFLAGS) 		
		@echo		

finale:
		@echo
		@echo ---------------------- Finished! ---------------------------------
		@echo ==================================================================		
     
removeObjects:
		@echo =============================================
		@echo            Cleaning Directory
		@echo =============================================		
		@echo		
		@rm -f ./*.o ./*.mod				
		@rm -f $(SRC)/*.o $(SRC)/*.mod		
		@rm -f $(MODS)/*.o $(MODS)/*.mod			
		@rm -f $(SRC)/*/*.o $(SRC)/*/*.mod
		@rm -f $(OBJ)/*.o		

removeExec:
		@rm -f $(CUR_DIR)/generalNucleation
		@rm -f $(CUR_DIR)/generalNucleation.exe            




# ====================================
#        Dependencies
# ====================================
$(CUR_DIR)/neighbortable.mod: $(SRC)/ETableFunctions.f90 $(SRC)/ForceFieldFunctions.f90 $(SRC)/Common.f90
$(CUR_DIR)/simparameters.mod: $(SRC)/Common.f90 $(SRC)/ForceFieldFunctions.f90
$(CUR_DIR)/cbmc_module.mod: $(CBMC)/CBMC.f90
$(SRC_SWAP): $(CUR_DIR)/neighbortable.mod
$(CUR_DIR)/avbmc_module.mod: $(SRC_SWAP) $(SRC)/ETableFunctions.f90
$(CUR_DIR)/simplemcmoves_module.mod: $(SRC)/BasicMovement.f90
$(CUR_DIR)/movetypemodule.mod: $(SRC)/MCMove_Module.f90 $(CUR_DIR)/avbmc_module.mod $(CUR_DIR)/simplemcmoves_module.mod $(CUR_DIR)/cbmc_module.mod
$(CUR_DIR)/energypointers.mod: $(ESUB)/EnergyPointers.f90 $(ESUB)/EnergyInterfaceFunctions.f90 $(CUR_DIR)/e_interface.mod
		
