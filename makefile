# Makefile for GPUPBTE #
EXECUTABLE=GPUPBTE.x
BUILD_DIR := build
BIN_DIR := bin
RELEASE_DIR := release
DEBUG_DIR := debug
CC := nvfortran
LL := nvfortran
NVINC := -L./external/spglib-1.7.4-intel16/lib
SRCS1ST := src/Cal_ACC_V3.f90 src/Cal_ACC_V4.f90 src/Data.f90 src/Disp.f90 src/IFC.f90 src/Job_Task.f90 src/Math.f90		\
src/Quantum.f90 src/Smearing.f90 src/Symm.f90 src/Symmetry.f90 src/WriteRead.f90 
SRCS2ND := src/Cal_DOS_3ph.f90 src/Cal_DOS_4ph.f90 src/Cal_DOS.f90 src/Cal_GrpVel.f90 src/Cal_RT_ISO.f90 src/Cal_RTA.f90	\
src/Cal_RTGPU_PH3.f90 src/Cal_RTGPU_PH4.f90 src/Cal_TC.f90 src/Cdiagh.f90 src/Gen_QPoint.f90 src/GPUPBTE.f90 				\
src/Read_Input.f90 src/RT.f90 src/Symm_Operation.f90 
MODS := accrelate_v3.mod accrelate_v4.mod crystal.mod mesh.mod task.mod tmpspace.mod disp.mod ifc.mod job_task.mod 			\
vectorfunctions.mod quantum.mod smearingfunctions.mod symm.mod symmetry.mod writeread.mod
SRCS := $(SRCS1ST) $(SRCS2ND)
OBJS1ST := $(subst src/,$(BUILD_DIR)/$(RELEASE_DIR)/,$(addsuffix .o,$(basename $(SRCS1ST))))
OBJS2ND := $(subst src/,$(BUILD_DIR)/$(RELEASE_DIR)/,$(addsuffix .o,$(basename $(SRCS2ND))))
OBJS := $(OBJS1ST) $(OBJS2ND)
NVCCFLAGS := -Mmkl -Mcuda -lsymspg -traceback -acc=gpu -gpu=cuda10.2 -Wall
NVCCFLAGS_RELEASE= -fast -DNDEBUG -module $(BUILD_DIR)/$(RELEASE_DIR)
NVCCFLAGS_DEBUG := -DDEBUG -g -module $(BUILD_DIR)/$(DEBUG_DIR)
all: release
release: $(BIN_DIR)/$(RELEASE_DIR)/$(EXECUTABLE)
debug: $(BIN_DIR)/$(DEBUG_DIR)/$(EXECUTABLE)
# Release
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_ACC_V3.o : src/Cal_ACC_V3.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_ACC_V3.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_ACC_V4.o : $(BUILD_DIR)/$(RELEASE_DIR)/IFC.o $(BUILD_DIR)/$(RELEASE_DIR)/Data.o 			\
src/Cal_ACC_V4.f90 
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_ACC_V4.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Data.o :src/Data.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Data.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Disp.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/IFC.o src/Disp.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Disp.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/IFC.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o src/IFC.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/IFC.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Job_Task.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
src/Job_Task.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Job_Task.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Math.o :src/Math.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Math.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o : $(BUILD_DIR)/$(RELEASE_DIR)/Data.o src/Quantum.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Quantum.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o src/Smearing.f90
	@mkdir -p $(dir $@) 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Smearing.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Symm.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o src/Symm.f90 
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Symm.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Symmetry.o :src/Symmetry.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Symmetry.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/WriteRead.o :$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
src/WriteRead.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/WriteRead.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_DOS_3ph.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 			\
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o src/Cal_DOS_3ph.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_DOS_3ph.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_DOS_4ph.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 			\
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o src/Cal_DOS_4ph.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_DOS_4ph.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_DOS.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o src/Cal_DOS.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_DOS.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_GrpVel.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
$(BUILD_DIR)/$(RELEASE_DIR)/Disp.o src/Cal_GrpVel.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_GrpVel.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_RT_ISO.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o				\
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o $(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o src/Cal_RT_ISO.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_RT_ISO.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_RTA.o: $(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
$(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o src/Cal_RTA.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_RTA.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_RTGPU_PH3.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 			\
$(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o $(BUILD_DIR)/$(RELEASE_DIR)/IFC.o $(BUILD_DIR)/$(RELEASE_DIR)/Cal_ACC_V3.o 			\
src/Cal_RTGPU_PH3.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_RTGPU_PH3.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_RTGPU_PH4.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 			\
$(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o $(BUILD_DIR)/$(RELEASE_DIR)/IFC.o $(BUILD_DIR)/$(RELEASE_DIR)/Cal_ACC_V4.o 			\
src/Cal_RTGPU_PH4.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_RTGPU_PH4.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cal_TC.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 					\
$(BUILD_DIR)/$(RELEASE_DIR)/Quantum.o src/Cal_TC.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cal_TC.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Cdiagh.o: src/Cdiagh.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Cdiagh.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Gen_QPoint.o : $(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 			\
src/Gen_QPoint.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Gen_QPoint.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/GPUPBTE.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/IFC.o 					\
$(BUILD_DIR)/$(RELEASE_DIR)/Disp.o $(BUILD_DIR)/$(RELEASE_DIR)/WriteRead.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 				\
$(BUILD_DIR)/$(RELEASE_DIR)/Job_Task.o src/GPUPBTE.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/GPUPBTE.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/Read_Input.o: $(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Math.o 			\
$(BUILD_DIR)/$(RELEASE_DIR)/Smearing.o src/Read_Input.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Read_Input.f90 -o $@
$(BUILD_DIR)/$(RELEASE_DIR)/RT.o:$(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o src/RT.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/RT.f90 -o $@ 
$(BUILD_DIR)/$(RELEASE_DIR)/Symm_Operation.o: $(BUILD_DIR)/$(RELEASE_DIR)/Data.o $(BUILD_DIR)/$(RELEASE_DIR)/Symm.o 		\
$(BUILD_DIR)/$(RELEASE_DIR)/Symmetry.o src/Symm_Operation.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -c src/Symm_Operation.f90  -o $@	
$(BIN_DIR)/$(RELEASE_DIR)/$(EXECUTABLE) : $(OBJS)
	@mkdir -p $(dir $@)
	$(LL) $^ $(NVCCFLAGS) $(NVCCFLAGS_RELEASE) $(NVINC) -o $@
	@chmod 777 $(BIN_DIR)/$(RELEASE_DIR)/$(EXECUTABLE)
# Debug
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_ACC_V3.o : src/Cal_ACC_V3.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_ACC_V3.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_ACC_V4.o : $(BUILD_DIR)/$(DEBUG_DIR)/IFC.o $(BUILD_DIR)/$(DEBUG_DIR)/Data.o 					\
src/Cal_ACC_V4.f90 
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_ACC_V4.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Data.o :src/Data.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Data.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Disp.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/IFC.o src/Disp.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Disp.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/IFC.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o src/IFC.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/IFC.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Job_Task.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
src/Job_Task.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Job_Task.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Math.o :src/Math.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Math.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o : $(BUILD_DIR)/$(DEBUG_DIR)/Data.o src/Quantum.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Quantum.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o src/Smearing.f90
	@mkdir -p $(dir $@) 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Smearing.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Symm.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o src/Symm.f90 
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Symm.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Symmetry.o :src/Symmetry.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Symmetry.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/WriteRead.o :$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
src/WriteRead.f90
	@mkdir -p $(dir $@)
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/WriteRead.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_DOS_3ph.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o src/Cal_DOS_3ph.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_DOS_3ph.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_DOS_4ph.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o src/Cal_DOS_4ph.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_DOS_4ph.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_DOS.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 						\
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o src/Cal_DOS.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_DOS.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_GrpVel.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
$(BUILD_DIR)/$(DEBUG_DIR)/Disp.o src/Cal_GrpVel.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_GrpVel.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_RT_ISO.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o					\
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o $(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o src/Cal_RT_ISO.f90
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_RT_ISO.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_RTA.o: $(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 						\
$(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o src/Cal_RTA.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_RTA.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_RTGPU_PH3.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 				\
$(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o $(BUILD_DIR)/$(DEBUG_DIR)/IFC.o $(BUILD_DIR)/$(DEBUG_DIR)/Cal_ACC_V3.o 					\
src/Cal_RTGPU_PH3.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_RTGPU_PH3.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_RTGPU_PH4.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 				\
$(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o $(BUILD_DIR)/$(DEBUG_DIR)/IFC.o $(BUILD_DIR)/$(DEBUG_DIR)/Cal_ACC_V4.o 					\
src/Cal_RTGPU_PH4.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_RTGPU_PH4.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cal_TC.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 						\
$(BUILD_DIR)/$(DEBUG_DIR)/Quantum.o src/Cal_TC.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cal_TC.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Cdiagh.o: src/Cdiagh.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Cdiagh.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Gen_QPoint.o : $(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
src/Gen_QPoint.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Gen_QPoint.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/GPUPBTE.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/IFC.o 						\
$(BUILD_DIR)/$(DEBUG_DIR)/Disp.o $(BUILD_DIR)/$(DEBUG_DIR)/WriteRead.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 					\
$(BUILD_DIR)/$(DEBUG_DIR)/Job_Task.o src/GPUPBTE.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/GPUPBTE.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/Read_Input.o: $(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Math.o 					\
$(BUILD_DIR)/$(DEBUG_DIR)/Smearing.o src/Read_Input.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Read_Input.f90 -o $@
$(BUILD_DIR)/$(DEBUG_DIR)/RT.o:$(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o src/RT.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/RT.f90 -o $@ 
$(BUILD_DIR)/$(DEBUG_DIR)/Symm_Operation.o: $(BUILD_DIR)/$(DEBUG_DIR)/Data.o $(BUILD_DIR)/$(DEBUG_DIR)/Symm.o 				\
$(BUILD_DIR)/$(DEBUG_DIR)/Symmetry.o src/Symm_Operation.f90 
	$(CC) $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -c src/Symm_Operation.f90  -o $@
$(BIN_DIR)/$(DEBUG_DIR)/$(EXECUTABLE) : $(subst $(RELEASE_DIR),$(DEBUG_DIR),$(OBJS))
	@mkdir -p $(dir $@)
	$(LL) $^ $(NVCCFLAGS) $(NVCCFLAGS_DEBUG) $(NVINC) -o $@
	@chmod 777 $(BIN_DIR)/$(DEBUG_DIR)/$(EXECUTABLE)
.PHONY : all release debug clean help 
# Clean generated files.
clean:
	@echo "Clean"
	@rm -rf $(BUILD_DIR) 
	@rm -rf $(BIN_DIR)
# Provide usage instructions
help: 

	@echo "                                    Welcome to GPUPBTE !                             Version 1.0"       
	@echo "   Usage:"
	@echo "     	make help       	Shows this help documentation"
	@echo "     	make all        	Build the default configuraiton (release)"
	@echo "     	make release    	Build the release executable ($(BIN_DIR)/$(RELEASE_DIR)/$(EXECUTABLE)) "
	@echo "     	make debug      	Build the release executable ($(BIN_DIR)/$(DEBUG_DIR)/$(EXECUTABLE)) "
	@echo "     	make clean      	Clean the build and bin directories"