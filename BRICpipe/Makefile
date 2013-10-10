# Makefile that generates co-registered T1W, T2W, FLAIR and T2*W volumes,
# tissue masks, and basal ganglia masks for ID segmentation.
# Author: Andreas Glatz <a.glatz@sms.ed.ac.uk>
#


# General variables ###########################################################
MISSING_SFX = .missing
COST_SFX = .cost
KEEP_FLAG ?= 1
SCRIPT_DIR ?= 
#$(info SCRIPT_DIR = $(SCRIPT_DIR))


# Fsl vars ####################################################################
# We keep this one flexible since there is a difference
# between the different Linux distributions.
FSL_PFX = $(shell $(SCRIPT_DIR)setup_environment.sh FSL_PFX)
#$(info FSL_PFX = $(FSL_PFX))
FSL_BIN_PFX = $(FSL_PFX)/bin/
FSL_DATA = $(FSL_PFX)/data/
# .nii.gz prefered so we can easily check if target was
# created or not... see 'check_target'.
FSL_EXT_VOL =
ifeq ($(FSLOUTPUTTYPE),NIFTI_GZ)
	FSL_EXT_VOL = .nii.gz
endif
ifeq ($(FSLOUTPUTTYPE),NIFTI)
	FSL_EXT_VOL = .nii
endif
ifeq ($(FSL_EXT_VOL),)
	$(error FSLOUTPUTTYPE not supported!))
endif
FSL_EXT_VOL_SED = sed 's+^.*\\/++g' | sed 's+\\([^.]*\\)\\..*+\\1+g' | uniq
FSL_EXT_MAT = .mat
FSL_EXT_MAT_SED = sed s+\\.mat$$++g


# Slicer vars #################################################################
SLC_PFX = $(shell $(SCRIPT_DIR)setup_environment.sh SLC_PFX)
#$(info SLC_PFX = $(SLC_PFX))


# Input data ##################################################################
INPUT_CSV = Input.csv
OUTPUT_CSV = Output.csv

# Local files #################################################################
# Template for registered T1W volume
L_T1W = T1W
# Registered T1W volume
L_T1W_VOL = $(L_T1W)$(FSL_EXT_VOL)
L_T1W1_VOL = $(L_T1W)1$(FSL_EXT_VOL)
# Template for registered, brain-extracted T1W
L_T1W_BRAIN = $(L_T1W)_brain
# Registered, brain-extracted T1W
L_T1W_BRAIN_VOL = $(L_T1W_BRAIN)$(FSL_EXT_VOL)
# Template for registerd, brain-extracted, bias-corrected T1W
L_T1W_BRAIN_RESTORE = $(L_T1W_BRAIN)_restore
# Registerd, brain-extracted, bias-corrected T1W
L_T1W_BRAIN_RESTORE_VOL = $(L_T1W_BRAIN_RESTORE)$(FSL_EXT_VOL)

# Template for unmodified (unregistered) T1W
L_T1W_NOTREG = $(L_T1W)_notreg
# Unmodified (unregistered) T1W volume
L_T1W_NOTREG_VOL = $(L_T1W_NOTREG)$(FSL_EXT_VOL)
L_T1W_NOTREG_TRUNC_VOL = $(L_T1W_NOTREG)_trunc$(FSL_EXT_VOL)
L_T1W_NOTREG_MAT = $(L_T1W_NOTREG)$(FSL_EXT_MAT)
L_T1W_NOTREG_NOISE_MASK = $(L_T1W_NOTREG)_noise_mask$(FSL_EXT_VOL)
# Template for unregistered, brain-extracted T1W (1st try)
L_T1W_NOTREG_BRAIN1 = $(L_T1W_NOTREG)_brain1
# Brain mask of unregistered, brain-extracted T1W (1st try)
L_T1W_NOTREG_BRAIN1_MASK = $(L_T1W_NOTREG_BRAIN1)_mask$(FSL_EXT_VOL)
# Unregistered, brain-extracted T1W (1st try)
L_T1W_NOTREG_BRAIN1_VOL = $(L_T1W_NOTREG_BRAIN1)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted
# T1W (1st try) to unregistered, brain-extracted GRE
L_T1W_NOTREG_BRAIN1_MAT = $(L_T1W_NOTREG_BRAIN1)$(FSL_EXT_MAT)
# Inverse of above matrix
L_T1W_NOTREG_BRAIN1_INV_MAT = $(L_T1W_NOTREG_BRAIN1)_inv$(FSL_EXT_MAT)
# Template for unregistered, brain-extracted T1W (final)
L_T1W_NOTREG_BRAIN = $(L_T1W_NOTREG)_brain
# Unregistered, brain-extracted T1W (final)
L_T1W_NOTREG_BRAIN_VOL = $(L_T1W_NOTREG_BRAIN)$(FSL_EXT_VOL)
L_T1W_NOTREG_BRAIN_MAT = $(L_T1W_NOTREG_BRAIN)$(FSL_EXT_MAT)
# Brain mask of unregistered, brain-extracted T1W (final)
L_T1W_NOTREG_BRAIN_MASK = $(L_T1W_NOTREG_BRAIN)_mask$(FSL_EXT_VOL)
# Template for unregistered, brain-extracted, bias_corrected T1W
L_T1W_NOTREG_BRAIN_RESTORE = $(L_T1W_NOTREG_BRAIN)_restore
# Unregistered, brain-extracted, bias_corrected T1W
L_T1W_NOTREG_BRAIN_RESTORE_VOL = $(L_T1W_NOTREG_BRAIN_RESTORE)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted,
# bias_corrected T1W to unregistered, brain-extracted, bias-corrected GRE
L_T1W_NOTREG_BRAIN_RESTORE_MAT = $(L_T1W_NOTREG_BRAIN_RESTORE)$(FSL_EXT_MAT)
L_T1W_NOTREG_BRAIN_RESTORE_INV_MAT = $(L_T1W_NOTREG_BRAIN_RESTORE)_inv$(FSL_EXT_MAT)
L_T1W_NOTREG_BRAIN_RESTORE_TOFLAIR_MAT = $(L_T1W_NOTREG_BRAIN_RESTORE)_toflair$(FSL_EXT_MAT)
L_T1W_NOTREG_BRAIN_RESTORE_TOT2W_MAT = $(L_T1W_NOTREG_BRAIN_RESTORE)_tot2w$(FSL_EXT_MAT)
L_T1W_NOTREG_BRAIN_RESTORE_SUSAN_VOL = $(L_T1W_NOTREG_BRAIN_RESTORE)_susan$(FSL_EXT_VOL)

# Template for registered FLAIR volume
L_FLAIR = FLAIR
# Registered FLAIR volume
L_FLAIR_VOL = $(L_FLAIR)$(FSL_EXT_VOL)
# Template for registered, brain-extracted FLAIR
L_FLAIR_BRAIN = $(L_FLAIR)_brain
# Registered, brain-extracted FLAIR
L_FLAIR_BRAIN_VOL = $(L_FLAIR_BRAIN)$(FSL_EXT_VOL)
# Template for registerd, brain-extracted, bias-corrected FLAIR
L_FLAIR_BRAIN_RESTORE = $(L_FLAIR_BRAIN)_restore
# Registerd, brain-extracted, bias-corrected FLAIR
L_FLAIR_BRAIN_RESTORE_VOL = $(L_FLAIR_BRAIN_RESTORE)$(FSL_EXT_VOL)
# Template for unmodified (unregistered) FLAIR
L_FLAIR_NOTREG = $(L_FLAIR)_notreg
# Unmodified (unregistered) FLAIR volume
L_FLAIR_NOTREG_VOL = $(L_FLAIR_NOTREG)$(FSL_EXT_VOL)
L_FLAIR_NOTREG_MAT = $(L_FLAIR_NOTREG)$(FSL_EXT_MAT)
# Template for unregistered, brain-extracted FLAIR (1st try)
L_FLAIR_NOTREG_BRAIN1 = $(L_FLAIR_NOTREG)_brain1
# Brain mask of unregistered, brain-extracted FLAIR (1st try)
L_FLAIR_NOTREG_BRAIN1_MASK = $(L_FLAIR_NOTREG_BRAIN1)_mask$(FSL_EXT_VOL)
# Unregistered, brain-extracted FLAIR (1st try)
L_FLAIR_NOTREG_BRAIN1_VOL = $(L_FLAIR_NOTREG_BRAIN1)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted
# FLAIR (1st try) to unregistered, brain-extracted GRE
L_FLAIR_NOTREG_BRAIN1_MAT = $(L_FLAIR_NOTREG_BRAIN1)$(FSL_EXT_MAT)
L_FLAIR_NOTREG_BRAIN1_INV_MAT = $(L_FLAIR_NOTREG_BRAIN1)_inv$(FSL_EXT_MAT)
# Template for unregistered, brain-extracted FLAIR (final)
L_FLAIR_NOTREG_BRAIN = $(L_FLAIR_NOTREG)_brain
# Unregistered, brain-extracted FLAIR (final)
L_FLAIR_NOTREG_BRAIN_VOL = $(L_FLAIR_NOTREG_BRAIN)$(FSL_EXT_VOL)
L_FLAIR_NOTREG_BRAIN_MAT = $(L_FLAIR_NOTREG_BRAIN)$(FSL_EXT_MAT)
# Brain mask of unregistered, brain-extracted FLAIR (final)
L_FLAIR_NOTREG_BRAIN_MASK = $(L_FLAIR_NOTREG_BRAIN)_mask$(FSL_EXT_VOL)
# Template for unregistered, brain-extracted, bias_corrected FLAIR
L_FLAIR_NOTREG_BRAIN_RESTORE = $(L_FLAIR_NOTREG_BRAIN)_restore
# Unregistered, brain-extracted, bias_corrected FLAIR
L_FLAIR_NOTREG_BRAIN_RESTORE_VOL = $(L_FLAIR_NOTREG_BRAIN_RESTORE)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted,
# bias_corrected FLAIR to unregistered, brain-extracted, bias-corrected GRE
L_FLAIR_NOTREG_BRAIN_RESTORE_MAT = $(L_FLAIR_NOTREG_BRAIN_RESTORE)$(FSL_EXT_MAT)
L_FLAIR_NOTREG_BRAIN_RESTORE_INV_MAT = $(L_FLAIR_NOTREG_BRAIN_RESTORE)_inv$(FSL_EXT_MAT)

# Template for registered T2W volume
L_T2W = T2W
# Registered T2W volume
L_T2W_VOL = $(L_T2W)$(FSL_EXT_VOL)
# Template for registered, brain-extracted T2W
L_T2W_BRAIN = $(L_T2W)_brain
# Registered, brain-extracted T2W
L_T2W_BRAIN_VOL = $(L_T2W_BRAIN)$(FSL_EXT_VOL)
# Template for registerd, brain-extracted, bias-corrected T2W
L_T2W_BRAIN_RESTORE = $(L_T2W_BRAIN)_restore
# Registerd, brain-extracted, bias-corrected T2W
L_T2W_BRAIN_RESTORE_VOL = $(L_T2W_BRAIN_RESTORE)$(FSL_EXT_VOL)
# Template for unmodified (unregistered) T2W
L_T2W_NOTREG = $(L_T2W)_notreg
# Unmodified (unregistered) T2W volume
L_T2W_NOTREG_VOL = $(L_T2W_NOTREG)$(FSL_EXT_VOL)
L_T2W_NOTREG_MAT = $(L_T2W_NOTREG)$(FSL_EXT_MAT)
# Template for unregistered, brain-extracted T2W (1st try)
L_T2W_NOTREG_BRAIN1 = $(L_T2W_NOTREG)_brain1
# Brain mask of unregistered, brain-extracted T2W (1st try)
L_T2W_NOTREG_BRAIN1_MASK = $(L_T2W_NOTREG_BRAIN1)_mask$(FSL_EXT_VOL)
# Unregistered, brain-extracted T2W (1st try)
L_T2W_NOTREG_BRAIN1_VOL = $(L_T2W_NOTREG_BRAIN1)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted
# T2W (1st try) to unregistered, brain-extracted GRE
L_T2W_NOTREG_BRAIN1_MAT = $(L_T2W_NOTREG_BRAIN1)$(FSL_EXT_MAT)
L_T2W_NOTREG_BRAIN1_INV_MAT = $(L_T2W_NOTREG_BRAIN1)_inv$(FSL_EXT_MAT)
# Template for unregistered, brain-extracted T2W (final)
L_T2W_NOTREG_BRAIN = $(L_T2W_NOTREG)_brain
# Unregistered, brain-extracted T2W (final)
L_T2W_NOTREG_BRAIN_VOL = $(L_T2W_NOTREG_BRAIN)$(FSL_EXT_VOL)
L_T2W_NOTREG_BRAIN_MAT = $(L_T2W_NOTREG_BRAIN)$(FSL_EXT_MAT)
# Brain mask of unregistered, brain-extracted T2W (final)
L_T2W_NOTREG_BRAIN_MASK = $(L_T2W_NOTREG_BRAIN)_mask$(FSL_EXT_VOL)
# Template for unregistered, brain-extracted, bias_corrected T2W
L_T2W_NOTREG_BRAIN_RESTORE = $(L_T2W_NOTREG_BRAIN)_restore
# Unregistered, brain-extracted, bias_corrected T2W
L_T2W_NOTREG_BRAIN_RESTORE_VOL = $(L_T2W_NOTREG_BRAIN_RESTORE)$(FSL_EXT_VOL)
# Transformation matrix of registration of unregistered, brain-extracted,
# bias_corrected T2W to unregistered, brain-extracted, bias-corrected GRE
L_T2W_NOTREG_BRAIN_RESTORE_MAT = $(L_T2W_NOTREG_BRAIN_RESTORE)$(FSL_EXT_MAT)
L_T2W_NOTREG_BRAIN_RESTORE_INV_MAT = $(L_T2W_NOTREG_BRAIN_RESTORE)_inv$(FSL_EXT_MAT)

# Template for unmodified GRE volume
L_GRE = GRE
# Unmodified GRE volume
L_GRE_VOL = $(L_GRE)$(FSL_EXT_VOL)
# Template for brain-extracted GRE
L_GRE_BRAIN = $(L_GRE)_brain
# Brain-extracted GRE
#L_GRE_BRAIN_BET_VOL = $(L_GRE_BRAIN)_bet$(FSL_EXT_VOL)
L_GRE_BRAIN_VOL = $(L_GRE_BRAIN)$(FSL_EXT_VOL)
# Brain mask of brain-extracted GRE
L_GRE_BRAIN_BET_MASK = $(L_GRE_BRAIN)_bet_mask$(FSL_EXT_VOL)
L_GRE_BRAIN_MASK = $(L_GRE_BRAIN)_mask$(FSL_EXT_VOL)
L_GRE_BRAIN_MASK1 = $(L_GRE_BRAIN)_mask1$(FSL_EXT_VOL)
# Template for brain-extracted, bias-corrected GRE
L_GRE_BRAIN_RESTORE = $(L_GRE_BRAIN)_restore
# Brain-extracted, bias-corrected GRE
L_GRE_BRAIN_RESTORE_VOL = $(L_GRE_BRAIN_RESTORE)$(FSL_EXT_VOL)

# Tissue masks
L_SE_MASK = SE_mask$(FSL_EXT_VOL)
L_SE_LORD = SE_mask.txt
L_SE_HARD_MASK = SE_hard_mask$(FSL_EXT_VOL)
L_SE_HARD_LORD = SE_hard_mask.txt
L_WM_MASK = WM_mask$(FSL_EXT_VOL)
L_GM_MASK = GM_mask$(FSL_EXT_VOL)
L_CS_MASK = CS_mask$(FSL_EXT_VOL)
L_AR_MASK = AR_mask$(FSL_EXT_VOL)
L_AR_ROI_MASK = AR_roi_mask$(FSL_EXT_VOL)

# Basal ganglia mask from registered T1W volume
L_BG_MASK = BG_mask$(FSL_EXT_VOL)
L_BG_NOTREG_MASK = BG_notreg_mask$(FSL_EXT_VOL)
L_BG_NOTREG_MAT = BG_notreg$(FSL_EXT_MAT)
L_BG_NOTREG_MNI_VOL = $(FSL_PFX)/data/standard/MNI152_T1_1mm_brain$(FSL_EXT_VOL)
L_RO_MASK = RO_mask$(FSL_EXT_VOL)

# Reference masks
L_FE_REF_MASK = FE_ref_mask$(FSL_EXT_VOL)
L_FE_REF_NOTREG_MASK = FE_ref_notreg_mask$(FSL_EXT_VOL)
L_FE_ROI = T2swHypoT1wiso
L_FE_ROI_MASK = $(L_FE_ROI)_mask$(FSL_EXT_VOL)
L_FE_ROI_BIN_MASK = $(L_FE_ROI)_bin_mask$(FSL_EXT_VOL)
L_FE_ROI_BIN_NOTREG_MASK = $(L_FE_ROI)_bin_notreg_mask$(FSL_EXT_VOL)
L_FE_ROI_BIN_NOTREG_MNI_MASK = $(L_FE_ROI)_bin_notreg_mni_mask$(FSL_EXT_VOL)
L_IC_REF_MASK = IC_ref_mask$(FSL_EXT_VOL)
L_WM_REF_MASK = WM_ref_mask$(FSL_EXT_VOL)

# Rsync stats files
EXT_RSYNC = .rsync
L_T1W_RSYNC = $(L_T1W_NOTREG)$(EXT_RSYNC)
L_GRE_RSYNC = $(L_GRE)$(EXT_RSYNC)
L_FLAIR_RSYNC = $(L_FLAIR)$(EXT_RSYNC)
L_T2W_RSYNC = $(L_T2W)$(EXT_RSYNC)
L_FEM_RSYNC = FE_mask$(EXT_RSYNC)
L_RIM_RSYNC = REF_ICV_mask$(EXT_RSYNC)
L_RWM_RSYNC = REF_WM_mask$(EXT_RSYNC)

# MNI registration
L_T1W_MNI_REF_VOL = ../asps057/T1W_notreg_trunc_mni2mm.nii.gz
L_T1W_BRAIN_MNI_REF_VOL = ../asps057/T1W_notreg_brain_restore_mni2mm.nii.gz
L_T1W_MNI_REF_FNIRT_CFG = ../node_scripts/T1_2_MNI152_2mm_ge65.cnf
L_T1W_NOTREG_TRUNC_MNI = $(L_T1W_NOTREG)_trunc_mni
L_T1W_NOTREG_TRUNC_MNI_VOL = $(L_T1W_NOTREG_TRUNC_MNI)$(FSL_EXT_VOL)
L_T1W_NOTREG_TRUNC_MNI_WARP = $(L_T1W_NOTREG_TRUNC_MNI)_warp$(FSL_EXT_VOL)
L_T1W_NOTREG_TRUNC_MNI_MAT = $(L_T1W_NOTREG_TRUNC_MNI)$(FSL_EXT_MAT)

# Copy and convert remote files rules #########################################
.PHONY: all
all:	$(L_T1W_BRAIN_RESTORE_VOL) $(L_RO_MASK) \
		$(L_GRE_BRAIN_RESTORE_VOL) $(L_WM_MASK) $(L_GM_MASK) \
		$(L_CS_MASK) $(L_AR_MASK) $(L_AR_ROI_MASK) \
		$(L_FE_REF_MASK) $(L_FE_ROI_MASK) \
		$(INPUT_CSV) $(OUTPUT_CSV) \
		$(L_T2W_BRAIN_RESTORE_VOL) \
		$(L_T1W_NOTREG_TRUNC_MNI_VOL) $(L_FE_ROI_BIN_NOTREG_MNI_MASK)
	@OUT_FILES=`echo "$+" | sed s+$(INPUT_CSV)++g | sed s+$(OUTPUT_CSV)++g`; \
	 DST_PATH=`cat $(OUTPUT_CSV)`; \
	 $(copy_result)
	@$(remove_intermediate)

define copy_result
if [ ! -z $${DST_PATH} ]; \
then \
	mkdir -p $${DST_PATH}; \
	$(SCRIPT_DIR)copy_result.sh $${OUT_FILES} $${DST_PATH}; \
fi;
endef

define remove_intermediate
if [ $(KEEP_FLAG) -eq 0 ]; \
then \
	ALL_FILES=output_all.rmi; \
	EXC_FILES=output_exc.rmi; \
	rm -f $${ALL_FILES} $${EXC_FILES}; \
	ls | sort > $${ALL_FILES}; \
	EXC_FILES_TMP=tmp_$${EXC_FILES}; \
	rm -f $${EXC_FILES_TMP}; \
	for FILE in $+ `ls *$(MISSING_SFX) | xargs`; \
	do \
		echo $${FILE} >> $${EXC_FILES_TMP}; \
	done; \
	for FILE in `ls *$(EXT_RSYNC) | xargs`; \
	do \
		TMP=`cat $${FILE} | cut -d ' ' -f 2`; \
		ls $${TMP}.* >> $${EXC_FILES_TMP}; \
	done; \
	cat $${EXC_FILES_TMP} | sort > $${EXC_FILES}; \
	REM_FILES=`comm -23 $${ALL_FILES} $${EXC_FILES} | xargs`; \
	echo "Removing $${REM_FILES}..."; \
	rm -rf $${REM_FILES} $${EXC_FILES_TMP}; \
	echo "done."; \
fi;
endef

# Required targets ############################################################
$(L_T1W_NOTREG_VOL): $(L_T1W_RSYNC)
	@$(convert_target)
	@$(check_target)

$(L_GRE_VOL): $(L_GRE_RSYNC)
	@$(convert_target)
	@$(check_target)

$(L_FLAIR_NOTREG_VOL): $(L_FLAIR_RSYNC)
	@$(convert_target)
	@$(check_target)

$(L_T2W_NOTREG_VOL): $(L_T2W_RSYNC)
	@$(convert_target)
	@$(check_target)

# Optional targets ############################################################
$(L_FE_REF_MASK): $(L_FEM_RSYNC) $(L_GRE_VOL)
	@$(convert_target)
	@$(check_target)

$(L_IC_REF_MASK): $(L_RIM_RSYNC) $(L_GRE_VOL)
	@$(convert_target)
	@$(check_target)

$(L_WM_REF_MASK): $(L_RWM_RSYNC) $(L_GRE_VOL)
	@$(convert_target)
	@$(check_target)

# Here's the tricky part: we need the following 'dummy' targets so that 
# rsync gets called _each_ time we wanna build one of the above targets.
# This way we ensure 1) that the files are in sync and 2) that the
# dependent targets (e.g. fast) are just rebuild (e.g. 'fast' is invocated)
# if there is an actual change to the file associated with the above targets
# (e.g. T1W.nii.gz).
.PHONY: $(L_T1W_RSYNC) $(L_GRE_RSYNC) $(L_FEM_RSYNC) \
		$(L_RIM_RSYNC) $(L_RWM_RSYNC)
$(L_T1W_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 1`; \
	 $(copy_target)
	@$(check_target)

$(L_GRE_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 2`; \
	 $(copy_target)
	@$(check_target)

$(L_FLAIR_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 3`; \
	 $(copy_target)
	@$(check_target)

$(L_T2W_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 4`; \
	 $(copy_target)
	@$(check_target)

$(L_FEM_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 5`; \
	 $(copy_target)
	@$(check_target)

$(L_RIM_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 6`; \
	 $(copy_target)
	@$(check_target)

$(L_RWM_RSYNC): $(INPUT_CSV)
	@IN=`cat $< | cut -d ',' -f 7`; \
	 $(copy_target)
	@$(check_target)

# Sleep 1 is important to avoid a denial of service from the ssh server! The
# ssh server has a limit of the number of parallel requests. The sleep 1
# ensures that not too many requests are made when executing this Makefile
# with the -j <#CPUs> option.
define copy_target
rm -f $@; \
echo -n "$(CUR_DIR)/$@: synchronizing $${IN}..."; \
if [ "$${IN}" != "" ]; \
then \
	perl -le 'use Time::HiRes;$$i=sprintf("%0.1f", rand(30)/10+0.1);print($$i);Time::HiRes::sleep($$i)'; \
	sh -c "rsync -azi $${IN} . > $@.tmp"; \
	NAME=`cat $@.tmp | cut -d ' ' -f 2- | $(FSL_EXT_VOL_SED)`; \
	if [ "$${NAME}" != "" ]; \
	then \
		echo "u $${NAME}" > $@; \
	else \
		L_IN=`echo $${IN} | sed 's/^.*\///g'`; \
		if [ "$${L_IN}" != "" ]; \
		then \
			TMP=`ls $${L_IN} 2> /dev/zero`; \
			NAME=`echo $${TMP} | $(FSL_EXT_VOL_SED)`; \
			echo $$NAME; \
			if [ "$${NAME}" != "" ]; \
			then \
				echo "s $${NAME}" > $@; \
			else \
				echo "Fatal error: Cannot find file $${L_IN}!"; \
			fi; \
		else \
			echo "Fatal error: No remote file specified!"; \
		fi; \
	fi; \
	rm -f $@.tmp; \
else \
	echo "Fatal error: No remote path specified!"; \
fi; \
echo "done."
endef

define convert_target
FILE=`echo "$+" | cut -d ' ' -f 1`; \
TEMP=`echo "$+" | cut -d ' ' -f 2`; \
if [ -f $${FILE} ]; \
then \
	MODE=`cat $${FILE} | cut -d ' ' -f 1`; \
	NAME=`cat $${FILE} | cut -d ' ' -f 2`; \
	if [ "$${MODE}" = "u" -o ! -f $@ ]; \
	then \
		echo -n "$(CUR_DIR)/$@: updating..."; \
		ls $$NAME.* &> /dev/zero; \
		if [ $$? -eq 0 ]; \
		then \
			$(SCRIPT_DIR)update_target.sh $(SCRIPT_DIR) "$$NAME" $@; \
		fi; \
		echo "done."; \
	fi; \
fi; \
if [ ! -f $@ -a -f $$TEMP ]; \
then \
	$(FSL_BIN_PFX)fslmaths $$TEMP -mul 0 $@ -odt char; \
fi;
endef

define check_target
MISSING=0; \
if [ `basename $@ .gz` != `basename $@` ]; \
then \
	gunzip -t $@ &> /dev/zero || MISSING=1; \
else \
	test -f $@ || MISSING=1; \
fi; \
if [ $$MISSING -eq 1 ]; \
then \
	echo $(CUR_DIR)/$@ > $@$(MISSING_SFX); \
else \
	rm -f $@$(MISSING_SFX); \
fi;
endef


# Brain extract rules #########################################################
$(L_GRE_BRAIN_BET_MASK): $(L_GRE_VOL)
	@BET_COEF=`$(SCRIPT_DIR)setup_environment.sh BET_GRE`; \
	 $(brain_mask)
	@$(check_target)

$(L_GRE_BRAIN_MASK1): $(L_GRE_BRAIN_BET_MASK) $(L_IC_REF_MASK)
	@VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $(L_IC_REF_MASK)`; \
	 if [ $$? -eq 0 -a $$VOL -gt 0 ]; \
	 then \
		$(FSL_BIN_PFX)fslmaths $(L_IC_REF_MASK) $@; \
	 else \
		$(FSL_BIN_PFX)fslmaths $(L_GRE_BRAIN_BET_MASK) $@; \
	 fi;
	@$(check_target)

#$(L_GRE_BRAIN_BET_VOL): $(L_GRE_VOL) $(L_GRE_BRAIN_BET_MASK)
#	@$(extract_brain)
#	@$(check_target)

$(L_GRE_BRAIN_VOL): $(L_GRE_VOL) $(L_GRE_BRAIN_MASK1)
	@$(extract_brain)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)prepareinput.mt $@ tmp_$@; \
	 mv tmp_$@ $@;
	@$(check_target)

$(L_GRE_BRAIN_MASK): $(L_GRE_BRAIN_VOL)
	@fslmaths $< -bin $@ -odt char;
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_MASK): $(L_T1W_NOTREG_VOL)
	@BET_COEF=`$(SCRIPT_DIR)setup_environment.sh BET_T1W`; \
	 $(brain_mask)
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN_MASK): $(L_FLAIR_NOTREG_VOL)
	@BET_COEF=`$(SCRIPT_DIR)setup_environment.sh BET_FLAIR`; \
	 $(brain_mask)
	@$(check_target)

$(L_T2W_NOTREG_BRAIN_MASK): $(L_T2W_NOTREG_VOL)
	@BET_COEF=`$(SCRIPT_DIR)setup_environment.sh BET_T2W`; \
	 $(brain_mask)
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_VOL): $(L_T1W_NOTREG_VOL) $(L_T1W_NOTREG_BRAIN_MASK)
	@$(extract_brain)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)prepareinput.mt $@ tmp_$@; \
	 mv tmp_$@ $@;
	@$(check_target)	

$(L_T1W_NOTREG_TRUNC_VOL): $(L_T1W_NOTREG_VOL) $(L_T1W_NOTREG_BRAIN_VOL)
	@$(SCRIPT_DIR)exec_matlab.sh	$(SCRIPT_DIR)truncvol.mt $< \
									tmp_$(L_T1W_NOTREG_BRAIN_VOL).txt \
									$@
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN_VOL): $(L_FLAIR_NOTREG_VOL) $(L_FLAIR_NOTREG_BRAIN_MASK)
	@$(extract_brain)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)prepareinput.mt $@ tmp_$@; \
	 mv tmp_$@ $@;
	@$(check_target)

$(L_T2W_NOTREG_BRAIN_VOL): $(L_T2W_NOTREG_VOL) $(L_T2W_NOTREG_BRAIN_MASK)
	@$(extract_brain)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)prepareinput.mt $@ tmp_$@; \
	 mv tmp_$@ $@;
	@$(check_target)

#$(L_T1W_NOTREG_BRAIN_MASK): $(L_GRE_BRAIN_MASK) $(L_T1W_NOTREG_BRAIN1_VOL) \
#							$(L_T1W_NOTREG_BRAIN1_INV_MAT)
#	@$(apply_transform_mask)
#	@$(check_target)

#$(L_FLAIR_NOTREG_BRAIN_MASK): $(L_GRE_BRAIN_MASK) $(L_FLAIR_NOTREG_BRAIN1_VOL) \
#							  $(L_FLAIR_NOTREG_BRAIN1_INV_MAT)
#	@$(apply_transform_mask)
#	@$(check_target)

#$(L_T2W_NOTREG_BRAIN_MASK): $(L_GRE_BRAIN_MASK) $(L_T2W_NOTREG_BRAIN1_VOL) \
#							$(L_T2W_NOTREG_BRAIN1_INV_MAT)
#	@$(apply_transform_mask)
#	@$(check_target)

#$(L_T1W_NOTREG_BRAIN_VOL): $(L_T1W_NOTREG_BRAIN1_VOL) $(L_T1W_NOTREG_BRAIN_MASK)
#	@$(extract_brain)
#	@$(check_target)

#$(L_FLAIR_NOTREG_BRAIN_VOL): $(L_FLAIR_NOTREG_VOL) $(L_FLAIR_NOTREG_BRAIN_MASK)
#	@$(extract_brain)
#	@$(check_target)

#$(L_T2W_NOTREG_BRAIN_VOL): $(L_T2W_NOTREG_VOL) $(L_T2W_NOTREG_BRAIN_MASK)
#	@$(extract_brain)
#	@$(check_target)

$(L_T1W_BRAIN_VOL): $(L_T1W_VOL) $(L_GRE_BRAIN_MASK)
	@$(extract_brain)
	@$(check_target)

$(L_FLAIR_BRAIN_VOL): $(L_FLAIR_VOL) $(L_GRE_BRAIN_MASK)
	@$(extract_brain)
	@$(check_target)

$(L_T2W_BRAIN_VOL): $(L_T2W_VOL) $(L_GRE_BRAIN_MASK)
	@$(extract_brain)
	@$(check_target)

define brain_mask
VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $<`; \
if [ $$? -eq 0 -a $$VOL -gt 0 ]; \
then \
	NAME=`echo $< | $(FSL_EXT_VOL_SED)`; \
	BET_COEF="-m -n $$BET_COEF"; \
	echo -n "$(CUR_DIR)/$@: extracting brain using bet "; \
	echo -n "with $$BET_COEF..."; \
	$(FSL_BIN_PFX)bet $$NAME $${NAME}_tmp $$BET_COEF; \
	mv $${NAME}_tmp_mask$(FSL_EXT_VOL) $@; \
	rm -f $${NAME}_tmp*; \
	echo "done."; \
else \
	echo -n "$(CUR_DIR)/$@: Info: Input is blank. Copying..."; \
	$(FSL_BIN_PFX)fslmaths $< $@; \
	echo "done."; \
fi
endef

define extract_brain
MASK=`echo "$+" | cut -d ' ' -f 2`; \
VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $<`; \
if [ $$? -eq 0 -a $$VOL -gt 0 ]; \
then \
	echo -n "$(CUR_DIR)/$@: applying brain mask..."; \
	$(FSL_BIN_PFX)fslmaths $< -mul $$MASK $@; \
	echo "done."; \
else \
	echo -n "$(CUR_DIR)/$@: Info: Input is blank. Copying..."; \
	$(FSL_BIN_PFX)fslmaths $< $@; \
	echo "done."; \
fi
endef


# Bias field correction rules #################################################
$(L_T1W_NOTREG_BRAIN_RESTORE_VOL):	$(L_T1W_NOTREG_BRAIN_VOL) \
									$(L_T1W_NOTREG_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN_RESTORE_VOL):$(L_FLAIR_NOTREG_BRAIN_VOL) \
									$(L_FLAIR_NOTREG_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_T2W_NOTREG_BRAIN_RESTORE_VOL):$(L_T2W_NOTREG_BRAIN_VOL) \
								  $(L_T2W_NOTREG_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_GRE_BRAIN_RESTORE_VOL): $(L_GRE_BRAIN_VOL) $(L_GRE_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_T1W_BRAIN_RESTORE_VOL): $(L_T1W_BRAIN_VOL) $(L_GRE_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_FLAIR_BRAIN_RESTORE_VOL): $(L_FLAIR_BRAIN_VOL) $(L_GRE_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

$(L_T2W_BRAIN_RESTORE_VOL): $(L_T2W_BRAIN_VOL) $(L_GRE_BRAIN_MASK)
	@$(bias_correct)
	@$(check_target)

define bias_correct
NAME=`echo "$+" | cut -d ' ' -f 1`; \
MASK=`echo "$+" | cut -d ' ' -f 2`; \
LOG=$@.log; \
echo -n "$(CUR_DIR)/$@: removing bias field using N4 ($${LOG})..."; \
LIB_PATH1=$(SLC_PFX)/lib/Slicer3/Plugins/; \
LIB_PATH2=$(SLC_PFX)/lib/InsightToolkit/; \
env LD_LIBRARY_PATH=$${LIB_PATH1}:$${LIB_PATH2} \
$${LIB_PATH1}/N4ITKBiasFieldCorrection \
	--outputimage $@ --outputbiasfield bf_$@ \
	--maskimage $${MASK} --inputimage $${NAME} > $${LOG}; \
echo "done."
endef


# Co-registering rules ########################################################
#$(L_T1W_NOTREG_MAT): $(L_T1W_NOTREG_BRAIN1_VOL) $(L_GRE_BRAIN_VOL)
#	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS1`; \
#	 $(register)
#	@$(check_target)

#$(L_T1W1_VOL):	$(L_T1W_NOTREG_BRAIN1_VOL) $(L_GRE_BRAIN_RESTORE_VOL) \
#				$(L_T1W_NOTREG_BRAIN1_MAT)
#	@FLIRT_COEF="-interp sinc"; \
#	 $(apply_transform)
#	@$(check_target)

$(L_T1W_NOTREG_BRAIN1_MAT):	$(L_T1W_NOTREG_BRAIN1_VOL) $(L_GRE_BRAIN_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS1`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN1_MAT): $(L_T1W_NOTREG_BRAIN_VOL) \
							  $(L_FLAIR_NOTREG_BRAIN1_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS1`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_T2W_NOTREG_BRAIN1_MAT): $(L_T1W_NOTREG_BRAIN_VOL) \
							$(L_T2W_NOTREG_BRAIN1_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS1`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_T1W_NOTREG_BRAIN1_INV_MAT):	$(L_T1W_NOTREG_BRAIN1_MAT)
	@$(invert_matrix)
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN1_INV_MAT): $(L_FLAIR_NOTREG_BRAIN1_MAT) \
								  $(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@$(invert_concat_matrix)
	@$(check_target)

$(L_T2W_NOTREG_BRAIN1_INV_MAT): $(L_T2W_NOTREG_BRAIN1_MAT) \
								$(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@$(invert_concat_matrix)
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_MAT):	$(L_T1W_NOTREG_BRAIN_VOL) $(L_GRE_BRAIN_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS2`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_RESTORE_MAT): $(L_T1W_NOTREG_BRAIN_RESTORE_VOL) \
								   $(L_GRE_BRAIN_RESTORE_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS2`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_RESTORE_TOFLAIR_MAT): $(L_T1W_NOTREG_BRAIN_RESTORE_VOL) \
										   $(L_FLAIR_NOTREG_BRAIN_RESTORE_VOL)
	@FLIRT_COEF=`$(SCRIPT_DIR)setup_environment.sh FLIRT_PASS2`; \
	 REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_FLAIR_NOTREG_BRAIN_RESTORE_MAT): $(L_T1W_NOTREG_BRAIN_RESTORE_MAT) \
									 $(L_T1W_NOTREG_BRAIN_RESTORE_TOFLAIR_MAT)
	@$(invert_concat_matrix)
	@$(check_target)

$(L_T2W_NOTREG_BRAIN_RESTORE_MAT): $(L_T2W_NOTREG_BRAIN_RESTORE_VOL) \
								   $(L_T1W_NOTREG_BRAIN_RESTORE_VOL) \
								   $(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@mv $@ tmp_$@; \
	 convert_xfm -omat $@ -concat $(L_T1W_NOTREG_BRAIN_RESTORE_MAT) tmp_$@
	@$(check_target)

$(L_T1W_VOL): $(L_T1W_NOTREG_BRAIN_VOL) $(L_GRE_BRAIN_RESTORE_VOL) \
			  $(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@FLIRT_COEF="-interp sinc"; \
	 $(apply_transform)
	@$(check_target)

$(L_FLAIR_VOL): $(L_FLAIR_NOTREG_BRAIN_VOL) $(L_GRE_BRAIN_RESTORE_VOL) \
			    $(L_FLAIR_NOTREG_BRAIN_RESTORE_MAT)
	@FLIRT_COEF="-interp sinc"; \
	 $(apply_transform)
	@$(check_target)

$(L_T2W_VOL): $(L_T2W_NOTREG_BRAIN_VOL) $(L_GRE_BRAIN_RESTORE_VOL) \
			  $(L_T2W_NOTREG_BRAIN_RESTORE_MAT)
	@FLIRT_COEF="-interp sinc"; \
	 $(apply_transform)
	@$(check_target)

define invert_matrix
echo -n "$(CUR_DIR)/$@: inverting matrix..."; \
convert_xfm -omat $@ -inverse $<; \
echo "done."
endef

define invert_concat_matrix
BTOC=`echo $+ | cut -d ' ' -f 1`; \
BTOA=`echo $+ | cut -d ' ' -f 2`; \
ATOB=tmp_$@; \
convert_xfm -omat $$ATOB -inverse $$BTOA; \
convert_xfm -omat $@ -concat $$BTOC $$ATOB; \
rm -f $$ATOB;
endef

define apply_transform
IN=`echo "$+" | cut -d ' ' -f 1`; \
REF=`echo "$+" | cut -d ' ' -f 2`; \
MAT=`echo "$+" | cut -d ' ' -f 3`; \
echo -n "$(CUR_DIR)/$@: applying transform using flirt with $$FLIRT_COEF..."; \
$(FSL_BIN_PFX)flirt $$FLIRT_COEF -in $$IN -ref $$REF -init $$MAT -out $@ -applyxfm;\
echo "done."
endef

define apply_transform_mask
IN=`echo "$+" | cut -d ' ' -f 1`; \
IN_VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $$IN`; \
IN_STAT=$$?; \
REF=`echo "$+" | cut -d ' ' -f 2`; \
REF_VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $$REF`; \
REF_STAT=$$?; \
MAT=`echo "$+" | cut -d ' ' -f 3`; \
if [ $$IN_STAT -eq 0 -a $$REF_STAT -eq 0 -a \
     $$IN_VOL -gt 0 -a $$REF_VOL -gt 0 -a -f $$MAT ]; \
then \
	if [ x$${THR} = x ]; \
	then \
		THR=80; \
	fi; \
	echo -n "$(CUR_DIR)/$@: applying transform mask using "; \
	echo -n "flirt with $$FLIRT_COEF..."; \
	$(FSL_BIN_PFX)fslmaths $${IN} -bin -mul 100 tmp_$${IN}; \
	$(FSL_BIN_PFX)flirt $$FLIRT_COEF -in tmp_$${IN} -ref $$REF \
						-init $$MAT -out $@ -applyxfm; \
	$(FSL_BIN_PFX)fslmaths $@ -thr $${THR} -bin $@ -odt char; \
	rm -f tmp_$${IN}; \
	echo "done."; \
else \
	echo -n "$(CUR_DIR)/$@: Info: No transformation needed. Copying..."; \
	$(FSL_BIN_PFX)fslmaths $$IN $@; \
	echo "done."; \
fi
endef

define register
IN=`echo "$+" | cut -d ' ' -f 1`; \
IN_VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $$IN`; \
IN_STAT=$$?; \
REF=`echo "$+" | cut -d ' ' -f 2`; \
REF_VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $$REF`; \
REF_STAT=$$?; \
if [ $$IN_STAT -eq 0 -a $$REF_STAT -eq 0 -a \
	 $$IN_VOL -gt 0 -a $$REF_VOL -gt 0 ]; \
then \
	echo -n "$(CUR_DIR)/$@: registering using flirt with $$FLIRT_COEF, $$REF_COEF..."; \
	$(FSL_BIN_PFX)fslmaths $$IN -bin -kernel sphere 3 -ero -bin \
							tmp_ero_in_mask -odt char; \
	$(FSL_BIN_PFX)fslmaths $$REF $$REF_COEF -bin \
							tmp_ero_ref_mask -odt char; \
	$(FSL_BIN_PFX)flirt $$FLIRT_COEF \
						-in $$IN -ref $$REF \
						-inweight tmp_ero_in_mask -refweight tmp_ero_ref_mask \
						-omat $@; \
	rm -f tmp_ero_in_mask.* tmp_ero_ref_mask.*; \
	echo "done."; \
else \
	echo -n "$(CUR_DIR)/$@: Info: Input is blank. Creating unity matrix..."; \
	$(create_unitymat); \
	echo "done."; \
fi
endef

define create_unitymat
echo "1  0  0  0  " > $@; \
echo "0  1  0  0  " >> $@; \
echo "0  0  1  0  " >> $@; \
echo "0  0  0  1  " >> $@
endef


# Basal ganglia mask rules ####################################################
$(L_T1W_NOTREG_NOISE_MASK): $(L_T1W_NOTREG_VOL)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)estnoise.mt $< 0.2 $@
	@$(check_target)

$(L_T1W_NOTREG_BRAIN_RESTORE_SUSAN_VOL): $(L_T1W_NOTREG_BRAIN_RESTORE_VOL) \
										 $(L_T1W_NOTREG_NOISE_MASK)
	@N=`cat $(L_T1W_NOTREG_NOISE_MASK).txt`; \
	 $(FSL_BIN_PFX)susan $< $$((N*4)) 3 3 0 0 $@
	@$(check_target)

#	 $(FSL_BIN_PFX)fslmaths $< -bin -kernel sphere 3 -ero -bin tmp_$< -odt char; \
#	 $(FSL_BIN_PFX)first_flirt $< $$OUT -b -inweight tmp_$< -strucweight $(SCRIPT_DIR)MNI_BG_mask; \

$(L_BG_NOTREG_MAT): $(L_T1W_NOTREG_BRAIN_RESTORE_SUSAN_VOL)
	@echo -n "$(CUR_DIR)/$@: first_flirt ..."; \
	 OUT=`echo $@ | $(FSL_EXT_VOL_SED)`; \
	 $(FSL_BIN_PFX)first_flirt $< $$OUT -b; \
	 mv $${OUT}.mat $@; rm -f tmp_$<; \
	 echo "done.";
	@$(check_target)

$(L_BG_NOTREG_MASK): $(L_T1W_NOTREG_BRAIN_RESTORE_VOL) $(L_BG_NOTREG_MAT)
	@$(run_first)
	@$(check_target)

$(L_BG_MASK): $(L_BG_NOTREG_MASK) $(L_GRE_BRAIN_VOL) $(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@FLIRT_COEF="-interp nearestneighbour"; \
	 $(apply_transform)
	@$(check_target)

define run_first
IN=`echo $< | $(FSL_EXT_VOL_SED)`; \
MAT=`echo $+ | cut -d ' ' -f 2`; \
OUT=$${IN}_first; \
FIRST_COEF="-s L_Thal,L_Caud,L_Puta,L_Pall,R_Thal,R_Caud,R_Puta,R_Pall"; \
echo -n "$(CUR_DIR)/$@: fitting model using first with $$FIRST_COEF..."; \
$(FSL_BIN_PFX)run_first_all -a $$MAT -b $$FIRST_COEF -m none -i $${IN} -o $$OUT; \
mv $${OUT}_all_none_firstseg$(FSL_EXT_VOL) $@; \
rm -rf $${OUT}* tmp_$${IN}; \
echo "done."
endef


# Normal-appearing white matter mask rules ####################################
$(L_SE_MASK): $(L_T1W_BRAIN_RESTORE_VOL) $(L_GRE_BRAIN_RESTORE_VOL)
	@echo -n "$(CUR_DIR)/$@: segmenting tissue..."; \
	 T1W=`echo "$+" | cut -d ' ' -f 1`; T1W=`echo $$T1W | $(FSL_EXT_VOL_SED)`; \
	 GRE=`echo "$+" | cut -d ' ' -f 2`; GRE=`echo $$GRE | $(FSL_EXT_VOL_SED)`; \
	 OUT=T1W_GRE_brain_restore; \
	 $(FSL_BIN_PFX)fast -n 4 -S 2 -N -I 0 -o $$OUT $$T1W $$GRE; \
	 mv $${OUT}_mixeltype$(FSL_EXT_VOL) $@; \
	 mv $${OUT}_seg$(FSL_EXT_VOL) tmp_$@; \
	 echo "done.";
	@$(check_target)

#$(L_SE_LORD): $(L_GRE_BRAIN_RESTORE_VOL) $(L_SE_MASK)
#	@echo -n "$(CUR_DIR)/$@: determine label order..."; \
#	 GRE=`echo "$+" | cut -d ' ' -f 1`; \
#	 MASK=`echo "$+" | cut -d ' ' -f 2`; \
#	 $(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)order_labels.mt $$GRE $$MASK $@; \
#	 echo "done."
#	@$(check_target)

$(L_SE_HARD_MASK): $(L_SE_MASK)
	@cp tmp_$< $@;
	@$(check_target)

$(L_SE_HARD_LORD): $(L_T1W_BRAIN_RESTORE_VOL) $(L_SE_HARD_MASK)
	@echo -n "$(CUR_DIR)/$@: determine label order..."; \
	 GRE=`echo "$+" | cut -d ' ' -f 1`; \
	 MASK=`echo "$+" | cut -d ' ' -f 2`; \
	 $(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)order_labels.mt $$GRE $$MASK $@; \
	 echo "done."
	@$(check_target)

$(L_WM_MASK): $(L_SE_HARD_MASK) $(L_WM_REF_MASK) $(L_SE_HARD_LORD)
	@VOL=`$(SCRIPT_DIR)check_empty.sh $(FSL_BIN_PFX) $(L_WM_REF_MASK)`; \
	 if [ $$? -eq 0 -a $$VOL -gt 0 ]; \
	 then \
		$(FSL_BIN_PFX)fslmaths $(L_WM_REF_MASK) $@; \
	 else \
		IDX=`$(SCRIPT_DIR)setup_environment.sh WM_IDX`; \
		TISS=`cat $(L_SE_HARD_LORD) | cut -d ' ' -f $$IDX`; \
		$(thr_mask) \
	 fi;
	@$(check_target)

$(L_GM_MASK): $(L_SE_HARD_MASK) $(L_SE_HARD_LORD)
	@IDX=`$(SCRIPT_DIR)setup_environment.sh GM_IDX`; \
	 TISS=`cat $(L_SE_HARD_LORD) | cut -d ' ' -f $$IDX`; \
	 $(thr_mask)
	@$(check_target)

$(L_CS_MASK): $(L_SE_HARD_MASK) $(L_SE_HARD_LORD)
	@IDX=`$(SCRIPT_DIR)setup_environment.sh CS_IDX`; \
	 TISS=`cat $(L_SE_HARD_LORD) | cut -d ' ' -f $$IDX`; \
	 $(thr_mask)
	@$(check_target)

$(L_AR_MASK): $(L_SE_HARD_MASK) $(L_SE_HARD_LORD)
	@IDX=`$(SCRIPT_DIR)setup_environment.sh AR_IDX`; \
	 TISS=`cat $(L_SE_HARD_LORD) | cut -d ' ' -f $$IDX`; \
	 $(thr_mask)
	@$(check_target)

tmp_$(L_AR_ROI_MASK): $(L_AR_MASK) $(L_RO_MASK) $(L_GRE_BRAIN_RESTORE_VOL)
	@DIM=2; VOLFAC=0.5; $(mask2maskroi)
	@$(check_target)

$(L_AR_ROI_MASK): tmp_$(L_AR_ROI_MASK) $(L_RO_MASK) $(L_GRE_BRAIN_RESTORE_VOL)
	@DIM=3; VOLFAC=0.5; $(mask2maskroi)
	@$(check_target)

define thr_mask
echo -n "$(CUR_DIR)/$@: generating tissue mask $@..."; \
$(FSL_BIN_PFX)fslmaths $< -mul 0 $@; \
for tiss in $$TISS; \
do \
	$(FSL_BIN_PFX)fslmaths $< -thr $$tiss -uthr $$tiss \
							  -bin -add $@ -bin $@ -odt char; \
done; \
echo "done.";
endef


# Iron mask rules ##############################################################
$(L_RO_MASK): $(L_BG_MASK)
	@$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)bg2roi.mt $< $@
	@$(check_target)

$(L_FE_ROI_MASK): $(L_FE_REF_MASK) $(L_RO_MASK) $(L_GRE_BRAIN_RESTORE_VOL)
	@DIM=3; VOLFAC=0.5; $(mask2maskroi)
	@$(check_target)

define mask2maskroi
MASK=`echo "$+" | cut -d ' ' -f 1`; \
ROI=`echo "$+" | cut -d ' ' -f 2`; \
WEIGHT=`echo "$+" | cut -d ' ' -f 3`; \
echo -n "$(CUR_DIR)/$@: generating reduced mask $@..."; \
$(SCRIPT_DIR)exec_matlab.sh $(SCRIPT_DIR)mask2maskroi.mt $$MASK $$DIM $$ROI $$WEIGHT $$VOLFAC $@; \
echo "done.";
endef

$(L_T1W_NOTREG_BRAIN_RESTORE_INV_MAT):	$(L_T1W_NOTREG_BRAIN_RESTORE_MAT)
	@$(invert_matrix)
	@$(check_target)

$(L_FE_ROI_BIN_MASK): $(L_FE_ROI_MASK)
	@$(FSL_BIN_PFX)fslmaths $< -bin $@ -odt char
	@$(check_target)

$(L_FE_ROI_BIN_NOTREG_MASK):	$(L_FE_ROI_BIN_MASK) $(L_T1W_NOTREG_BRAIN_VOL) \
								$(L_T1W_NOTREG_BRAIN_RESTORE_INV_MAT)
	@THR=30;$(apply_transform_mask)
	@$(check_target)

$(L_FE_REF_NOTREG_MASK): $(L_FE_REF_MASK) $(L_T1W_NOTREG_BRAIN_VOL) \
						 $(L_T1W_NOTREG_BRAIN_RESTORE_INV_MAT)
	@THR=30;$(apply_transform_mask)
	@$(check_target)


# Nonlinear registration rules ################################################
$(L_T1W_NOTREG_TRUNC_MNI_MAT):	$(L_T1W_NOTREG_BRAIN_RESTORE_VOL) \
								$(L_T1W_BRAIN_MNI_REF_VOL)
	@REF_COEF="-bin -kernel sphere 3 -ero"; \
	 $(register)
	@$(check_target)

$(L_T1W_NOTREG_TRUNC_MNI_WARP): $(L_T1W_NOTREG_TRUNC_VOL) \
								$(L_T1W_NOTREG_TRUNC_MNI_MAT) \
								$(L_T1W_MNI_REF_FNIRT_CFG)
	@$(nregister)
	@$(check_target)

define nregister
IN=`echo "$+" | cut -d ' ' -f 1 | $(FSL_EXT_VOL_SED)`; \
AFF=`echo "$+" | cut -d ' ' -f 2`; \
CFG=`echo "$+" | cut -d ' ' -f 3`; \
echo -n "$(CUR_DIR)/$@: registering using fnirt with ..."; \
$(FSL_BIN_PFX)fnirt --in=$$IN --aff=$$AFF --cout=$@ --config=$$CFG; \
echo "done."
endef

$(L_T1W_NOTREG_TRUNC_MNI_VOL):	$(L_T1W_NOTREG_TRUNC_VOL) \
								$(L_T1W_BRAIN_MNI_REF_VOL) \
								$(L_T1W_NOTREG_TRUNC_MNI_WARP)
	@$(apply_warp)
	@$(check_target)

$(L_FE_ROI_BIN_NOTREG_MNI_MASK):$(L_FE_ROI_BIN_NOTREG_MASK) \
								$(L_T1W_BRAIN_MNI_REF_VOL) \
								$(L_T1W_NOTREG_TRUNC_MNI_WARP)
	@$(apply_warp_mask)
	@$(check_target)

define apply_warp
IN=`echo "$+" | cut -d ' ' -f 1`; \
REF=`echo "$+" | cut -d ' ' -f 2`; \
VOL=`echo "$+" | cut -d ' ' -f 3`; \
echo -n "$(CUR_DIR)/$@: applying warp using applywarp ..."; \
$(FSL_BIN_PFX)applywarp --ref=$$REF --in=$$IN --warp=$$VOL --out=$@; \
echo "done."
endef

define apply_warp_mask
IN=`echo "$+" | cut -d ' ' -f 1`; \
REF=`echo "$+" | cut -d ' ' -f 2`; \
VOL=`echo "$+" | cut -d ' ' -f 3`; \
echo -n "$(CUR_DIR)/$@: applying warp to mask using applywarp ..."; \
$(FSL_BIN_PFX)fslmaths $$IN -bin -mul 100 tmp_$${IN}; \
$(FSL_BIN_PFX)applywarp --ref=$$REF --in=tmp_$${IN} --warp=$$VOL --out=$@; \
rm -f tmp_$${IN}; \
$(FSL_BIN_PFX)fslmaths $@ -thr 60 -bin $@ -odt char; \
echo "done."
endef


# Missing targets rules #######################################################
.PHONY: clean-mtg print-mtg
clean-mtg:
	@rm -f *$(MISSING_SFX)

print-mtg:
	@MISSING_LIST=`ls *$(MISSING_SFX) 2> /dev/zero`; \
	 if [ $$? -eq 0 ]; \
	 then \
		cat $$MISSING_LIST; \
	 fi;

