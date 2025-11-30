# ------------------- Required for MSVC nmake ---------------------------------
# This file should be included at the top of a MAKEFILE as follows:


CPU = AMD64

MODEL  = LQR_min_obs
TARGET = cgxe
MODULE_SRCS 	= m_ez6b5k7jBQ3Ohnvv6n8cu.c
MODEL_SRC	= LQR_min_obs_cgxe.c
MODEL_REG = LQR_min_obs_cgxe_registry.c
MAKEFILE    = LQR_min_obs_cgxe.mak
MATLAB_ROOT	= D:\Program Files
BUILDARGS   =

#--------------------------- Tool Specifications ------------------------------
#
#
MSVC_ROOT1 = $(MSDEVDIR:SharedIDE=vc)
MSVC_ROOT2 = $(MSVC_ROOT1:SHAREDIDE=vc)
MSVC_ROOT  = $(MSVC_ROOT2:sharedide=vc)

# Compiler tool locations, CC, LD, LIBCMD:
CC     = cl.exe
LD     = link.exe
LIBCMD = lib.exe
#------------------------------ Include/Lib Path ------------------------------


USER_INCLUDES   =  /I "F:\IITM\sem 5\ME4010\Ball_bot" /I "F:\IITM\sem 5\ME4010\Ball_bot\slprj\_cprj"

MLSL_INCLUDES     = \
    /I "D:\Program Files\extern\include" \
    /I "D:\Program Files\simulink\include" \
    /I "D:\Program Files\rtw\c\src"
COMPILER_INCLUDES = /I "$(MSVC_ROOT)\include"

THIRD_PARTY_INCLUDES   =  /I "F:\IITM\sem 5\ME4010\Ball_bot\slprj\_cgxe\LQR_min_obs\src"
INCLUDE_PATH = $(MLSL_INCLUDES) $(USER_INCLUDES) $(THIRD_PARTY_INCLUDES)
LIB_PATH     = "$(MSVC_ROOT)\lib"
CFLAGS = /c /Zp8 /GR /w /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMX_COMPAT_64 /DMATLAB_MEXCMD_RELEASE=R2018a /DMATLAB_MEX_FILE /nologo /MD   
LDFLAGS = /nologo /dll /MANIFEST /OPT:NOREF /export:mexFunction /export:mexfilerequiredapiversion  
#----------------------------- Source Files -----------------------------------

USER_OBJS =

AUX_SRCS = D:\Program Files\extern\version\c_mexapi_version.c  

REQ_SRCS  = $(MODEL_SRC) $(MODEL_REG) $(MODULE_SRCS) $(AUX_SRCS)
REQ_OBJS = $(REQ_SRCS:.cpp=.obj)
REQ_OBJS2 = $(REQ_OBJS:.c=.obj)
OBJS = $(REQ_OBJS2) $(USER_OBJS) $(AUX_ABS_OBJS)
OBJLIST_FILE = LQR_min_obs_cgxe.mol
TMWLIB = "D:\Program Files\extern\lib\win64\microsoft\libmx.lib" "D:\Program Files\extern\lib\win64\microsoft\libmex.lib" "D:\Program Files\extern\lib\win64\microsoft\libmat.lib" "D:\Program Files\extern\lib\win64\microsoft\libfixedpoint.lib" "D:\Program Files\extern\lib\win64\microsoft\libut.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwmathutil.lib" "D:\Program Files\extern\lib\win64\microsoft\libemlrt.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwcgxert.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwcgxeooprt.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwslexec_simbridge.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwslccrt.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwstringutil.lib" "D:\Program Files\extern\lib\win64\microsoft\libcovrt.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwsl_sfcn_cov_bridge.lib" "D:\Program Files\extern\lib\win64\microsoft\libmwdsp_halidesim.lib" 
THIRD_PARTY_LIBS = 

#--------------------------------- Rules --------------------------------------

MEX_FILE_NAME_WO_EXT = $(MODEL)_$(TARGET)
MEX_FILE_NAME = $(MEX_FILE_NAME_WO_EXT).mexw64
all : $(MEX_FILE_NAME) 


$(MEX_FILE_NAME) : $(MAKEFILE) $(OBJS)
	@echo ### Linking ...
	$(LD) $(LDFLAGS) /OUT:$(MEX_FILE_NAME) /map:"$(MEX_FILE_NAME_WO_EXT).map" $(TMWLIB) $(THIRD_PARTY_LIBS) @$(OBJLIST_FILE)
	@echo ### Created $@

.c.obj :
	@echo ### Compiling "$<"
	$(CC) $(CFLAGS) $(INCLUDE_PATH) "$<"

.cpp.obj :
	@echo ### Compiling "$<"
	$(CC) $(CFLAGS) $(INCLUDE_PATH) "$<"

