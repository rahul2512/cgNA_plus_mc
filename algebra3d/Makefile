# See the Makefile.config for the required settings
# (like the choice of compilers and compilation flags)
include Makefile.config

# Project directories
INC_DIR         = include
SRC_DIR         = src
BUILD_DIR       = build
BIN_DIR         = bin
LIB_DIR         = lib
TEST_DIR        = test

# Name of the library
LIB_NAME        = algebra3d

# All required headers 
# (additional can be set in Makefile.config)
INC             = -I$(INC_DIR) $(ADDITIONAL_INC)
INC_TEST        = $(ADDITIONAL_INC_TEST)

# All required linking flags
# (additional can be set in Makefile.config)
LDFLAGS         += -lm $(ADDITIONAL_LDFLAGS)
LDFLAGS_TEST    += $(ADDITIONAL_LDFLAGS_TEST)

vpath   %.cpp   $(SRC_DIR)
vpath   %.h     $(INC_DIR)/$(LIB_NAME)

###############################################################################

OBJS            = $(BUILD_DIR)/global.o \
                  $(BUILD_DIR)/Matrix3.o $(BUILD_DIR)/Vector3.o \
                  $(BUILD_DIR)/Quaternion.o \

TARGETS         = $(OBJS) \
                  $(LIB_DIR)/lib$(LIB_NAME).so \
                  $(LIB_DIR)/lib$(LIB_NAME).a
                  
TARGETS_TIMING  = $(BIN_DIR)/timingD \
                  $(BIN_DIR)/timingF

all: $(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) $(TARGETS)

timing: $(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) $(TARGETS) $(TARGETS_TIMING) 

###############################################################################
# Only for tests - unused in the default make rule

OBJS_TEST_D     = $(BUILD_DIR)/GlobalTestD.o \
                  $(BUILD_DIR)/Matrix3TestD.o \
                  $(BUILD_DIR)/Vector3TestD.o \
                  $(BUILD_DIR)/QuaternionTestD.o \
                  $(BUILD_DIR)/RotationTestD.o \

OBJS_TEST_F     = $(BUILD_DIR)/GlobalTestF.o \
                  $(BUILD_DIR)/Matrix3TestF.o \
                  $(BUILD_DIR)/Vector3TestF.o \
                  $(BUILD_DIR)/QuaternionTestF.o \
                  $(BUILD_DIR)/RotationTestF.o \

TARGETS_TEST    = $(BIN_DIR)/TesterD \
                  $(BIN_DIR)/TesterF \
                  $(BIN_DIR)/TesterStaticD \
                  $(BIN_DIR)/TesterStaticF

all_test: $(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) $(TARGETS) $(TARGETS_TEST)

###############################################################################

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

$(BIN_DIR):
	mkdir $(BIN_DIR)

$(LIB_DIR):
	mkdir $(LIB_DIR)

doc: $(INC_DIR)/$(LIB_NAME)/*
	doxygen

$(BUILD_DIR)/%.o:	%.cpp %.h global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ \
	    $<

$(LIB_DIR)/lib$(LIB_NAME).so:	$(OBJS)
	$(CXX) \
	    -shared \
	    -o "$(PWD)/$@" \
	    $^

$(LIB_DIR)/lib$(LIB_NAME).a:	$(OBJS)
	ar -cvq $@ $(OBJS)

$(BIN_DIR)/%D:   %.cpp $(LIB_DIR)/lib$(LIB_NAME).a
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -o $@ \
	    $< \
	    $(LIB_DIR)/lib$(LIB_NAME).a \
	    $(LDFLAGS)

$(BIN_DIR)/%F:   %.cpp $(LIB_DIR)/lib$(LIB_NAME).a
	$(CXX) \
	    $(CXXFLAGS) \
	    -D FPTYPE_FLOAT \
	    $(INC) \
	    -o $@ \
	    $< \
	    $(LIB_DIR)/lib$(LIB_NAME).a \
	    $(LDFLAGS)

###############################################################################
# Only for tests - unused in the default make rule
$(BIN_DIR)/TesterD: $(TEST_DIR)/Tester.cpp $(LIB_DIR)/lib$(LIB_NAME).so \
                    $(OBJS_TEST_D) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -o $@ \
		$< \
	    $(OBJS_TEST_D) \
	    $(LDFLAGS) -lcppunit -ldl \
	    -L$(LIB_DIR) -l$(LIB_NAME)

$(BIN_DIR)/TesterF: $(TEST_DIR)/Tester.cpp $(LIB_DIR)/lib$(LIB_NAME).so \
                    $(OBJS_TEST_F) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    -D FPTYPE_FLOAT \
	    $(INC) \
	    -o $@ \
		$< \
	    $(OBJS_TEST_F) \
	    $(LDFLAGS) $(LDFLAGS_TEST) \
	    -L$(LIB_DIR) -l$(LIB_NAME)

$(BIN_DIR)/TesterStaticD:   $(TEST_DIR)/Tester.cpp $(LIB_DIR)/lib$(LIB_NAME).a \
                            $(OBJS_TEST_D) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -o $@ \
		$< \
	    $(OBJS_TEST_D) \
		$(LIB_DIR)/lib$(LIB_NAME).a \
	    $(LDFLAGS) $(LDFLAGS_TEST)

$(BIN_DIR)/TesterStaticF:   $(TEST_DIR)/Tester.cpp $(LIB_DIR)/lib$(LIB_NAME).a \
                            $(OBJS_TEST_F) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    -D FPTYPE_FLOAT \
	    $(INC) \
	    -o $@ \
		$< \
	    $(OBJS_TEST_F) \
		$(LIB_DIR)/lib$(LIB_NAME).a \
	    $(LDFLAGS) $(LDFLAGS_TEST)

$(BUILD_DIR)/%D.o:	$(TEST_DIR)/%.cpp $(OBJS) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ \
	    $<

$(BUILD_DIR)/%F.o:	$(TEST_DIR)/%.cpp $(OBJS) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    -D FPTYPE_FLOAT \
	    $(INC) \
	    -c -o $@ \
	    $<

$(BUILD_DIR)/%.o:	$(TEST_DIR)/%.cpp $(OBJS) global.h
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ \
	    $<
###############################################################################

clean:
	@rm -fr $(BUILD_DIR)/* $(BIN_DIR)/* $(LIB_DIR)/*
	@echo
	@echo "Cleaned " `date`

.PHONY: all timing all_test clean
