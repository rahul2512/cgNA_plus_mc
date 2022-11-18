# See the Makefile.config for the required settings
# (like the choice of compilers and compilation flags)
include Makefile.config

# All platform-independent/known required headers
# (additional can be set in Makefile.config)
INC      = -I. $(ADDITIONAL_INC)

# All platform-independent/known required linking flags
# (additional can be set in Makefile.config)
LDFLAGS += $(ADDITIONAL_LDFLAGS)

OBJS     = rand.o cgDNAutils.o helpers.o

BINARIES = build_seq_data  run_cgDNAmc_from_seq
 
all: ${OBJS} $(BINARIES)

%: %.cpp $(OBJS)
	$(CXX) \
	    $(CXXFLAGS) \
	    $(INC) \
	    -o $@ \
	    $< \
	    ${OBJS} \
	    $(LDFLAGS)

%.o: %.c %.h
	${CXX} \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ $<

%.o: %.c
	${CXX} \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ $<

%.o: %.cpp %.h
	${CXX} \
	    $(CXXFLAGS) \
	    $(INC) \
	    -c -o $@ $<

clean:
	rm -f *.o $(BINARIES)

.PHONY: all, clean
