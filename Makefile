# Set the Vrui installation directory:
VRUIDIR = /opt/local/Vrui-1.0
# Specify the main target
TARGET = VRProteinShop
# Default build type
TYPE = debug
# Which directories contain source files
DIRS = source/GL source/Protein source
# Which libraries are linked
LIBS = GLU
# Dynamic libraries
DLIBS = -lglu
# Frameworks for MAC
FRAMEWORKS = 

# The next blocks change some variables depending on the build type
ifeq ($(TYPE),debug)
  # Build debug version of the applications, using the debug version of Vrui:
  include $(VRUIDIR)/etc/Vrui.debug.makeinclude
  CFLAGS += -g2 -O0 -DCLIENT_SIDE -D__LITTLE_ENDIAN=1234 -D__BIG_ENDIAN=4321 -DVRPROTEINSHOP_CONFIGFILENAME='"ProteinShop.cfg"' -DVRPROTEINSHOP_STANDARDSDIR='"Standards"'
  LFLAGS += 
endif

ifeq ($(TYPE), release)
  # Build release version of the applications, using the release version of Vrui:
  include $(VRUIDIR)/etc/Vrui.makeinclude
  CFLAGS += -g0 -O3 -DCLIENT_SIDE -D__LITTLE_ENDIAN=1234 -D__BIG_ENDIAN=4321 -DVRPROTEINSHOP_CONFIGFILENAME='"ProteinShop.cfg"' -DVRPROTEINSHOP_STANDARDSDIR='"Standards"'
  LFLAGS += 
endif

# Add directories to the include and library paths
INCPATH = $(DIRS)
LIBPATH = 

# Which files to add to backups, apart from the source code
EXTRA_FILES = Makefile
# The compiler
C++ = g++

# Store TARGET file.
EXECDIR = bin
# Store source files.
SRCDIR = source
# Store object (.o) and dependancy (.d) files.
OBJDIR = objects
# List of source (.cpp) files.
SOURCE := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cpp))
# List of header (.h) files.
HEADERS := $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
# List of object (.o) files.
OBJECTS := $(addprefix $(OBJDIR)/, $(SOURCE:.cpp=.o))
# List of dependancy (.d) files.
DFILES := $(addprefix $(OBJDIR)/,$(SOURCE:.cpp=.d))

# Specify phony rules. These are rules that are not real files.
.PHONY: clean backup dirs all

ALL = $(TARGET)

all: $(ALL)

# Main target. The @ in front of a command prevents make from displaying
# it to the standard output.
VRProteinShop: dirs $(OBJECTS)
		@echo Linking $(EXECDIR)/$(TARGET).
		@$(C++) -o $(EXECDIR)/$(TARGET) $(OBJECTS) $(VRUI_LINKFLAGS) $(LFLAGS) $(foreach LIBRARY, \
			$(LIBS),-l$(LIBRARY)) $(foreach LIB,$(LIBPATH),-L$(LIB)) $(foreach FRAMEWORK,$(FRAMEWORKS),-framework $(FRAMEWORK))

# Rule for creating object file and .d file, the sed magic is to add
# the object path at the start of the file because the files gcc
# outputs assume it will be in the same dir as the source file.
$(OBJDIR)/%.o: %.cpp
		@echo Creating object file for $*.
		@$(C++) -Wp,-MMD,$(OBJDIR)/$*.dd $(VRUI_CFLAGS) $(CFLAGS) $(foreach INC,$(INCPATH),-I$(INC))\
				$(foreach MACRO,$(MACROS),-D$(MACRO)) -c $< -o $@
		@sed -e '1s/^\(.*\)$$/$(subst /,\/,$(dir $@))\1/' $(OBJDIR)/$*.dd > $(OBJDIR)/$*.d
		@rm -f $(OBJDIR)/$*.dd

# Empty rule to prevent problems when a header is deleted.
%.hpp: ;

# Cleans up the objects, .d files and executables.
clean:
		@echo Making clean.
		@-rm -f $(foreach DIR,$(DIRS),$(OBJDIR)/$(DIR)/*.d $(OBJDIR)/$(DIR)/*.o)
		@-rm -f $(EXECDIR)/$(TARGET)

# Backup the source files.
backup:
		@-if [ ! -e .backup ]; then mkdir .backup; fi;
		@tar cvf .backup/backup_`date +%d-%m-%y_%H.%M`.tar $(SOURCE) $(HEADERS) $(EXTRA_FILES)

# Create necessary directories
dirs:
		@-if [ ! -e $(OBJDIR) ]; then mkdir $(OBJDIR); fi;
		@-if [ ! -e $(EXECDIR) ]; then mkdir $(EXECDIR); fi;
		@-$(foreach DIR,$(DIRS), if [ ! -e $(OBJDIR)/$(DIR) ]; \
		then mkdir $(OBJDIR)/$(DIR); fi; )

# Includes the .d files so it knows the exact dependencies for every
# source.
-include $(DFILES)
