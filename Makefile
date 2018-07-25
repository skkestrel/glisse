.PHONY: all clean

NVCC = /usr/local/cuda-9.2/bin/nvcc
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj
TARGETS_DIR = targets
DOCOPT_DIR = docopt
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES)) $(OBJ_DIR)/docopt/docopt.o
DEP_FILES = $(OBJ_FILES:%.o=%.d)

WFLAGS = -Wcast-align -Wshadow -Wcast-qual -Wconversion -Wdisabled-optimization \
-Wfloat-equal -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k \
-Wimport  -Winit-self -Winvalid-pch -Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn -Wpointer-arith -Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused \
-Wunused-parameter -Wvariadic-macros -Wwrite-strings

CPPFLAGS = ${WFLAGS} -g --std=c++11 -Wall -Wextra -Wpedantic ${WFLAGS} -O3 -DNO_CUDA # -fsanitize=address

LDFLAGS = -pthread # -lasan

sr:
	@mkdir -p $(BIN_DIR)
	@$(NVCC) $(TARGETS_DIR)/main.cpp $(DOCOPT_DIR)/docopt.cpp $(SRC_FILES) $(SRC_DIR)/*.cu -lineinfo -g -maxrregcount 64 -arch=sm_35 --std=c++11 --compiler-options "-Wall -Wextra ${WFLAGS} -fstack-protector" -o $(BIN_DIR)/sr -O0

clean:
	rm -r $(OBJ_DIR)/* 
	rm -r $(BIN_DIR)/*

$(OBJ_DIR)/docopt/docopt.o: $(DOCOPT_DIR)/docopt.cpp
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/docopt
	$(info CPP $(OBJ_DIR)/docopt/docopt.o)
	@g++ $(DOCOPT_DIR)/docopt.cpp $(CPPFLAGS) -c -o $(OBJ_DIR)/docopt/docopt.o

-include $(DEP_FILES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p obj
	$(info CPP $@)
	@g++ $(CPPFLAGS) -MMD -c -o $@ $<

define make-target =
@mkdir -p $(BIN_DIR)
@mkdir -p $(OBJ_DIR)/targets
$(info CPP $(OBJ_DIR)/targets/$1.o)
@g++ $(TARGETS_DIR)/$1.cpp $(CPPFLAGS) -c -o $(OBJ_DIR)/targets/$1.o
$(info LD $(BIN_DIR)/$2)
@g++ $(OBJ_DIR)/targets/$1.o $(OBJ_FILES) $(LDFLAGS) -o $(BIN_DIR)/$2
endef

cpu: $(OBJ_FILES)
	$(call make-target,main,sr_cpu)

convert-state: $(OBJ_FILES)
	$(call make-target,convert_state,convert-state)

filter-state: $(OBJ_FILES)
	$(call make-target,filter_state,filter-state)

prune-track: $(OBJ_FILES)
	$(call make-target,prune_track,prune-track)

find-max-e: $(OBJ_FILES)
	$(call make-target,find_max_e,find-max-e)
