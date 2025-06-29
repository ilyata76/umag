
# ======================================================================
#  UMAG / Makefile  —  Poetry-driven build & dev helpers
# ======================================================================

.PHONY:   deps venv install update build-debug build-release \
          test lint format clean clean-all shell run

.DEFAULT_GOAL := build-all

# =====================================
#  Paths & execs
# =====================================

PROJECT_ROOT := $(realpath $(dir $(firstword $(MAKEFILE_LIST))))

POETRY := $(shell which poetry)

ifndef POETRY
	$(error Poetry executable not found in PATH)
endif

PYTHON        := $(POETRY) run python
PIP           := $(POETRY) run pip
BLACK         := $(POETRY) run black
CLANG_FORMAT  := $(shell which clang-format)

# =====================================
#  Builds
# =====================================

build-all: # e.g. : make build-spindynapy FLAGS=-v
	ASAN_OPTIONS=detect_leaks=0 DEBUG=0 $(POETRY) install $(FLAGS)

build-spindynapy:  # e.g. : make build-spindynapy FLAGS=-v
	ASAN_OPTIONS=detect_leaks=0 DEBUG=1 $(PIP) install $(PROJECT_ROOT)/spindynapy $(FLAGS)

build-release:  # Релизная сборка с DEBUG=FALSE
	DEBUG=0 ASAN_OPTIONS=detect_leaks=0 $(PIP) install $(PROJECT_ROOT)/spindynapy $(FLAGS)

venv:
	$(shell $(POETRY) env activate)

format:
ifdef CLANG_FORMAT
	find $(PROJECT_ROOT)/spindynapy/core_src -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec $(CLANG_FORMAT) -style=file -i {} \;
endif
	$(BLACK) $(shell find $(PROJECT_ROOT) -maxdepth 1 -name '*.py') \
	                $(shell find $(PROJECT_ROOT)/spindynapy -maxdepth 1 -name '*.py') \
	                $(shell find $(PROJECT_ROOT)/examples -maxdepth 1 -name '*.py');

# =====================================
#  Cleaning
# =====================================

clean-spindynapy:
	rm -rf $(PROJECT_ROOT)/spindynapy/lib/*
	rm -rf $(PROJECT_ROOT)/spindynapy/*.so
	rm -rf $(PROJECT_ROOT)/spindynapy/core/*
	rm -rf $(PROJECT_ROOT)/spindynapy/build/*
	rm -rf $(PROJECT_ROOT)/spindynapy/__pycache__
	rm -rf $(PROJECT_ROOT)/spindynapy/spindynapy.egg-info
	rm -rf $(PROJECT_ROOT)/spindynapy/compile_commands.json

clean-all: clean-spindynapy
	rm -rf $(PROJECT_ROOT)/spindynapy/build-deps/*
	rm -rf $(PROJECT_ROOT)/.pytest_cache/
	rm -rf $(PROJECT_ROOT)/.mypy_cache
	rm -rf $(PROJECT_ROOT)/__pycache__
	rm -rf $(PROJECT_ROOT)/.cache

# =====================================
#  Run main script
# =====================================

run:
	ASAN_OPTIONS=detect_leaks=0 $(PYTHON) $(PROJECT_ROOT)/main.py
