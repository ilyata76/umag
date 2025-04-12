
# Пути
PROJECT_ROOT := $(PWD)

# Бинарные файлы
CMAKE := cmake
POETRY := poetry
PIP = $(POETRY) run pip
CLANG_FORMAT = clang-format

build-all:
	$(POETRY) install

build-spindynapy:  # e.g. : make build-spindynapy FLAGS=-v
	$(PIP) install $(PROJECT_ROOT)/spindynapy $(FLAGS)

format:
	find ./spindynapy/core_src -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec $(CLANG_FORMAT) -style=file -i {} \;

clean-spindynapy:
	rm -rf $(PROJECT_ROOT)/spindynapy/lib/*
	rm -rf $(PROJECT_ROOT)/spindynapy/*.so
	rm -rf $(PROJECT_ROOT)/spindynapy/core/*
	rm -rf $(PROJECT_ROOT)/spindynapy/build/*
	rm -rf $(PROJECT_ROOT)/spindynapy/__pycache__
	rm -rf $(PROJECT_ROOT)/spindynapy/spindynapy.egg-info

clean-all: clean-spindynapy
	rm -rf $(PROJECT_ROOT)/spindynapy/build-deps/*
	rm -rf $(PROJECT_ROOT)/.pytest_cache/
	rm -rf $(PROJECT_ROOT)/.mypy_cache
	rm -rf $(PROJECT_ROOT)/__pycache__

run:
	$(POETRY) run python main.py
