
# Пути
PROJECT_ROOT := $(PWD)

# Бинарные файлы
CMAKE := cmake
POETRY := poetry
PIP = $(POETRY) run pip
CLANG_FORMAT = clang-format

build-all:
	$(POETRY) install

build-spindynapy:
	$(PIP) install $(PROJECT_ROOT)/spindynapy

format:
	find ./spindynapy/core_src -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec $(CLANG_FORMAT) -style=file -i {} \;

clean-spindynapy:
	rm -rf $(PROJECT_ROOT)/spindynapy/lib/*
	rm -rf $(PROJECT_ROOT)/spindynapy/*.so
	rm -rf $(PROJECT_ROOT)/spindynapy/core/*
	rm -rf $(PROJECT_ROOT)/spindynapy/build/*
	rm -rf $(PROJECT_ROOT)/spindynapy/__pycache__
	rm -rf $(PROJECT_ROOT)/spindynapy/spindynapy.egg-info

run:
	$(POETRY) run python main.py
