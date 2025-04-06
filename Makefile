
# Пути
PROJECT_ROOT := $(PWD)

# Бинарные файлы
CMAKE := cmake
POETRY := poetry
PIP = $(POETRY) run pip

build-all:
	$(POETRY) install

build-spindynapy:
	$(PIP) install $(PROJECT_ROOT)/spindynapy

clean-spindynapy:
	rm -rf $(PROJECT_ROOT)/spindynapy/lib/*
	rm -rf $(PROJECT_ROOT)/spindynapy/build/*
	rm -rf $(PROJECT_ROOT)/spindynapy/__pycache__
	rm -rf $(PROJECT_ROOT)/spindynapy/spindynapy.egg-info

run:
	$(POETRY) run python main.py
