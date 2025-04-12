# umag
Требования к системе:
- C++ compiler, CMake
- Python3.12 headers (`apt install python3-dev`)
- make utility + Poetry Python PM
- clang-format

Делится на два: python-lib `spindynapy` & `main.py` для конечного применения. Можно и `numpy` использовать (!!!)

## Run main.py
```sh
make build-all  # build library & configure poetry venv
make run  # run sample
```
