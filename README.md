# umag
Требования к системе:
- C++ compiler, CMake C++20 `ninja-build` `g++13`
- Python3.12 headers (`apt install python3-dev` + `which pyton` `setuptools`) 
- make utility + Poetry Python PM
- clang-format

Делится на два: python-lib `spindynapy` & `main.py` для конечного применения. Можно и `numpy` использовать (!!!)

## Run main.py
```sh
make build-all  # build library & configure poetry venv
make run  # run sample
```

```sh
pyenv global 3.12
pipx install poetry # sudo apt install -y libffi-dev tk-dev tcl-dev liblzma-dev
# apt install make
poetry config virtualenvs.in-project true
```

```sh
git submodule update --init --recursive
sudo apt install cmake make
sudo apt install cmake-doc ninja-build cmake-format
sudo apt install -y libstdc++6 libstdc++-12-dev
sudo apt install --reinstall g++-12 gcc-12 libstdc++-12-dev

sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install -y g++-13

sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100
```

OPENMP SIMPLE PARALLELING

simulate_many_steps releases GIL when simulate_one NOT

: 1745750733:0;sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 20