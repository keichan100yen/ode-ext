# [ode-ext](https://github.com/keichan100yen/ode-ext)

## A fast Ordinary Differential Equations solver using expression template for C++

ode-ext is a reference implementation of [my study](https://waseda.repo.nii.ac.jp/?action=repository_uri&item_id=18718&file_id=20&file_no=1) 
using [kv library](https://github.com/mskashi/kv).

## <font color="Red">Build and execute</font>

```bash
git clone https://github.com/keichan100yen/ode-ext
cd ode-ext
git submodule init; git submodule update
cd example; mkdir build; cd build; cmake ..; make
./ode-ext
```

## <font color="Red">Known issues</font>

* Nonlinear functions are not implemented using expression template.
