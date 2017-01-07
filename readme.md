# [ode-ext](https://github.com/keichan100yen/ode-ext)

## A fast Ordinary Differential Equations solver using expression template for C++

ode-ext is a reference implementation of [my study](http://dspace.wul.waseda.ac.jp/dspace/bitstream/2065/40295/1/Honbun-6338.pdf) 
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
