Build steps

1. Clone repo
```bash
git clone https://gitlab.com/Yzoni/bter-graph-generator && \
git clone https://github.com/google/googletest bter-graph-generator/test/lib
```
2. Create build directory
```bash
mkdir out && cd out
```

3. Make
```bash
cmake ../bter-graph-generator && make
```

4. Run
```bash
./bter_run
```