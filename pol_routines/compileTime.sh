set -euo pipefail

g++ -std=c++17 main2.cpp helperF.cpp integerMath.cpp interpAlgo.cpp polyMath.cpp int128g.cpp -I../header_FILES -march=native -O2 -finline-functions -funroll-loops -flto -DNDEBUG  -o time
