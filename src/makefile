# PT_Poissonプロジェクト makefile
# cソースコードによる.dll用

# 製作者 : 棚橋秀斗(shuuto570@outlook.com)

CC = clang
# CC = gcc

# Windows環境向け
ALL: multistep.dll

multistep.dll : multisteplib_poisson.c
	$(CC) -shared -O2 -o multistep.dll multisteplib_poisson.c
