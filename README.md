# ポアソン分布を想定した変化点検出

乱数列がポアソン分布に従うことを想定し、AICを用いて複数個の変化点を推定するプログラム

制作者: 棚橋秀斗[tanabashi@s.okayama-u.ac.jp]

## 注意事項

本プログラムはWindows環境での実行を想定して作成しています。
他OSでは実行できない可能性が高いです。

本プログラムを実行あるいは流用したことにより生じる一切の事象について製作者は責任を負いません。 自己責任でお願いいたします。

## 実行方法

まず`src`ディレクトリ内で`multisteplib_poisson.c`をDLLファイルにコンパイルします。
コンパイラには[`Clang`](https://clang.llvm.org)または[`GCC`](https://gcc.gnu.org)を推奨します。  
PowerShell（またはコマンドプロンプト）で`src`ディレクトリに移動したうえで以下のコマンドを実行してください。
[`GNU Make`](https://www.gnu.org/software/make/)向けに`src/makefile`を用意しています。

```powershell
# GNU Make が利用できる場合
make
# Clangでコンパイルする場合
clang -shared -O2 -o multistep.dll multisteplib_poisson.c
# GCCでコンパイルする場合
gcc -shared -O2 -o multistep.dll multisteplib_poisson.c
```

`src/parametertracking.py`で定義される`ParamTrackPoisson`クラスがポアソン分布に対する変化点検出法です。
以下のように`parametertracking.py`を実行することで、`src/testdata/testdata.xlsx`に記載されたデータを読み取って変化点検出を実行します。

```powershell
python .\parametertracking.py .\testdata\testdata.xlsx
```

変化点検出での解析結果が`testdata_analyzed.xlsx`として出力されます。
