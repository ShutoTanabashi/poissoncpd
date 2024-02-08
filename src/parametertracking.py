"""
parametertracking
=================

変化点検出プログラム
ポアソン分布版

class
_____
ParamTrackPoisson

function
________
gen_outpath
"""
# 製作者
# 棚橋秀斗(shuuto570@outlook.com)

import sys
import os
import numpy as np
import pandas as pd
import ctypes

# dllの呼び出し
_multistep = ctypes.WinDLL("./multistep.dll")
_multistep.detect_change.restype = None
_multistep.detect_change.argtypes = [
    ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=1, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=1, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.intc, ndim=2, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C'),
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C')]

class ParamTrackPoisson:
    """
    変数列に対するポアソン分布の強度変化の追跡

    Atributes
    _________
    data : pandas.DataFrame
        ポアソン変数列を"人数"列に保持したデータフレーム
    T : int
        データ個数
    X : pandas.Series
        ポアソン変数列
    model : pandas.DataFrame
        ポアソン分布に基づいた確率モデル
            λ_i : 最尤強度
            τ_i : 変化点
    mle : pandas.Series
        最尤推定量
    AIC : pandas.Series
        赤池情報量基準
    L : pandas.Series
        対数尤度
    Bk : pandas.Series
        パラメータ数
    """
    
    def __init__(self, data):
        """
        Parameters
        __________
        data : pandas.DataFrame
            ポアソン変数列を"人数"列に保持したデータフレーム
        """
        self._data = data
        self._index = data.index
        self._T = ctypes.c_int(len(data))
        self._t1 = self.T + 1
        self._X = np.insert(data.loc[:, "人数"].to_numpy(dtype=np.intc), 0, 0)
        self._mle = np.full(self._t1, np.nan, dtype=np.double)
        self._tau = np.full((self.T, self._t1), -1, dtype=np.intc)
        self._AIC = np.full(self.T, np.nan, dtype=np.double)
        self._L = np.full(self.T, np.nan, dtype=np.double)
        self._D = np.full((self.T, self._t1), np.nan, dtype=np.double)
        self._Bk = np.full(self.T, 0, dtype=np.intc)
        self._model = self._detect()

    @property
    def data(self):
        df_data = pd.DataFrame([self.X, self.mle]).T
        name_X = self.X.name
        return df_data.astype({name_X: pd.Int64Dtype()})

    @property
    def T(self):
        return self._T.value
    
    @property
    def X(self):
        return pd.Series(np.delete(self._X, 0), 
                         index=self._index, name="人数")

    @property
    def model(self):
        return self._model
    
    @property
    def mle(self):
        return pd.Series(np.delete(self._mle, 0),
                         index=self._index, name="mle")
        
    @property
    def AIC(self):
        return pd.Series(self._AIC,
                         index=np.arange(0, self.T), name="AIC")
    
    @property
    def L(self):
        return pd.Series(self._L,
                         index=np.arange(0, self.T), name="対数尤度")
    
    @property
    def Bk(self):
        return pd.Series(self._Bk,
                         index=np.arange(0, self.T), name="パラメータ数")

    def _detect(self):
        """
        ポアソン分布に基づく母数変化の検出

        Returns
        _______
        model : pandas.DataFrame
            ポアソン分布に基づいた確率モデル
            λ_i : 最尤強度
            τ_i : 変化点

        Notes
        _____
        c言語により実装された変化点検出プログラムを利用
        以下の変数をポインタで受け渡した後書き換える
            self._Bk, self._tau, self._mle, self._L, self._AIC, self._D
        """
        # 最尤推定量計算
        _multistep.detect_change(self._T, ctypes.c_int(self._t1),
                                self._X, self._Bk, self._tau,
                                self._mle, self._L, self._AIC, self._D)
        # 変化点の検出
        # ステップ数
        lis_tau = []
        # 強度の最尤推定量
        lis_ll = []
        ll = np.nan
        for i, m in zip(reversed(list(self._index)), reversed(self._mle)):
            if ll != m:
                lis_tau.insert(0, i)
                lis_ll.insert(0, m)
                ll = m
        model = pd.DataFrame({'λ_i': lis_ll, 'τ_i': lis_tau},
                             index=np.arange(0, len(lis_tau), 1))
        model.index.name = "i"

        return model

    def to_excel(self, path):
        """
        Excelファイルに出力する

        Parameters
        __________
        path : str
            出力するExcelファイルのパス

        Notes
        _____
        MLE，model，AICの3つのシートを出力する
            MLEシート : 各ステップにおける強度の最尤推定量
            modelシート : 推定強度と変化点
            AICシート : 変化点個数毎のAIC
        """
        df_MLE = self._data.copy()
        df_MLE.insert(1, "最尤推定量", self.mle)
        df_AIC = pd.DataFrame([self.L, self.Bk, self.AIC]).T

        with pd.ExcelWriter(path) as writer:    # pylint: disable=abstract-class-instantiated
            # (上のコメントはpylint向けコマンドなので無視してよい)
            df_MLE.to_excel(writer, sheet_name="MLE")
            self.model.to_excel(writer, sheet_name="model")
            df_AIC.to_excel(writer, sheet_name="AIC")

    def internal_to_csv(self, path):
        """
        次の内部計算結果をcsvで出力
            - 対数尤度の動的計画法部分
            - 変化点推移の全データ
        
        Parameters
        __________
        path : str
            出力するcsvファイルのパス
            要拡張子
        
        Notes
        _____
        出力されるcsvファイルのファイル名は次の通り
            対数尤度の動的計画法部分 -> ファイルパスに"_D"が追加
            変化点推移の全データ -> ファイルパスに"_tau"が追加
        """
        [path_file, _] = os.path.splitext(path)

        df_D = pd.DataFrame(self._D,
                            index=np.arange(0, self.T),
                            columns=np.arange(0, self._t1))
        df_D.index.name = r"変化点個数\変数番号"
        df_D.to_csv(path_file + "_D.csv", encoding="shift-jis")

        df_tau = pd.DataFrame(self._tau,
                              index=np.arange(0, self.T),
                              columns=np.arange(0, self._t1))
        df_tau.index.name = r"変化点個数\変数番号"
        df_tau.to_csv(path_file + "_tau.csv", encoding="shift-jis")
    
    @staticmethod
    def read_excel(path):
        """
        Excelファイルの"乱数列"シートからParamTrackPoissonクラスのインスタンス生成

        Parameters
        __________
        path : str
            読み取るExcelのファイルパス
        
        Return
        ______
        ptp : ParamTrackPoisson
            Excel内のデータを用いたParamTrackPoissonインスタンス
        """
        # データ呼び出し
        data = pd.read_excel(path, sheet_name="乱数列", header=0, index_col=0)
        if ("人数" not in data.columns.values):
            raise IndexError("ポアソン変数列の見出しは 人数 にしてください")
        return ParamTrackPoisson(data)


def gen_outpath(inpath, folda="Output", option="_analyzed", ext="xlsx"):
    """
    pathを元にした出力ファイルのパス生成

    Parameters
    __________
    inpath : str
        入力ファイルパス
    folda : str, default "Output"
        出力先のフォルダ名
    option : str, default "_analyzed"
        ファイル名に追加する文字列
    ext : str, default "xlsx"
        出力ファイルの拡張子

    Returns
    _______
    outpath : str
        出力ファイルパス
    """
    [file, _] = os.path.splitext(os.path.basename(inpath))
    outpath = os.path.join(folda, file + option + "." + ext)
    return outpath

if __name__ == "__main__":
    # コマンドライン引数
    input = sys.argv[1]
    output = gen_outpath(input, folda=".")
    ptp = ParamTrackPoisson.read_excel(input)
    ptp.to_excel(output)
    print("Computation completed.")
