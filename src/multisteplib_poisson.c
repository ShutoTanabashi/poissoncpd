/**********************************/
//
// 変化点検出モデル（ポアソン分布）
//
/**********************************/
// Encording utf-8

/* ヘッダ・ファイルの読み込み */
#include "math.h"
#include "float.h"

/* プロトタイプ宣言 */
double cal_mle(int [], int, int);
double cal_Dk(int [], int , int);
double cal_ll_poisson(int, int , int size, int [], int [][size],
                      double [][size]);
double log_fac_x(int, int []);
int cal_AIC(int, int size, double, int [], 
            double [], double [], double [][size]);
void set_mle(int, int size, int, int [], int [][size], double []);            
void detect_change(int , int size, int [], int [], int [][size],
                   double [], double [], double [], double [][size]);


double cal_mle(int X[], int start, int end){
    /*
	 * 最尤推定量計算
	 * 
	 *  X[size] ポアソン変数列, start 計算区間始点, end 計算区間終点
	 */
    double mle; // 最尤推定量
    int i;

    mle = 0;
    for (i = end; i > start; i--)
    {
        mle += X[i];
    }
    mle = mle / (double)(end - start);

    return mle;
}


double cal_Dk(int X[], int start, int end){
    /*
     * 対数尤度の動的計画法部分 $D_k(\tau_k)$ の計算
     * 
     * X[size] ポアソン変数列, start 計算区間始点, end 計算区間終点
     */
    
    double Dk; // 対数尤度の動的計画法部分
    double mle; //最尤推定量
    int i;

    Dk = 0.0;
    mle = cal_mle(X, start, end);

    for(i = end; i > start; i--){
        if(X[i] == 0){
            Dk -= mle;
        } else
        {
            Dk += X[i] * log(mle) - mle;
        }
    }
    
    return Dk;
}


double cal_ll_poisson(int k, int tau_k, int size, int X[], 
                      int TAU[][size], double D[][size]){
    /*
	 * 最尤推定量計算
	 * 
	 * k 変化点個数,
     * tau_k 変化点,
     * size 配列数,
     * X[size] ポアソン変数列, 
	 * TAU[size][size] 変化点履歴,
     * D[size][size] 対数尤度,
     * 
     * NOTE:D[][]はnanで初期化されているものとする
     *      動的計画法にはメモ化再帰関数を利用
     *      この関数ではD(メモ)のほかTAUを変更します．
	 */

    // test
    // printf("k:tau_k = %d:%d\n", k, tau_k);

    if (!isnan(D[k][tau_k]))
    {
        return D[k][tau_k];
    } else if (k <= 0)
    {
        D[k][tau_k] = cal_Dk(X, 0, tau_k);
        TAU[k][tau_k] = 0;
        
        return D[k][tau_k];
    }
    

    double ell1, ell2; // 対数尤度
    double mle;        // 最尤推定量
    double logmle;     //最尤推定についての対数

    // 初期化
    ell1 = -DBL_MAX;
    ell2 = -DBL_MAX;

    // 追加データにて強度が変化しているシナリオ
    if (k < tau_k)
    {
        ell1 = cal_ll_poisson(k - 1, tau_k - 1, size, X, TAU, D)
                + cal_Dk(X, tau_k - 1, tau_k);
    }

    // 追加データが以前の強度から変化していないシナリオ
    if (k < tau_k - 1 && !(k == 1 && tau_k < 2))
    {
        if (TAU[k][tau_k - 1] < 0)
        {
            (void)cal_ll_poisson(k, tau_k - 1, size, X, TAU, D);
        }

        ell2 = cal_ll_poisson(k - 1, TAU[k][tau_k - 1], size, X, TAU, D)
                + cal_Dk(X, TAU[k][tau_k - 1], tau_k);
    }

    // 最適な方をメモ
    if (ell1 >= ell2)
    {
        D[k][tau_k] = ell1;
        TAU[k][tau_k] = tau_k - 1;
    }
    else
    {
        D[k][tau_k] = ell2;
        TAU[k][tau_k] = TAU[k][tau_k - 1];
    }

    return D[k][tau_k];
}


double log_fac_x(int T, int X[]){
    /*
     * 対数尤度について各モデルで共通となる定数部分を計算
     * 
     * T 解析対象データの個数, X[T+1] ポアソン変数列
     */

    double C;   // 対数尤度の定数部
    int i, j; 

    C = 0.0;
    for (i = 1; i <= T; i++)
    {
        if (X[i] != 0)
        {
            for (j = X[i]; j > 1; j--)
            {
                C -= log(j);
            }
        }
    }
    
    return C;
}


int cal_AIC(int T, int size, double C, int Bk[], double L[], double AIC[], double D[][size]){
    /*
     * モデル毎の赤池情報量基準の計算し，それが最小値となるindexを返す
     * T 解析対象データの個数, size 配列数, C 対数尤度の定数部
	 * L[T] 対数尤度, Bk[T] パラメータ数, AIC[T] 赤池情報量基準 
	 * D[T][size] 対数尤度の動的計画法部分, 
     * 
     * NOTE:この関数ではL,Bk,AICの値を変更します．
     */

    double minAIC;  // 最小のAIC
    int K;          // minAICとなるモデルの番号
    int i;

    minAIC = DBL_MAX;
    for (i = 0; i < T; i++)
    {
        L[i] = D[i][T] + C;
        Bk[i] = 2 * i + 1;
        AIC[i] = -2.0 * L[i] + 2.0 * Bk[i];
        // printf("AIC[ %d ]:%lf\n", i, AIC[i]); //経過の出力

        if (minAIC > AIC[i])
        {
            minAIC = AIC[i];
            K = i;
        }
    }

    return K;
}


void set_mle(int T, int size, int K, 
             int X[], int TAU[][size], double mle[]){
    /*
     * 変数列それぞれの位置に対応した強度の最尤推定量を配列に書き込む
     * 
     * T 解析対象データの個数, size 配列数, K モデルの変化点数,
	 * X[size] ポアソン変数列, TAU[T][size] 変化点履歴, 
	 * mle[size] 強度の最尤推定量
     * 
     * NOTE:この関数ではmleの値を変更します．
     */

    double mle_ij;  // 強度の最尤推定量
    int i, j, k;

    j = T;
    // printf("%d	\n%d	",  K, j); //経過の出力
    for (i = K; i >= 0; i--)
    {
        // printf(u8"%d_%d to %d\n",i, j, TAU[i][j]); //経過の出力
        mle_ij = cal_mle(X, TAU[i][j], j);
        for (k = j; k > TAU[i][j]; k--)
        {
            mle[k] = mle_ij;
        }
        j = TAU[i][j];
        // printf("%d	", j); //経過の出力
    }
}


void detect_change(int T, int size, int X[], int Bk[], 
                   int TAU[][size],double mle[], double L[],
                   double AIC[], double D[][size]){
    /* 
	 * 変化点検出
	 * 
	 * T 解析対象データの個数, size 配列数,
	 * X[size] ポアソン変数列, TAU[T][size] 変化点履歴, 
	 * D[T][size] 対数尤度の動的計画法部分, mle[size] 強度の最尤推定量,
	 * L[T] 対数尤度, Bk[T] パラメータ数, AIC[T] 赤池情報量基準
	 */
    
    /* 変数宣言 */
    int K; // 最大尤度となる変化点個数
    double C; // 対数尤度の定数部

    /* カウント用変数 */
    int i, j, k;

    /* 変化点検出*/
    for (i = 0; i < T; i++)
    {
        // printf("Change num [%d]\n",i);
        (void)cal_ll_poisson(i, T, size, X, TAU, D);
    }

    /* 状態変化数の推定 */
    C = log_fac_x(T,X);
    K = cal_AIC(T, size, C, Bk, L, AIC, D);

    /* 最尤推定量計算 */
    set_mle(T, size, K, X, TAU, mle);
}
