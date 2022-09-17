#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>

#define		NUM_CITY		100	//都市の総数
#define		NUM_GENE		200		//遺伝子の総数

//関数宣言
void main_GA();
void flg_op(double v, long x[]);//前世代の最適経路
void flg_op1(long x[]);
void mutate_gene();
void mutate_gene_each(long i, long j);
//------localで使用する関数-------
void local_gene();//ローカルサーチ関数
void local_each(long p[]);
void local_sort(long i, long j, long p[]);
void output_gene(long i);
void select_gene();
void select_gene1();
//------selectで使用する関数------
void kakuno(long wk[]);
long hikaku(long a, long b);
//------------------------------
void kousa_gene();
//------kousaで使用する関数------
void ex_kousa(long p1[], long p2[], long mask[], long x);
void kousa_gene_pair(long p1[], long p2[]);
long num_where(long p1[], long p2);
//------------------------------
void sort_gene();
//------sortで使用する関数-------
void sort();
double full_kyori(long i);
double kyori(long i, long j);
void ex(long i, long j);
//-------------------------------
void mk_city();
void init_gene();
void init_gene_each(long i, long j);
long check_same(long i, long j);
void final_process();

//======================外部変数の定義====================================
//-----メイン関数で設定するパラメータ-----
long Num_select, Num_kousa; 
double P_mut, P_mut1, P_mut2;
long G_fin;
long G_disp_GA;
long Count;
long Limit;
long Top[NUM_CITY];
long F_local;

//-----main_GAの外部変数-----
long Gene[NUM_GENE][NUM_CITY];	//都市を回る経路

double City[NUM_CITY][2];		//都市の座標

//-----̧ファイルポインタ-----
FILE *Fp_city;
FILE *Fp_deb;

//メイン関数
int main(void){
	if((Fp_deb = fopen("local.dat", "w") ) == NULL){ printf("err for fopen\n"); exit(1); }//ファイルオープン
	F_local = 0;
//******************************************* parameter設定*****************************************************
	//==============GAの制御=============
	Num_select  = NUM_GENE/2;//選択によって上位Num_select個の遺伝子を残す(残りの遺伝子は,上位Num_selec

	
	Num_kousa	= NUM_GENE/4; //交叉を行う遺伝子の数
	P_mut  		= P_mut1;//突然変異を起こす確率(初期値)
    P_mut1      = 5.e-3;
    P_mut2      = P_mut1 * 10; //
	G_fin	    = 50000;			//繰り返し世代数
	Count       = 0; //突然変異確率操作用
	Limit = 100;

	//=============出力制御=============
	G_disp_GA		= G_fin/100;		//main_GAで世代数を画面に表示する確率
	
	//=============実行=================
	mk_city();			//都市の作成
	main_GA();			//GAの実行
}


void mk_city(){
	long i;
    long j;
	char line[128];

	if ( (Fp_city = fopen("city.dat", "r") ) == NULL ){ 
		printf("error in fopen\n"); 
		exit(1);
	}    
	fgets(line, 128, Fp_city);//1行目をlineに格納
	i = 0;
    while(fscanf(Fp_city, "%ld %lf %lf", &i, &City[i][0], &City[i][1]) != EOF)//都市の座標を格納
	i++;
    fclose(Fp_city);
    for (i = 0; i < NUM_CITY; i++){
        printf("i=%ld,x=%lf,y=%lf\n", i,City[i][0], City[i][1]);//読み込みチェック
    }
}



void main_GA(){
	long g;
	long i, j;
	double path, top_path = 0;
    long wk[NUM_CITY];
    
	init_gene();
	
	
	for (g = 0; ;g++){
		if (g % G_disp_GA == 0) printf("g = %ld in main_GA\n", g);
		if (g == G_fin){ final_process(); return; } 
	
		sort_gene();
     
        if(g != 0) flg_op(path, wk);//当然変異の確率の操作と、最適経路の総距離の比較のための関数
        
    
        for(i = 0; i< NUM_CITY; i++) wk[i] = Gene[0][i]; //最適経路

		path = full_kyori(0);//前世代の最適経路の総距離の保存

		if(g % 100 == 0) printf("%lf\n",full_kyori(0));//デバック用
        select_gene();
		if(g >= G_fin / 2) {
            select_gene1();//エリート選択
            Limit = 1000;//突然変異確率の操作count数
        }
		if(g >= G_fin-300) P_mut = P_mut1;
		kousa_gene();	
		mutate_gene();
		if(g != 0 && g % 100 == 0) local_gene();

 	
    }
}



void init_gene(){
	long i, j;
	for (i = 0; i < NUM_GENE; i++){
		for (j = 0; j < NUM_CITY; j++){
			init_gene_each(i, j);	
		}
	}
}


void init_gene_each(long i, long j){
	long n_big = 10000;
	long k;
	for (k = 0; k < n_big; k++){
		Gene[i][j] = rand()%NUM_CITY;		
		if (check_same(i, j) == 0) return;	
	}
	printf("err in init-gene-each\n"); exit(1);
}

long check_same(long i, long j){
	long k;
	for (k = 0; k < j; k++){
		if (Gene[i][j] == Gene[i][k]) return 1;
	}
	return 0;
}

void flg_op(double v, long x[]){
	double pre_path;
	double path;
	path = full_kyori(0);//現在の最適経路の総距離
	pre_path = v;//前世代の最適経路の総距離
	

	if(pre_path < path) {
		flg_op1(x);//前世代の総距離の方が短い場合、flg_op1を呼び出す
		Count = 1;//初期化
	}
    else if(pre_path == path) { //同じ総距離の場合
		
		Count++;//前回の解と一致している場合、countをインクリメントする.countが100以上の時、P_mutは高い値にする
    }

    else {//前回の解と一致してない場合、P_mutをデフォに戻す
        P_mut = P_mut1;//突然変異確率を元に戻す
        Count = 1;//初期化
    }
	if(Count == Limit) P_mut = P_mut2;//突然変異確率を上げる

}

void flg_op1(long x[]){
    long i, a;
	double pre_path;
	double path;

	a = rand() % NUM_GENE;
    for(i = 0; i < NUM_CITY; i++) Gene[0][i] = x[i]; //前世代の最適経路をランダムの場所に格納
}


void sort_gene(){
	double wk1, wk2;
    long i, j;
    for(i = 0; i < NUM_GENE - 1; i++){//i番目の遺伝子
        for(j = i + 1; j < NUM_GENE; j++){//一個次の遺伝子
            wk1 = full_kyori(i);//i番目の遺伝子の総距離
            wk2 = full_kyori(j);//j番目の遺伝子の総距離
            if(wk1 > wk2){//総距離の比較
                ex(i ,j);//交換する関数の呼び出し
            }
        }
    }
	
}

double full_kyori(long i){
    long j, n1, n2;
    double wk, sum;
    sum = 0.e0; //初期化
    for(j = 0; j < NUM_CITY ; j++){
    	n1 = Gene[i][j];
    	if(j != NUM_CITY-1) n2 = Gene[i][j + 1];
    	else n2 = Gene[i][0];
    	wk = kyori(n1, n2);//i番目の遺伝子のCity間の距離
    	sum += wk;//合計
    }
    return sum;
}

double kyori(long i, long j){
    double wk1, wk2;
    wk1 = City[i][0] - City[j][0];
    wk2 = City[i][1] - City[j][1];
    return sqrt(wk1 * wk1 + wk2 * wk2);//都市間の距離
}

void ex(long i, long j){
    double wk;
    long k;
    for(k = 0; k < NUM_CITY; k++){
        wk = Gene[i][k];//一時的に格納
    	Gene[i][k] = Gene[j][k];//i番目の遺伝子にｊ番目の遺伝子を格納
    	Gene[j][k] = wk;//j番目の遺伝子にi番目の遺伝子を格納
    }
}



void select_gene(){
	long i, j;
	long wk[300];
	long a, b;
	i =0;
	while(i < Num_select){
		a = rand() % NUM_GENE;
		b = rand() % NUM_GENE;
		if(a != b){
			wk[i] = hikaku(a, b); //遺伝子情報を格納
			i++;
		}
	 }
	kakuno(wk);
}

long hikaku(long a, long b){
	if(full_kyori(a) < full_kyori(b)) return a;
	else return b;
}

void kakuno(long wk[]){
	long wk1[300][NUM_CITY];
	long i, j;
	for(i = 0; i < Num_select; i++){
		for(j = 0; j < NUM_CITY; j++){
			wk1[i][j] = Gene[wk[i]][j];
		}
	}
	for(i = 0; i < NUM_GENE; i++){
		for(j = 0; j < NUM_CITY; j++){
			Gene[i][j] = wk1[i % Num_select][j];
		}
	}
}

void select_gene1(){
	long i, j;
    for(i = 0; i < NUM_GENE; i++){
        for(j = 0; j < NUM_CITY; j++) Gene[i][j] = Gene[i % (NUM_GENE / 2)][j]; //上位半分をコピー
    }
}


void kousa_gene(){
	long p1[NUM_CITY], p2[NUM_CITY];
    long i, j;
    if (Num_kousa % 2 != 0){ printf("err: kousa-gene\n"); exit(1); } //Num_kousa が偶数でなければエラー出力
    for(i = 0; i < Num_kousa/2; i++){
        for(j = 0; j < NUM_CITY; j++){
            p1[j] = Gene[2*i][j]; //遺伝子情報をp1にコピー
            p2[j] = Gene[2*i + 1][j]; //遺伝子情報をp2にコピー
		}
	
    	kousa_gene_pair(p1, p2); //p1 と p2 の間で実際に交叉を行う

    	for (j = 0; j < NUM_CITY; j++){
        	Gene[2 * i][j] = p1[j]; //交叉後の遺伝子情報をGeneにコピー
        	Gene[2*i + 1][j] = p2[j];
    	}
	}
}
    

void kousa_gene_pair(long p1[], long p2[]){
    long mask[NUM_CITY];
	long num_mask;
    long i, j, n = 0;
	num_mask = 0;
	mask[0] = 300;//nになり得ない適当な値を初期値として格納
    for(i = 0;; i++){
        n = num_where(p1, p2[n]); //P2[n]の都市がp1のどの要素番号か返す関数num_where
        if(mask[0] == n) break; //既出の要素番号が出たら終了
        mask[num_mask] = n; //mask配列に格納
        num_mask++;//mask配列の要素をインクリメント
    }
    ex_kousa(p1, p2, mask, num_mask);//交叉用の交換関数
}

long num_where(long p1[], long p2){
    long i;
    for(i = 0; i < NUM_CITY; i++){
		if(p1[i] == p2) return i;//都市の一致が終了条件
    }
    printf("errr");
	exit(1);
}

void ex_kousa(long p1[], long p2[], long mask[], long x){
    long wk, i;
    for(i = 0; i < x; i++){ //mask配列の部分だけ交換
    	wk = p1[mask[i]];
    	p1[mask[i]] = p2[mask[i]];
    	p2[mask[i]] = wk;
    }
}




void mutate_gene(){
	double u;
	long i, j;
	for (i = 0; i < NUM_GENE; i++){
		for (j = 0; j < NUM_CITY; j++){
			u = 1.e0 * rand() / RAND_MAX;
			if (u < P_mut) mutate_gene_each(i, j); //確率 P_mut で,Gene[i][j]の突然変異を行う
		}
	}
}

void mutate_gene_each(long i, long j){
	long n, m, k;
	n = Gene[i][j];
	for(;;){
		m = rand() % NUM_CITY;
		if (m != n) break;
	}

	Gene[i][j] = m;
	
	for(k = 0; k < NUM_CITY; k++){
		if(Gene[i][k] == m) break;
	}
	Gene[i][k] = n;
}


void final_process(){	
	long i;
	FILE *fp;
	if((fp = fopen("out_path.dat", "w") ) == NULL){ printf("err for fopen\n"); exit(1); }//ファイルオープン
	fprintf(fp, "Number_of_city, x_of_city, y_of_city\n");//第 1 行目の記入
	//最適解の経路を表示する
	
	sort_gene();//総距離が短い順にsortする
	for (i = 0; i < NUM_CITY; i++){ 
		fprintf(fp, "%ld %lf %lf\n",Gene[0][i],City[Gene[0][i]][0], City[Gene[0][i]][1]);//100mainGAで一番良い遺伝子を表示
	}
	fprintf(fp, "%lf\n",full_kyori(0));

}

void local_gene(){
    long i,j;
    long p[NUM_CITY];
    long p_bef[NUM_CITY];

    for(i = 0; i < NUM_GENE; i++){
        for(j = 0; j < NUM_CITY; j++){
            p[j] = Gene[i][j];	//遺伝子情報を1次元配列pに一時的に格納
            p_bef[j] = Gene[i][j];
        }

        local_each(p);

        for(j = 0; j < NUM_CITY; j++){
        	Gene[i][j] = p[j];//一次元配列pをGeneに戻す
        }
    if (F_local == 1){
   	    fprintf(Fp_deb, "before local_sort in local_each\n"); 
   	    for (j = 0; j < NUM_CITY; j++) fprintf(Fp_deb, "%ld", p_bef[j]);
   	    fprintf(Fp_deb, "\n");
   	    
		fprintf(Fp_deb, "after local_sort in local_each\n"); 
		output_gene(i);
		
	}
	
	
	}
	
}

void local_each(long p[]){
    long i, j, g, k;
    double vii, v_ij, vij, vjj;
	

    
	for(i = 0; i < NUM_CITY-2; i++){
        for(j = i+2; j < NUM_CITY; j++){
            
            vij = kyori(p[i], p[j]);
            v_ij =kyori(p[i + 1], p[(j + 1) % NUM_CITY]);
            vii = kyori(p[i], p[i + 1]);
            vjj = kyori(p[j], p[(j + 1) % NUM_CITY]);
            
            if((vij + v_ij) < (vii + vjj)) {
            	F_local = 1;
				local_sort( i+1, j, p);//localサーチの条件	
		    }
				
		}
	}		
}

void output_gene(long i){
	long j;
	for (j = 0; j < NUM_CITY; j++){
		fprintf(Fp_deb, "%ld", Gene[i][j]);
	}
}


void local_sort(long i, long j, long p[]){
    long q[NUM_CITY];
    long k, s;

    for(k = 0; k < i; k++) q[k] = p[k];

    s = 0;
    for(k = i; k < j+1; k++) {
        q[k] = p[j-s];
        s++;
    }

    for(k = j + 1; k < NUM_CITY; k++) q[k] = p[k];
    //qを作成
    //qをpに格納
	for(k = 0; k < NUM_CITY; k++) p[k] = q[k];//qをpに格納
}


